#ifndef SDINIT_H_
#define SDINIT_H_

#include <string>
#include <vector>
#include <random>

#include <AMReX_ParmParse.H>
#include <AMReX_Enum.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Random.H>
#include <AMReX_GpuContainers.H>

#include "ERF_Constants.H"
#include "ERF_MaterialProperties.H"

AMREX_ENUM(SDInitShape,
    uniform,
    bubble,
    null
);

/*! \brief Distribution type enum for GPU compatibility */
AMREX_ENUM(SDDistributionType,
    mass_constant,
    mass_exponential,
    radius_log_normal,
    radius_lognormal_autorange
);

/*! \brief List of super-droplet initializations */
namespace SupDropInit
{
    /*! Maximum number of vapour/condensate species */
    const int num_species_max = 10;

    /*! Maximum number of aerosols */
    const int num_aerosols_max = 8;
}

AMREX_ENUM(SDMultiplicityType,
    constant, sampled
);

/*! \brief GPU-compatible structure holding distribution parameters */
struct SDDistributionParams {
    SDDistributionType dist_type;  /*!< Type of distribution */
    amrex::Real mass_min;          /*!< Minimum mass */
    amrex::Real mass_max;          /*!< Maximum mass */
    amrex::Real mass_mean;         /*!< Mean mass */
    amrex::Real radius_min;        /*!< Minimum radius */
    amrex::Real radius_max;        /*!< Maximum radius */
    amrex::Real radius_mean;       /*!< Mean radius (mu for log-normal) */
    amrex::Real radius_gstd;       /*!< Geometric std dev for radius */
    amrex::Real density;           /*!< Material density */
    amrex::Real numdens;           /*!< Number density for multiplicity */
    amrex::Real cell_volume;       /*!< Cell volume for multiplicity scaling */
    // Pre-computed values for log-normal CDF inversion (truncated distribution)
    amrex::Real cdf_min;           /*!< CDF at radius_min */
    amrex::Real cdf_max;           /*!< CDF at radius_max */
    amrex::Real sigma;             /*!< log(radius_gstd) */
    amrex::Real lnrng;             /*!< log(radius_max) - log(radius_min) */
    amrex::Real lnmin;             /*!< log(radius_min) or log(mass_min) */
    amrex::Real delta;             /*!< mass_mean - mass_min for exponential */
    int sampled_mult;              /*!< 1 if sampled multiplicity, 0 if constant */
};

/*! \brief Inverse error function approximation for GPU
 *  \param[in] x Input value in (-1, 1)
 *  \return Approximate inverse error function value
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real SD_erfinv_gpu(const amrex::Real x) {
    amrex::Real a = amrex::Real(0.147);
    amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real term = std::log(amrex::Real(1) - x * x + eps);
    amrex::Real p1 = amrex::Real(2) / (PI * a) + term / amrex::Real(2);
    amrex::Real p2 = term / a;
    amrex::Real sign = (x >= 0) ? amrex::Real(1) : amrex::Real(-1);
    return sign * std::sqrt(std::sqrt(p1 * p1 - p2) - p1);
}

/*! \brief Generate mass from distribution parameters on GPU
 *  \param[in] params Distribution parameters
 *  \param[in] engine Random engine
 *  \param[out] mult_contribution Contribution to multiplicity (added, not set)
 *  \return Generated mass value
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real SD_sample_mass_gpu(const SDDistributionParams& params,
                                const amrex::RandomEngine& engine,
                                amrex::Real& mult_contribution)
{
    using amrex::Random;
    amrex::Real mass = 0.0;
    amrex::Real u = Random(engine);

    switch (params.dist_type) {
        case SDDistributionType::mass_constant:
            mass = params.mass_mean;
            mult_contribution = u;  // Random [0,1] for later rescaling
            break;

        case SDDistributionType::mass_exponential: {
            if (params.sampled_mult) {
                // Sampled multiplicity: sample uniformly in log-space, weight by exponential PDF
                amrex::Real lnval = params.lnmin + u * params.lnrng;
                mass = std::exp(lnval);
                mult_contribution = (params.numdens * params.cell_volume) * std::exp(-mass / params.delta);
            } else {
                // Constant multiplicity: true exponential distribution via inverse transform
                // If U ~ Uniform(0,1), then -delta*ln(U) ~ Exponential(1/delta)
                mass = -params.delta * std::log(u) + params.mass_min;
                mult_contribution = 0.0;  // Not used for constant multiplicity
            }
            break;
        }

        case SDDistributionType::radius_log_normal:
        case SDDistributionType::radius_lognormal_autorange: {
            amrex::Real dry_r;
            if (params.sampled_mult) {
                // Sampled multiplicity: sample uniformly in log-space, weight by log-normal PDF
                amrex::Real lnval = params.lnmin + u * params.lnrng;
                dry_r = std::exp(lnval);
                // Log-normal PDF contribution for multiplicity weighting
                amrex::Real term = std::exp(-std::log(dry_r/params.radius_mean)*std::log(dry_r/params.radius_mean)
                                            /(amrex::Real(2)*params.sigma*params.sigma));
                mult_contribution = (params.numdens * params.cell_volume) / (params.sigma * std::sqrt(amrex::Real(2)*PI)) * term;
            } else {
                // Constant multiplicity: true log-normal via inverse CDF
                // Map u from [0,1] to [cdf_min, cdf_max] for truncated distribution
                amrex::Real u_trunc = params.cdf_min + u * (params.cdf_max - params.cdf_min);
                // Inverse CDF: r = mu * exp(sigma * sqrt(2) * erfinv(2*u - 1))
                amrex::Real z = SD_erfinv_gpu(amrex::Real(2) * u_trunc - amrex::Real(1)) * std::sqrt(amrex::Real(2));
                dry_r = params.radius_mean * std::exp(params.sigma * z);
                mult_contribution = 0.0;  // Not used for constant multiplicity
            }
            mass = four_thirds_pi * dry_r * dry_r * dry_r * params.density;
            break;
        }
    }
    return mass;
}

/*! \brief Super-droplets initial properties */
class SDInitProperties
{
    public:

        using MatVec = std::vector<std::unique_ptr<MaterialProperties>>;

        /*!< Initial number of super-droplets per cell */
        int m_ppc = 1;

        /*! Initial distribution type */
        SDInitShape m_type = SDInitShape::uniform;

        /*!< Initial number density (m^{-3}) of physical particles */
        amrex::Real m_numdens = -1;
        /*!< Maximum multiplicity */
        amrex::Real m_max_multiplicity = 1000000;

        /*! Minimum species mass */
        std::vector<amrex::Real> m_mass_species_min;
        /*! Maximum species mass */
        std::vector<amrex::Real> m_mass_species_max;
        /*! Mean species mass */
        std::vector<amrex::Real> m_mass_species_mean;
        /*! Minimum species dry radius */
        std::vector<amrex::Real> m_radius_species_min;
        /*! Maximum species dry radius */
        std::vector<amrex::Real> m_radius_species_max;
        /*! Mean species dry radius */
        std::vector<amrex::Real> m_radius_species_mean;
        /*! Standard deviation of species dry radius */
        std::vector<amrex::Real> m_radius_species_geom_std;
        /*! Initial distribution type for species */
        std::vector<SDDistributionType> m_species_init_type;

        /*! Minimum aerosol mass */
        std::vector<amrex::Real> m_mass_aerosol_min;
        /*! Maximum aerosol mass */
        std::vector<amrex::Real> m_mass_aerosol_max;
        /*! Mean aerosol mass */
        std::vector<amrex::Real> m_mass_aerosol_mean;
        /*! Minimum aerosol dry radius */
        std::vector<amrex::Real> m_radius_aerosol_min;
        /*! Maximum aerosol dry radius */
        std::vector<amrex::Real> m_radius_aerosol_max;
        /*! Mean aerosol dry radius */
        std::vector<amrex::Real> m_radius_aerosol_mean;
        /*! Standard deviation of aerosol dry radius */
        std::vector<amrex::Real> m_radius_aerosol_geom_std;
        /*! Initial distribution type for aerosol */
        std::vector<SDDistributionType> m_aerosol_init_type;

        /*! Number of species */
        int m_num_species = 0;
        /*! Number of aerosols */
        int m_num_aerosols = 0;

        /*! box within which to place particles*/
        amrex::RealBox m_particle_domain;
        /*! Particle distribution shape parameters
        Box shape, p1 -> lo end, p2 -> hi end
        Bubble shape, p1 -> center, p2 -> radius */
        amrex::Vector<amrex::Real> m_init_particle_p1;
        amrex::Vector<amrex::Real> m_init_particle_p2;

        /*! multiplicity type */
        SDMultiplicityType m_mult_type;

        /*! \brief Set default values for initialization parameters
         * \param[in] a_geom Simulation geometry information
         * \param[in] a_species_mat Vector of species material properties
         * \param[in] a_aerosol_mat Vector of aerosol material properties
         */
        virtual void setDefaults ( const amrex::Geometry& a_geom,
                                   const MatVec& a_species_mat,
                                   const MatVec& a_aerosol_mat );

        /*! \brief Read super-droplets initialization parameters from input file
         * \param[in] a_prefix Prefix for parameter parser
         * \param[in] a_key Key identifier for initialization parameters
         * \param[in] a_geom Simulation geometry information
         * \param[in] a_species_mat Vector of species material properties
         * \param[in] a_aerosol_mat Vector of aerosol material properties
         */
        virtual void readInputs ( const std::string& a_prefix,
                                  const std::string& a_key,
                                  const amrex::Geometry& a_geom,
                                  const MatVec& a_species_mat,
                                  const MatVec& a_aerosol_mat );

        /*! \brief Print super-droplets initialization parameters to screen
         * \param[in] a_species_mat Vector of species material properties
         * \param[in] a_aerosol_mat Vector of aerosol material properties
         */
        virtual void printParameters ( const MatVec& a_species_mat,
                                       const MatVec& a_aerosol_mat ) const;

        /*! \brief Get a distribution with constant multiplicity
         * \param[out] a_mass Output vector of particle masses
         * \param[in] a_np Number of particles
         * \param[in] a_density Density of the particle material
         * \param[in] a_init_type Type of initialization distribution
         * \param[in] a_mass_min Minimum mass value
         * \param[in] a_mass_max Maximum mass value
         * \param[in] a_mass_mean Mean mass value
         * \param[in] a_radius_min Minimum radius value
         * \param[in] a_radius_max Maximum radius value
         * \param[in] a_radius_mean Mean radius value
         * \param[in] a_radius_gstd Geometric standard deviation for radius
         * \param[in,out] a_rng Random number generator
         */
        void getDistribution( amrex::Vector<amrex::Real>& a_mass,
                              int a_np,
                              amrex::Real a_density,
                              SDDistributionType a_init_type,
                              amrex::Real a_mass_min,
                              amrex::Real a_mass_max,
                              amrex::Real a_mass_mean,
                              amrex::Real a_radius_min,
                              amrex::Real a_radius_max,
                              amrex::Real a_radius_mean,
                              amrex::Real a_radius_gstd,
                              std::mt19937& a_rng ) const;

        /*! \brief Get a distribution with sampled multiplicity
         * \param[out] a_mass Output vector of particle masses
         * \param[out] a_mult Output vector of particle multiplicities
         * \param[in] a_dV Cell volume
         * \param[in] a_np Number of particles
         * \param[in] a_density Density of the particle material
         * \param[in] a_init_type Type of initialization distribution
         * \param[in] a_mass_min Minimum mass value
         * \param[in] a_mass_max Maximum mass value
         * \param[in] a_mass_mean Mean mass value
         * \param[in] a_radius_min Minimum radius value
         * \param[in] a_radius_max Maximum radius value
         * \param[in] a_radius_mean Mean radius value
         * \param[in] a_radius_gstd Geometric standard deviation for radius
         * \param[in,out] a_rng Random number generator
         */
        void getDistribution( amrex::Vector<amrex::Real>& a_mass,
                              amrex::Vector<amrex::Real>& a_mult,
                              amrex::Real a_dV,
                              int a_np,
                              amrex::Real a_density,
                              SDDistributionType a_init_type,
                              amrex::Real a_mass_min,
                              amrex::Real a_mass_max,
                              amrex::Real a_mass_mean,
                              amrex::Real a_radius_min,
                              amrex::Real a_radius_max,
                              amrex::Real a_radius_mean,
                              amrex::Real a_radius_gstd,
                              std::mt19937& a_rng ) const;

        /*! \brief Compute the aerosol mass distribution
         * \param[out] a_mass Output vector of aerosol masses
         * \param[in] a_idx Aerosol species index
         * \param[in] a_np Number of particles
         * \param[in] a_density Density of the aerosol material
         * \param[in,out] a_rng Random number generator
         */
        void getAerosolDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                      const int a_idx,
                                      const int a_np,
                                      const amrex::Real a_density,
                                      std::mt19937& a_rng ) const
        {
            getDistribution( a_mass,
                             a_np,
                             a_density,
                             m_aerosol_init_type[a_idx],
                             m_mass_aerosol_min[a_idx],
                             m_mass_aerosol_max[a_idx],
                             m_mass_aerosol_mean[a_idx],
                             m_radius_aerosol_min[a_idx],
                             m_radius_aerosol_max[a_idx],
                             m_radius_aerosol_mean[a_idx],
                             m_radius_aerosol_geom_std[a_idx],
                             a_rng);
        }

        /*! \brief Compute the aerosol mass distribution with sampled multiplicity
         * \param[out] a_mass Output vector of aerosol masses
         * \param[out] a_mult Output vector of particle multiplicities
         * \param[in] a_dV Cell volume for scaling multiplicity
         * \param[in] a_idx Aerosol species index
         * \param[in] a_np Number of particles to generate
         * \param[in] a_density Density of the aerosol material
         * \param[in,out] a_rng Random number generator
         */
        void getAerosolDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                      amrex::Vector<amrex::Real>& a_mult,
                                      amrex::Real a_dV,
                                      int a_idx,
                                      int a_np,
                                      amrex::Real a_density,
                                      std::mt19937& a_rng ) const
        {
            getDistribution( a_mass,
                             a_mult,
                             a_dV,
                             a_np,
                             a_density,
                             m_aerosol_init_type[a_idx],
                             m_mass_aerosol_min[a_idx],
                             m_mass_aerosol_max[a_idx],
                             m_mass_aerosol_mean[a_idx],
                             m_radius_aerosol_min[a_idx],
                             m_radius_aerosol_max[a_idx],
                             m_radius_aerosol_mean[a_idx],
                             m_radius_aerosol_geom_std[a_idx],
                             a_rng);
        }

        /*! \brief Compute the species mass distribution with constant multiplicity
         * \param[out] a_mass Output vector of species masses
         * \param[in] a_idx Species index
         * \param[in] a_np Number of particles to generate
         * \param[in] a_density Density of the species material
         * \param[in,out] a_rng Random number generator
         */
        void getSpeciesDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                      const int a_idx,
                                      const int a_np,
                                      const amrex::Real a_density,
                                      std::mt19937& a_rng ) const
        {
            getDistribution( a_mass,
                             a_np,
                             a_density,
                             m_species_init_type[a_idx],
                             m_mass_species_min[a_idx],
                             m_mass_species_max[a_idx],
                             m_mass_species_mean[a_idx],
                             m_radius_species_min[a_idx],
                             m_radius_species_max[a_idx],
                             m_radius_species_mean[a_idx],
                             m_radius_species_geom_std[a_idx],
                             a_rng);
        }

        /*! \brief Compute the species mass distribution with sampled multiplicity
         * \param[out] a_mass Output vector of species masses
         * \param[out] a_mult Output vector of particle multiplicities
         * \param[in] a_dV Cell volume for scaling multiplicity
         * \param[in] a_idx Species index
         * \param[in] a_np Number of particles to generate
         * \param[in] a_density Density of the species material
         * \param[in,out] a_rng Random number generator
         */
        void getSpeciesDistribution ( amrex::Vector<amrex::Real>& a_mass,
                                      amrex::Vector<amrex::Real>& a_mult,
                                      amrex::Real a_dV,
                                      int a_idx,
                                      int a_np,
                                      amrex::Real a_density,
                                      std::mt19937& a_rng ) const
        {
            getDistribution( a_mass,
                             a_mult,
                             a_dV,
                             a_np,
                             a_density,
                             m_species_init_type[a_idx],
                             m_mass_species_min[a_idx],
                             m_mass_species_max[a_idx],
                             m_mass_species_mean[a_idx],
                             m_radius_species_min[a_idx],
                             m_radius_species_max[a_idx],
                             m_radius_species_mean[a_idx],
                             m_radius_species_geom_std[a_idx],
                             a_rng);
        }

        /*! \brief Get GPU-compatible distribution parameters for a species
         * \param[in] a_idx Species index
         * \param[in] a_density Density of the species material
         * \param[in] a_cell_volume Cell volume for multiplicity scaling
         * \param[in] a_sampled_mult Whether using sampled multiplicity mode
         * \return SDDistributionParams structure for GPU use
         */
        [[nodiscard]] SDDistributionParams getSpeciesDistParams(
            int a_idx,
            amrex::Real a_density,
            amrex::Real a_cell_volume,
            bool a_sampled_mult) const
        {
            return makeDistributionParams(
                m_species_init_type[a_idx],
                m_mass_species_min[a_idx],
                m_mass_species_max[a_idx],
                m_mass_species_mean[a_idx],
                m_radius_species_min[a_idx],
                m_radius_species_max[a_idx],
                m_radius_species_mean[a_idx],
                m_radius_species_geom_std[a_idx],
                a_density,
                a_cell_volume,
                a_sampled_mult);
        }

        /*! \brief Get GPU-compatible distribution parameters for an aerosol
         * \param[in] a_idx Aerosol species index
         * \param[in] a_density Density of the aerosol material
         * \param[in] a_cell_volume Cell volume for multiplicity scaling
         * \param[in] a_sampled_mult Whether using sampled multiplicity mode
         * \return SDDistributionParams structure for GPU use
         */
        [[nodiscard]] SDDistributionParams getAerosolDistParams(
            int a_idx,
            amrex::Real a_density,
            amrex::Real a_cell_volume,
            bool a_sampled_mult) const
        {
            return makeDistributionParams(
                m_aerosol_init_type[a_idx],
                m_mass_aerosol_min[a_idx],
                m_mass_aerosol_max[a_idx],
                m_mass_aerosol_mean[a_idx],
                m_radius_aerosol_min[a_idx],
                m_radius_aerosol_max[a_idx],
                m_radius_aerosol_mean[a_idx],
                m_radius_aerosol_geom_std[a_idx],
                a_density,
                a_cell_volume,
                a_sampled_mult);
        }

        /*! \brief Create GPU-compatible distribution parameters structure
         * \param[in] a_init_type Distribution type
         * \param[in] a_mass_min Minimum mass
         * \param[in] a_mass_max Maximum mass
         * \param[in] a_mass_mean Mean mass
         * \param[in] a_radius_min Minimum radius
         * \param[in] a_radius_max Maximum radius
         * \param[in] a_radius_mean Mean radius
         * \param[in] a_radius_gstd Geometric standard deviation
         * \param[in] a_density Material density
         * \param[in] a_cell_volume Cell volume
         * \param[in] a_sampled_mult Whether using sampled multiplicity mode
         * \return SDDistributionParams structure
         */
        [[nodiscard]] SDDistributionParams makeDistributionParams(
            SDDistributionType a_init_type,
            amrex::Real a_mass_min,
            amrex::Real a_mass_max,
            amrex::Real a_mass_mean,
            amrex::Real a_radius_min,
            amrex::Real a_radius_max,
            amrex::Real a_radius_mean,
            amrex::Real a_radius_gstd,
            amrex::Real a_density,
            amrex::Real a_cell_volume,
            bool a_sampled_mult) const;

        /*! \brief Determine whether multiplicity is sampled or constant
         * \return True if multiplicity is sampled, false if constant
         */
        [[nodiscard]] inline bool sampledMultiplicity() const
        {
            return (m_mult_type == SDMultiplicityType::sampled);
        }

        /*! \brief Calculate the volume of the particle domain
         *
         * For uniform distribution type, returns the box volume.
         * For bubble distribution type, returns the volume of the bubble.
         *
         * \return Volume of the particle domain in cubic meters
         */
        [[nodiscard]] inline amrex::Real volume() const
        {
            amrex::Real vol = zero;
            if (m_type == SDInitShape::uniform) {
                vol = m_particle_domain.volume();
            } else if (m_type == SDInitShape::bubble) {
                const auto& radius = m_particle_domain.hi();
                vol = four_thirds_pi*radius[0]*radius[1]*radius[2];
            }
            return vol;
        }

        virtual int numSDPerCell (const amrex::Real) const = 0;
        virtual amrex::Real numParticlesPerCell (const amrex::Real) const = 0;
};

/*! \brief Super-droplets initialization structure */
class SDInjection : public SDInitProperties
{
    public:

        virtual ~SDInjection() = default;

        /*!< Injection rate (number of physical particles per meter^3 second) */
        amrex::Real m_inj_rate = 0;
        /*!< Injection rate (number of SDs per meter^3 second) */
        amrex::Real m_sd_inj_rate = 0;

        /*!< Initial number density of super-droplets */
        amrex::Real m_numdens_sd = -1;

        /*!< Injection domain velocity */
        amrex::Vector<amrex::Real> m_domain_vel = {zero,zero,zero};

        /*!< Start time for injection */
        amrex::Real m_tstart = zero;

        /*!< Stop time for injection */
        amrex::Real m_tstop = amrex::Real(1.0e17);

        /*! \brief Update time-dependent quantities for particle injection
         *
         * This function updates:
         * 1. Number density based on injection rate and timestep
         * 2. Super-droplet number density based on SD injection rate
         * 3. Position of the injection domain based on domain velocity
         *
         * \param[in] a_dt Timestep size in seconds
         */
        inline void updateDt (const amrex::Real a_dt)
        {
            this->m_numdens = m_inj_rate * a_dt;
            m_numdens_sd = (m_sd_inj_rate>0 ? std::max(m_sd_inj_rate*a_dt, one) : -1);

            if (this->m_type == SDInitShape::uniform) {
                amrex::Vector<amrex::Real> lo = {zero, zero, zero};
                amrex::Vector<amrex::Real> hi = {zero, zero, zero};
                for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                    lo[dir] = this->m_particle_domain.lo(dir) + m_domain_vel[dir] * a_dt;
                    hi[dir] = this->m_particle_domain.hi(dir) + m_domain_vel[dir] * a_dt;
                }
                this->m_particle_domain.setLo(lo);
                this->m_particle_domain.setHi(hi);
            } else if (m_type == SDInitShape::bubble) {
                amrex::Vector<amrex::Real> center = {zero, zero, zero};
                for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                    center[dir] = this->m_particle_domain.lo(dir) + m_domain_vel[dir] * a_dt;
                }
                this->m_particle_domain.setLo(center);
            }
        }

        /*! read super-droplets injection parameters */
        void readInputs ( const std::string&,
                          const std::string&,
                          const amrex::Geometry&,
                          const MatVec&,
                          const MatVec& ) override
        {
            amrex::Abort("SDInjection::readInputs(): Do not use this interface");
        }

        /*! \brief Read super-droplet injection parameters from input file
         *
         * This function reads injection-specific parameters like injection rate,
         * timing parameters, and domain velocity.
         *
         * \param[in] a_prefix Prefix for parameter parser
         * \param[in] a_geom Simulation geometry information
         * \param[in] a_species_mat Vector of species material properties
         * \param[in] a_aerosol_mat Vector of aerosol material properties
         * \param[in] a_dt Current timestep size
         */
        void readInputs ( const std::string& a_prefix,
                          const amrex::Geometry& a_geom,
                          const MatVec& a_species_mat,
                          const MatVec& a_aerosol_mat,
                          const amrex::Real a_dt );

        /*! print super-droplets injection parameters to screen */
        void printParameters ( const MatVec&, const MatVec& ) const override;

        /*! \brief Compute number of super-droplets to inject per grid cell
         *
         * Determines the number of super-droplets to inject based on either:
         * 1. The super-droplet number density and cell volume, or
         * 2. A constant number per cell (m_ppc) if number density is not set
         *
         * \param[in] a_dv Volume of the grid cell
         * \return Number of super-droplets to inject in this cell
         */
        [[nodiscard]] inline int numSDPerCell (const amrex::Real a_dv) const override
        {
            int num_sd_per_cell = 0;
            if (m_numdens_sd >= 0) {
                num_sd_per_cell = static_cast<int>(std::ceil(m_numdens_sd*a_dv));
            } else {
                num_sd_per_cell = this->m_ppc;
            }
            return num_sd_per_cell;
        }

        /*! Compute number of physical particles per grid cell */
        [[nodiscard]] inline amrex::Real numParticlesPerCell (const amrex::Real a_dv) const override
        {
            amrex::Real num_par_per_cell = zero;
            if (this->m_numdens >= 0) {
                num_par_per_cell = std::ceil(this->m_numdens*a_dv);
            } else {
                num_par_per_cell = 1;
            }
            return num_par_per_cell;
        }

};

/*! \brief Super-droplets initialization structure */
class SDInitialization : public SDInitProperties
{
    public:

        virtual ~SDInitialization() = default;

        /*!< Initial number density of super-droplets */
        amrex::Real m_numdens_sd_init = -1;

        /*! read super-droplets initialization parameters */
        void readInputs ( const std::string&,
                          const std::string&,
                          const amrex::Geometry&,
                          const MatVec&,
                          const MatVec& ) override
        {
            amrex::Abort("SDInjection::readInputs(): Do not use this interface");
        }

        /*! \brief Read super-droplet initialization parameters from input file
         *
         * This function reads initialization-specific parameters including
         * distribution types, particle counts, and sizing parameters.
         *
         * \param[in] a_prefix Prefix for parameter parser
         * \param[in] a_geom Simulation geometry information
         * \param[in] a_species_mat Vector of species material properties
         * \param[in] a_aerosol_mat Vector of aerosol material properties
         */
        void readInputs ( const std::string& a_prefix,
                          const amrex::Geometry& a_geom,
                          const MatVec& a_species_mat,
                          const MatVec& a_aerosol_mat );

        /*! print super-droplets initialization parameters to screen */
        void printParameters ( const MatVec&, const MatVec& ) const override;

        /*! Compute number of super-droplets per grid cell */
        [[nodiscard]] inline int numSDPerCell (const amrex::Real a_dv) const override
        {
            int num_sd_per_cell = 0;
            if (m_numdens_sd_init >= 0) {
                num_sd_per_cell = static_cast<int>(std::ceil(m_numdens_sd_init*a_dv));
            } else {
                num_sd_per_cell = this->m_ppc;
            }
            return num_sd_per_cell;
        }

        /*! Compute number of physical particles per grid cell */
        [[nodiscard]] inline amrex::Real numParticlesPerCell (const amrex::Real a_dv) const override
        {
            amrex::Real num_par_per_cell = zero;
            if (this->m_numdens >= 0) {
                num_par_per_cell = std::ceil(this->m_numdens*a_dv);
            } else {
                num_par_per_cell = 1;
            }
            return num_par_per_cell;
        }

};

#endif
