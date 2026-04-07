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

#include "ERF_Constants.H"
#include "ERF_MaterialProperties.H"

AMREX_ENUM(SDInitShape,
    uniform,
    bubble,
    null
);

/*! \brief List of super-droplet initializations */
namespace SupDropInit
{
    /*! Maximum number of vapour/condensate species */
    const int num_species_max = 10;

    /*! Maximum number of aerosols */
    const int num_aerosols_max = 8;

    const std::string attrib_init_const = "mass_constant";
    const std::string attrib_init_exp = "mass_exponential";
    const std::string attrib_init_lnr = "radius_log_normal";
    const std::string attrib_init_lnr_auto = "radius_lognormal_autorange";
}

AMREX_ENUM(SDMultiplicityType,
    constant, sampled
);

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
        std::vector<std::string> m_species_init_type;

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
        std::vector<std::string> m_aerosol_init_type;

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
                              const int a_np,
                              const amrex::Real a_density,
                              const std::string& a_init_type,
                              const amrex::Real a_mass_min,
                              const amrex::Real a_mass_max,
                              const amrex::Real a_mass_mean,
                              const amrex::Real a_radius_min,
                              const amrex::Real a_radius_max,
                              const amrex::Real a_radius_mean,
                              const amrex::Real a_radius_gstd,
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
                              const amrex::Real a_dV,
                              const int a_np,
                              const amrex::Real a_density,
                              const std::string& a_init_type,
                              const amrex::Real a_mass_min,
                              const amrex::Real a_mass_max,
                              const amrex::Real a_mass_mean,
                              const amrex::Real a_radius_min,
                              const amrex::Real a_radius_max,
                              const amrex::Real a_radius_mean,
                              const amrex::Real a_radius_gstd,
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
                vol = (amrex::Real(4.0)/three)*PI*radius[0]*radius[1]*radius[2];
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
         * one Number density based on injection rate and timestep
         * two Super-droplet number density based on SD injection rate
         * three Position of the injection domain based on domain velocity
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
         * one The super-droplet number density and cell volume, or
         * two A constant number per cell (m_ppc) if number density is not set
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
