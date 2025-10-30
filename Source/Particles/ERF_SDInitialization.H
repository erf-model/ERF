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

using MatVec = std::vector<std::unique_ptr<MaterialProperties>>;

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

        /*! set default values */
        virtual void setDefaults ( const amrex::Geometry&,
                                   const MatVec&,
                                   const MatVec& );

        /*! read super-droplets initialization parameters */
        virtual void readInputs ( const std::string&,
                                  const std::string&,
                                  const amrex::Geometry&,
                                  const MatVec&,
                                  const MatVec& );

        /*! print super-droplets initialization parameters to screen */
        virtual void printParameters ( const MatVec&, const MatVec& ) const;

        /*! Get a distribution with constant multiplicity */
        void getDistribution( amrex::Vector<amrex::Real>&,
                              const int,
                              const amrex::Real,
                              const std::string&,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              std::mt19937& ) const;

        /*! Get a distribution with sampled multiplicity */
        void getDistribution( amrex::Vector<amrex::Real>&,
                              amrex::Vector<amrex::Real>&,
                              const amrex::Real,
                              const int,
                              const amrex::Real,
                              const std::string&,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              const amrex::Real,
                              std::mt19937& ) const;

        /*! Compute the mass distribution */
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

        /*! Compute the mass distribution (sampled multiplicity) */
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

        /*! Compute the mass distribution */
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

        /*! Compute the mass distribution (sampled multiplicity) */
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

        /* return whether multiplicity is sampled or constant */
        [[nodiscard]] inline bool sampledMultiplicity() const
        {
            return (m_mult_type == SDMultiplicityType::sampled);
        }

        /* return volume of particle domain */
        [[nodiscard]] inline auto volume() const
        {
            amrex::Real vol = 0.0;
            if (m_type == SDInitShape::uniform) {
                vol = m_particle_domain.volume();
            } else if (m_type == SDInitShape::bubble) {
                const auto& radius = m_particle_domain.hi();
                vol = (4.0/3.0)*PI*radius[0]*radius[1]*radius[2];
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

        /*!< Injection rate (number of physical particles per meter^3 second) */
        amrex::Real m_inj_rate = 0;
        /*!< Injection rate (number of SDs per meter^3 second) */
        amrex::Real m_sd_inj_rate = 0;

        /*!< Initial number density of super-droplets */
        amrex::Real m_numdens_sd = -1;

        /*!< Injection domain velocity */
        amrex::Vector<amrex::Real> m_domain_vel = {0.0,0.0,0.0};

        /*!< Start time for injection */
        amrex::Real m_tstart = 0.0;

        /*!< Stop time for injection */
        amrex::Real m_tstop = 1.0e17;

        /*! update dt-based quantities */
        inline void updateDt (const amrex::Real a_dt)
        {
            this->m_numdens = m_inj_rate * a_dt;
            m_numdens_sd = (m_sd_inj_rate>0 ? std::max(m_sd_inj_rate*a_dt, 1.0) : -1);

            if (this->m_type == SDInitShape::uniform) {
                amrex::Vector<amrex::Real> lo = {0.0, 0.0, 0.0};;
                amrex::Vector<amrex::Real> hi = {0.0, 0.0, 0.0};;
                for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
                    lo[dir] = this->m_particle_domain.lo(dir) + m_domain_vel[dir] * a_dt;
                    hi[dir] = this->m_particle_domain.hi(dir) + m_domain_vel[dir] * a_dt;
                }
                this->m_particle_domain.setLo(lo);
                this->m_particle_domain.setHi(hi);
            } else if (m_type == SDInitShape::bubble) {
                amrex::Vector<amrex::Real> center = {0.0, 0.0, 0.0};;
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

        /*! read super-droplets injection parameters */
        void readInputs ( const std::string&,
                          const amrex::Geometry&,
                          const MatVec&,
                          const MatVec&,
                          const amrex::Real );

        /*! print super-droplets injection parameters to screen */
        void printParameters ( const MatVec&, const MatVec& ) const override;

        /*! Compute number of super-droplets per grid cell */
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
            amrex::Real num_par_per_cell = 0.0;
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

        /*! read super-droplets initialization parameters */
        void readInputs ( const std::string&,
                          const amrex::Geometry&,
                          const MatVec&,
                          const MatVec& );

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
            amrex::Real num_par_per_cell = 0.0;
            if (this->m_numdens >= 0) {
                num_par_per_cell = std::ceil(this->m_numdens*a_dv);
            } else {
                num_par_per_cell = 1;
            }
            return num_par_per_cell;
        }

};

#endif
