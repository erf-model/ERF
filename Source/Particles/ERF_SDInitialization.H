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

#include "ERF_MaterialProperties.H"

/*! \brief List of super-droplet initializations */
namespace SupDropInit
{
    /*! Uniform distribution throughout domain */
    const std::string init_uniform = "uniform";
    /*! Uniform distribution inside a bubble */
    const std::string init_bubble = "bubble";
    /*! No particles are initialized; nucleation generates them */
    const std::string init_null = "null";

    /*! Maximum number of aerosols */
    const int num_aerosols_max = 8;

    const std::string attrib_init_const = "mass_constant";
    const std::string attrib_init_exp = "mass_exponential";
    const std::string attrib_init_lnr = "radius_log_normal";
}

AMREX_ENUM(SDMultiplicityType,
    constant, sampled
);

/*! \brief Super-droplets initialization structure */
struct SDInitialization
{

    /*! Initial distribution type */
    std::string m_type = SupDropInit::init_uniform;

    /*!< Maximum multiplicity */
    amrex::Real m_max_multiplicity = 1000000;

    /*!< Initial number density (m^{-3}) of physical particles */
    amrex::Real m_numdens_init = -1;
    /*!< Initial number density of super-droplets */
    amrex::Real m_numdens_sd_init = -1;
    /*!< Initial number of super-droplets per cell */
    int m_ppc_init = 1;

    /*! type of condensate initial distribution */
    std::string m_condensate_init_type = SupDropInit::attrib_init_const;
    /*! Initial mean condensate mass */
    amrex::Real m_mass_condensate_mean = 4.1887902e-42; // => radius = 1e-15 [m];
    /*! Initial minimum condensate mass */
    amrex::Real m_mass_condensate_min = 4.1887902e-45; // => radius = 1e-16 [m];
    /*! Initial minimum condensate radius */
    amrex::Real m_radius_condensate_min = 1e-15;
    /*! Initial maximum condensate radius */
    amrex::Real m_radius_condensate_max = 1e-15;

    /*!< Number of seed particles per cell */
    int m_ppc_seed = 0;
    /*!< Condensate mass of seed */
    amrex::Real m_seed_mass = 0.0;

    /*! Minimum aerosol mass */
    std::vector<amrex::Real> m_mass_aerosol_min;
    /*! Maximum aerosol mass */
    std::vector<amrex::Real> m_mass_aerosol_max;
    /*! Mean aerosol mass */
    std::vector<amrex::Real> m_mass_aerosol_mean;
    /*! Minimum aerosol dry radius */
    std::vector<amrex::Real> m_radius_aerosol_min;
    /*! Mean aerosol dry radius */
    std::vector<amrex::Real> m_radius_aerosol_mean;
    /*! Maximum aerosol dry radius */
    std::vector<amrex::Real> m_radius_aerosol_max;
    /*! Standard deviation of aerosol dry radius */
    std::vector<amrex::Real> m_radius_aerosol_geom_std;
    /*! Initial distribution type for aerosol */
    std::vector<std::string> m_aerosol_init_type;

    /*! Number of aerosols */
    int m_num_aerosols = 0;

    /*! box within which to place particles*/
    amrex::RealBox m_init_particle_box;

    /*! multiplicity type */
    SDMultiplicityType m_mult_type;

    /*! set default values */
    void setDefaults ( const std::vector<std::unique_ptr<MaterialProperties>>& );

    /*! read super-droplets initialization parameters */
    void readInputs ( const std::string&,
                      const amrex::Geometry&,
                      const std::vector<std::unique_ptr<MaterialProperties>>& );

    /*! print super-droplets initialization parameters to screen */
    void printParameters ( const std::vector<std::unique_ptr<MaterialProperties>>& ) const;

    /*! Compute the aerosol mass distribution */
    void getAerosolDistribution ( amrex::Vector<amrex::Real>&, int, int, amrex::Real, std::mt19937& ) const;

    /*! Compute the aerosol mass distribution (sampled multiplicity) */
    void getAerosolDistribution ( amrex::Vector<amrex::Real>&, amrex::Vector<amrex::Real>&, int, int, amrex::Real, std::mt19937& ) const;

    /*! Compute the condensate mass distribution */
    void getCondensateDistribution ( amrex::Vector<amrex::Real>&, int, amrex::Real, std::mt19937& ) const;

    /*! Compute number of super-droplets per grid cell */
    [[nodiscard]] inline int numSDPerCell (const amrex::Real a_dv) const
    {
        int num_sd_per_cell = 0;
        if (m_numdens_sd_init >= 0) {
            num_sd_per_cell = static_cast<int>(std::ceil(m_numdens_sd_init*a_dv));
        } else {
            num_sd_per_cell = m_ppc_init;
        }
        return num_sd_per_cell;
    }

    /*! Compute number of physical particles per grid cell */
    [[nodiscard]] inline amrex::Real numParticlesPerCell (const amrex::Real a_dv) const
    {
        amrex::Real num_par_per_cell = 0.0;
        if (m_numdens_init >= 0) {
            num_par_per_cell = std::ceil(m_numdens_init*a_dv);
        } else {
            num_par_per_cell = 1;
        }
        return num_par_per_cell;
    }

    /* return whether multiplicity is sampled or constant */
    [[nodiscard]] inline bool sampledMultiplicity() const
    {
        return (m_mult_type == SDMultiplicityType::sampled);
    }

};

#endif
