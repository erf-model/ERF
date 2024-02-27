#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include "MaterialProperties.H"
#include "ERFPC.H"

/*! \brief Real-type array-of-struct attributes of a super-droplet */
using SuperDropletsRealIdxAoS = ERFParticlesRealIdxAoS;

/*! \brief Int-type array-of-struct attributes of a super-droplet */
using SuperDropletsIntIdxAoS = ERFParticlesIntIdxAoS;

/*! \brief Real-type struct-of-array attributes of a super-droplet */
using SuperDropletsRealIdxSoA = ERFParticlesRealIdxSoA;

/*! \brief Real-type struct-of-array attributes added during runtime */
struct SuperDropletsRealIdxSoA_RT
{
    enum {
        radius = 0,     /*!< Radius of physical particles */
        multiplicity,   /*!< Number of droplets that this super-droplet represents
                             (Real type to handle large values) */
        sol_mass,       /*!< Solute mass */
        sd_mass,        /*!< Superdroplet mass */
        ncomps
    };
};

/*! \brief Int-type struct-of-array attributes of a super-droplet */
using SuperDropletsIntIdxSoA = ERFParticlesIntIdxSoA;

/*! \brief List of super-droplet initializations */
namespace SuperDropletInitializations
{
    /*! Uniform distribution throughout domain */
    const std::string init_uniform = "uniform";
    /*! No particles are initialized; nucleation generates them */
    const std::string init_null = "null";

}

class SuperDropletPC : public ERFPC
{
    public:

        /*! \brief Constructor */
        SuperDropletPC (  amrex::ParGDBBase* a_gdb,
                          const std::string& a_vap_mat,
                          const std::vector<std::string>& a_aerosol_mat,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_gdb, a_name)
        {
            define(a_vap_mat, a_aerosol_mat);
        }

        /*! \brief Constructor */
        SuperDropletPC (  const amrex::Geometry& a_geom,
                          const amrex::DistributionMapping& a_dmap,
                          const amrex::BoxArray& a_ba,
                          const std::string& a_vap_mat,
                          const std::vector<std::string>& a_aerosol_mat,
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_geom, a_dmap, a_ba, a_name)
        {
            define(a_vap_mat, a_aerosol_mat);
        }

        /*! \brief Set vapour material */
        virtual inline void setVapourMaterial ( const std::string& a_name )
        {
            m_vapour_mat = std::make_unique<MaterialProperties>(a_name);
        }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::string& a_name )
        {
            m_aerosol_mat.push_back(std::make_unique<MaterialProperties>(a_name));
        }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::vector<std::string>& a_names )
        {
            for (int i = 0; i < a_names.size(); i++) {
                setAerosolMaterial( a_names[i] );
            }
        }

        /*! \brief Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("ERFPCPC::varNames()");
            return {    AMREX_D_DECL("xvel","yvel","zvel"),
                        "particle_mass",
                        "radius",
                        "multiplicity",
                        "aerosol_mass",
                        "superdroplet_mass"};
        }

        /*! \brief Initialize super-droplets in domain */
        inline virtual void InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr)
        {
            InitializeParticles( m_initialization_type, a_ptr );
        }

        /*! \brief Initialize super-droplets in domain */
        virtual void InitializeParticles (  const std::string&,
                                            const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! \brief Particle initialization - uniform distribution */
        void initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);
        /*! \brief Set super-droplets attributes from a given condensate mass density */
        virtual void SetAttributes (  amrex::MultiFab& );

        /*! \brief Evolve particles for one time step */
        virtual void EvolveParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

        /*! \brief Advect particles for one time step */
        virtual void AdvectParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>&,
                                         const bool, const bool );

        /*! \brief Return the number of super-droplets
            This returns the number of AMReX particles being used in the simulation. */
        inline virtual amrex::Long NumSuperDroplets ()
        {
            return ERFPC::TotalNumberOfParticles();
        }

        /*! \brief Return the total number of physical particles */
        virtual amrex::Real TotalNumberOfParticles ();

        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity( amrex::MultiFab&, const int& a_comp=0 ) const;

        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity( amrex::MultiFab&, const int& a_comp=0 ) const;


    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        amrex::Real m_max_multiplicity;     /*!< Maximum multiplicity */
        amrex::Real m_numdens_init;         /*!< Initial number density of physical particles */
        amrex::Real m_numdens_sd_init;      /*!< Initial number density of super-droplets */
        amrex::Real m_mass_condensate_init; /*!< Initial condensate mass of each physical particle */
        amrex::Real m_mass_aerosol_init;    /*!< Initial aerosol mass of each physical particle */

        int m_num_sd_per_cell;  /*!< Number of super-droplets per cell */

        bool m_nucleate_particles;  /*!< nucleate new superdroplets from vapour */

        std::unique_ptr<MaterialProperties> m_vapour_mat; /*!< Vapour material */
        std::vector<std::unique_ptr<MaterialProperties>> m_aerosol_mat; /*!< Aerosol materials */

        /*! \brief read inputs from file */
        virtual void readInputs ();

        /*! \brief Particle initialization - null (no particles are initialized) */
        void initializeParticlesNull (const std::unique_ptr<amrex::MultiFab>&) { }

    private:

        /*! \brief define super-droplets */
        inline void define ( const std::string& a_vap_mat,
                             const std::vector<std::string>& a_aerosol_mat)
        {
            m_aerosol_mat.clear();

            add_superdroplet_attributes();
            readInputs();

            setVapourMaterial( a_vap_mat );
            setAerosolMaterial( a_aerosol_mat );
        }

        /*! \brief add super-droplet method-specific particle attributes */
        void add_superdroplet_attributes();

};

#endif
#endif
