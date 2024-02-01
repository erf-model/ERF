#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include "ERFPC.H"

/*! \brief Real-type array-of-struct attributes of a super-droplet */
using SuperDropletsRealIdxAoS = ERFParticlesRealIdx;

/*! \brief Real-type struct-of-array attributes of a super-droplet */
struct SuperDropletsRealIdxSoA
{
    enum {
        radius = 0, /*!< Droplet radius */
        sol_mass,   /*!< Solute mass */
        ncomps
    };
};

/*! \brief Int-type array-of-struct attributes of a super-droplet */
using SuperDropletsIntIdxAoS = ERFParticlesIntIdx;

/*! \brief Int-type struct-of-array attributes of a super-droplet */
struct SuperDropletsIntIdxSoA
{
    enum {
        multiplicity = 0, /*!< Number of droplets that this super-droplet represents */
        ncomps
    };
};

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
                          const std::string& a_name = "super_droplets" )
                      : ERFPC (a_gdb, a_name)
        {
            add_superdroplet_attributes();
            readInputs();
        }

        /*! \brief Constructor */
        SuperDropletPC (  const amrex::Geometry&            a_geom,
                          const amrex::DistributionMapping& a_dmap,
                          const amrex::BoxArray&            a_ba,
                          const std::string&                a_name = "super_droplets" )
                      : ERFPC (a_geom, a_dmap, a_ba, a_name)
        {
            add_superdroplet_attributes();
            readInputs();
        }

        /*! \brief Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("ERFPCPC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass","radius","aerosol_mass"};
        }

        /*! \brief Initialize super-droplets in domain */
        inline virtual void InitializeParticles (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr)
        {
            InitializeParticles( m_initialization_type, a_ptr );
        }

        /*! \brief Initialize super-droplets in domain */
        virtual void InitializeParticles (  const std::string&,
                                            const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! \brief Set super-droplets attributes from a given condensate mass density */
        virtual void SetAttributes (  amrex::MultiFab& );

        /*! \brief Evolve particles for one time step */
        virtual void EvolveParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

        /*! \brief Return the number of super-droplets
            This returns the number of AMReX particles being used in the simulation. */
        inline virtual amrex::Long NumSuperDroplets ()
        {
            return ERFPC::TotalNumberOfParticles();
        }

        /*! \brief Return the total number of particles */
        virtual amrex::Long TotalNumberOfParticles ();

        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity( amrex::MultiFab&, const int& a_comp=0 ) const;

        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity( amrex::MultiFab&, const int& a_comp=0 ) const;

    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        int m_max_multiplicity;      /*!< Maximum multiplicity */
        amrex::Real m_numdens_init;  /*!< Initial number density of particles */
        amrex::Real m_radius_init;   /*!< Initial radius */

        /*! \brief read inputs from file */
        virtual void readInputs ();

        /*! \brief Particle initialization - null (no particles are initialized) */
        void initializeParticlesNull (const std::unique_ptr<amrex::MultiFab>&) { }
        /*! \brief Particle initialization - uniform distribution */
        void initializeParticlesUniformDistribution (const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        bool m_nucleate_particles;  /*!< nucleate new superdroplets from vapour */

    private:

        /*! \brief add super-droplet method-specific particle attributes */
        void add_superdroplet_attributes();
};

#endif
#endif
