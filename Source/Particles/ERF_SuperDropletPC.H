#ifndef SUPERDROPLET_PC_H_
#define SUPERDROPLET_PC_H_

#ifdef ERF_USE_PARTICLES

#include <random>

#include "ERF_IndexDefines.H"
#include "ERF_MaterialProperties.H"
#include "ERFPC.H"
#include "ERF_SDInitialization.H"

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
        term_vel,       /*!< Terminal velocity */
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        cond_tendency,  /*!< Condensation/evaporation tendency */
#endif
        uid,            /*!< unique ID */
        ncomps
    };
};

/*! \brief Int-type struct-of-array attributes of a super-droplet */
using SuperDropletsIntIdxSoA = ERFParticlesIntIdxSoA;

/*! \brief Types of terminal velocity models */
enum struct SDTerminalVelocityType {
    RogersYau,
    AtlasUlbrich,
    CloudRainShima
};

/*! \brief Types of coalescence kernels */
enum struct SDCoalescenceKernelType {
    golovin,
    sedimentation,
    Longs,
    Halls
};

/*! \brief Types of phase change ODE time integrator */
enum struct SDMassChangeTIMethod {
    RK3BS, RK4, BE, CN, DIRK2
};

using SDAerosolMassArr = amrex::GpuArray<amrex::ParticleReal*,SupDropInit::num_aerosols_max>;
using MFPtr = std::unique_ptr<amrex::MultiFab>;
using BCTypeArr = amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>;

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
            define( a_vap_mat,
                    a_aerosol_mat,
                    a_gdb->ParticleBoxArray(m_lev),
                    a_gdb->ParticleDistributionMap(m_lev) );
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
            define(a_vap_mat, a_aerosol_mat, a_ba, a_dmap);
        }

        /*! \brief Destructor */
        ~SuperDropletPC()
        {
            if (m_mass_change_logging) {
                fclose(m_mass_change_log);
            }
        }

        /*! \brief Set vapour material */
        virtual inline void setVapourMaterial ( const std::string& a_name )
        {
            m_vapour_mat = std::make_unique<MaterialProperties>(a_name);
        }

        /*! \brief Get vapour material */
        [[nodiscard]] virtual inline MaterialProperties& getVapourMaterial() const
        {
            return (*m_vapour_mat);
        }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::string& a_name )
        {
            m_aerosol_mat.push_back(std::make_unique<MaterialProperties>(a_name));
        }

        /*! \brief Set aerosol material */
        virtual inline void setAerosolMaterial ( const std::vector<std::string>& a_names )
        {
            for (unsigned long i = 0; i < a_names.size(); i++) {
                setAerosolMaterial( a_names[i] );
            }
        }

        /*! \brief Get real-type particle attribute names */
        [[nodiscard]] virtual amrex::Vector<std::string> varNames () const override;

        /*! \brief Get Eulerian plot variable names */
        [[nodiscard]] virtual amrex::Vector<std::string> meshPlotVarNames () const override;

        /*! \brief Compute mesh variable from particles */
        virtual void computeMeshVar( const std::string&,
                                     amrex::MultiFab&,
                                     const int ) const override;

        /*! \brief Initialize super-droplets in domain */
        virtual void InitializeParticles (const MFPtr& a_ptr) override;

        /*! \brief Initialize particles for a given set of initialization parameters */
        void initializeParticles ( const MFPtr&, const SDInitialization& );

        /*! \brief Set initial number of superdroplets per cell as a box with uniform density */
        void setNumSDBoxDistribution( amrex::iMultiFab&,
                                      const int,
                                      const MFPtr&,
                                      const amrex::RealBox& );

        /*! \brief Set initial number of superdroplets per cell as a bubble with uniform density */
        void setNumSDBubbleDistribution( amrex::iMultiFab&,
                                      const int,
                                      const MFPtr&,
                                      const amrex::RealBox& );

        /*! \brief Set super-droplets attributes from a given condensate mass density */
        virtual void SetAttributes (  amrex::MultiFab& );

        /*! \brief Scale the SDM number density with air density  */
        virtual void DensityScaling ( const amrex::MultiFab& );

        /*! \brief Evolve particles for one time step */
        virtual void EvolveParticles ( int,
                                       amrex::Real,
                                       amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                       const amrex::Vector<MFPtr>& ) override
        {
            amrex::Abort("SuperDropletPC::EvolveParticles() is intentionally disabled.");
        }

        /*! \brief Advect particles for one time step */
        virtual void AdvectParticles (   int,
                                         amrex::Real,
                                         amrex::Real,
                                         const amrex::MultiFab* const,
                                         const amrex::MultiFab&,
                                         const amrex::MultiFab&,
                                         const amrex::MultiFab&,
                                         const amrex::Vector<MFPtr>&,
                                         const BCTypeArr& );

        /*! \brief Grow/shrink particles for one time step */
        virtual void MassChange (   int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab&,
                                    const amrex::Vector<MFPtr>& );

        /*! \brief Coalescence of super-droplets */
        virtual void Coalescence (  int,
                                    amrex::Real,
                                    const amrex::MultiFab&,
                                    const amrex::MultiFab& );

        /*! \brief Recycle super-droplets */
        virtual void Recycle ( const int, const amrex::Vector<MFPtr>& );

        /*! \brief Return the number of super-droplets
            This returns the number of AMReX particles being used in the simulation. */
        [[nodiscard]] inline virtual amrex::Long NumSuperDroplets ()
        {
            return ERFPC::TotalNumberOfParticles();
        }

        /*! \brief Return the total number of physical particles */
        [[nodiscard]] virtual amrex::Real TotalNumberOfParticles ();

        /*! \brief Return the number of deactivated super-droplets */
        [[nodiscard]] virtual amrex::Long NumSDDeactivated ();

        /*! \brief Compute and print diagnostics */
        virtual void Diagnostics (const int&, const amrex::Real&, const bool);
        /*! \brief Compute mass density distribution */
        virtual void ComputeDistributions (  const int&,
                                             const amrex::ParticleReal&,
                                             const amrex::ParticleReal& );
#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        /*! \brief Compute binned size distributions */
        virtual void ComputeBinnedDistributions (  const int&);
        /*! \brief Compute cell-wise binned size distributions */
        virtual void ComputeBinnedDistributionsCell (  const int&, const amrex::Real&);
#endif

        /*! \brief Compute the particle number density on a mesh */
        virtual void numberDensity ( amrex::MultiFab&, const int& a_comp=0 ) const;
        /*! \brief Compute the particle mass density on a mesh */
        virtual void massDensity ( amrex::MultiFab&, const int& a_comp=0 ) const;
        /*! \brief Compute the condensate mass density on a mesh */
        virtual void massDensityCondensate ( amrex::MultiFab& a_mf ) const
        {
            massDensityCondensate(a_mf, 0.0, DBL_MAX);
        }
        /*! \brief Compute the condensate mass density on a mesh */
        virtual void massDensityCondensate ( amrex::MultiFab&,
                                             const amrex::Real,
                                             const amrex::Real,
                                             const int& a_comp = 0) const;
        /*! \brief Compute the mass flux component on a mesh */
        virtual void massFlux ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the condensate mass flux component on a mesh */
        virtual void massFluxCondensate ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the aerosol mass density on a mesh */
        virtual void aerosolMassDensity ( amrex::MultiFab&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the aerosol mass flux component on a mesh */
        virtual void aerosolMassFlux ( amrex::MultiFab&, const int&, const int&, const int& a_comp=0 ) const;
        /*! \brief Compute the particle effective radius on a mesh */
        virtual void effectiveRadius ( amrex::MultiFab&, const int& a_comp=0 ) const;

        /*! \brief Applies boundary treatment to the particle container */
        virtual void applyBoundaryTreatment ( int, const amrex::Vector<MFPtr>&, const BCTypeArr& );

    protected:

        const int m_lev = 0; /*!< AMR Level: always 0 for super-droplets PC */

        SDCoalescenceKernelType m_coalescence_kernel; /*!< Coalescence kernel */
        bool m_include_brownian_coalescence; /*!< Include Brownian coalescence? */

        SDTerminalVelocityType m_term_vel_type; /*!< Terminal velocity model */

        int m_num_sd_per_cell;  /*!< Number of super-droplets per cell */
        bool m_density_scaling; /*!< Scale initial number density with air density */
        bool m_nucleate_particles;  /*!< nucleate new superdroplets from vapour */

        int m_num_aerosols; /*!< Number of aerosols */
        std::unique_ptr<MaterialProperties> m_vapour_mat; /*!< Vapour material */
        std::vector<std::unique_ptr<MaterialProperties>> m_aerosol_mat; /*!< Aerosol materials */

        /* Mass change equation Newton solver parameters */
        amrex::Real m_newton_rtol;  /*!< Newton solver - relative tolerance */
        amrex::Real m_newton_atol;  /*!< Newton solver - absolute tolerance */
        amrex::Real m_newton_stol;  /*!< Newton solver - step tolerance */
        int m_newton_maxits;        /*!< Newton solver - maximum iterations */

        /* Mass change equation (un)convergence logs */
        bool m_mass_change_logging; /*!< Whether to log unconverged particles */
        FILE* m_mass_change_log; /*!< File handle for log with unconverged info */
        std::string m_mass_change_log_fname; /*!< Name of unconverged info log file */

        amrex::Real m_mass_change_cfl; /*!< CFL for phase change equation */
        SDMassChangeTIMethod m_mass_change_ti; /*!< time integrator for mass change ODE */
        long m_num_unconverged_particles; /*!< total number of unconverged particles */

        amrex::IntVect m_coalescence_bin_size; /*!< bin size for coalescence */

        int m_distribution_grid_size; /*!< Size of 1D grid to compute distributions on */

#ifdef ERF_USE_ML_UPHYS_DIAGNOSTICS
        amrex::Real m_bindist_rmin; /*!< min radius for binned distribution */
        amrex::Real m_bindist_rmax; /*!< max radius for binned distribution */

        /*! Cell-wise mass distribution */
        amrex::MultiFab m_mass_ln_R_mf;
        /*! Cell-wise number distribution */
        amrex::MultiFab m_num_ln_R_mf;
#endif

        /*! maximum coalescence time scale */
        amrex::Real m_t_coalescence;

        /*! Initialization parameters objects */
        int m_num_initializations = 1;
        /*! vector of initializations */
        std::vector< std::unique_ptr<SDInitialization> > m_initializations;

        /*! sigma_0 parameter for mass distribution calculation */
        amrex::Real m_sigma0;

        /*! place the initial SDs randomly within cells */
        bool m_place_randomly_in_cells;

        /*! random engine */
        std::mt19937 m_rndeng;

        /*! \brief read inputs from file */
        virtual void readInputs () override;

        /*! \brief Particle initialization - null (no particles are initialized) */
        void initializeParticlesNull (const MFPtr&) { }

    private:

        /*! \brief define super-droplets */
        void define ( const std::string&,
                      const std::vector<std::string>&,
                      const amrex::BoxArray&,
                      const amrex::DistributionMapping& );

        /*! \brief add super-droplet method-specific particle attributes */
        void add_superdroplet_attributes();

        /* do not use */
        using ERFPC::massDensity;

};

#endif
#endif
