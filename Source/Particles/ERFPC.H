#ifndef ERF_PC_H_
#define ERF_PC_H_

#ifdef ERF_USE_PARTICLES

#include <string>
#include <AMReX_Particles.H>
#include <AMReX_ParticleUtil.H>
#include <ERF_Constants.H>

/** Particle cell-binning functor for DenseBins.  Maps a particle to its
 *  cell-tile index via amrex::getParticleCell + amrex::getTileIndex. */
struct GetParticleBinERF
{
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> plo;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxi;
    amrex::Box                                   domain;
    amrex::IntVect                               bin_size;
    amrex::Box                                   box;
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (const P& p) const noexcept
    {
        amrex::Box tbx;
        auto iv  = amrex::getParticleCell(p, plo, dxi, domain);
        auto tid = amrex::getTileIndex(iv, box, true, bin_size, tbx);
        return static_cast<unsigned int>(tid);
    }
};

/**
 * @brief Index mapping for integer particle attributes.
 */
struct ERFParticlesIntIdx
{
    enum {
        ncomps = 0
    };
};

/**
 * @brief Index mapping for real-valued particle attributes.
 */
struct ERFParticlesRealIdx
{
    enum {
        vx = 0,
        vy,
        vz,
        mass,
        temperature,
        ncomps
    };
};

namespace ERFParticleInitializations
{
    /* list of particle initializations */
    const std::string init_box_uniform  = "box";
}

namespace ERFParticleNames
{
    const std::string tracers = "tracer_particles";
}

/**
 * @brief Particle container for ERF particles.
 */
class ERFPC : public amrex::ParticleContainer<  0,                              // AoS real attributes
                                                0,                              // AoS integer attributes
                                                ERFParticlesRealIdx::ncomps,    // SoA real attributes
                                                ERFParticlesIntIdx::ncomps,     // SoA integer attributes
                                                amrex::DefaultAllocator >
{
    public:

        /*! Constructor */
        ERFPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< 0,                              // AoS real attributes
                                        0,                              // AoS integer attributes
                                        ERFParticlesRealIdx::ncomps,    // SoA real attributes
                                        ERFParticlesIntIdx::ncomps,     // SoA integer attributes
                                        amrex::DefaultAllocator> (a_gdb)
        {
            BL_PROFILE("ERFPCPC::ERFPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Constructor */
        ERFPC ( const amrex::Geometry&            a_geom,
                const amrex::DistributionMapping& a_dmap,
                const amrex::BoxArray&            a_ba,
                const std::string&                a_name = "particles" )
            : amrex::ParticleContainer< 0,                              // AoS real attributes
                                        0,                              // AoS integer attributes
                                        ERFParticlesRealIdx::ncomps,    // SoA real attributes
                                        ERFParticlesIntIdx::ncomps,     // SoA integer attributes
                                        amrex::DefaultAllocator> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("ERFPCPC::ERFPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (const double time, const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! Evolve particles for one time step */
        virtual void EvolveParticles (int                                                    a_lev,
                                      double                                                 a_dt_lev,
                                      amrex::Vector<amrex::Vector<amrex::MultiFab>>&         a_flow_vars,
                                      const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd );

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> varNames () const
        {
            BL_PROFILE("ERFPCPC::varNames()");
            return {AMREX_D_DECL("xvel","yvel","zvel"),"mass","temperature"};
        }

        /*! Get real-type particle attribute names */
        virtual amrex::Vector<std::string> meshPlotVarNames () const
        {
            BL_PROFILE("ERFPCPC::varNames()");
            return {"mass_density"};
        }

        /*! Uses midpoint method to advance particles using flow velocity. */
        virtual void AdvectWithFlow (amrex::MultiFab*                        a_umac,
                                     int                                     a_lev,
                                     double                                  a_dt,
                                     const std::unique_ptr<amrex::MultiFab>& a_z_height);

        /*! Uses midpoint method to advance particles falling under gravity. */
        virtual void AdvectWithGravity (int                                     a_lev,
                                        double                                  a_dt,
                                        const std::unique_ptr<amrex::MultiFab>& a_z_height);

        /*! Interpolates flow temperature to particles */
        virtual void ComputeTemperature (const amrex::MultiFab&                  a_ucons,
                                         int                                     a_lev,
                                         double                                /*a_dt*/,
                                         const std::unique_ptr<amrex::MultiFab>& a_z_height);

        /*! Compute mass density */
        virtual void massDensity ( amrex::MultiFab&,
                                   const amrex::MultiFab& a_z_phys_nd,
                                   const int&, const int& a_comp = 0) const;

        /*! Compute mesh variable from particles */
        virtual void computeMeshVar(    const std::string&      a_var_name,
                                        amrex::MultiFab&        a_mf,
                                        const amrex::MultiFab&  a_z_phys_nd,
                                        const int               a_lev) const
        {
            if (a_var_name == "mass_density") {
                massDensity( a_mf, a_z_phys_nd, a_lev );
            } else {
                a_mf.setVal(0.0);
            }
        }

        /*! Specify if particles should advect with flow */
        inline void setAdvectWithFlow (bool a_flag)
        {
            BL_PROFILE("ERFPCPC::setAdvectWithFlow()");
            m_advect_w_flow = a_flag;
        }
        /*! Specify if particles fall under gravity */
        inline void setAdvectWithGravity (bool a_flag)
        {
            BL_PROFILE("ERFPCPC::setAdvectWithGravity()");
            m_advect_w_gravity = a_flag;
        }

        // the following functions should ideally be private or protected, but need to be
        // public due to CUDA extended lambda capture rules

        /*! \brief Helper for ParticleToMesh operations on terrain-following grids.
         *
         * Generic CIC deposition from particles to a MultiFab. The caller provides
         * a device lambda `value_func(ptd, i)` that returns the per-particle scalar
         * to deposit (e.g. mass, number, flux). Implementation in ERFPCParticleToMesh.H.
         */
        template<typename ValueFunc>
        void ERFPCParticleToMesh(amrex::MultiFab& a_mf,
                                 const amrex::MultiFab& a_z_phys_nd,
                                 int a_lev, int a_comp,
                                 ValueFunc&& value_func) const;

        /*! Default particle initialization */
        void initializeParticlesUniformDistributionInBox (const std::unique_ptr<amrex::MultiFab>& a_ptr,
                                                          const amrex::RealBox& particle_box);

        /*! Replace pos(2) (zeta) with the corresponding physical z, in place,
         *  on every level.  Used around plotfile writes so the on-disk
         *  particles are in physical coordinates. */
        void ConvertZetaToZ (const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd);

        /*! Inverse of ConvertZetaToZ: replace pos(2) (physical z) with zeta. */
        void ConvertZToZeta (const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd);

        /*! Route fine-level particles whose position is outside the level's
         *  BA down to coarser levels, so per-level Redistribute on the fine
         *  level does not assert. */
        void ExtractAndRouteOORParticles ( int a_lev );

        /*! Diagnostic: count particles per level and halo (for AMR) */
        void CountParticlesPerLevelAndHalo ( int a_finest_level );

        /*! Split coarse-level particles in newly-refined cells after AMR regrid
         *  to restore the per-cell super-droplet count.  Handles cascading
         *  refinement across multiple new levels in a single call (see the
         *  SuperDropletPC override for details).  Default no-op. */
        virtual void SplitParticlesForRefinement ( int /*a_finest_level*/ ) {}

        /*! Split new entrants and merge departees at AMR level boundaries.
         *  Called every Redistribute(z_phys_nd) to maintain uniform per-cell
         *  super-droplet density across levels.  Default no-op. */
        virtual void SplitMergeAtLevelBoundary () {}

        /*! Regrid-time de-refinement merge.  For cells that were in the old
         *  fine BoxArray at level `a_lev` but are no longer in the current
         *  one, reduces the coarse-level super-droplet count from `np_bin`
         *  to `np_bin / merge_factor`.  Default no-op. */
        virtual void MergeParticlesAtDerefinement ( int /*a_lev*/,
                                                    const amrex::BoxArray& /*a_old_fine_ba*/,
                                                    const amrex::IntVect& /*a_ref_ratio*/ ) {}

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */
        bool m_advect_w_gravity;            /*!< advect under gravitational force */

        amrex::RealBox m_particle_box;      /*!< box within which to place particles */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        double m_inject_start_time;    /*!< particles will only be initialized/injected if time >= m_inject_start_time */

        bool m_initialized = false;         /*!< have the particles been initialized */

        bool m_stable_redistribute;         /*!< use stable redistribute for deterministic simulations */

        /*! read inputs from file */
        virtual void readInputs ();

    private:

        bool place_randomly_in_cells; /*!< place particles at random positions? */
};

#endif
#endif
