#ifndef ERF_PC_H_
#define ERF_PC_H_

#ifdef ERF_USE_PARTICLES

#include <string>
#include <AMReX_Particles.H>
#include <ERF_Constants.H>

struct ERFParticlesIntIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct ERFParticlesRealIdxAoS
{
    enum {
        ncomps = 0
    };
};

struct ERFParticlesIntIdxSoA
{
    enum {
        ncomps = 0
    };
};

struct ERFParticlesRealIdxSoA
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

/** Particle binner using AMReX's getParticleCell for DenseBins. */
struct GetParticleBinERF
{
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> plo;
    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> dxi;
    amrex::Box domain;
    amrex::IntVect bin_size;
    amrex::Box box;

    template <typename P>
    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (const P& p) const noexcept
    {
        amrex::Box tbx;
        auto iv = amrex::getParticleCell(p, plo, dxi, domain);
        auto tid = amrex::getTileIndex(iv, box, true, bin_size, tbx);
        return static_cast<unsigned int>(tid);
    }
};

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_location_idata (  P& /*a_p*/,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& /*a_plo*/,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& /*a_dxi*/,
                              const amrex::Array4<amrex::Real const>&  /*a_height_arr*/ )
{
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int compute_k_from_z (amrex::Real z,
                      amrex::Real plo_z,
                      amrex::Real dxi_z,
                      int k_max,
                      const amrex::Real* zlevels,
                      int nz_lev0,
                      int ref_ratio)
{
    if (zlevels == nullptr || nz_lev0 <= 0) {
        int k = static_cast<int>(amrex::Math::floor((z - plo_z) * dxi_z));
        return amrex::max(0, amrex::min(k, k_max));
    }

    int k = static_cast<int>(amrex::Math::floor((z - plo_z) * dxi_z));
    k = amrex::max(0, amrex::min(k, k_max));

    auto fine_z = [=] AMREX_GPU_HOST_DEVICE (int k_f) -> amrex::Real {
        int kc = k_f / ref_ratio;
        int sub = k_f % ref_ratio;
        // kc is clamped to valid range for safety
        kc = amrex::max(0, amrex::min(kc, nz_lev0));
        if (kc >= nz_lev0) { return zlevels[nz_lev0]; }
        return zlevels[kc]
             + static_cast<amrex::Real>(sub)
             * (zlevels[kc+1] - zlevels[kc])
             / static_cast<amrex::Real>(ref_ratio);
    };

    for (int iter = 0; iter <= k_max; iter++) {
        amrex::Real z_lo = fine_z(k);
        amrex::Real z_hi = fine_z(k + 1);
        if (z >= z_hi && k < k_max) {
            k++;
        } else if (z < z_lo && k > 0) {
            k--;
        } else {
            break;
        }
    }
    return k;
}

class ERFPC : public amrex::ParticleContainer<  ERFParticlesRealIdxAoS::ncomps,  // AoS real attributes
                                                ERFParticlesIntIdxAoS::ncomps,   // AoS integer attributes
                                                ERFParticlesRealIdxSoA::ncomps,  // SoA real attributes
                                                ERFParticlesIntIdxSoA::ncomps,   // SoA integer attributes
                                                amrex::DefaultAllocator >
{
    public:

        /*! Constructor */
        ERFPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< ERFParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        ERFParticlesIntIdxAoS::ncomps,    // AoS integer attributes
                                        ERFParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        ERFParticlesIntIdxSoA::ncomps,    // SoA integer attributes
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
            : amrex::ParticleContainer< ERFParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        ERFParticlesIntIdxAoS::ncomps,    // AoS real attributes
                                        ERFParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        ERFParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator> ( a_geom, a_dmap, a_ba )
        {
            BL_PROFILE("ERFPCPC::ERFPC()");
            m_name = a_name;
            readInputs();
        }

        /*! Initialize particles in domain */
        virtual void InitializeParticles (const amrex::Real time, const std::unique_ptr<amrex::MultiFab>& a_ptr = nullptr);

        /*! Evolve particles for one time step */
        virtual void EvolveParticles (   int,
                                         amrex::Real,
                                         amrex::Vector<amrex::Vector<amrex::MultiFab>>&,
                                         const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& );

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
        virtual void AdvectWithFlow (  amrex::MultiFab*,
                                       int,
                                       amrex::Real,
                                       const std::unique_ptr<amrex::MultiFab>& );

        /*! Uses midpoint method to advance particles falling under gravity. */
        virtual void AdvectWithGravity ( int,
                                         amrex::Real,
                                         const std::unique_ptr<amrex::MultiFab>& );

        /*! Interpolates flow temperature to particles */
        virtual void ComputeTemperature ( const amrex::MultiFab&,
                                          int,
                                          amrex::Real,
                                          const std::unique_ptr<amrex::MultiFab>& );

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

        /*! Fix k-indices for all particles after AMR regrid/removal */
        void FixKIndexAMR (const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd);

        /*! Diagnostic: count particles per level and halo (for AMR) */
        void CountParticlesPerLevelAndHalo ( int a_finest_level );

    protected:

        bool m_advect_w_flow;               /*!< advect with flow velocity */
        bool m_advect_w_gravity;            /*!< advect under gravitational force */

        amrex::RealBox m_particle_box;      /*!< box within which to place particles */

        std::string m_name;                 /*!< name of this particle species */

        std::string m_initialization_type;  /*!< initial particle distribution type */
        int m_ppc_init;                     /*!< initial number of particles per cell */

        amrex::Real m_inject_start_time;    /*!< particles will only be initialized/injected if time >= m_inject_start_time */

        bool m_initialized = false;         /*!< have the particles been initialized */

        bool m_stable_redistribute;         /*!< use stable redistribute for deterministic simulations */

        amrex::Gpu::DeviceVector<amrex::Real> m_zlevels_d; /*!< Level-0 z-coordinate node positions for non-uniform vertical grids */

        /*! read inputs from file */
        virtual void readInputs ();

    private:

        bool place_randomly_in_cells; /*!< place particles at random positions? */
};

#endif
#endif
