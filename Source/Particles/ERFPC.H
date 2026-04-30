#ifndef ERF_PC_H_
#define ERF_PC_H_

#ifdef ERF_USE_PARTICLES

#include <string>
#include <AMReX_Particles.H>
#include <ERF_Constants.H>

struct ERFParticlesIntIdxAoS
{
    enum {
        k = 0,        // cell index at the particle's current AMR level
                      // — read directly by AMReX cic_/mac_interpolate_mapped_z
        k_finest,     // cell index at the FINEST AMR level — used only by
                      // ERFParticlesAssignor; lets multi-level Redistribute
                      // place particles correctly across levels via scaling
        ncomps
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

/*! \brief Custom particle Assignor for AMR-aware z-cell placement.
 *
 *  The Assignor returns the (i,j,k) IntVect at the queried level for a
 *  particle.  In z, it converts idata(k_finest) — which always holds the
 *  cell index at the FINEST AMR level — to the queried-level cell index
 *  by dividing by the cumulative refinement ratio between the queried and
 *  finest levels.  The ratio is derived from the queried-level dxi[2] and
 *  the static s_dxi_z_finest, set by ERFPC::RefreshAssignorState() after
 *  every regrid.  Single-level mode (sentinel = -1) falls back to using
 *  idata(k) directly.
 */
struct ERFParticlesAssignor
{
    static AMREX_GPU_MANAGED amrex::Real s_dxi_z_finest;

    template <typename P>
    AMREX_GPU_HOST_DEVICE
    amrex::IntVect operator() (  P const& p,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                 amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                 const amrex::Box& domain ) const noexcept
    {
        int k_idx;
        if (s_dxi_z_finest > amrex::Real(0)) {
            int ratio = static_cast<int>(amrex::Math::round(s_dxi_z_finest / dxi[AMREX_SPACEDIM-1]));
            if (ratio < 1) { ratio = 1; }
            k_idx = p.idata(ERFParticlesIntIdxAoS::k_finest) / ratio;
        } else {
            k_idx = p.idata(ERFParticlesIntIdxAoS::k);
        }
        amrex::IntVect iv(
            AMREX_D_DECL(   int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                            int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                            k_idx ) );
        iv[0] += domain.smallEnd()[0];
        iv[1] += domain.smallEnd()[1];
        iv[2] = amrex::max(domain.smallEnd()[2],
                           amrex::min(iv[2], domain.bigEnd()[2]));
        return iv;
    }
};

/*! \brief Convert a finest-level cell index to a level cell index. */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int finestK_to_levelK (int k_finest, int ref_to_finest) noexcept
{
    return (ref_to_finest > 1) ? (k_finest / ref_to_finest) : k_finest;
}

/*! \brief Convert a level cell index to a finest-level cell index (bottom-aligned). */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int levelK_to_finestK (int k_level, int ref_to_finest) noexcept
{
    return (ref_to_finest > 1) ? (k_level * ref_to_finest) : k_level;
}

/** Particle binner using ERFParticlesAssignor (k-index for z) for DenseBins. */
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
        auto iv = ERFParticlesAssignor()(p, plo, dxi, domain);
        auto tid = amrex::getTileIndex(iv, box, true, bin_size, tbx);
        return static_cast<unsigned int>(tid);
    }
};

/*! \brief Refine the particle's idata(k) and idata(k_finest) after a position
 *  update.  Walks the level's terrain heights starting from the current
 *  idata(k) to find the cell containing pos(2); writes the result to
 *  idata(k) and to idata(k_finest) = level_k * a_ref_to_finest.
 */
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void update_location_idata (  P& a_p,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& a_dxi,
                              const amrex::Array4<amrex::Real const>&  a_height_arr,
                              int a_ref_to_finest = 1 )
{
    if (a_height_arr) {
        // Compute horizontal indices from particle position
        int ix = int(amrex::Math::floor((a_p.pos(0)-a_plo[0])*a_dxi[0]));
        int iy = int(amrex::Math::floor((a_p.pos(1)-a_plo[1])*a_dxi[1]));

        // Check if horizontal indices are within the valid array bounds.
        // We need ix, ix+1, iy, iy+1 to be valid for bilinear interpolation.
        // If particle has moved outside the local tile region (after advection
        // but before redistribution), skip the k-index update. The particle
        // will be redistributed to the correct tile and can update then.
        if (ix   < a_height_arr.begin[0] || ix+1 >= a_height_arr.end[0] ||
            iy   < a_height_arr.begin[1] || iy+1 >= a_height_arr.end[1]) {
            return;
        }

        amrex::Real lx = (a_p.pos(0)-a_plo[0])*a_dxi[0] - static_cast<amrex::Real>(ix);
        amrex::Real ly = (a_p.pos(1)-a_plo[1])*a_dxi[1] - static_cast<amrex::Real>(iy);

        // Clamp k to the valid range of the height array so the iteration
        // can proceed even with partial-z AMR tiles.
        int klo = a_height_arr.begin[2];
        int khi = a_height_arr.end[2] - 2; // need iz+1 < end[2]
        if (khi < klo) { return; }
        int kcur = a_p.idata(ERFParticlesIntIdxAoS::k);
        kcur = amrex::max(klo, amrex::min(kcur, khi));

        while (1) {
            if (kcur < klo || kcur > khi) {
                break;
            }

            auto zlo = a_height_arr(ix  ,iy  ,kcur  ) * (one-lx) * (one-ly) +
                       a_height_arr(ix+1,iy  ,kcur  ) *      lx  * (one-ly) +
                       a_height_arr(ix  ,iy+1,kcur  ) * (one-lx) * ly +
                       a_height_arr(ix+1,iy+1,kcur  ) *      lx  * ly;
            auto zhi = a_height_arr(ix  ,iy  ,kcur+1) * (one-lx) * (one-ly) +
                       a_height_arr(ix+1,iy  ,kcur+1) *      lx  * (one-ly) +
                       a_height_arr(ix  ,iy+1,kcur+1) * (one-lx) * ly +
                       a_height_arr(ix+1,iy+1,kcur+1) *      lx  * ly;

            if (a_p.pos(2) > zhi && kcur < khi) {
                kcur += 1;
            } else if (a_p.pos(2) <= zlo && kcur > klo) {
                kcur -= 1;
            } else {
                break;
            }
        }

        a_p.idata(ERFParticlesIntIdxAoS::k) = kcur;
        a_p.idata(ERFParticlesIntIdxAoS::k_finest) =
            levelK_to_finestK(kcur, a_ref_to_finest);
    }
}

/*! \brief Compute vertical cell index k from particle z-position.
 *
 * Handles both uniform and non-uniform (stretched) vertical grids.
 * For non-uniform grids, uses the level-0 z-levels array with the
 * cumulative refinement ratio to compute fine-level cell boundaries.
 * Starts from an initial guess and iterates (like update_location_idata).
 *
 * \param[in] z         Particle z-position (physical coordinates)
 * \param[in] plo_z     Domain low z from Geometry
 * \param[in] dxi_z     Inverse cell size in z from Geometry
 * \param[in] k_max     Maximum valid k index (domain.bigEnd(2))
 * \param[in] zlevels   Level-0 z-level node positions (nullptr if uniform)
 * \param[in] nz_lev0   Number of level-0 cells (nz_levels - 1)
 * \param[in] ref_ratio Cumulative refinement ratio from level 0 to target level
 * \return              Cell index k in [0, k_max]
 */
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
        // Uniform grid: direct computation
        int k = static_cast<int>(amrex::Math::floor((z - plo_z) * dxi_z));
        return amrex::max(0, amrex::min(k, k_max));
    }

    // Non-uniform z-levels: initial guess from uniform formula, then iterate.
    // Fine-level cell boundary positions are linearly interpolated from
    // the level-0 z-levels: for fine cell k_f, the coarse cell is k_f/R
    // and the sub-cell position is (k_f%R)/R within that coarse cell.
    int k = static_cast<int>(amrex::Math::floor((z - plo_z) * dxi_z));
    k = amrex::max(0, amrex::min(k, k_max));

    // Helper: compute z-position of fine cell boundary k_f
    // fine_z(k_f) = zlevels[k_f/R] + (k_f%R) * (zlevels[k_f/R+1] - zlevels[k_f/R]) / R
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
                                                amrex::DefaultAllocator,
                                                ERFParticlesAssignor >
{
    public:

        /*! Constructor */
        ERFPC ( amrex::ParGDBBase* a_gdb,
                const std::string& a_name = "particles" )
            : amrex::ParticleContainer< ERFParticlesRealIdxAoS::ncomps,   // AoS real attributes
                                        ERFParticlesIntIdxAoS::ncomps,    // AoS integer attributes
                                        ERFParticlesRealIdxSoA::ncomps,   // SoA real attributes
                                        ERFParticlesIntIdxSoA::ncomps,    // SoA integer attributes
                                        amrex::DefaultAllocator,
                                        ERFParticlesAssignor> (a_gdb)
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
                                        amrex::DefaultAllocator,
                                        ERFParticlesAssignor> ( a_geom, a_dmap, a_ba )
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

        /*! Refresh the static Assignor state used to convert idata(k_finest)
         *  to a queried-level cell index.  Must be called once after construction
         *  and after every regrid (before any Assignor-using operation). */
        void RefreshAssignorState ();

        /*! Extract and route out-of-range particles to proper levels */
        void ExtractAndRouteOORParticles ( int a_lev,
                                           const amrex::Vector<std::unique_ptr<amrex::MultiFab>>& a_z_phys_nd );

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
