#ifndef NOAHMP_H
#define NOAHMP_H

#include <string>
#include <vector>
#include <memory>
#include <limits>
#include <functional>

#include <AMReX_ParmParse.H>
#include <AMReX_Interpolater.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>
#include <AMReX_MFIter.H>
#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

#include <ERF_NullSurf.H>
#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_DataStruct.H>
#include <ERF_SurfacePrecipitation.H>

#include <NoahmpIO.H>
#include <ERF_NOAHMP_Fields.H>

// Small POD helpers used by the Advance_With_State
namespace erf_noahmp {
    // Typed precip accumulator views collected once per Advance: per slot, the
    // borrowed native accumulator, its native->kg/m^2 factor, and whether it is filled.
    struct PrecipSlots {
        const amrex::MultiFab* accum [NoahmpPrecipSlot::NumSlots] = {};
        amrex::Real            factor[NoahmpPrecipSlot::NumSlots] = {};
        bool                   has   [NoahmpPrecipSlot::NumSlots] = {};
        bool                   have_any = false;
    };

    // Cells the precip guard clamped this call (per rank); reported after the box loop.
    struct ClampedPrecipCell { int i; int j; amrex::Real raw_mm; };

    // Cells where frozen components violated MP_SNOW+MP_GRAUP <= MP_RAINNC and were
    // rescaled (per rank); reported after the loop.
    struct InvariantPrecipCell { int i; int j; amrex::Real froz_mm; amrex::Real total_mm; };

    // Per-surface-box handles co-indexed by idb: the NoahmpIO block plus its
    // pinned input/output staging buffers.
    struct NoahmpBlockViews {
        NoahmpIO_type*             io;      // &noahmpio_vect[idb]
        amrex::Array4<amrex::Real> input;   // noahmp_input_tmp[idb]->array()
        amrex::Array4<amrex::Real> output;  // noahmp_output_tmp[idb]->array()
    };
}

class NOAHMP : public NullSurf {

    using FabPtr = std::shared_ptr<amrex::MultiFab>;

public:
    // Constructor
    NOAHMP () {}

    // Destructor
    virtual ~NOAHMP () = default;

    // Set thermo and grid properties
    void
    Define (SolverChoice& /*sc*/) override
    {
        // Parse solver constants from sc here if needed.
    }

    // Initialize data structures
    void
    Init (const int& lev,
          const amrex::MultiFab& cons_in,
          const amrex::Geometry& geom,
          const amrex::Geometry& geom0,
          amrex::Vector<amrex::BCRec>& domain_bcs_type,
          amrex::IntVect& refRatio,
          const amrex::Real& dt,
          amrex::Vector<amrex::Vector<std::string>>& nc_init_file) override;

    // Wrapper to do all the updating
    void
    Advance_With_State (const int& lev,
                        amrex::MultiFab& cons_in,
                        amrex::MultiFab& xvel_in,
                        amrex::MultiFab& yvel_in,
                        amrex::MultiFab* hfx3_out,
                        amrex::MultiFab* qfx3_out,
                        const SurfacePrecipAccumulationSources& precip_sources,
                        const amrex::Real& elapsed_time,
                        const amrex::Real& dt,
                        const int& nstep,
                        const bool updated_lev0) override;

    void
    Plot_Landfile (const int& nstep) override;

    // Get state vars from lsm class
    amrex::MultiFab*
    Lsm_Data_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmDataMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_data_size && lsmIdx>=0);
        return lsm_fab_data[lsmIdx].get();
    }

    // Get flux vars from lsm class
    amrex::MultiFab*
    Lsm_Flux_Ptr (const int& varIdx) override
    {
        int lsmIdx = LsmFluxMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_flux_size && lsmIdx>=0);
        return lsm_fab_flux[lsmIdx].get();
    }

    // Get lsm geometry
    amrex::Geometry
    Lsm_Geom ( ) override { return m_lsm_geom; }

    // Get number of vars. Parent queries this before Init (to size lsm_data), so
    // resolve the soil-layer count here on first use rather than relying on Init.
    int
    Lsm_Data_Size () override { m_ensure_nsoil_resolved(); return NOAHMP::m_lsm_data_size; }

    // Get number of fluxes lsm class contains
    int
    Lsm_Flux_Size () override { return NOAHMP::m_lsm_flux_size; }

    // Get variable names
    std::string
    Lsm_DataName (const int& varIdx) override
    {
        int lsmIdx = LsmDataMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_data_size && lsmIdx>=0);
        return LsmDataName[lsmIdx];
    }

    // Get variable index from name
    int
    Lsm_DataIndex (std::string varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<m_lsm_data_size; ++idx) {
            if (lc_varname == amrex::toLower(LsmDataName[idx])) {
                varIdx = idx;
            }
        }
        return varIdx;
    }

    // Get flux variable names
    std::string
    Lsm_FluxName (const int& varIdx) override
    {
        int lsmIdx = LsmFluxMap[varIdx];
        AMREX_ALWAYS_ASSERT(lsmIdx < NOAHMP::m_lsm_flux_size && lsmIdx>=0);
        return LsmFluxName[lsmIdx];
    }

    // Get flux variable index from name
    int
    Lsm_FluxIndex (std::string varname) override
    {
        int varIdx = -1;
        std::string lc_varname = amrex::toLower(varname);
        for (int idx(0); idx<LsmFlux_NOAHMP::NumVars; ++idx) {
            if (lc_varname == amrex::toLower(LsmFluxName[idx])) {
                varIdx = idx;
            }
        }
        return varIdx;
    }

    void
    Lsm_Set_Lev0_Data_Ptr (const int& varIdx,
                           amrex::MultiFab* lev0_data_ptr) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < m_lsm_data_size);
        lsm_lev0_data[varIdx] = lev0_data_ptr;
    }

    void
    Lsm_Set_Lev0_Flux_Ptr (const int& varIdx,
                           amrex::MultiFab* lev0_flux_ptr) override
    {
        AMREX_ALWAYS_ASSERT(varIdx >= 0 && varIdx < int(LsmFlux_NOAHMP::NumVars));
        lsm_lev0_flux[varIdx] = lev0_flux_ptr;
    }

    bool
    Get_Update_Status () override { return m_updated; }

    // Checkpoint/restore the NoahMP substep counter via m_itimestep (valid on every
    // rank, unlike noahmpio_vect[0] which is absent on land-free ranks).
    int
    Get_LSM_Step () const override
    {
        return m_itimestep;
    }

    void
    Set_LSM_Step (int step) override
    {
        m_itimestep = step;
        for (auto& noahmpio : noahmpio_vect) {
            noahmpio.itimestep = step;
        }
    }

    // Write/read the full NoahMP prognostic state to/from a NetCDF restart dir so a
    // restart reproduces a cold-start trajectory bitwise (#3255). ERF_NOAHMP_IO.cpp.
    void
    Write_Lsm_Restart (const std::string& dir) const override;

    void
    Read_Lsm_Restart (const std::string& dir) override;

    // --- Advance_With_State pipeline helpers (ERF_NOAHMP_Advance.cpp) ---------

    // Fine-level (lev>0) coupling data/flux fill by interpolation from level 0,
    // taken when there is no per-level NetCDF land file. Sets m_updated.
    void
    interp_from_lev0 (const int& lev,
                      amrex::MultiFab& cons_in,
                      const int& nstep);

    // Subcycling gate: true if the model fires this ERF step (advances the lockstep
    // counter + sets m_updated). Identical on every rank.
    bool
    time_to_fire (const amrex::Real& elapsed_time);

    // Per surface box: stage ERF forcing into the pinned input buffer and copy into
    // the NoahmpIO arrays (steps 1-3). blk bundles the block + its staging buffers.
    void
    stage_forcing (const amrex::MFIter& mfi,
                   const erf_noahmp::NoahmpBlockViews& blk,
                   const amrex::Box& bx,
                   const int& klo,
                   const int& lev,
                   const bool is_moist,
                   amrex::MultiFab& cons_in,
                   amrex::MultiFab& xvel_in,
                   amrex::MultiFab& yvel_in,
                   const erf_noahmp::PrecipSlots& precip,
                   amrex::Vector<erf_noahmp::ClampedPrecipCell>& clamped_cells,
                   amrex::Vector<erf_noahmp::InvariantPrecipCell>& invariant_cells);

    // Per surface box: read NoahmpIO results into the pinned output buffer and scatter
    // them into the ERF coupling fields (steps 5-6), applying the -9999 fill guard.
    void
    read_results (const amrex::MFIter& mfi,
                  const erf_noahmp::NoahmpBlockViews& blk,
                  const amrex::Box& bx,
                  const amrex::Box& gbx,
                  amrex::MultiFab& cons_in);

    // --- Precip bookkeeping (ERF_NOAHMP_Precip.cpp) --------------------------

    // Collect the typed precip accumulator views for this Advance.
    erf_noahmp::PrecipSlots
    collect_precip_sources (const SurfacePrecipAccumulationSources& precip_sources);

    // Lazily allocate + seed the per-level "previous cumulative accum" snapshots
    // (cold start seeds from current; restart consumes m_precip_accum_restored).
    void
    prepare_precip_snapshots (const int& lev,
                              const erf_noahmp::PrecipSlots& precip);

    // Advance the snapshots to the current accumulation after the land call.
    void
    advance_precip_snapshots (const int& lev,
                              const erf_noahmp::PrecipSlots& precip);

    // Report (per rank, never silent) any cells clamped or rescaled by the guard.
    void
    report_precip_diagnostics (const int& lev,
                               const amrex::Vector<erf_noahmp::ClampedPrecipCell>& clamped_cells,
                               const amrex::Vector<erf_noahmp::InvariantPrecipCell>& invariant_cells);

    // 3 soil variable groups: smois, sh2o, tslb. All 3*nsoil buffer-size occurrences use this.
    static constexpr int m_num_soil_groups = 3;

    // Index of soil field (g in {0:smois,1:sh2o,2:tslb}, k in [0,nsoil)) within
    // lsm_fab_data. Output-buffer index is NoahmpOutputComp::NumComps + g*m_nsoil + k.
    int soil_data_idx (int g, int k) const {
        return LsmData_NOAHMP::NumVars + g*m_nsoil + k;
    }

    // Resolve the soil-layer count and total LSM data size exactly once, before Init
    // (the parent queries Lsm_Data_Size() first). Init later asserts vs namelist NSOIL.
    bool m_nsoil_resolved = false;
    void m_ensure_nsoil_resolved () {
        if (m_nsoil_resolved) { return; }
        amrex::ParmParse pp_nsoil("erf");
        pp_nsoil.queryAdd("lsm_nsoil", m_nsoil);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nsoil >= 1, "erf.lsm_nsoil must be >= 1");
        m_lsm_data_size = LsmData_NOAHMP::NumVars + m_num_soil_groups*m_nsoil;
        m_nsoil_resolved = true;
    }

    // VisMF file prefix for the level-l, slot-s previous-accumulation snapshot.
    static std::string
    m_precip_snapshot_name (const std::string& dir, int l, int s) {
        return dir + "/precip_accum_prev_Level_" + std::to_string(l)
                   + "_slot_" + std::to_string(s);
    }

private:

    // Runs fn(noahmpio) per local block on a land-owning MPI subcommunicator, so a
    // rank with zero local boxes doesn't join the collective nf90_create in
    // WriteLand/WriteRestart/ReadRestart. ERF_NOAHMP_IO.cpp.
    void with_land_comm (const std::function<void(NoahmpIO_type&)>& fn) const;

    // lsm data var count = fixed 2D fields + runtime soil profile (3 x m_nsoil).
    // Finalized once m_nsoil is known; the fixed-field count is a safe lower bound.
    int m_lsm_data_size = LsmData_NOAHMP::NumVars;

    // number of lsm flux variables
    int m_lsm_flux_size = LsmFlux_NOAHMP::NumVars;

    // LsmData map (state indices -> LsmVar enum)
    amrex::Vector<int> LsmDataMap;

    // LsmFlux map (state indices -> LsmFlux enum)
    amrex::Vector<int> LsmFluxMap;

    // Lsm data names
    amrex::Vector<std::string> LsmDataName;

    // Lsm flux names
    amrex::Vector<std::string> LsmFluxName;

    // geometry for atmosphere
    amrex::Geometry m_geom;

    // geometry for atmosphere at lev==0
    amrex::Geometry m_geom0;

    // geometry for lsm
    amrex::Geometry m_lsm_geom;

    // refinement ratio
    amrex::IntVect m_refRatio;

    // domain bc types
    amrex::Vector<amrex::BCRec> m_domain_bcs_type;

    // timestep
    amrex::Real m_dt;

    // Authoritative NoahMP substep counter, kept on the class so the firing decision
    // is identical on every rank (including land-free ranks with empty noahmpio_vect).
    int m_itimestep = 0;

    // NoahMP timestep (DTBL), broadcast to all ranks in Init so the firing
    // decision is well-defined on land-free ranks too.
    amrex::Real m_dtbl = std::numeric_limits<amrex::Real>::max();

    // Default -1; not written at runtime (single-slab LSM has no z extent to track).
    int khi_lsm = -1;

    // Number of grid points in z
    int m_nz_lsm = 1;

    // Size of grid spacing in z
    amrex::Real m_dz_lsm = one;

    // LSM data variables, sized at runtime to NumVars + 3*m_nsoil (fixed 2D fields +
    // the SMOIS/SH2O/TSLB soil profile). Vector, since m_nsoil is only known in Init.
    amrex::Vector<FabPtr> lsm_fab_data;

    // LSM flux variables
    amrex::Array<FabPtr, LsmFlux_NOAHMP::NumVars> lsm_fab_flux;

    // LSM data variables at level 0 (runtime-sized, see lsm_fab_data)
    amrex::Vector<amrex::MultiFab*> lsm_lev0_data;

    // LSM flux variables at level 0
    amrex::Array<amrex::MultiFab*, LsmFlux_NOAHMP::NumVars> lsm_lev0_flux {};

    // Number of Noah-MP soil layers, set from erf.lsm_nsoil (identical on every rank,
    // so no broadcast needed for the collective writes). Asserted vs namelist NSOIL.
    int m_nsoil = 4;

    // Vector to store NoahmpIO objects at a level
    NoahmpIO_vector noahmpio_vect;

    // Pinned memory buffers for host/device coupling to Noah-MP
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_input_tmp;
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_output_tmp;

    // Previous cumulative precip accumulator (native units) per slot per level;
    // interval precip is (current - this) * native_to_kg_m2. Only provided slots exist.
    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, NoahmpPrecipSlot::NumSlots>>
        m_precip_accum_prev;

    // Checkpointed previous-accumulation snapshots staged by Read_Lsm_Restart, consumed
    // (reset) by the first Advance. Empty on a cold start or legacy checkpoint.
    amrex::Array<std::unique_ptr<amrex::MultiFab>, NoahmpPrecipSlot::NumSlots> m_precip_accum_restored;

    // This instance's AMR level (set in Init); keys the per-level snapshot filenames.
    int m_lev = 0;

    // Update status
    bool m_updated = false;

    // Do we have a nc init file?
    bool m_has_nc_file = false;
};
#endif
