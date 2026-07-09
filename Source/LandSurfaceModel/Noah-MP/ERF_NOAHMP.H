#ifndef NOAHMP_H
#define NOAHMP_H

#include <string>
#include <vector>
#include <memory>
#include <limits>

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

// External include from the noahmp library
#include <NoahmpIO.H>

namespace LsmData_NOAHMP {
   enum {
         // 2D RRTMGP coupling vars
         t_sfc = 0,
         sfc_emis,
         sfc_alb_dir_vis,
         sfc_alb_dir_nir,
         sfc_alb_dif_vis,
         sfc_alb_dif_nir,
         cos_zenith_angle,
         sw_flux_dn,
         sw_flux_dn_dir_vis,
         sw_flux_dn_dir_nir,
         sw_flux_dn_dif_vis,
         sw_flux_dn_dif_nir,
         lw_flux_dn,
         // 2D land return-term diagnostics (read from Noah-MP, output only -- NOT
         // fed back into ERF dynamics; matches how WRF carries these). Auto-written
         // to the LSM plotfile + checkpoint via Plot_Lsm_Data / lsm_data checkpoint.
         grdflx,          // ground/snow heat flux [W m-2]
         fira,            // total net longwave radiation [W m-2] (FIRAXY)
         sav,             // solar radiation absorbed by vegetation [W m-2] (SAVXY)
         sag,             // solar radiation absorbed by ground [W m-2] (SAGXY)
         albedo,          // total grid surface albedo [-] (ALBEDO)
         sfcrunoff,       // accumulated surface runoff [m] (SFCRUNOFF)
         udrunoff,        // accumulated sub-surface runoff [m] (UDRUNOFF)
         smstav,          // soil moisture availability [-] (SMSTAV; NOT computed by this
                          // Noah-MP core -> emitted as lsm_undefined placeholder, not 0)
         smstot,          // total soil moisture [mm] (SMSTOT; NOT computed by this
                          // Noah-MP core -> emitted as lsm_undefined placeholder, not 0)
         NumVars
         // NOTE: the 3D soil-profile state (SMOIS/SH2O/TSLB) is NOT enumerated here.
         // It is appended at RUNTIME, one 2D field per soil layer, in three
         // contiguous groups of m_nsoil starting at index NumVars:
         //   [NumVars           .. NumVars +   nsoil) -> smois_1..smois_nsoil
         //   [NumVars +   nsoil .. NumVars + 2*nsoil) -> sh2o_1 ..sh2o_nsoil
         //   [NumVars + 2*nsoil .. NumVars + 3*nsoil) -> tslb_1 ..tslb_nsoil
         // so the diagnostics work for any Noah-MP soil-layer count, not just 4.
         // See NOAHMP::soil_data_idx / NOAHMP::soil_out_idx.
  };
}

namespace LsmFlux_NOAHMP {
   enum {
         // 2D SurfaceLayer flux vars
         t_flux = 0,
         q_flux,
         tau13,
         tau23,
         NumVars
  };
}


// Fixed index for each typed surface-precip source slot, in the field order of
// SurfacePrecipAccumulationSources. Used to key the per-slot previous-accumulation
// snapshots and the water-equivalent interval deltas fed to Noah-MP.
namespace NoahmpPrecipSlot {
    enum {
        total = 0,   // total non-convective precip accumulator (if the scheme exposes one)
        rain,        // explicit rain accumulator (schemes without a total slot)
        snow,        // snow (+ice) accumulator
        graupel,     // graupel accumulator
        hail,        // hail accumulator: intentionally folded into MP_GRAUP (this Noah-MP
                     // core has no hail input; not exact WRF MP_HAIL parity -- see ERF_NOAHMP.cpp)
        NumSlots
    };
}

namespace NoahmpInputComp {
    enum {
        u_phy = 0,
        v_phy,
        t_phy,
        qv_curr,
        p8w,
        swdown,
        glw,
        coszen,
        rainbl,      // accumulated precip into land model over the interval [mm]
        sr_in,       // frozen fraction of precip [-]
        mp_rainnc,   // total non-convective precip this interval [mm] (WRF MP_RAINNC)
        mp_snow,     // snow (+ice) subset of mp_rainnc [mm] (WRF MP_SNOW)
        mp_graup,    // graupel subset of mp_rainnc [mm] (WRF MP_GRAUP)
        NumComps
    };
}

namespace NoahmpOutputComp {
    enum {
        hfx = 0,
        lh,
        tau_ew,
        tau_ns,
        tsk,
        emiss,
        albsfcdir_vis,
        albsfcdir_nir,
        albsfcdif_vis,
        albsfcdif_nir,
        // land return-term diagnostics (output only)
        o_grdflx,
        o_fira,
        o_sav,
        o_sag,
        o_albedo,
        o_sfcrunoff,
        o_udrunoff,
        o_smstav,
        o_smstot,
        NumComps
        // NOTE: the per-layer soil outputs (SMOIS/SH2O/TSLB) are appended at
        // RUNTIME after NumComps, mirroring the LsmData_NOAHMP layout above:
        //   [NumComps           .. NumComps +   nsoil) -> smois layers
        //   [NumComps +   nsoil .. NumComps + 2*nsoil) -> sh2o  layers
        //   [NumComps + 2*nsoil .. NumComps + 3*nsoil) -> tslb  layers
        // The pinned host output buffer is sized NumComps + 3*nsoil to match.
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
        // NOTE: We should parse constants from sc here if needed,
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

    // Get number of vars lsm class contains. The parent queries this BEFORE
    // NOAHMP::Init runs (to size its lsm_data array in MakeNewLevel), so resolve
    // the soil-layer count here on first use rather than relying on Init.
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
        lsm_lev0_data[varIdx] = lev0_data_ptr;
    }

    void
    Lsm_Set_Lev0_Flux_Ptr (const int& varIdx,
                           amrex::MultiFab* lev0_flux_ptr) override
    {
        lsm_lev0_flux[varIdx] = lev0_flux_ptr;
    }

    bool
    Get_Update_Status () override { return m_updated; }

    // Checkpoint/restore the NoahMP substep counter. Uses the class member
    // m_itimestep (valid on every rank) rather than noahmpio_vect[0] (absent on
    // land-free ranks).
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

    // Write/read the full NoahMP prognostic state (soil, snow incl. layers,
    // canopy, aquifer, albedo history, ...) to/from a NetCDF restart directory
    // so a restart reproduces a cold-start trajectory bitwise (issue #3255).
    // Each local block writes its tile into the global-domain file collectively.
    void
    Write_Lsm_Restart (const std::string& dir) const override
    {
        for (const NoahmpIO_type& noahmpio : noahmpio_vect) {
            const_cast<NoahmpIO_type&>(noahmpio).WriteRestart(dir);
        }

        // Also save the precip snapshots (the cumulative precip Noah-MP last consumed).
        // Restoring these on restart lets the next land call difference against the same
        // baseline, so no precip is dropped or double-counted across the restart.
        if (m_lev < int(m_precip_accum_prev.size())) {
            for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
                if (m_precip_accum_prev[m_lev][s]) {
                    amrex::VisMF::Write(*m_precip_accum_prev[m_lev][s],
                                        m_precip_snapshot_name(dir, m_lev, s));
                }
            }
        }
    }

    void
    Read_Lsm_Restart (const std::string& dir) override
    {
        for (NoahmpIO_type& noahmpio : noahmpio_vect) {
            noahmpio.ReadRestart(dir);
        }

        // Read the saved precip snapshots into a staging store (the live snapshots are
        // not allocated until the first Advance, which then copies these in). Legacy
        // checkpoints have no such files -> store stays empty -> Advance cold-seeds.
        for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
            const std::string name = m_precip_snapshot_name(dir, m_lev, s);
            if (amrex::FileExists(name + "_H")) {
                auto mf = std::make_unique<amrex::MultiFab>();
                amrex::VisMF::Read(*mf, name);
                m_precip_accum_restored[s] = std::move(mf);
            }
        }
    }

private:
    // number of lsm data variables = fixed 2D fields (LsmData_NOAHMP::NumVars) plus
    // the runtime soil profile (3 groups x m_nsoil layers). Finalized in Init once
    // m_nsoil is known; the fixed-field count is a safe lower bound before then.
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

    // Authoritative NoahMP substep counter, kept on the class so the firing
    // decision in Advance_With_State is identical on every rank (including
    // land-free ranks, whose noahmpio_vect is empty).
    int m_itimestep = 0;

    // NoahMP timestep (DTBL), broadcast to all ranks in Init so the firing
    // decision is well-defined on land-free ranks too.
    amrex::Real m_dtbl = std::numeric_limits<amrex::Real>::max();

    // domain klo-1 or lsm khi
    int khi_lsm;

    // Number of grid points in z
    int m_nz_lsm = 1;

    // Size of grid spacing in z
    amrex::Real m_dz_lsm = one;

    // LSM data variables. Sized at runtime to LsmData_NOAHMP::NumVars + 3*m_nsoil
    // (the fixed 2D fields plus the per-layer SMOIS/SH2O/TSLB soil profile), so the
    // soil diagnostics work for any Noah-MP soil-layer count. A Vector (not a
    // fixed-size Array) is required because m_nsoil is only known after the
    // namelist is read in Init.
    amrex::Vector<FabPtr> lsm_fab_data;

    // LSM flux variables
    amrex::Array<FabPtr, LsmFlux_NOAHMP::NumVars> lsm_fab_flux;

    // LSM data variables at level 0 (runtime-sized, see lsm_fab_data)
    amrex::Vector<amrex::MultiFab*> lsm_lev0_data;

    // LSM flux variables at level 0
    amrex::Array<amrex::MultiFab*, LsmFlux_NOAHMP::NumVars> lsm_lev0_flux;

    // Number of Noah-MP soil layers (NSOIL). Set in Init from erf.lsm_nsoil, which
    // is identical on every rank (including land-free ranks) so all agree on the LSM
    // data size without a broadcast -- required because the plotfile/checkpoint
    // writes are collective. Asserted against the namelist NSOIL on land ranks.
    // Default 4 = the standard Noah-MP configuration.
    int m_nsoil = 4;

    // Index of soil field (group g in {0:smois,1:sh2o,2:tslb}, layer k in [0,nsoil))
    // within lsm_fab_data / LsmDataName. The matching pinned-output-buffer index is
    // NoahmpOutputComp::NumComps + g*m_nsoil + k (computed inline in the scatter
    // kernel, since that runs on the device where this host member is not callable).
    int soil_data_idx (int g, int k) const {
        return LsmData_NOAHMP::NumVars + g*m_nsoil + k;
    }

    // Resolve the soil-layer count (erf.lsm_nsoil, default 4) and the total LSM data
    // size exactly once, independent of Init. Must be callable before Init because
    // the parent (ERF_MakeNewLevel) queries Lsm_Data_Size() to size its lsm_data
    // array first; Init later asserts m_nsoil matches the namelist NSOIL.
    bool m_nsoil_resolved = false;
    void m_ensure_nsoil_resolved () {
        if (m_nsoil_resolved) { return; }
        amrex::ParmParse pp_nsoil("erf");
        pp_nsoil.query("lsm_nsoil", m_nsoil);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_nsoil >= 1, "erf.lsm_nsoil must be >= 1");
        m_lsm_data_size = LsmData_NOAHMP::NumVars + 3*m_nsoil;
        m_nsoil_resolved = true;
    }

    // Vector to store NoahmpIO objects at a level
    NoahmpIO_vector noahmpio_vect;

    // Pinned memory buffers for host/device coupling to Noah-MP
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_input_tmp;
    amrex::Vector<std::unique_ptr<amrex::FArrayBox>> noahmp_output_tmp;

    // Previous cumulative surface precip accumulator (scheme-native units) at the last
    // land-model call, one snapshot per typed source slot (NoahmpPrecipSlot), per level.
    // The interval precip fed to Noah-MP is (current - this) * native_to_kg_m2, i.e. the
    // native->kg/m^2 (water-equivalent mm) factor is applied when the delta is formed.
    // Only slots the active scheme provides are allocated; the rest stay null. Synced to
    // current on the first call and on restart so a restored accumulation does not
    // produce a spurious precip spike.
    amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, NoahmpPrecipSlot::NumSlots>>
        m_precip_accum_prev;

    // Checkpointed previous-accumulation snapshots staged by Read_Lsm_Restart, consumed
    // (reset) by the first Advance. Empty on a cold start or legacy checkpoint.
    amrex::Array<std::unique_ptr<amrex::MultiFab>, NoahmpPrecipSlot::NumSlots>
        m_precip_accum_restored;

    // This instance's AMR level (set in Init); keys the per-level snapshot filenames.
    int m_lev = 0;

    // VisMF file prefix for the level-l, slot-s previous-accumulation snapshot.
    static std::string
    m_precip_snapshot_name (const std::string& dir, int l, int s) {
        return dir + "/precip_accum_prev_Level_" + std::to_string(l)
                   + "_slot_" + std::to_string(s);
    }

    // Plot frequency
    int m_plot_int_1 = -1;

    // Update status
    bool m_updated = false;

    // Do we have a nc init file?
    bool m_has_nc_file = false;
};
#endif
