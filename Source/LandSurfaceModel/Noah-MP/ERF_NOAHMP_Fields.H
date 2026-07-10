#ifndef ERF_NOAHMP_FIELDS_H
#define ERF_NOAHMP_FIELDS_H

// ===========================================================================
//  ERF <-> Noah-MP coupling contract
// ---------------------------------------------------------------------------
//  This header is the single place that declares WHAT crosses the ERF/Noah-MP
//  boundary. Four (well, five) enums describe TWO different boundaries; see
//  dev/spec-noahmp-api.md and dev/spec-noahmp-reorg.md.
//
//    - LsmData_NOAHMP / LsmFlux_NOAHMP : ERF-internal coupling fields, read
//      and written by the rest of ERF (radiation, surface layer). Persisted to
//      the LSM plotfile / checkpoint BY NAME.
//    - NoahmpInputComp / NoahmpOutputComp : transient per-step staging-buffer
//      components used only inside Advance_With_State. Never persisted.
//    - NoahmpPrecipSlot : the typed precip accumulator slots.
//
//  Mechanical (pure 1:1) coupled fields are declared ONCE below as X-macro
//  lists (NOAHMP_*_FIELDS). Expanding a list with different X definitions
//  generates the enum entries, the name strings, and the host<->NoahmpIO copy
//  lines from one source of truth, so they cannot drift. Fields that carry
//  per-field math (winds, EOS, precip, banded albedos, fluxes, soil) are NOT
//  table-driven and stay explicit in the .cpp -- see the "NOT table-driven"
//  note in dev/spec-noahmp-reorg.md §4.
// ===========================================================================

// ---------------------------------------------------------------------------
//  Field registries (X-macro lists). Each entry is a coupled field.
// ---------------------------------------------------------------------------

// ERF-internal state fields (LsmData_NOAHMP), fixed 2D. One token per field;
// the enum entry and the plotfile/checkpoint name string are both generated
// from this list, so the enum order and the names can never disagree.
//   t_sfc,sfc_emis,albedos,cos_zenith,downwelling fluxes  : radiation coupling
//   grdflx..smstot                                        : land return-term
//                                                           diagnostics (output)
// The 3D soil profile (SMOIS/SH2O/TSLB) is appended at RUNTIME after NumVars,
// 3 groups of m_nsoil layers -- see NOAHMP::soil_data_idx. Not listed here.
#define NOAHMP_LSMDATA_FIELDS(X)                                               \
    X(t_sfc) X(sfc_emis)                                                       \
    X(sfc_alb_dir_vis) X(sfc_alb_dir_nir)                                      \
    X(sfc_alb_dif_vis) X(sfc_alb_dif_nir)                                      \
    X(cos_zenith_angle) X(sw_flux_dn)                                          \
    X(sw_flux_dn_dir_vis) X(sw_flux_dn_dir_nir)                                \
    X(sw_flux_dn_dif_vis) X(sw_flux_dn_dif_nir)                                \
    X(lw_flux_dn)                                                              \
    X(grdflx) X(fira) X(sav) X(sag) X(albedo)                                  \
    X(sfcrunoff) X(udrunoff) X(smstav) X(smstot)

// Forcing ERF -> Noah-MP, mechanical 1:1 copies. (comp-name, NoahmpIO member)
//   3D atmospheric members use the k/j transpose: noahmpio->MEMBER(i,1,j).
//   2D surface members use noahmpio->MEMBER(i,j).
// The computed forcing (u_phy/v_phy via face-average, t_phy/p8w via EOS) is
// staged into the input buffer explicitly in the device kernel; its host copy
// to the 3D members is still mechanical and IS generated from NOAHMP_INPUT_3D.
#define NOAHMP_INPUT_3D(X)                                                     \
    X(u_phy,U_PHY) X(v_phy,V_PHY) X(t_phy,T_PHY)                               \
    X(qv_curr,QV_CURR) X(p8w,P8W)
#define NOAHMP_INPUT_2D(X)                                                     \
    X(swdown,SWDOWN) X(glw,GLW) X(coszen,COSZEN)

// Results Noah-MP -> ERF, mechanical 1:1 copies (member -> output buffer).
// Split HEAD/TAIL only so the enum can interleave the banded albedos between
// them in their historical positions; the copy loop uses HEAD+TAIL together
// (copy order is irrelevant -- each assignment is independent).
#define NOAHMP_OUTPUT_2D_HEAD(X)                                               \
    X(hfx,HFX) X(lh,LH) X(tau_ew,TAU_EW) X(tau_ns,TAU_NS)                      \
    X(tsk,TSK) X(emiss,EMISS)
#define NOAHMP_OUTPUT_2D_TAIL(X)                                               \
    X(o_grdflx,GRDFLX) X(o_fira,FIRAXY) X(o_sav,SAVXY) X(o_sag,SAGXY)          \
    X(o_albedo,ALBEDO) X(o_sfcrunoff,SFCRUNOFF) X(o_udrunoff,UDRUNOFF)
#define NOAHMP_OUTPUT_2D(X) NOAHMP_OUTPUT_2D_HEAD(X) NOAHMP_OUTPUT_2D_TAIL(X)

// ---------------------------------------------------------------------------
//  Enums generated from the registries (plus the non-mechanical entries that
//  stay in their fixed positions). Component order is preserved exactly.
// ---------------------------------------------------------------------------

namespace LsmData_NOAHMP {
   enum {
#define NOAHMP_X(name) name,
         NOAHMP_LSMDATA_FIELDS(NOAHMP_X)
#undef NOAHMP_X
         NumVars
         // NOTE: the 3D soil-profile state (SMOIS/SH2O/TSLB) is NOT enumerated
         // here. It is appended at RUNTIME, one 2D field per soil layer, in three
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
        hail,        // hail accumulator: folded into MP_GRAUP (coupled MP_HAIL set to 0) -- see ERF_NOAHMP_Advance.cpp
        NumSlots
    };
}

namespace NoahmpInputComp {
    enum {
#define NOAHMP_X(name,member) name,
        NOAHMP_INPUT_3D(NOAHMP_X)
        NOAHMP_INPUT_2D(NOAHMP_X)
#undef NOAHMP_X
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
#define NOAHMP_X(name,member) name,
        NOAHMP_OUTPUT_2D_HEAD(NOAHMP_X)
#undef NOAHMP_X
        // banded albedos (NOT table-driven: read at band 1=VIS / 2=NIR)
        albsfcdir_vis,
        albsfcdir_nir,
        albsfcdif_vis,
        albsfcdif_nir,
#define NOAHMP_X(name,member) name,
        NOAHMP_OUTPUT_2D_TAIL(NOAHMP_X)
#undef NOAHMP_X
        // sentinels: not computed by this Noah-MP core -> emitted as lsm_undefined
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

#endif /* ERF_NOAHMP_FIELDS_H */
