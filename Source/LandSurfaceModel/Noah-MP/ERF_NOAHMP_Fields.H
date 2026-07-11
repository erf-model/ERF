#ifndef ERF_NOAHMP_FIELDS_H
#define ERF_NOAHMP_FIELDS_H

// ===========================================================================
//  ERF <-> Noah-MP coupling contract. Single place declaring WHAT crosses the
//  boundary. See dev/spec-noahmp-api.md and dev/spec-noahmp-reorg.md.
//
//    - LsmData_NOAHMP / LsmFlux_NOAHMP : ERF-internal coupling fields (radiation,
//      surface layer). Persisted to the LSM plotfile / checkpoint BY NAME.
//    - NoahmpInputComp / NoahmpOutputComp : transient per-step staging-buffer
//      components used only inside Advance_With_State. Never persisted.
//    - NoahmpPrecipSlot : the typed precip accumulator slots.
//
//  Mechanical (1:1) fields live in the NOAHMP_*_FIELDS registries below. Each is
//  expanded with an "emitter" macro to generate, from one source of truth, the
//  enum entries, name strings, and host<->NoahmpIO copies -- so they cannot
//  drift. Fields with per-field math stay explicit in the .cpp.
//
//  To add a mechanical field: add one X(...) row to the relevant registry,
//  keeping its historical enum position (order is an invariant -- reordering
//  desyncs the pinned buffers). Nothing else needs to change.
// ===========================================================================

// ---------------------------------------------------------------------------
//  Field registries. Each X(...) row is one coupled field:
//    X(comp)          -- LsmData: ERF-internal, no NoahmpIO member (kept by name)
//    X(comp, member)  -- boundary: `comp` enum entry, `member` NoahmpIO array
// ---------------------------------------------------------------------------

// ERF-internal state fields (LsmData_NOAHMP), fixed 2D:
//   t_sfc..lw_flux_dn : radiation coupling
//   grdflx..smstot    : land return-term diagnostics
// The 3D soil profile (SMOIS/SH2O/TSLB) is appended at RUNTIME after NumVars,
// 3 groups of m_nsoil layers -- see NOAHMP::soil_data_idx.
#define NOAHMP_LSMDATA_FIELDS(X)                                               \
    /* radiation coupling */                                                  \
    X(t_sfc)             X(sfc_emis)                                          \
    X(sfc_alb_dir_vis)   X(sfc_alb_dir_nir)                                   \
    X(sfc_alb_dif_vis)   X(sfc_alb_dif_nir)                                   \
    X(cos_zenith_angle)  X(sw_flux_dn)                                        \
    X(sw_flux_dn_dir_vis) X(sw_flux_dn_dir_nir)                               \
    X(sw_flux_dn_dif_vis) X(sw_flux_dn_dif_nir)                               \
    X(lw_flux_dn)                                                             \
    /* land return-term diagnostics */                                       \
    X(grdflx) X(fira) X(sav) X(sag) X(albedo)                                 \
    X(sfcrunoff) X(udrunoff) X(smstav) X(smstot)

// Forcing ERF -> Noah-MP, mechanical 1:1 copies:  X( ERF comp, NoahmpIO member )
//   3D members use the k/j transpose noahmpio->MEMBER(i,1,j); 2D use MEMBER(i,j).
// Computed forcing (u_phy/v_phy, t_phy/p8w) is staged into the buffer in the
// device kernel; only the host copy to NoahmpIO is mechanical and generated here.
#define NOAHMP_INPUT_3D_FIELDS(X)                                             \
    /*    ERF comp      NoahmpIO member  */                                   \
    X(    u_phy,        U_PHY   )                                             \
    X(    v_phy,        V_PHY   )                                             \
    X(    t_phy,        T_PHY   )                                             \
    X(    qv_curr,      QV_CURR )                                             \
    X(    p8w,          P8W     )

#define NOAHMP_INPUT_2D_FIELDS(X)                                             \
    /*    ERF comp      NoahmpIO member  */                                   \
    X(    swdown,       SWDOWN  )                                             \
    X(    glw,          GLW     )                                             \
    X(    coszen,       COSZEN  )

// Results Noah-MP -> ERF, mechanical 1:1 copies:  X( ERF comp, NoahmpIO member )
// Split HEAD/TAIL only so the enum can interleave the non-mechanical banded
// albedos between them (see NoahmpOutputComp); the copy loop uses both together.
#define NOAHMP_OUTPUT_2D_FIELDS_HEAD(X)                                       \
    /*    ERF comp      NoahmpIO member  */                                   \
    X(    hfx,          HFX     )                                             \
    X(    lh,           LH      )                                             \
    X(    tau_ew,       TAU_EW  )                                             \
    X(    tau_ns,       TAU_NS  )                                             \
    X(    tsk,          TSK     )                                             \
    X(    emiss,        EMISS   )

#define NOAHMP_OUTPUT_2D_FIELDS_TAIL(X)                                       \
    /*    ERF comp      NoahmpIO member  */                                   \
    X(    o_grdflx,     GRDFLX    )                                           \
    X(    o_fira,       FIRAXY    )                                           \
    X(    o_sav,        SAVXY     )                                           \
    X(    o_sag,        SAGXY     )                                           \
    X(    o_albedo,     ALBEDO    )                                           \
    X(    o_sfcrunoff,  SFCRUNOFF )                                           \
    X(    o_udrunoff,   UDRUNOFF  )

#define NOAHMP_OUTPUT_2D_FIELDS(X) NOAHMP_OUTPUT_2D_FIELDS_HEAD(X)            \
                                   NOAHMP_OUTPUT_2D_FIELDS_TAIL(X)

// Noah-MP result -> ERF, mechanical 1:1. Each row ties the three index-aligned
// sites in the scatter kernel:  X( ERF Array4 alias, LsmData comp, Output comp ).
//   - alias   : local Array4<Real> bound from lsm_fab_data
//   - LsmData : persisted LsmData_NOAHMP field it writes into
//   - Output  : NoahmpOutputComp buffer component it reads from
// Drives the bind, the valid-cell copy, and the -9999 fill (lsm_undefined) branch.
// Banded albedos are listed for bind/copy/sentinel, but their NoahmpIO READ stays
// explicit; SMSTAV/SMSTOT are sentinels and are not listed.
#define NOAHMP_RESULT_FIELDS(X)                                               \
    /*   ERF Array4 alias  LsmData comp       Output comp    */               \
    X(   TSK,              t_sfc,             tsk           )                  \
    X(   EMISS,            sfc_emis,          emiss         )                  \
    X(   ALBSFCDIR_VIS,    sfc_alb_dir_vis,   albsfcdir_vis )                  \
    X(   ALBSFCDIR_NIR,    sfc_alb_dir_nir,   albsfcdir_nir )                  \
    X(   ALBSFCDIF_VIS,    sfc_alb_dif_vis,   albsfcdif_vis )                  \
    X(   ALBSFCDIF_NIR,    sfc_alb_dif_nir,   albsfcdif_nir )                  \
    X(   GRDFLX_o,         grdflx,            o_grdflx      )                  \
    X(   FIRA_o,           fira,              o_fira        )                  \
    X(   SAV_o,            sav,               o_sav         )                  \
    X(   SAG_o,            sag,               o_sag         )                  \
    X(   ALBEDO_o,         albedo,            o_albedo      )                  \
    X(   SFCRUNOFF_o,      sfcrunoff,         o_sfcrunoff   )                  \
    X(   UDRUNOFF_o,       udrunoff,          o_udrunoff    )

// ---------------------------------------------------------------------------
//  Reusable field emitters -- the per-row "action" applied when a registry is
//  expanded, e.g. NOAHMP_LSMDATA_FIELDS(NOAHMP_QUOTE) -> "t_sfc", "sfc_emis", ...
//  Copy-loop emitters reference local buffers and live next to each loop in the .cpp.
#define NOAHMP_ENUM(comp)          comp,    // (comp)         -> enum entry
#define NOAHMP_QUOTE(comp)         #comp,   // (comp)         -> "comp" literal
#define NOAHMP_COMP(comp, member)  comp,    // (comp, member) -> enum entry

// ---------------------------------------------------------------------------
//  Enums generated from the registries (plus fixed non-mechanical entries).
//  Component order is preserved exactly.
// ---------------------------------------------------------------------------

namespace LsmData_NOAHMP {
   enum {
         NOAHMP_LSMDATA_FIELDS(NOAHMP_ENUM)
         NumVars
         // The 3D soil profile (SMOIS/SH2O/TSLB) is appended at RUNTIME after
         // NumVars, three contiguous groups of m_nsoil layers (any NSOIL):
         //   [NumVars           .. NumVars +   nsoil) -> smois
         //   [NumVars +   nsoil .. NumVars + 2*nsoil) -> sh2o
         //   [NumVars + 2*nsoil .. NumVars + 3*nsoil) -> tslb
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
// SurfacePrecipAccumulationSources. Keys the per-slot previous-accumulation
// snapshots and the water-equivalent interval deltas fed to Noah-MP.
namespace NoahmpPrecipSlot {
    enum {
        total = 0,   // total non-convective precip (if the scheme exposes one)
        rain,        // explicit rain (schemes without a total slot)
        snow,        // snow (+ice)
        graupel,     // graupel
        hail,        // folded into MP_GRAUP (coupled MP_HAIL set to 0)
        NumSlots
    };
}

namespace NoahmpInputComp {
    enum {
        NOAHMP_INPUT_3D_FIELDS(NOAHMP_COMP)
        NOAHMP_INPUT_2D_FIELDS(NOAHMP_COMP)
        // computed precip staging (explicit, not table-driven)
        rainbl,      // accumulated precip over the interval [mm]
        sr_in,       // frozen fraction [-]
        mp_rainnc,   // total non-convective precip this interval [mm]
        mp_snow,     // snow (+ice) subset of mp_rainnc [mm]
        mp_graup,    // graupel subset of mp_rainnc [mm]
        NumComps
    };
}

namespace NoahmpOutputComp {
    enum {
        NOAHMP_OUTPUT_2D_FIELDS_HEAD(NOAHMP_COMP)
        // banded albedos (explicit: read at band 1=VIS / 2=NIR)
        albsfcdir_vis,
        albsfcdir_nir,
        albsfcdif_vis,
        albsfcdif_nir,
        NOAHMP_OUTPUT_2D_FIELDS_TAIL(NOAHMP_COMP)
        // sentinels: not computed by this core -> emitted as lsm_undefined
        o_smstav,
        o_smstot,
        NumComps
        // Per-layer soil outputs (SMOIS/SH2O/TSLB) are appended at RUNTIME after
        // NumComps, mirroring the LsmData_NOAHMP layout; the pinned output buffer
        // is sized NumComps + 3*nsoil to match.
    };
}

#endif /* ERF_NOAHMP_FIELDS_H */
