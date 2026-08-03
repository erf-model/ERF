#ifndef ERF_NOAHMP_FIELDS_H
#define ERF_NOAHMP_FIELDS_H

// Append X(...) rows to add coupling fields; enum order is an invariant -- do not reorder.

// ERF-internal state fields (LsmData_NOAHMP), fixed 2D. The 3D soil profile
// (SMOIS/SH2O/TSLB) is appended at RUNTIME after NumVars -- see soil_data_idx.
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
    X(sfcrunoff) X(udrunoff) X(smstav) X(smstot)                              \
    /* native Noah-MP 2-m bundle */                                           \
    X(noahmp_temperature_2m_vegetated)                                       \
    X(noahmp_temperature_2m_bare)                                            \
    X(noahmp_water_vapor_mixing_ratio_2m_vegetated)                           \
    X(noahmp_water_vapor_mixing_ratio_2m_bare)                                \
    X(noahmp_vegetation_fraction)

// Forcing ERF -> Noah-MP, mechanical 1:1 copies: X( ERF comp, NoahmpIO member ).
// 3D members use the k/j transpose MEMBER(i,1,j); 2D use MEMBER(i,j).
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

// Results Noah-MP -> ERF, mechanical 1:1 copies: X( ERF comp, NoahmpIO member ).
// Split HEAD/TAIL so the enum can interleave the banded albedos; copy loop uses both.
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
    X(    o_udrunoff,   UDRUNOFF  )                                           \
    X(    o_noahmp_t2m_veg,  T2MVXY )                                         \
    X(    o_noahmp_t2m_bare, T2MBXY )                                         \
    X(    o_noahmp_q2m_veg,  Q2MVXY )                                         \
    X(    o_noahmp_q2m_bare, Q2MBXY )                                         \
    X(    o_noahmp_fveg,     FVEGXY )

#define NOAHMP_OUTPUT_2D_FIELDS(X) NOAHMP_OUTPUT_2D_FIELDS_HEAD(X)            \
                                   NOAHMP_OUTPUT_2D_FIELDS_TAIL(X)

// Noah-MP result -> ERF, mechanical 1:1: X( ERF Array4 alias, LsmData comp, Output
// comp ). Banded albedos are listed but read explicitly; SMSTAV/SMSTOT are not listed.
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
    X(   UDRUNOFF_o,       udrunoff,          o_udrunoff    )                  \
    X(   T2MVXY_o,         noahmp_temperature_2m_vegetated, o_noahmp_t2m_veg ) \
    X(   T2MBXY_o,         noahmp_temperature_2m_bare,      o_noahmp_t2m_bare) \
    X(   Q2MVXY_o,         noahmp_water_vapor_mixing_ratio_2m_vegetated, o_noahmp_q2m_veg ) \
    X(   Q2MBXY_o,         noahmp_water_vapor_mixing_ratio_2m_bare,      o_noahmp_q2m_bare) \
    X(   FVEGXY_o,         noahmp_vegetation_fraction, o_noahmp_fveg       )

// ---------------------------------------------------------------------------
//  Reusable field emitters -- the per-row action applied when a registry is
//  expanded, e.g. NOAHMP_LSMDATA_FIELDS(NOAHMP_QUOTE) -> "t_sfc", "sfc_emis", ...
#define NOAHMP_ENUM(comp)          comp,    // (comp)         -> enum entry
#define NOAHMP_QUOTE(comp)         #comp,   // (comp)         -> "comp" literal
#define NOAHMP_COMP(comp, member)  comp,    // (comp, member) -> enum entry

// ---------------------------------------------------------------------------
//  Enums generated from the registries (plus fixed entries); order is preserved.
// ---------------------------------------------------------------------------

namespace LsmData_NOAHMP {
   enum {
         NOAHMP_LSMDATA_FIELDS(NOAHMP_ENUM)
         NumVars
         // Soil profile (SMOIS/SH2O/TSLB) appended at RUNTIME after NumVars as three
         // contiguous groups of m_nsoil layers. See NOAHMP::soil_data_idx.
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

// Fixed index per typed surface-precip source slot (field order of
// SurfacePrecipAccumulationSources), keying the snapshots and interval deltas.
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
        // Per-layer soil outputs appended at RUNTIME after NumComps; the output
        // buffer is sized NumComps + 3*nsoil to match.
    };
}

#endif /* ERF_NOAHMP_FIELDS_H */
