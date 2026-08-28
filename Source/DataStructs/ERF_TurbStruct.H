#ifndef ERF_TURB_STRUCT_H_
#define ERF_TURB_STRUCT_H_

#include <ERF_MYNNStruct.H>

/**
 * @brief Large-eddy simulation closure type.
 */
AMREX_ENUM(LESType, None, Smagorinsky, Smagorinsky2D, Deardorff);

/**
 * @brief Reynolds-averaged turbulence closure type.
 */
AMREX_ENUM(RANSType, None, kEqn);

/**
 * @brief Planetary boundary-layer closure type.
 */
AMREX_ENUM(PBLType, None, MYJ, MYNN25, MYNNEDMF, YSU, YSUNew, MRF, SHOC, EAMXX_SHOC, NATIVE_SHOC);

/**
 * @brief Thermodynamic variable used for stability stratification.
 */
AMREX_ENUM(StratType, theta, thetav, thetal);

/**
 * @brief Query a scalar or per-level input value.
 * @tparam T Value type accepted by ParmParse::query.
 * @param pp ParmParse object to query.
 * @param query_string Input key to query.
 * @param query_var Destination value for the requested level.
 * @param lev AMR level index.
 * @param maxlev Maximum AMR level configured for the run.
 * @return ParmParse query result, or zero if the key is absent.
 */
template <typename T>
int
query_one_or_per_level (
  const amrex::ParmParse& pp,
  const char* query_string,
  T& query_var,
  const int lev,
  const int maxlev)
{
  int count = pp.countval(query_string);
  if (count == 0) {
    return 0; // nothing to do
  } else if (count == 1) {
    // In this case, we assume the same value is used at every level
    return pp.query(query_string, query_var);
  } else if (count >= maxlev + 1) {
    // In this case, there may be more values than levels so we
    // will just read the number of values we need and ignore the
    // extra levels
    return pp.query(query_string, query_var, lev);
  } else {
    // In this case, there is more than one value AND there
    // are more levels than values; as an example, if max_level = 2
    // but count = 2 so we have three levels but only two values --
    // it is not clear how to interpret this so we abort
    amrex::Error(
      "For parmparse variable " + pp.prefixedName(query_string) +
      ": if specified, specify once total or at least once for each level");
    return 0; // avoid compiler warning
  }
}

/**
 * @brief Query a scalar or per-level enum input value using case-insensitive matching.
 * @tparam T Enum type accepted by ParmParse::query_enum_case_insensitive.
 * @param pp ParmParse object to query.
 * @param query_string Input key to query.
 * @param query_var Destination enum value for the requested level.
 * @param lev AMR level index.
 * @param maxlev Maximum AMR level configured for the run.
 * @return ParmParse query result, or zero if the key is absent.
 */
template <typename T>
int
query_one_or_per_level_enum_case_insensitive (
  const amrex::ParmParse& pp,
  const char* query_string,
  T& query_var,
  const int lev,
  const int maxlev)
{
  int count = pp.countval(query_string);
  if (count == 0) {
    return 0; // nothing to do
  } else if (count == 1) {
    // In this case, we assume the same value is used at every level
    return pp.query_enum_case_insensitive(query_string, query_var);
  } else if (count >= maxlev + 1) {
    // In this case, there may be more values than levels so we
    // will just read the number of values we need and ignore the
    // extra levels
    return pp.query_enum_case_insensitive(query_string, query_var, lev);
  } else {
    // In this case, there is more than one value AND there
    // are more levels than values; as an example, if max_level = 2
    // but count = 2 so we have three levels but only two values --
    // it is not clear how to interpret this so we abort
    amrex::Error(
      "For parmparse variable " + pp.prefixedName(query_string) +
      ": if specified, specify once total or at least once for each level");
    return 0; // avoid compiler warning
  }
}

/**
 * Container holding quantities related to turbulence parametrizations
 */
struct TurbChoice
{
public:
  /**
   * @brief Read turbulence options for one AMR level from the input parameter database.
   * @param lev AMR level index.
   * @param max_level Maximum AMR level configured for the run.
   * @param pp_prefix ParmParse prefix for the ERF input namespace.
   */
  void init_params (int lev, int max_level, std::string pp_prefix)
  {
    amrex::ParmParse pp(pp_prefix);

    // Which LES closure?
    query_one_or_per_level(pp, "les_type", les_type, lev, max_level);

    // Handle 2-D Smag
    if (les_type == LESType::Smagorinsky2D) {
        les_type = LESType::Smagorinsky;
        smag2d = true;
    }

    // Which RANS closure?
    query_one_or_per_level(pp, "rans_type", rans_type, lev, max_level);

    if ((rans_type != RANSType::None) && (les_type != LESType::None)) {
        amrex::Error("Hybrid RANS-LES not implemented");
    }

    // Which PBL Closure
    query_one_or_per_level_enum_case_insensitive(pp, "pbl_type", pbl_type, lev, max_level);
    if (pbl_type == PBLType::SHOC) {
      static bool warned_legacy_shoc = false;
      if (!warned_legacy_shoc) {
        amrex::Warning("erf.pbl_type = SHOC is deprecated; use erf.pbl_type = EAMXX_SHOC");
        warned_legacy_shoc = true;
      }
      pbl_type = PBLType::EAMXX_SHOC;
    }

    // Do some more stuff for PBL Modeling
    if (pbl_type != PBLType::None) {
      // Check for compatibility between PBL, LES, Molec Transport
      if (les_type != LESType::None) {
        amrex::Print() << "Selected a PBL model and an LES model: "
                       << "Using PBL for vertical transport, LES for horizontal"
                       << std::endl;
      }
      if (les_type == LESType::Smagorinsky) {
        if (!smag2d)
          amrex::Error("If using Smagorinsky with a PBL model, the 2-D "
                       "formulation should be used");
      } else if (les_type == LESType::Deardorff) {
        amrex::Error(
          "It is not recommended to use Deardorff LES and a PBL model");
      }

      if (pbl_type == PBLType::MYNN25 || pbl_type == PBLType::MYNNEDMF) {
        query_one_or_per_level(pp, "pbl_mynn_A1", pbl_mynn.A1, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_A2", pbl_mynn.A2, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_B1", pbl_mynn.B1, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_B2", pbl_mynn.B2, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_C1", pbl_mynn.C1, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_C2", pbl_mynn.C2, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_C3", pbl_mynn.C3, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_C4", pbl_mynn.C4, lev, max_level);
        query_one_or_per_level(pp, "pbl_mynn_C5", pbl_mynn.C5, lev, max_level);
        pbl_mynn_level2.init_coeffs(
          pbl_mynn.A1, pbl_mynn.A2, pbl_mynn.B1, pbl_mynn.B2, pbl_mynn.C1,
          pbl_mynn.C2, pbl_mynn.C3, pbl_mynn.C4, pbl_mynn.C5);
        query_one_or_per_level(
          pp, "pbl_mynn_diffuse_moistvars", pbl_mynn.diffuse_moistvars, lev,
          max_level);
        query_one_or_per_level(
          pp, "pbl_mynn_SMmin", pbl_mynn.SMmin, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mynn_SMmax", pbl_mynn.SMmax, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mynn_SHmin", pbl_mynn.SHmin, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mynn_SHmax", pbl_mynn.SHmax, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mynn_SQfactor", pbl_mynn.SQfac, lev, max_level);

        // SQ = SQfac*SM, so its bound is tied to SHmax and SQfac; the default member
        // initializer was evaluated before either could be overridden here.
        pbl_mynn.SQmax = pbl_mynn.SQfac * pbl_mynn.SHmax;

      } else if (pbl_type == PBLType::YSU || pbl_type == PBLType::YSUNew) {
        query_one_or_per_level(
          pp, "pbl_ysu_coriolis_freq", pbl_ysu_coriolis_freq, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_ysu_use_consistent_coriolis",
          pbl_ysu_use_consistent_coriolis, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_ysu_force_over_water", pbl_ysu_force_over_water, lev,
          max_level);
        query_one_or_per_level(
          pp, "pbl_ysu_land_Ribcr", pbl_ysu_land_Ribcr, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_ysu_unst_Ribcr", pbl_ysu_unst_Ribcr, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_liquid_theta", enable_ysu_liquid_theta, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_countergradient", enable_ysu_countergradient, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_terrain_pblh_floor", enable_ysu_terrain_pblh_floor, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_sat_limiter", enable_ysu_sat_limiter, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_topdown", enable_ysu_topdown, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_entrainment", enable_ysu_entrainment, lev, max_level);
        query_one_or_per_level(
          pp, "enable_ysu_cloud_pblh", enable_ysu_cloud_pblh, lev, max_level);
        query_one_or_per_level(
          pp, "enable_mrf_unbounded_vpert", enable_mrf_unbounded_vpert, lev, max_level);
        query_one_or_per_level(
          pp, "ysu_qcloud_threshold", ysu_qcloud_threshold, lev, max_level);
        query_one_or_per_level(
          pp, "ysu_moistvars", ysu_moistvars, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_ysunew_highres_bounds", pbl_ysunew_highres_bounds, lev, max_level);
      } else if (pbl_type == PBLType::MRF) {
        query_one_or_per_level(
          pp, "pbl_mrf_coriolis_freq", pbl_mrf_coriolis_freq, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mrf_Ribcr", pbl_mrf_Ribcr, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mrf_const_b", pbl_mrf_const_b, lev, max_level);
        query_one_or_per_level(pp, "pbl_mrf_sf", pbl_mrf_sf, lev, max_level);
        query_one_or_per_level(
          pp, "mrf_moistvars", mrf_moistvars, lev, max_level);
        query_one_or_per_level(
          pp, "enable_mrf_countergradient", enable_mrf_countergradient, lev, max_level);
        query_one_or_per_level(
          pp, "enable_mrf_cloud_adjustment", enable_mrf_cloud_adjustment, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mrf_highres_bounds", pbl_mrf_highres_bounds, lev, max_level);
        query_one_or_per_level(
          pp, "enable_mrf_unbounded_vpert", enable_mrf_unbounded_vpert, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_mrf_use_zero_ri_extent", pbl_mrf_use_zero_ri_extent, lev, max_level);
        // Scale-aware PBL-LES blending parameters.
        query_one_or_per_level(
          pp, "pbl_blend_length", pbl_blend_length, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_blend_cs", pbl_blend_cs, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_blend_c_max", pbl_blend_c_max, lev, max_level);
        query_one_or_per_level(
          pp, "pbl_blend_use_smag", pbl_blend_use_smag, lev, max_level);
      } else if (pbl_type == PBLType::SHOC) {
#ifndef ERF_USE_SHOC
          amrex::Abort("You set use_shoc to true but didn't build with SHOC; you must rebuild the executable");
      }
#endif
      if (uses_eamxx_shoc()) {
#ifndef ERF_USE_EAMXX_SHOC
          amrex::Abort("PBLType::EAMXX_SHOC requested, but ERF was not built with ERF_ENABLE_EAMXX_SHOC=ON");
#endif
      }

      if (pbl_requires_surface_layer()) {
          std::string zlo_bc = "none";
          amrex::ParmParse pp_bc("zlo");
          pp_bc.get("type",zlo_bc);
          if (amrex::toLower(zlo_bc) != "surface_layer") {
               amrex::Abort("You must use the surface_layer BC at zlo with the selected PBL.");
          }
      }
    }

    // Flags for QKE/TKE equation
    if (pbl_type == PBLType::MYJ || pbl_type == PBLType::MYNN25 || pbl_type == PBLType::MYNNEDMF) {
        // Add sources/sinks to QKE/TKE? (MYJ does this inline)
        if (pbl_type == PBLType::MYNN25 || pbl_type == PBLType::MYNNEDMF) {
            use_pbl_tke = true;
        }
        // Advect QKE/TKE?
        query_one_or_per_level(pp, "advect_tke"    , advect_tke    , lev, max_level);
        // Apply numerical diffusion to QKE/TKE?
        query_one_or_per_level(pp, "diffuse_tke_3D", diffuse_tke_3D, lev, max_level);
    }

    // There is a default value of 1e-6 but the user can override the value here,
    //   and in addition can add perturbational values.  and in addition can add perturbational values.
    pp.queryAdd("tke_min",tke_min);

    // LES constants...
    query_one_or_per_level(pp, "Cs", Cs, lev, max_level);

    query_one_or_per_level(pp, "Pr_t", Pr_t, lev, max_level);
    query_one_or_per_level(pp, "Sc_t", Sc_t, lev, max_level);

    // Compute relevant forms of diffusion parameters
    Pr_t_inv = one / Pr_t;
    Sc_t_inv = one / Sc_t;

    if (les_type == LESType::Deardorff) {
        query_one_or_per_level(pp, "Ck", Ck, lev, max_level);
        query_one_or_per_level(pp, "Ce", Ce, lev, max_level);
        query_one_or_per_level(pp, "Ce_wall", Ce_wall, lev, max_level);
    }

    // To quantify atmospheric stability for subgrid modeling
    query_one_or_per_level(pp, "thermal_stratification", strat_type, lev, max_level);
    if (strat_type == StratType::theta) {
        amrex::Print() << "Thermal stratification based on gradient of potential temperature" << std::endl;
    } else if (strat_type == StratType::thetav) {
        amrex::Print() << "Thermal stratification based on gradient of virtual potential temperature" << std::endl;
    } else if (strat_type == StratType::thetal) {
        amrex::Print() << "Thermal stratification based on gradient of linearized liquid-water potential temperature" << std::endl;
    }

    // k-eqn constants
    query_one_or_per_level(pp, "Cmu0", Cmu0, lev, max_level);
    query_one_or_per_level(pp, "Cb", Cb, lev, max_level);
    query_one_or_per_level(pp, "Rt_crit", Rt_crit, lev, max_level);
    query_one_or_per_level(pp, "Rt_min", Rt_min, lev, max_level);
    query_one_or_per_level(pp, "max_geom_lscale", l_g_max, lev, max_level);
    query_one_or_per_level(pp, "dirichlet_k", dirichlet_k, lev, max_level);

    // Common inputs (LES or RANS)
    if (!query_one_or_per_level(pp, "sigma_k", sigma_k, lev, max_level) && rans_type == RANSType::kEqn) {
        amrex::Print() << "Overriding default sigma_k for k-eqn RANS" << std::endl;
        sigma_k = one;
    };
    query_one_or_per_level(pp, "theta_ref", theta_ref, lev, max_level);

    query_one_or_per_level(pp, "mix_isotropic", mix_isotropic, lev, max_level);
    query_one_or_per_level(pp, "use_Ri_correction", use_Ri_correction, lev, max_level);
    query_one_or_per_level(pp, "Ri_crit", Ri_crit, lev, max_level);

    // Set common flags
    use_kturb =
      ((les_type != LESType::None) || (rans_type != RANSType::None) ||
       (pbl_type != PBLType::None));
    use_keqn =
      ((les_type == LESType::Deardorff) || (rans_type == RANSType::kEqn));
    use_tke =
      ((les_type == LESType::Deardorff) || (rans_type == RANSType::kEqn) ||
       (pbl_type == PBLType::MYJ)       || (pbl_type == PBLType::MYNN25) ||
       (pbl_type == PBLType::MYNNEDMF)  || uses_shoc_family());

    if (use_tke) {
        query_one_or_per_level(pp, "init_tke_from_ustar", init_tke_from_ustar, lev, max_level);
    }

    // Validate inputs
    if (les_type == LESType::Smagorinsky) {
      if (Cs == 0) {
        amrex::Error("Need to specify Cs for Smagorsinky LES");
      }
      if (smag2d && mix_isotropic) {
        amrex::Print() << "Turning off mix_isotropic for 2-D Smagorinsky" << std::endl;
        mix_isotropic = false;
      }
    }
  }

  /**
   * @brief Validate turbulence options against physical boundary conditions.
   * @param phys_bc_type Physical boundary-condition types.
   */
  void check_params (amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>& phys_bc_type)
  {
      // BC compatibility
      if (pbl_requires_surface_layer() &&
          phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::surface_layer ) {
          amrex::Abort("The selected PBL model requires MOST at lower boundary");
      }
      if ( (les_type == LESType::Deardorff) && (Ce_wall > 0) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::surface_layer) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::slip_wall) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::no_slip_wall) )
      {
          amrex::Warning("Deardorff LES assumes wall at zlo when applying Ce_wall");
      }
  }

  /**
   * @brief Print turbulence settings for one AMR level.
   * @param lev AMR level index.
   */
  void display (int lev)
  {
    amrex::Print() << "Turbulence Settings at level " << lev << std::endl;

    if (
      les_type == LESType::None && rans_type == RANSType::None &&
      pbl_type == PBLType::None) {
      amrex::Print() << "    Using DNS model at level " << lev << std::endl;
    } else if (les_type == LESType::Smagorinsky) {
      if (smag2d) {
        amrex::Print() << "    Using 2D Smagorinsky LES model at level " << lev << std::endl;
      } else {
        amrex::Print() << "    Using Smagorinsky LES model at level " << lev << std::endl;
      }
      if (use_Ri_correction) {
        amrex::Print() << "    Smagorinsky uses Richardson number correction with Ri_crit = "
                       << Ri_crit << std::endl;
      }
    } else if (les_type == LESType::Deardorff) {
      amrex::Print() << "    Using Deardorff LES model at level " << lev << std::endl;
    } else if (rans_type == RANSType::kEqn) {
      amrex::Print()
        << "    Using Axell & Liungman one-equation RANS k model at level " << lev
        << std::endl;
    } else if (pbl_type == PBLType::MYJ) {
      amrex::Print() << "    Using MYJ PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MYNN25) {
      amrex::Print() << "    Using MYNN2.5 PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MYNNEDMF) {
      amrex::Print() << "    Using MYNNEDMF PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::YSU) {
      amrex::Print() << "    Using YSU PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::YSUNew) {
      amrex::Print() << "    Using YSU PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MRF) {
      amrex::Print() << "    Using MRF PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::EAMXX_SHOC) {
      amrex::Print() << "    Using EAMxx SHOC PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::NATIVE_SHOC) {
      amrex::Print() << "    Using native SHOC PBL model at level " << lev << std::endl;
    } else {
      amrex::Error("Unknown turbulence model");
    }

    if (les_type != LESType::None) {
      if (les_type == LESType::Smagorinsky) {
        amrex::Print() << "    Cs                      : " << Cs << std::endl;
      }
      if (les_type == LESType::Deardorff) {
        amrex::Print() << "    Ce                      : " << Ce      << std::endl;
        amrex::Print() << "    Ce at wall              : " << Ce_wall << std::endl;
        amrex::Print() << "    Ck                      : " << Ck      << std::endl;
        amrex::Print() << "    sigma_k                 : " << sigma_k << std::endl;

        // Sullivan et al 1994, Eqn 14
        amrex::Real Cs_equiv = std::sqrt(Ck * std::sqrt(Ck / Ce));
        amrex::Print() << "    equivalent Cs           : " << Cs_equiv
                       << std::endl;
      }
      amrex::Print() << "    isotropic mixing        : " << mix_isotropic
                     << std::endl;
    }

    if (rans_type != RANSType::None) {
      if (rans_type == RANSType::kEqn) {
        amrex::Print() << "Cmu0                        : " << Cmu0 << std::endl;
        amrex::Print() << "sigma_k                     : " << sigma_k << std::endl;
        amrex::Print() << "Cb                          : " << Cb << std::endl;
        amrex::Print() << "Rt_crit                     : " << Rt_crit << std::endl;
        amrex::Print() << "Rt_min                      : " << Rt_min << std::endl;
        amrex::Print() << "max_geom_lscale             : " << l_g_max << std::endl;
      }
    }

    if ((les_type  == LESType::Deardorff) ||
        (rans_type == RANSType::kEqn)) {
      if (theta_ref > 0) {
        amrex::Print() << "    reference theta         : " << theta_ref << std::endl;
      } else {
        amrex::Print() << "    reference theta         : n/a" << std::endl;
      }
    }

    if ((les_type != LESType::None) || (rans_type != RANSType::None)) {
      amrex::Print() << "    Pr_t                    : " << Pr_t << std::endl;
      amrex::Print() << "    Sc_t                    : " << Sc_t << std::endl;
    }

    if (pbl_type == PBLType::MYNN25 || pbl_type == PBLType::MYNNEDMF) {
      amrex::Print() << "    pbl_mynn_A1 : " << pbl_mynn.A1 << std::endl;
      amrex::Print() << "    pbl_mynn_A2 : " << pbl_mynn.A2 << std::endl;
      amrex::Print() << "    pbl_mynn_B1 : " << pbl_mynn.B1 << std::endl;
      amrex::Print() << "    pbl_mynn_B2 : " << pbl_mynn.B2 << std::endl;
      amrex::Print() << "    pbl_mynn_C1 : " << pbl_mynn.C1 << std::endl;
      amrex::Print() << "    pbl_mynn_C2 : " << pbl_mynn.C2 << std::endl;
      amrex::Print() << "    pbl_mynn_C3 : " << pbl_mynn.C3 << std::endl;
      amrex::Print() << "    pbl_mynn_C4 : " << pbl_mynn.C4 << std::endl;
      amrex::Print() << "    pbl_mynn_C5 : " << pbl_mynn.C5 << std::endl;
    } else if (pbl_type == PBLType::YSU || pbl_type == PBLType::YSUNew) {
      amrex::Print() << "    pbl_ysu_coriolis_freq              : "
                     << pbl_ysu_coriolis_freq << std::endl;
      amrex::Print() << "    pbl_ysu_use_consistent_coriolis    : "
                     << pbl_ysu_use_consistent_coriolis << std::endl;
      amrex::Print() << "    pbl_ysu_force_over_water           : "
                     << pbl_ysu_force_over_water << std::endl;
      amrex::Print() << "    pbl_ysu_land_Ribcr                 : "
                     << pbl_ysu_land_Ribcr << std::endl;
      amrex::Print() << "    pbl_ysu_unst_Ribcr                 : "
                     << pbl_ysu_unst_Ribcr << std::endl;
      amrex::Print() << "    enable_ysu_liquid_theta            : "
                     << enable_ysu_liquid_theta << std::endl;
      amrex::Print() << "    enable_ysu_countergradient         : "
                     << enable_ysu_countergradient << std::endl;
      amrex::Print() << "    enable_ysu_terrain_pblh_floor      : "
                     << enable_ysu_terrain_pblh_floor << std::endl;
      amrex::Print() << "    enable_ysu_sat_limiter             : "
                     << enable_ysu_sat_limiter << std::endl;
      amrex::Print() << "    enable_ysu_topdown                 : "
                     << enable_ysu_topdown << std::endl;
      amrex::Print() << "    enable_ysu_entrainment             : "
                     << enable_ysu_entrainment << std::endl;
      amrex::Print() << "    enable_ysu_cloud_pblh              : "
                     << enable_ysu_cloud_pblh << std::endl;
      amrex::Print() << "    ysu_qcloud_threshold               : "
                     << ysu_qcloud_threshold << std::endl;
      amrex::Print() << "    ysu_moistvars                      : "
                     << ysu_moistvars << std::endl;
      amrex::Print() << "    pbl_ysunew_highres_bounds          : "
                     << pbl_ysunew_highres_bounds << std::endl;
      amrex::Print() << "    enable_mrf_unbounded_vpert         : "
                     << enable_mrf_unbounded_vpert << std::endl;
    } else if (pbl_type == PBLType::MRF) {
      amrex::Print() << "    pbl_mrf_coriolis_freq : " << pbl_mrf_coriolis_freq
                     << std::endl;
      amrex::Print() << "    pbl_mrf_Ribcr         : " << pbl_mrf_Ribcr
                     << std::endl;
      amrex::Print() << "    pbl_mrf_const_b       : " << pbl_mrf_const_b
                     << std::endl;
      amrex::Print() << "    pbl_mrf_sf            : " << pbl_mrf_sf
                     << std::endl;
      amrex::Print() << "    mrf_moistvars         : " << mrf_moistvars
                     << std::endl;
      amrex::Print() << "    enable_mrf_countergradient : " << enable_mrf_countergradient
                     << std::endl;
      amrex::Print() << "    enable_mrf_cloud_adjustment : " << enable_mrf_cloud_adjustment
                     << std::endl;
      amrex::Print() << "    pbl_mrf_highres_bounds     : " << pbl_mrf_highres_bounds
                     << std::endl;
      amrex::Print() << "    enable_mrf_unbounded_vpert : " << enable_mrf_unbounded_vpert
                     << std::endl;
      amrex::Print() << "    pbl_mrf_use_zero_ri_extent : " << pbl_mrf_use_zero_ri_extent
                     << std::endl;
      amrex::Print() << "    pbl_blend_length   : " << pbl_blend_length
                     << std::endl;
      amrex::Print() << "    pbl_blend_cs       : " << pbl_blend_cs
                     << std::endl;
      amrex::Print() << "    pbl_blend_c_max    : " << pbl_blend_c_max
                     << std::endl;
      amrex::Print() << "    pbl_blend_use_smag : " << pbl_blend_use_smag
                     << std::endl;
    }
  }

  // LES model
  LESType les_type = LESType::None;

  // Turbulent Prandtl number
  amrex::Real Pr_t = one / three;
  amrex::Real Pr_t_inv = three;

  // Turbulent Schmidt number
  amrex::Real Sc_t = one;
  amrex::Real Sc_t_inv = one;

  // Smagorinsky
  amrex::Real Cs = zero;
  bool smag2d = false;

  // Deardorff
  amrex::Real Ce = amrex::Real(0.93);
  amrex::Real Ce_wall = zero; // if > 0, then set Ce to this at k=0
  amrex::Real Ck = amrex::Real(0.1);

  // k-eqn RANS coefficients (Axell & Liungman 2001)
  amrex::Real Cmu0 = amrex::Real(0.5562);
  amrex::Real Cb = amrex::Real(0.35);
  amrex::Real Rt_crit = -one;
  amrex::Real Rt_min = -three;
  amrex::Real l_g_max = amrex::Real(30.0); // ~ kappa * (amrex::Real(0.1) * zi)

  // Deardorff or k-eqn RANS
  // - diffusivity of tke is 1/sigma_k times the eddy viscosity
  amrex::Real sigma_k = amrex::Real(0.5);
  // - reference potential temperature used to quantify the stratification in
  //   a stable region
  amrex::Real theta_ref = zero;
  // - how to quantify stratification effects
  StratType strat_type = StratType::theta;

  // Anisotropic length scales
  bool mix_isotropic = true;

  bool use_Ri_correction = true;
  amrex::Real Ri_crit = fourth;

  // RANS type
  RANSType rans_type = RANSType::None;

  bool dirichlet_k = false;

  // PBL model
  PBLType pbl_type = PBLType::None;

  /**
   * @brief Query whether this level uses the EAMxx SHOC PBL scheme.
   * @return True if PBLType::EAMXX_SHOC is selected.
   */
  bool uses_eamxx_shoc () const noexcept
  {
      return pbl_type == PBLType::EAMXX_SHOC;
  }

  /**
   * @brief Query whether this level uses the native SHOC PBL scheme.
   * @return True if PBLType::NATIVE_SHOC is selected.
   */
  bool uses_native_shoc () const noexcept
  {
      return pbl_type == PBLType::NATIVE_SHOC;
  }

  /**
   * @brief Query whether this level uses any SHOC-family PBL scheme.
   * @return True if native or EAMxx SHOC is selected.
   */
  bool uses_shoc_family () const noexcept
  {
      return uses_eamxx_shoc() || uses_native_shoc();
  }

  /**
   * @brief Query whether the selected PBL scheme requires surface-layer boundary conditions.
   * @return True if the PBL scheme requires a lower surface-layer boundary.
   */
  bool pbl_requires_surface_layer () const noexcept
  {
      return (pbl_type == PBLType::MYNN25)   ||
             (pbl_type == PBLType::MYNNEDMF) ||
             (pbl_type == PBLType::YSU)      ||
             (pbl_type == PBLType::YSUNew)      ||
             (pbl_type == PBLType::MRF)       ||
             uses_shoc_family();
  }

  /**
   * @brief Query whether the selected PBL scheme suppresses microphysics condensation.
   * @return True when the PBL scheme owns condensation handling.
   */
  bool pbl_suppresses_microphysics_condensation () const noexcept
  {
      return uses_shoc_family();
  }

  MYNNLevel25 pbl_mynn;
  MYNNLevel2 pbl_mynn_level2; // limiting in the decaying turbulence regime

  // Common Flags
  bool use_kturb = false; // Any turbulence modeling?
  bool use_keqn =
    false; // Any microscale turbulence modeling (LES, RANS) with TKE closure?
           // Then need to populate SmnSmn_lev for production term.
  bool use_pbl_tke =
    false; // Any mesoscale turbulence modeling (PBL) with TKE closure?
  bool use_tke = false; // Any TKE closure (meso or microscale)?

  // Initialize TKE/QKE with linear profiles whose surface values are a
  // function of the friction velocity calculated by the surface layer scheme
  // (e.g., as done in MYNN-EDMF)
  bool init_tke_from_ustar = false;

  // This is the value of tke_min in WRF 4.5
  amrex::Real tke_min = amrex::Real(1.e-6);

  // Model coefficients - YSU
  // TODO: Add parmparse for all of these above
  amrex::Real pbl_ysu_coriolis_freq =
    amrex::Real(1.0e-4); // 1e-4 is hardcoded in WRF, we let the user specify or take the
            // value from ERF coriolis forcing
  bool pbl_ysu_use_consistent_coriolis =
    false; // ignore input pbl_ysu_coriolis_freq, take value from ERF coriolis
           // forcing instead
  bool pbl_ysu_force_over_water =
    false; // Force YSU to act as if it is over water regardless of other inputs
           // (for testing)
  amrex::Real pbl_ysu_land_Ribcr =
    fourth; // Critical Bulk Richardson number of Land for stable conditions
  amrex::Real pbl_ysu_unst_Ribcr =
    zero; // Critical Bulk Richardson number for unstable conditions
  bool enable_ysu_liquid_theta =
    true; // if true, use liquid-water virtual potential temperature in YSU stability
  bool enable_ysu_countergradient =
    true; // if true, apply YSU-style HGAMT/HGAMQ/HGAMU/HGAMV countergradient correction
  bool enable_ysu_terrain_pblh_floor =
    true; // if true, apply terrain-following floor on PBL height
  bool enable_ysu_sat_limiter =
    false; // if true, apply saturation limiter to moisture countergradient
  bool enable_ysu_topdown =
    true; // if true, enable top-down mixing for cloud-topped boundary layers
  bool enable_ysu_entrainment =
    true; // if true, enable entrainment layer parameterization
  bool enable_ysu_cloud_pblh =
    true; // if true, enable cloud-based PBL height detection
  amrex::Real ysu_qcloud_threshold =
    amrex::Real(1.0e-4); // Cloud liquid water threshold for YSUNew (kg/kg)
  bool ysu_moistvars = false; // if true, adds turbulence to moisture
  bool pbl_ysunew_highres_bounds = false; // if true, apply high-resolution grid-dependent diffusivity bounds
  // Model coefficients - MRF
  amrex::Real pbl_mrf_coriolis_freq = amrex::Real(1.0e-4);
  amrex::Real pbl_mrf_Ribcr   = amrex::Real(0.5); // Critical Bulk Richardson number for MRF PBL model
  amrex::Real pbl_mrf_const_b = amrex::Real(7.8); // Constant b in MRF PBL model, used to compute the PBL height
  amrex::Real pbl_mrf_sf      = amrex::Real(0.1); // MRF surface flux, used to compute the PBL height
  // Scale-aware PBL-LES blending parameters.
  // Boutle et al. (2014): https://doi.org/10.1175/MWR-D-13-00229.1
  // Set pbl_blend_length = 0 to disable (default). Recovers original behaviour.
  amrex::Real pbl_blend_length   = amrex::Real(0.0);   ///< Boutle blending length L [m]. 0 = off.
  amrex::Real pbl_blend_cs       = amrex::Real(0.17);  ///< Smagorinsky coeff for K_h ceiling.
  amrex::Real pbl_blend_c_max    = amrex::Real(0.1);   ///< Power-law ceiling coeff [m^(2/3)/s].
  bool        pbl_blend_use_smag = true;               ///< Use Smagorinsky ceiling (else power-law).
  bool mrf_moistvars = false; // if true, adds turbulence to moisture
  bool enable_mrf_countergradient = false; // if true, apply WRF-style HGAMT/HGAMQ correction to PBL height
  bool enable_mrf_cloud_adjustment = false; // if true, apply cloud-aware stability function adjustments in MRF scheme
  bool pbl_mrf_highres_bounds = false; // if true, apply high-resolution grid-dependent diffusivity bounds instead of global forecast bounds
  bool enable_mrf_unbounded_vpert = false; // if true, does not limit VPERT to GAMCRT
  bool pbl_mrf_use_zero_ri_extent = false; // if true, use Ri=0 (pbli_zero_arr) for K-profile extent; default false uses Ri=0.5 corrector (pbli_arr) for WRF compatibility
  // MYNN2.5 PBL model
  // TKE/QKE stuff
  bool advect_tke = true; // if MYNN2.5 PBL is used default is turb transport in
                          // Z-direction only
  bool diffuse_tke_3D = true; // if numerical diffusion is turned on
};

/**
 * @var TurbChoice::les_type
 * @brief Selected LES closure.
 * @var TurbChoice::Pr_t
 * @brief Turbulent Prandtl number.
 * @var TurbChoice::Pr_t_inv
 * @brief Inverse turbulent Prandtl number.
 * @var TurbChoice::Sc_t
 * @brief Turbulent Schmidt number.
 * @var TurbChoice::Sc_t_inv
 * @brief Inverse turbulent Schmidt number.
 * @var TurbChoice::Cs
 * @brief Smagorinsky model coefficient.
 * @var TurbChoice::smag2d
 * @brief Whether the 2-D Smagorinsky formulation is used.
 * @var TurbChoice::Ce
 * @brief Deardorff dissipation coefficient.
 * @var TurbChoice::Ce_wall
 * @brief Wall value for the Deardorff dissipation coefficient.
 * @var TurbChoice::Ck
 * @brief Deardorff eddy-viscosity coefficient.
 * @var TurbChoice::Cmu0
 * @brief One-equation RANS Cmu0 coefficient.
 * @var TurbChoice::Cb
 * @brief One-equation RANS buoyancy coefficient.
 * @var TurbChoice::Rt_crit
 * @brief Critical turbulent Reynolds number.
 * @var TurbChoice::Rt_min
 * @brief Minimum turbulent Reynolds number.
 * @var TurbChoice::l_g_max
 * @brief Maximum geometric length scale.
 * @var TurbChoice::sigma_k
 * @brief TKE diffusivity coefficient denominator.
 * @var TurbChoice::theta_ref
 * @brief Reference potential temperature for stable stratification.
 * @var TurbChoice::strat_type
 * @brief Thermodynamic variable used for stability stratification.
 * @var TurbChoice::mix_isotropic
 * @brief Whether subgrid mixing uses isotropic length scales.
 * @var TurbChoice::use_Ri_correction
 * @brief Whether Richardson-number correction is applied.
 * @var TurbChoice::Ri_crit
 * @brief Critical Richardson number for stability correction.
 * @var TurbChoice::rans_type
 * @brief Selected RANS closure.
 * @var TurbChoice::dirichlet_k
 * @brief Whether TKE uses Dirichlet boundary treatment.
 * @var TurbChoice::pbl_type
 * @brief Selected PBL closure.
 * @var TurbChoice::pbl_mynn
 * @brief MYNN level-2.5 closure coefficients.
 * @var TurbChoice::pbl_mynn_level2
 * @brief MYNN level-2 closure coefficients for limiting.
 * @var TurbChoice::use_kturb
 * @brief Whether any turbulence model is active.
 * @var TurbChoice::use_keqn
 * @brief Whether a microscale TKE closure is active.
 * @var TurbChoice::use_pbl_tke
 * @brief Whether a mesoscale PBL TKE closure is active.
 * @var TurbChoice::use_tke
 * @brief Whether any TKE or QKE closure is active.
 * @var TurbChoice::init_tke_from_ustar
 * @brief Whether initial TKE/QKE profiles are based on surface friction velocity.
 * @var TurbChoice::tke_min
 * @brief Minimum TKE/QKE value.
 * @var TurbChoice::pbl_ysu_coriolis_freq
 * @brief Coriolis frequency used by YSU-family PBL schemes.
 * @var TurbChoice::pbl_ysu_use_consistent_coriolis
 * @brief Whether YSU uses the ERF Coriolis frequency.
 * @var TurbChoice::pbl_ysu_force_over_water
 * @brief Whether YSU is forced to use over-water behavior.
 * @var TurbChoice::pbl_ysu_land_Ribcr
 * @brief Critical bulk Richardson number over land for stable YSU conditions.
 * @var TurbChoice::pbl_ysu_unst_Ribcr
 * @brief Critical bulk Richardson number for unstable YSU conditions.
 * @var TurbChoice::enable_ysu_liquid_theta
 * @brief Whether YSU uses liquid-water virtual potential temperature for stability.
 * @var TurbChoice::enable_ysu_countergradient
 * @brief Whether YSU countergradient corrections are enabled.
 * @var TurbChoice::enable_ysu_terrain_pblh_floor
 * @brief Whether YSU applies a terrain-following PBL-height floor.
 * @var TurbChoice::enable_ysu_sat_limiter
 * @brief Whether YSU applies a saturation limiter to moisture countergradient terms.
 * @var TurbChoice::enable_ysu_topdown
 * @brief Whether YSU top-down mixing is enabled.
 * @var TurbChoice::enable_ysu_entrainment
 * @brief Whether YSU entrainment-layer parameterization is enabled.
 * @var TurbChoice::enable_ysu_cloud_pblh
 * @brief Whether YSU cloud-based PBL-height detection is enabled.
 * @var TurbChoice::ysu_qcloud_threshold
 * @brief Cloud liquid water threshold for YSUNew [kg/kg].
 * @var TurbChoice::ysu_moistvars
 * @brief Whether YSU applies turbulence to moisture variables.
 * @var TurbChoice::pbl_ysunew_highres_bounds
 * @brief Whether YSUNew applies high-resolution grid-dependent diffusivity bounds.
 * @var TurbChoice::pbl_mrf_coriolis_freq
 * @brief Coriolis frequency used by the MRF PBL scheme.
 * @var TurbChoice::pbl_mrf_Ribcr
 * @brief Critical bulk Richardson number for the MRF PBL scheme.
 * @var TurbChoice::pbl_mrf_const_b
 * @brief MRF constant used to compute PBL height.
 * @var TurbChoice::pbl_mrf_sf
 * @brief MRF surface flux value used to compute PBL height.
 * @var TurbChoice::mrf_moistvars
 * @brief Whether MRF applies turbulence to moisture variables.
 * @var TurbChoice::enable_mrf_countergradient
 * @brief Whether MRF countergradient corrections are enabled.
 * @var TurbChoice::enable_mrf_cloud_adjustment
 * @brief Whether MRF cloud-aware stability adjustments are enabled.
 * @var TurbChoice::pbl_mrf_highres_bounds
 * @brief Whether MRF applies high-resolution grid-dependent diffusivity bounds.
 * @var TurbChoice::enable_mrf_unbounded_vpert
 * @brief Whether MRF leaves VPERT unlimited by GAMCRT.
 * @var TurbChoice::pbl_mrf_use_zero_ri_extent
 * @brief Whether MRF uses the Ri=0 K-profile extent.
 * @var TurbChoice::advect_tke
 * @brief Whether TKE/QKE is advected.
 * @var TurbChoice::diffuse_tke_3D
 * @brief Whether three-dimensional numerical diffusion is applied to TKE/QKE.
 */
#endif
