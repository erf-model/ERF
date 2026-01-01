#ifndef ERF_TURB_STRUCT_H_
#define ERF_TURB_STRUCT_H_

#include <ERF_MYNNStruct.H>

AMREX_ENUM(LESType, None, Smagorinsky, Smagorinsky2D, Deardorff);

AMREX_ENUM(RANSType, None, kEqn);

AMREX_ENUM(PBLType, None, MYJ, MYNN25, MYNNEDMF, YSU, MRF, SHOC);

AMREX_ENUM(StratType, theta, thetav, thetal);

template <typename T>
void
query_one_or_per_level (
  const amrex::ParmParse& pp,
  const char* query_string,
  T& query_var,
  const int lev,
  const int maxlev)
{
  int count = pp.countval(query_string);
  if (count == 0) {
    return; // nothing to do
  } else if (count == 1) {
    pp.query(query_string, query_var);
  } else if (count == maxlev + 1) {
    pp.query(query_string, query_var, lev);
  } else {
    amrex::Error(
      "For parmparse variable " + pp.prefixedName(query_string) +
      ": if specified, specify once total or once for each level");
  }
}

/**
 * Container holding quantities related to turbulence parametrizations
 */
struct TurbChoice
{
public:
  void init_params (int lev, int max_level, std::string pp_prefix)
  {
    amrex::ParmParse pp(pp_prefix);

    // Which LES closure?
    std::string les_type_string = "None";
    query_one_or_per_level(pp, "les_type", les_type, lev, max_level);

    // Handle 2-D Smag
    if (les_type == LESType::Smagorinsky2D) {
        les_type = LESType::Smagorinsky;
        smag2d = true;
    }

    // Which RANS closure?
    std::string rans_type_string = "None";
    query_one_or_per_level(pp, "rans_type", rans_type, lev, max_level);

    if ((rans_type != RANSType::None) && (les_type != LESType::None)) {
        amrex::Error("Hybrid RANS-LES not implemented");
    }

    // Which PBL Closure
    static std::string pbl_type_string = "None";
    query_one_or_per_level(pp, "pbl_type", pbl_type, lev, max_level);

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
      } else if (pbl_type == PBLType::YSU) {
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
      } else if (pbl_type == PBLType::SHOC) {
#ifndef ERF_USE_SHOC
          amrex::Abort("You set use_shoc to true but didn't build with SHOC; you must rebuild the executable");
#endif
          std::string zlo_bc = "none";
          amrex::ParmParse pp_bc("zlo");
          pp_bc.get("type",zlo_bc);
          if (amrex::toLower(zlo_bc) != "surface_layer") {
               amrex::Abort("You must use the surface_layer BC at zlo with SHOC.");
          }
      }
    }

    // Flags for QKE/TKE euation
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

    // LES constants...
    query_one_or_per_level(pp, "Cs", Cs, lev, max_level);

    query_one_or_per_level(pp, "Pr_t", Pr_t, lev, max_level);
    query_one_or_per_level(pp, "Sc_t", Sc_t, lev, max_level);

    // Compute relevant forms of diffusion parameters
    Pr_t_inv = amrex::Real(1.0) / Pr_t;
    Sc_t_inv = amrex::Real(1.0) / Sc_t;

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

    // Common inputs (LES or RANS)
    query_one_or_per_level(pp, "sigma_k", sigma_k, lev, max_level);
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
       (pbl_type == PBLType::MYNNEDMF)  || (pbl_type == PBLType::SHOC));

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

  void check_params (amrex::GpuArray<ERF_BC, AMREX_SPACEDIM*2>& phys_bc_type)
  {
      // BC compatibility
      if ( ( (pbl_type == PBLType::MYNN25)   ||
             (pbl_type == PBLType::MYNNEDMF) ||
             (pbl_type == PBLType::YSU) ||
             (pbl_type == PBLType::MRF)
                 ) &&
          phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::surface_layer ) {
          amrex::Abort("MYNN2.5/MYNNEDMF/YSU/MRF PBL Model requires MOST at lower boundary");
      }
      if ( (les_type == LESType::Deardorff) && (Ce_wall > 0) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::surface_layer) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::slip_wall) &&
           (phys_bc_type[amrex::Orientation(amrex::Direction::z,amrex::Orientation::low)] != ERF_BC::no_slip_wall) )
      {
          amrex::Warning("Deardorff LES assumes wall at zlo when applying Ce_wall");
      }
  }

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
        << "    Using Axell & Liungman one-equation RANS k model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MYJ) {
      amrex::Print() << "    Using MYJ PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MYNN25) {
      amrex::Print() << "    Using MYNN2.5 PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MYNNEDMF) {
      amrex::Print() << "    Using MYNNEDMF PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::YSU) {
      amrex::Print() << "    Using YSU PBL model at level " << lev << std::endl;
    } else if (pbl_type == PBLType::MRF) {
      amrex::Print() << "    Using MRF PBL model at level " << lev << std::endl;
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
    } else if (pbl_type == PBLType::YSU) {
      amrex::Print() << "    pbl_ysu_coriolis_freq           : "
                     << pbl_ysu_coriolis_freq << std::endl;
      amrex::Print() << "    pbl_ysu_use_consistent_coriolis : "
                     << pbl_ysu_use_consistent_coriolis << std::endl;
      amrex::Print() << "    pbl_ysu_force_over_water        : "
                     << pbl_ysu_force_over_water << std::endl;
      amrex::Print() << "    pbl_ysu_land_Ribcr              : "
                     << pbl_ysu_land_Ribcr << std::endl;
      amrex::Print() << "    pbl_ysu_unst_Ribcr              : "
                     << pbl_ysu_unst_Ribcr << std::endl;
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
    }
  }

  // LES model
  LESType les_type;

  // Turbulent Prandtl number
  amrex::Real Pr_t = amrex::Real(1.0) / amrex::Real(3.0);
  amrex::Real Pr_t_inv = amrex::Real(3.0);

  // Turbulent Schmidt number
  amrex::Real Sc_t = 1.0;
  amrex::Real Sc_t_inv = 1.0;

  // Smagorinsky
  amrex::Real Cs = 0.0;
  bool smag2d = false;

  // Deardorff
  amrex::Real Ce = 0.93;
  amrex::Real Ce_wall = 0.0; // if > 0, then set Ce to this at k=0
  amrex::Real Ck = 0.1;

  // k-eqn RANS coefficients (Axell & Liungman 2001)
  amrex::Real Cmu0 = 0.5562;
  amrex::Real Cb = 0.35;
  amrex::Real Rt_crit = -1.0;
  amrex::Real Rt_min = -3.0;
  amrex::Real l_g_max = 1e99;

  // Deardorff or k-eqn RANS
  // - diffusivity of tke is 1/sigma_k times the eddy viscosity
  amrex::Real sigma_k = 0.5;
  // - reference potential temperature used to quantify the stratification in
  //   a stable region
  amrex::Real theta_ref = 0.0;
  // - how to quantify stratification effects
  StratType strat_type = StratType::theta;

  // Anistropic length scales
  bool mix_isotropic = true;

  bool use_Ri_correction = true;
  amrex::Real Ri_crit = 0.25;

  // RANS type
  RANSType rans_type;

  // PBL model
  PBLType pbl_type;

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

  // Model coefficients - YSU
  // TODO: Add parmparse for all of these above
  amrex::Real pbl_ysu_coriolis_freq =
    1.0e-4; // 1e-4 is hardcoded in WRF, we let the user specify or tske the
            // value from ERF coriolis forcing
  bool pbl_ysu_use_consistent_coriolis =
    false; // ignore input pbl_ysu_coriolis_freq, take value from ERF coriolis
           // forcing instead
  bool pbl_ysu_force_over_water =
    false; // Force YSU to act as if it is over water regardless of other inputs
           // (for testing)
  amrex::Real pbl_ysu_land_Ribcr =
    0.25; // Critical Bulk Richardson number of Land for stable conditions
  amrex::Real pbl_ysu_unst_Ribcr =
    0.0; // Critical Bulk Richardson number for unstable conditions
  // Model coefficients - MRF
  amrex::Real pbl_mrf_coriolis_freq = 1.0e-4;
  amrex::Real pbl_mrf_Ribcr =
    0.5; // Critical Bulk Richardson number for MRF PBL model
  amrex::Real pbl_mrf_const_b =
    7.8; // Constant b in MRF PBL model, used to compute the PBL height
  amrex::Real pbl_mrf_sf =
    0.1; // MRF surface flux, used to compute the PBL height
  bool mrf_moistvars = false; // if true, adds turbulence to moisture
  // MYNN2.5 PBL model
  // TKE/QKE stuff
  bool advect_tke = true; // if MYNN2.5 PBL is used default is turb transport in
                          // Z-direction only
  bool diffuse_tke_3D = true; // if numerical diffusion is turned on
};
#endif
