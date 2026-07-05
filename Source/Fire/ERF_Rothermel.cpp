#include <ERF_Rothermel.H>


using namespace amrex;

// Anderson FBFM13 fuel model database
// Reference: Anderson, H.E. (1982) "Aids to determining fuel models..."
// All values are standard published values for these fuel models

FuelModelParams get_anderson_fuel_params(int model_id)
{
    FuelModelParams fp;

    // All fuel loads in lb/ft² (standard Anderson 1982)
    // SAV in 1/ft (characteristic for 1-hour dead fuel; 10-hr = 109, 100-hr = 30 ft⁻¹ fixed)
    // Heat content in BTU/lb
    // Particle density in lb/ft³
    // Depth in ft
    // Reference: wildfire_levelset/src/fuel_database.H (Anderson 13 fuel models)

    switch (model_id) {
    case 1:  // Short grass (1 ft)
        fp.w_d1 = 0.034;   fp.w_d10 = 0.0;    fp.w_d100 = 0.0;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 3500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.0;    fp.Mx = 0.12;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 2:  // Timber (Grass and Understory)
        fp.w_d1 = 0.046;   fp.w_d10 = 0.023;  fp.w_d100 = 0.023;
        fp.w_lh = 0.023;   fp.w_lw = 0.0;
        fp.sigma_d1 = 3000.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.0;    fp.Mx = 0.15;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 3:  // Tall grass (2.5 ft)
        fp.w_d1 = 0.138;   fp.w_d10 = 0.0;    fp.w_d100 = 0.0;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 1500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 4:  // Chaparral (6 ft)
        fp.w_d1 = 0.230;   fp.w_d10 = 0.184;  fp.w_d100 = 0.092;
        fp.w_lh = 0.230;   fp.w_lw = 0.0;
        fp.sigma_d1 = 1739.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 6.0;    fp.Mx = 0.20;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 5:  // Brush (2 ft)
        fp.w_d1 = 0.046;   fp.w_d10 = 0.023;  fp.w_d100 = 0.0;
        fp.w_lh = 0.092;   fp.w_lw = 0.0;
        fp.sigma_d1 = 1739.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.0;    fp.Mx = 0.20;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 6:  // Dormant Brush, Hardwood Slash
        fp.w_d1 = 0.069;   fp.w_d10 = 0.115;  fp.w_d100 = 0.115;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 1750.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 7:  // Southern Rough
        fp.w_d1 = 0.052;   fp.w_d10 = 0.086;  fp.w_d100 = 0.069;
        fp.w_lh = 0.017;   fp.w_lw = 0.017;
        fp.sigma_d1 = 1550.0;  fp.sigma_lh = 1550.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.40;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 8:  // Closed Timber Litter
        fp.w_d1 = 0.069;   fp.w_d10 = 0.046;  fp.w_d100 = 0.115;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 2000.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 0.2;    fp.Mx = 0.30;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 9:  // Hardwood Litter
        fp.w_d1 = 0.134;   fp.w_d10 = 0.019;  fp.w_d100 = 0.007;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 2500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 0.2;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 10:  // Timber (Litter and Understory)
        fp.w_d1 = 0.138;   fp.w_d10 = 0.092;  fp.w_d100 = 0.230;
        fp.w_lh = 0.092;   fp.w_lw = 0.0;
        fp.sigma_d1 = 2000.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.0;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 11:  // Light Logging Slash
        fp.w_d1 = 0.069;   fp.w_d10 = 0.207;  fp.w_d100 = 0.253;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 1500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.0;    fp.Mx = 0.15;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 12:  // Medium Logging Slash
        fp.w_d1 = 0.184;   fp.w_d10 = 0.644;  fp.w_d100 = 0.759;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 1500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.3;    fp.Mx = 0.20;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 13:  // Heavy Logging Slash
        fp.w_d1 = 0.322;   fp.w_d10 = 1.058;  fp.w_d100 = 1.288;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 1500.0;  fp.sigma_lh = 1500.0;  fp.sigma_lw = 0.0;
        fp.delta = 3.0;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    default:
        // Default to model 1 if invalid ID
        amrex::Print() << "Warning: Unknown fuel model ID " << model_id
                       << ", using model 1" << std::endl;
        return get_anderson_fuel_params(1);
    }

    return fp;
}

RothermelComputed compute_rothermel_params(const FuelModelParams& fp,
                                           Real moisture_1hr,
                                           Real moisture_10hr,
                                           Real moisture_100hr)
{
    RothermelComputed rc;

    // Unit conversions
    const Real FT_MIN_TO_M_S = 0.00508;

    // Standard mineral contents (Rothermel 1972)
    const Real S_T = 0.0555;   // Total mineral content
    const Real S_e = 0.010;    // Effective mineral content

    // ===================================================================
    // 1. Total fuel load and weighted moisture (single-class path)
    // ===================================================================
    Real w_0 = fp.w_d1 + fp.w_d10 + fp.w_d100 + fp.w_lh + fp.w_lw;  // total oven-dry load [lb/ft²]
    
    // Compute weighted dead fuel moisture
    Real w_d = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real M_f = 0.0;  // weighted fuel moisture fraction
    if (w_d > 1.0e-6) {
        Real r_d1 = fp.w_d1 / w_d;
        Real r_d10 = fp.w_d10 / w_d;
        Real r_d100 = fp.w_d100 / w_d;
        M_f = r_d1 * moisture_1hr + r_d10 * moisture_10hr + r_d100 * moisture_100hr;
    }

    // ===================================================================
    // 2. Net fuel load (Eq. 24)
    // ===================================================================
    Real w_n = w_0 * (1.0 - S_T);

    // ===================================================================
    // 3. Bulk density and packing ratio (Eq. 24)
    // ===================================================================
    Real rho_b = w_0 / fp.delta;  // lb/ft³
    Real beta = rho_b / fp.rho_p;

    // ===================================================================
    // 4-7. Reaction velocity components (Eqs. 36-38)
    // ===================================================================
    // Use sigma_d1 as characteristic SAV
    Real sigma = fp.sigma_d1;
    sigma = amrex::max(sigma, 100.0);  // Minimum SAV guard

    Real beta_op = 3.348 * std::pow(sigma, -0.8189);           // Eq. 37: optimum packing ratio
    Real sigma_1p5 = std::pow(sigma, 1.5);
    Real Gamma_max = sigma_1p5 / (495.0 + 0.0594 * sigma_1p5); // Eq. 36: maximum reaction velocity
    Real A = 133.0 * std::pow(sigma, -0.7913);                 // Eq. 38: A coefficient
    Real beta_ratio = beta / beta_op;
    Real Gamma_prime = Gamma_max * std::pow(beta_ratio, A) * std::exp(A * (1.0 - beta_ratio)); // Eq. 38

    // ===================================================================
    // 8-10. Moisture damping (Eq. 29) and mineral damping (Eq. 30)
    // ===================================================================
    Real rm = amrex::min(M_f / fp.Mx, 1.0);
    Real eta_M = amrex::max(0.0, 1.0 - 2.59*rm + 5.11*rm*rm - 3.52*rm*rm*rm);  // Eq. 29
    Real eta_s = 0.174 * std::pow(S_e, -0.19);                                  // Eq. 30

    // ===================================================================
    // 11. Reaction intensity (Eq. 27)
    // ===================================================================
    Real I_R = Gamma_prime * w_n * fp.heat_content * eta_M * eta_s;
    I_R = amrex::max(I_R, 0.01);

    // ===================================================================
    // 12-15. Propagating flux ratio, heating number, heat of preignition, and R0 (Eqs. 1, 12, 14, 42)
    // ===================================================================
    Real xi = std::exp((0.792 + 0.681 * std::sqrt(sigma)) * (beta + 0.1)) 
              / (192.0 + 0.2595 * sigma);                                      // Eq. 42
    Real eps_h = std::exp(-138.0 / sigma);                                     // Eq. 14
    Real Q_ig = 250.0 + 1116.0 * M_f;                                          // Eq. 12: heat of preignition
    Real R0_ft_min = (I_R * xi) / (rho_b * eps_h * Q_ig);                      // Eq. 1: no-wind ROS [ft/min]

    // ===================================================================
    // 16-18. Wind factor coefficients (Eqs. 47-49)
    // ===================================================================
    Real C = 7.47 * std::exp(-0.133 * std::pow(sigma, 0.55));           // Eq. 47 (Rothermel 1972, raw sigma path)
    Real B = 0.02526 * std::pow(sigma, 0.54);                           // Eq. 48
    Real E = 0.715 * std::exp(-3.59e-4 * sigma);                        // Eq. 49

    // ===================================================================
    // 19. Packing ratio wind factor
    // ===================================================================
    Real beta_ratio_E = std::pow(beta_ratio, -E);

    // ===================================================================
    // 20. Slope factor coefficient (Eq. 51)
    // ===================================================================
    Real phi_s_const = 5.275 * std::pow(beta, -0.3);

    // ===================================================================
    // 21. MEWS wind speed cap (Andrews 2018 / Rothermel 1972)
    // ===================================================================
    // The formula phi_w_max = 0.9 * I_R mixes incompatible quantities:
    // phi_w (dimensionless) and I_R (BTU/ft²/min ~500 for fine fuels).
    // This produces a dimensionally incorrect cap of phi_w ~ 475 for FM1,
    // giving unrealistically high ROS (ROS ~ wind speed × 2-3).
    //
    // Instead, use fuel-type-based absolute cap on midflame wind speed,
    // consistent with published BEHAVE/BehavePlus validation tables:
    //   Fine fuels (sigma > 1000 ft⁻¹): cap at 300 ft/min (~1.5 m/s midflame)
    //   Coarse fuels (sigma <= 1000 ft⁻¹): cap at 500 ft/min (~2.5 m/s midflame)
    // These caps are bypassed when erf.fire.use_wind_limit = false.
    Real U_max_ftmin = (sigma > 1000.0) ? 300.0 : 500.0;

    // ===================================================================
    // Store results in RothermelComputed
    // ===================================================================
    rc.R0           = R0_ft_min * FT_MIN_TO_M_S;  // Convert ft/min → m/s
    rc.C            = C;
    rc.B            = B;
    rc.beta_ratio_E = beta_ratio_E;
    rc.beta         = beta;
    rc.phi_s_const  = phi_s_const;
    rc.U_max_ftmin  = U_max_ftmin;
    rc.wind_conv    = 196.85;    // m/s → ft/min
    rc.ros_conv     = 1.0;       // No double conversion: rc.R0 is already in m/s
    rc.I_R          = I_R;

    return rc;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
Real rothermel_ros_cell(
    Real ux_eff, Real uy_eff,
    Real sx, Real sy,
    const RothermelComputed& rc) noexcept
{
    // Effective wind speed [m/s]
    Real U_eff = std::sqrt(ux_eff*ux_eff + uy_eff*uy_eff);

    // Wind factor: φ_w = C * U_eff^B * (β/β_opt)^E
    // Convert to ft/min for Rothermel calculation
    Real U_eff_ftmin = U_eff * rc.wind_conv;

    // Apply MEWS cap
    U_eff_ftmin = amrex::min(U_eff_ftmin, rc.U_max_ftmin);

    // Slope magnitude and squared tangent (for Rothermel Eq. 51)
    Real tan_slope_sq = sx*sx + sy*sy;

    // Slope factor (Rothermel Eq. 51): φ_s = 5.275 * β^(-0.3) * tan²(φ)
    Real phi_s = rc.phi_s_const * tan_slope_sq;

    // Wind factor
    Real phi_w = rc.C * std::pow(U_eff_ftmin, rc.B) * rc.beta_ratio_E;

    // Rate of spread: R = R0 * (1 + φ_w + φ_s)
    Real ROS = rc.R0 * (1.0 + phi_w + phi_s);
    
    return amrex::max(ROS, 0.0);
}

void compute_ros_field(
    MultiFab& fire_ros,
    const MultiFab& fire_wind,
    const MultiFab& fire_slopes,
    const RothermelComputed& rc)
{
    for (MFIter mfi(fire_ros, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        Array4<Real> ros = fire_ros.array(mfi);
        Array4<const Real> wind = fire_wind.array(mfi);
        Array4<const Real> slopes = fire_slopes.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (const IntVect& iv) {
            int i = iv[0];
            int j = iv[1];
            int k = 0;

            Real ux = wind(i, j, k, 0);
            Real uy = wind(i, j, k, 1);
            Real sx = slopes(i, j, k, 0);
            Real sy = slopes(i, j, k, 1);

            ros(i, j, k) = rothermel_ros_cell(ux, uy, sx, sy, rc);
        });
    }
}
