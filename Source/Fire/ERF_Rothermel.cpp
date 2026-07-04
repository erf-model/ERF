#include <ERF_Rothermel.H>


using namespace amrex;

// Anderson FBFM13 fuel model database
// Reference: Anderson, H.E. (1982) "Aids to determining fuel models..."
// All values are standard published values for these fuel models

FuelModelParams get_anderson_fuel_params(int model_id)
{
    FuelModelParams fp;

    // All fuel loads in lb/ft² (will be converted if needed)
    // SAV in 1/ft
    // Heat content in BTU/lb
    // Particle density in lb/ft³
    // Depth in ft

    switch (model_id) {
    case 1:  // Short grass
        fp.w_d1 = 0.074;   fp.w_d10 = 0.0;    fp.w_d100 = 0.0;
        fp.w_lh = 0.019;   fp.w_lw = 0.0;
        fp.sigma_d1 = 8000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.0;    fp.Mx = 0.12;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 2:  // Timber-grass
        fp.w_d1 = 0.115;   fp.w_d10 = 0.045;  fp.w_d100 = 0.0;
        fp.w_lh = 0.05;    fp.w_lw = 0.0;
        fp.sigma_d1 = 8000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.0;    fp.Mx = 0.15;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 3:  // Tall grass
        fp.w_d1 = 0.325;   fp.w_d10 = 0.0;    fp.w_d100 = 0.0;
        fp.w_lh = 0.165;   fp.w_lw = 0.0;
        fp.sigma_d1 = 8000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 4:  // Chaparral
        fp.w_d1 = 0.184;   fp.w_d10 = 0.092;  fp.w_d100 = 0.023;
        fp.w_lh = 0.276;   fp.w_lw = 0.0;
        fp.sigma_d1 = 6800.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 6.0;    fp.Mx = 0.20;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 5:  // Shrub
        fp.w_d1 = 0.184;   fp.w_d10 = 0.138;  fp.w_d100 = 0.092;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 6500.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.20;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 6:  // Dormant brush
        fp.w_d1 = 0.092;   fp.w_d10 = 0.115;  fp.w_d100 = 0.092;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 6500.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 7:  // Southern rough
        fp.w_d1 = 0.230;   fp.w_d10 = 0.184;  fp.w_d100 = 0.138;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 6000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.5;    fp.Mx = 0.40;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 8:  // Closed timber litter
        fp.w_d1 = 0.276;   fp.w_d10 = 0.046;  fp.w_d100 = 0.023;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 7000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 0.2;    fp.Mx = 0.30;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 9:  // Hardwood litter
        fp.w_d1 = 0.138;   fp.w_d10 = 0.092;  fp.w_d100 = 0.115;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 6000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 0.2;    fp.Mx = 0.40;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 10:  // Timber litter-understory
        fp.w_d1 = 0.092;   fp.w_d10 = 0.138;  fp.w_d100 = 0.230;
        fp.w_lh = 0.046;   fp.w_lw = 0.023;
        fp.sigma_d1 = 6000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 1500.0;
        fp.delta = 1.0;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 11:  // Light logging slash
        fp.w_d1 = 0.230;   fp.w_d10 = 0.322;  fp.w_d100 = 0.046;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 6000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 1.5;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 12:  // Medium logging slash
        fp.w_d1 = 0.575;   fp.w_d10 = 0.575;  fp.w_d100 = 0.092;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 5500.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
        fp.delta = 2.0;    fp.Mx = 0.25;      fp.heat_content = 8000.0;
        fp.rho_p = 32.0;
        break;

    case 13:  // Heavy logging slash
        fp.w_d1 = 1.058;   fp.w_d10 = 1.268;  fp.w_d100 = 0.322;
        fp.w_lh = 0.0;     fp.w_lw = 0.0;
        fp.sigma_d1 = 5000.0;  fp.sigma_lh = 60.0;  fp.sigma_lw = 0.0;
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
    const Real LB_FT2_TO_KG_M2 = 4.88243;
    const Real FT_TO_M = 0.3048;
    const Real FT_MIN_TO_M_S = 0.00508;

    // Total fuel load [lb/ft²]
    Real w_d = fp.w_d1 + fp.w_d10 + fp.w_d100;
    Real w_l = fp.w_lh + fp.w_lw;
    Real w_0 = w_d + fp.w_lh + fp.w_lw;

    // Fuel load ratios
    Real r_d1 = (w_d > 1.0e-6) ? fp.w_d1 / w_d : 0.0;
    Real r_d10 = (w_d > 1.0e-6) ? fp.w_d10 / w_d : 0.0;
    Real r_d100 = (w_d > 1.0e-6) ? fp.w_d100 / w_d : 0.0;

    // Surface area to volume ratio (weighted)
    Real sigma = fp.sigma_d1 * r_d1 + fp.sigma_d1 * r_d10 + fp.sigma_d1 * r_d100;
    if (fp.w_lh > 1.0e-6) sigma += fp.sigma_lh * (fp.w_lh / w_0);
    sigma = amrex::max(sigma, 100.0);  // Minimum SAV

    // Packing ratio (dry bulk density / particle density)
    Real rho_b = w_0 / fp.delta;  // lb/ft³
    Real beta = rho_b / fp.rho_p;

    // Optimal packing ratio
    Real beta_opt = 3.348 / std::pow(sigma, 0.8189);

    // Fuel moisture (weighted average)
    Real M_x = moisture_1hr * r_d1 + moisture_10hr * r_d10 + moisture_100hr * r_d100;
    if (fp.w_lh > 1.0e-6) M_x += moisture_1hr * (fp.w_lh / w_0);

    // Moisture of extinction
    Real M_x_limit = fp.Mx;
    M_x = amrex::min(M_x, M_x_limit);

    // Reaction intensity (BTU/ft²/min)
    Real A = 133.0 / std::pow(sigma, 0.7913);
    Real E_s = 0.75 - 0.00023 * (M_x_limit - M_x * 100.0);
    Real gamma_x = (sigma / (460.0 + 25.9 * M_x)) * std::exp(fp.delta * (std::log(10.0) / 3.0));

    // Damping coefficients
    Real B_d = 0.02526 * std::pow(sigma, 0.54);
    Real C_d = 7.47 * std::exp(-0.8711 * std::pow(sigma, -0.55));
    Real E_dust = 0.595 * std::pow((M_x_limit - M_x * 100.0), -1.628);

    Real eta_M = 1.0 - 2.59 * (M_x / M_x_limit) + 5.11 * std::pow(M_x / M_x_limit, 2.0)
               - 3.861 * std::pow(M_x / M_x_limit, 3.0);
    eta_M = amrex::max(0.0, eta_M);

    Real eta_S = 0.40;  // Approximation for slope/aspect effect

    Real I_R = gamma_x * w_0 * fp.heat_content * eta_M * eta_S;
    I_R = amrex::max(I_R, 0.01);

    // No-wind, no-slope ROS [ft/min] (Eq. 3 from Rothermel 1972)
    Real R0_ft_min = A * std::pow(beta / beta_opt, B_d) *
                     std::exp((C_d * (1.0 - beta / beta_opt)));

    // Wind factor coefficients
    Real E_wind = 0.92 * std::pow(sigma / 1000.0, -0.3);
    Real B_wind = 0.02526 * std::pow(sigma, 0.54);

    // Convert to SI
    rc.R0 = R0_ft_min * FT_MIN_TO_M_S;
    rc.C = I_R / (173.0 + 60.0);  // Proportional constant (empirical)
    rc.B = 0.255;  // Standard exponent from Rothermel
    rc.beta_ratio_E = 1.0;  // (β/β_opt)^E, set E to give neutral effect
    rc.beta = beta;
    rc.phi_s_const = 0.0;  // Base slope factor (computed per-cell with terrain)
    rc.U_max_ftmin = amrex::max(R0_ft_min, 100.0);  // MEWS cap
    rc.wind_conv = 196.85;  // m/s to ft/min
    rc.ros_conv = FT_MIN_TO_M_S;
    rc.I_R = I_R;

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

    // Slope magnitude [dimensionless]
    Real slope_mag = std::sqrt(sx*sx + sy*sy);
    Real tan_slope = slope_mag;  // For small slopes, tan(θ) ≈ slope
    if (slope_mag > 0.1) tan_slope = std::atan(slope_mag);  // Exact for large slopes

    // Wind direction
    Real wind_dir = std::atan2(uy_eff, ux_eff);
    // Slope direction (maximum slope)
    Real slope_dir = std::atan2(sy, sx);
    // Angle between wind and slope
    Real angle = wind_dir - slope_dir;

    // Normalize to [-π, π]
    const Real pi = 3.14159265358979323846;
    while (angle > pi) angle -= 2.0*pi;
    while (angle < -pi) angle += 2.0*pi;

    // Slope factor: φ_s = tan(θ) * cos(angle)
    Real phi_s = tan_slope * std::cos(angle);

    // Wind factor
    Real phi_w = rc.C * std::pow(U_eff_ftmin, rc.B) * rc.beta_ratio_E;

    // Rate of spread: R = R0 * (1 + φ_w + φ_s)
    Real ROS = rc.R0 * (1.0 + phi_w + phi_s);
    ROS = amrex::max(ROS, 0.0);

    return ROS * rc.ros_conv;  // Convert ft/min to m/s
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
