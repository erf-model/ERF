// Fire behavior diagnostics: Byram fireline intensity and flame length.
//
// Byram (1959) fireline intensity:
//   I_B [kW/m] = h_fuel [kJ/kg] * w_consumed [kg/m²] * ros [m/s]
//   where w_consumed = fuel_load_initial * f_c
//
// Thomas (1963) flame length:
//   L [m] = 0.0775 * I_B^0.46    (I_B in kW/m)
//
// References:
//   Byram, G.M. (1959). Combustion of forest fuels. In: Forest Fire:
//     Control and Use, pp. 61-89.
//   Thomas, P.H. (1963). The size of flames from natural fires.
//     Ninth Int. Symposium on Combustion, pp. 844-859.

#ifndef ERF_FIRE_DIAGNOSTICS_H
#define ERF_FIRE_DIAGNOSTICS_H

#include <AMReX_MultiFab.H>
#include <cmath>
#include <string>

/**
 * @brief Compute Byram (1959) fireline intensity [kW/m].
 *
 * I_B = h_fuel_kJ_per_kg * (fuel_load_initial - fuel_load_remain) * ros
 *
 * The consumed fuel is the difference between initial and remaining loads.
 * Returns 0 for unburned cells (phi >= 0).
 *
 * @param[in] phi                  Level-set field value (burned if < 0)
 * @param[in] ros                  Rate of spread [m/s]
 * @param[in] fuel_load_initial    Initial fuel load [kg/m²]
 * @param[in] fuel_load_remain     Remaining fuel load [kg/m²]
 * @param[in] h_fuel_kJ_per_kg     Heat content [kJ/kg]
 * @return Byram fireline intensity [kW/m], non-negative
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_fireline_intensity_kW_per_m(
    amrex::Real phi,
    amrex::Real ros,
    amrex::Real fuel_load_initial,
    amrex::Real fuel_load_remain,
    amrex::Real h_fuel_kJ_per_kg) noexcept
{
    // Only burned cells contribute to fireline intensity
    if (phi >= 0.0) return 0.0;

    // Consumed fuel is the difference between initial and remaining
    amrex::Real w_consumed = fuel_load_initial - fuel_load_remain;

    // No intensity if consumed fuel or ROS is negligible
    if (w_consumed <= 0.0 || ros <= 0.0) {
        return 0.0;
    }

    // I_B = h * w_consumed * ros
    amrex::Real I_B = h_fuel_kJ_per_kg * w_consumed * ros;
    return amrex::max(0.0, I_B);
}

/**
 * @brief Compute Thomas (1963) flame length [m].
 *
 * L = 0.0775 * I_B^0.46  (I_B in kW/m)
 *
 * Returns 0 for I_B <= 0.
 *
 * @param[in] I_B_kW_per_m    Byram fireline intensity [kW/m]
 * @return Flame length [m], non-negative
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flame_length_m(amrex::Real I_B_kW_per_m) noexcept
{
    if (I_B_kW_per_m <= 0.0) return 0.0;

    // Thomas (1963): L = 0.0775 * I_B^0.46
    amrex::Real L = 0.0775 * std::pow(I_B_kW_per_m, 0.46);
    return amrex::max(0.0, L);
}

/**
 * @brief Fill fireline_intensity [kW/m] and flame_length [m] MultiFabs.
 *
 * Reads phi, ros, fuel_load_remain from fire layer fields.
 * Computes Byram fireline intensity and Thomas flame length for each cell.
 *
 * fuel_load_initial is the initial total fuel load for the fuel model [kg/m²].
 * h_fuel_kJ_per_kg = fp.heat_content [BTU/lb] * 2.326  (1 BTU/lb = 2.326 kJ/kg)
 *
 * This is a local operation: no MPI communication required.
 *
 * @param[out] fireline_intensity  Byram fireline intensity [kW/m]
 * @param[out] flame_length        Thomas flame length [m]
 * @param[in] fire_phi             Level-set field
 * @param[in] fire_ros             Rate of spread [m/s]
 * @param[in] fire_fuel_load       Remaining fuel load [kg/m²]
 * @param[in] fuel_load_initial_kg_m2  Initial fuel load [kg/m²]
 * @param[in] h_fuel_kJ_per_kg     Heat content [kJ/kg]
 */
inline void fill_fire_diagnostics(
    amrex::MultiFab& fireline_intensity,
    amrex::MultiFab& flame_length,
    const amrex::MultiFab& fire_phi,
    const amrex::MultiFab& fire_ros,
    const amrex::MultiFab& fire_fuel_load,
    amrex::Real fuel_load_initial_kg_m2,
    amrex::Real h_fuel_kJ_per_kg)
{
    using namespace amrex;

    for (MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const phi_arr = fire_phi.const_array(mfi);
        auto const ros_arr = fire_ros.const_array(mfi);
        auto const fuel_arr = fire_fuel_load.const_array(mfi);
        auto I_B_arr = fireline_intensity.array(mfi);
        auto L_arr = flame_length.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(const IntVect& iv) {
            int i = iv[0], j = iv[1], k = iv[2];

            Real phi_val = phi_arr(i, j, k);
            Real ros_val = ros_arr(i, j, k);
            Real w_remain = fuel_arr(i, j, k);

            // Compute Byram fireline intensity
            Real I_B = compute_fireline_intensity_kW_per_m(
                phi_val, ros_val, fuel_load_initial_kg_m2, w_remain, h_fuel_kJ_per_kg);
            I_B_arr(i, j, k) = I_B;

            // Compute Thomas flame length from Byram intensity
            Real L = compute_flame_length_m(I_B);
            L_arr(i, j, k) = L;
        });
    }
}

// Flame temperature estimation methods.
//
// Three methods are available, selected via fire_params.flame_temp_method:
//
// Method "byram_radiant":
//   Uses Byram's (1959) empirical relation between fireline intensity and
//   flame radiation temperature. Approximates the effective radiating
//   temperature of the flame column:
//     T_flame [K] = T_amb + 800.0 * (I_B / 1000.0)^0.25
//   where T_amb = 300 K (ambient), I_B is Byram fireline intensity [kW/m].
//   Reference: Byram, G.M. (1959). Combustion of forest fuels.
//     In Forest Fire: Control and Use, pp. 61-89.
//
// Method "mcalpine_heat":
//   Based on McAlpine & Xanthopoulos (1989) adiabatic flame temperature
//   from heat of combustion and stoichiometry:
//     T_flame [K] = T_amb + h_fuel [J/kg] * (1 - M_f) / (c_p_products * (1 + r_stoich))
//   where c_p_products = 1200 J/(kg·K) (specific heat of combustion products),
//   r_stoich = 6.0 (stoichiometric air-fuel ratio for cellulose), and
//   M_f is fuel moisture fraction. Clipped to T_amb + 100 < T_flame < 1500 K.
//   Reference: McAlpine, R.S. & Xanthopoulos, G. (1989). Predicted flame
//     temperature for laboratory fires burning under simulated slope and wind
//     conditions. *Int. J. Wildland Fire*, Proceedings.
//
// Method "nelson_emissivity":
//   Nelson (1980) method using flame emissivity and Stefan-Boltzmann law.
//   Relates observed radiation flux to an effective blackbody temperature:
//     Q_rad = sigma * epsilon * T_flame^4  =>  T_flame = (I_B * 1000 / (sigma * epsilon))^0.25
//   where sigma = 5.67e-8 W/(m²·K⁴) (Stefan-Boltzmann), epsilon = 0.9 (flame emissivity),
//   and I_B is fireline intensity [kW/m] converted to W/m.
//   Clipped to range [300, 1800] K.
//   Reference: Nelson, R.M. Jr. (1980). Flame characteristics for fires in
//     southern fuels. USDA Forest Service Research Paper SE-205.
//
// All methods return 0.0 at unburned cells (phi >= 0) and cells with I_B <= 0.

/**
 * @brief Byram (1959) radiant flame temperature
 *
 * T_flame [K] = T_amb + 800.0 * (I_B / 1000.0)^0.25
 *
 * @param[in] I_B_kWm    Byram fireline intensity [kW/m]
 * @param[in] T_amb      Ambient temperature [K]
 * @return Flame temperature [K], or 0 for I_B <= 0
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flame_temp_byram(amrex::Real I_B_kWm, amrex::Real T_amb) noexcept
{
    if (I_B_kWm <= 0.0) return 0.0;
    amrex::Real T_flame = T_amb + 800.0 * std::pow(I_B_kWm / 1000.0, 0.25);
    return amrex::max(0.0, T_flame);
}

/**
 * @brief McAlpine & Xanthopoulos (1989) adiabatic flame temperature
 *
 * T_flame [K] = T_amb + h_fuel [J/kg] * (1 - M_f) / (c_p_products * (1 + r_stoich))
 * where c_p_products = 1200 J/(kg·K), r_stoich = 6.0.
 *
 * @param[in] I_B_kWm       Byram fireline intensity [kW/m] (unused, for interface consistency)
 * @param[in] h_fuel_Jkg    Heat content of fuel [J/kg]
 * @param[in] M_f           Fuel moisture fraction
 * @param[in] T_amb         Ambient temperature [K]
 * @return Flame temperature [K], clipped to (T_amb+100, 1500)
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flame_temp_mcalpine(amrex::Real I_B_kWm, amrex::Real h_fuel_Jkg,
                                         amrex::Real M_f, amrex::Real T_amb) noexcept
{
    if (I_B_kWm <= 0.0) return 0.0;
    constexpr amrex::Real c_p_products = 1200.0;  // [J/(kg·K)]
    constexpr amrex::Real r_stoich = 6.0;         // stoichiometric air-fuel ratio
    amrex::Real rise = h_fuel_Jkg * (1.0 - M_f) / (c_p_products * (1.0 + r_stoich));
    amrex::Real T_flame = T_amb + rise;
    // Clip to physically reasonable range
    amrex::Real T_min = T_amb + 100.0;
    amrex::Real T_max = 1500.0;
    T_flame = amrex::max(T_min, amrex::min(T_flame, T_max));
    return T_flame;
}

/**
 * @brief Nelson (1980) emissivity-based flame temperature
 *
 * T_flame [K] = (I_B * 1000 / (sigma * epsilon))^0.25
 * where sigma = 5.67e-8 W/(m²·K⁴), epsilon = 0.9.
 *
 * @param[in] I_B_kWm    Byram fireline intensity [kW/m]
 * @return Flame temperature [K], clipped to [300, 1800]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flame_temp_nelson(amrex::Real I_B_kWm) noexcept
{
    if (I_B_kWm <= 0.0) return 0.0;
    constexpr amrex::Real sigma = 5.67e-8;  // Stefan-Boltzmann [W/(m²·K⁴)]
    constexpr amrex::Real epsilon = 0.9;    // flame emissivity
    amrex::Real Q_rad_Wm = I_B_kWm * 1000.0;  // convert kW/m to W/m
    amrex::Real T_flame = std::pow(Q_rad_Wm / (sigma * epsilon), 0.25);
    // Clip to reasonable range
    T_flame = amrex::max(300.0, amrex::min(T_flame, 1800.0));
    return T_flame;
}

/**
 * @brief Fill flame_temperature [K] MultiFab using the selected method
 *
 * Selects one of three flame temperature methods and fills the output field.
 * Returns 0 at unburned cells.
 *
 * @param[out] fire_flame_temp     Flame temperature field [K]
 * @param[in]  fireline_intensity  Byram intensity [kW/m]
 * @param[in]  fire_phi            Level-set field
 * @param[in]  method              Method selection: "byram_radiant", "mcalpine_heat", "nelson_emissivity"
 * @param[in]  h_fuel_Jkg          Heat content [J/kg] (for mcalpine_heat)
 * @param[in]  M_f                 Fuel moisture fraction (for mcalpine_heat)
 * @param[in]  T_amb               Ambient temperature [K]
 */
inline void fill_flame_temperature(amrex::MultiFab& fire_flame_temp,
                                    const amrex::MultiFab& fireline_intensity,
                                    const amrex::MultiFab& fire_phi,
                                    const std::string& method,
                                    amrex::Real h_fuel_Jkg,
                                    amrex::Real M_f,
                                    amrex::Real T_amb)
{
    using namespace amrex;

    int method_id = 0;
    if (method == "mcalpine_heat") {
        method_id = 1;
    } else if (method == "nelson_emissivity") {
        method_id = 2;
    }

    for (MFIter mfi(fire_phi); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const phi_arr = fire_phi.const_array(mfi);
        auto const I_B_arr = fireline_intensity.const_array(mfi);
        auto T_flame_arr = fire_flame_temp.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (phi_arr(i, j, k) >= 0.0) {
                T_flame_arr(i, j, k) = 0.0;
                return;
            }

            Real I_B = I_B_arr(i, j, k);
            Real T_flame = 0.0;

            if (method_id == 1) {
                T_flame = compute_flame_temp_mcalpine(I_B, h_fuel_Jkg, M_f, T_amb);
            } else if (method_id == 2) {
                T_flame = compute_flame_temp_nelson(I_B);
            } else {
                T_flame = compute_flame_temp_byram(I_B, T_amb);
            }

            T_flame_arr(i, j, k) = T_flame;
        });
    }
}

/**
 * @brief Compute a Balbi-style flame tilt angle diagnostic [degrees from vertical].
 *
 * Uses the local effective 2-D wind speed and a buoyancy velocity estimated
 * from Byram fireline intensity, Thomas flame length, and Stefan-Boltzmann
 * radiative flux. Returns 0 for unburned or zero-intensity cells.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_flame_tilt_angle_deg(
    amrex::Real I_B_kWm,
    amrex::Real wind_u,
    amrex::Real wind_v,
    amrex::Real rho_air,
    amrex::Real T_amb) noexcept
{
    if (I_B_kWm <= 0.0 || rho_air <= 0.0 || T_amb <= 0.0) {
        return 0.0;
    }

    constexpr amrex::Real g       = 9.81;
    constexpr amrex::Real PI      = 3.14159265358979323846;
    constexpr amrex::Real sigma   = 5.67e-8;
    constexpr amrex::Real Cp_air  = 1005.0;

    amrex::Real wind_speed = std::sqrt(wind_u * wind_u + wind_v * wind_v);
    if (wind_speed <= 0.0) {
        return 0.0;
    }

    amrex::Real flame_length = compute_flame_length_m(I_B_kWm);
    amrex::Real T_flame = compute_flame_temp_byram(I_B_kWm, T_amb);
    amrex::Real q_rad = sigma * (std::pow(T_flame, 4.0) - std::pow(T_amb, 4.0));
    q_rad = amrex::max(q_rad, 0.0);

    if (flame_length <= 0.0 || q_rad <= 0.0) {
        return 0.0;
    }

    amrex::Real buoyancy = std::sqrt(g * q_rad * flame_length / (rho_air * Cp_air * T_amb));
    if (buoyancy <= 0.0) {
        return 90.0;
    }

    amrex::Real tilt_rad = std::atan(wind_speed / buoyancy);
    return amrex::min(tilt_rad * 180.0 / PI, 90.0);
}

/**
 * @brief Fill flame tilt angle [degrees from vertical].
 *
 * Uses Byram fireline intensity and the local 2-component effective wind to
 * estimate a Balbi-style flame tilt angle diagnostic.
 */
inline void fill_flame_tilt_angle(
    amrex::MultiFab& fire_flame_tilt,
    const amrex::MultiFab& fireline_intensity,
    const amrex::MultiFab& fire_wind_eff,
    amrex::Real rho_air,
    amrex::Real T_amb)
{
    using namespace amrex;

    for (MFIter mfi(fire_flame_tilt); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();

        auto const I_B_arr = fireline_intensity.const_array(mfi);
        auto const wind_arr = fire_wind_eff.const_array(mfi);
        auto tilt_arr = fire_flame_tilt.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tilt_arr(i, j, k) = compute_flame_tilt_angle_deg(
                I_B_arr(i, j, k),
                wind_arr(i, j, k, 0),
                wind_arr(i, j, k, 1),
                rho_air,
                T_amb);
        });
    }
}

#endif
