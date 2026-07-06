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

#endif
