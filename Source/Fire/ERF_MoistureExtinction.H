// Fuel-dependent moisture of extinction (M_x).
//
// M_x varies with fuel surface-area-to-volume ratio (SAV) rather than
// using the constant 0.30 value in the base Rothermel (1972) model.
// Finer fuels (higher SAV) have lower M_x; coarser fuels (lower SAV)
// have higher M_x.
//
// Reference:
//   Rothermel, R.C. (1972). A mathematical model for predicting fire spread
//     in wildland fuels. USDA Forest Service Research Paper INT-115.
//   Anderson, H.E. (1970). Forest fuel ignitability. Fire Technology 6(4):312-319.
//
// Ported from hgopalan/wildfire_levelset src/moisture_of_extinction.H.

#ifndef ERF_MOISTURE_EXTINCTION_H
#define ERF_MOISTURE_EXTINCTION_H

#include <AMReX_REAL.H>
#include <AMReX_GpuDevice.H>
#include <cmath>

// Compute fuel-dependent moisture of extinction.
//
// M_x = M_x_min + 0.28 * (sigma_fuel / sigma_ref)^exponent
//
// Parameters:
//   sigma_fuel  — fuel SAV [ft^-1] (aggregate weighted SAV from fuel model)
//   sigma_ref   — reference SAV [ft^-1] (default 1739.0, 1-hr fuel)
//   M_x_min     — minimum M_x (default 0.12)
//   exponent    — SAV exponent (default -0.3)
//
// Returns M_x [fraction], clamped to [M_x_min, 0.45].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_moisture_of_extinction(
    amrex::Real sigma_fuel,
    amrex::Real sigma_ref,
    amrex::Real M_x_min,
    amrex::Real exponent) noexcept
{
    // Avoid division by zero
    if (sigma_ref <= 0.0 || sigma_fuel <= 0.0) {
        return 0.30;  // Default value
    }

    // Sigma ratio
    amrex::Real sigma_ratio = sigma_fuel / sigma_ref;

    // M_x = M_x_min + (M_x_max - M_x_min) × (σ_fuel / σ_ref)^exponent
    // Simplified: M_x = 0.12 + 0.28 × (σ_fuel / σ_ref)^(-0.3)
    amrex::Real M_x = M_x_min + 0.28 * std::pow(sigma_ratio, exponent);

    // Apply reasonable bounds
    M_x = amrex::max(M_x, M_x_min);
    M_x = amrex::min(M_x, 0.45);  // M_X_MAX

    return M_x;
}

// Convenience overload with Rothermel (1972) defaults.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_moisture_of_extinction(amrex::Real sigma_fuel) noexcept
{
    return compute_moisture_of_extinction(
        sigma_fuel,
        1739.0,   // SIGMA_REF_DEFAULT - reference SAV [ft^-1] (1-hr fuel)
        0.12,     // M_X_MIN_DEFAULT - minimum M_x
        -0.3      // EXPONENT_DEFAULT - SAV exponent
    );
}

// Compute the load-weighted effective dead fuel moisture for Rothermel (1972).
// M_f = (w_d1*M_1hr + w_d10*M_10hr + w_d100*M_100hr) / (w_d1 + w_d10 + w_d100)
// If total dead load is zero, returns moisture_1hr as a fallback.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_weighted_moisture(
    amrex::Real M_1hr, amrex::Real M_10hr, amrex::Real M_100hr,
    amrex::Real w_d1, amrex::Real w_d10, amrex::Real w_d100) noexcept
{
    const amrex::Real total_load = w_d1 + w_d10 + w_d100;

    if (total_load <= 0.0) {
        return M_1hr;  // Fallback if no dead fuel
    }

    amrex::Real M_f = (w_d1 * M_1hr + w_d10 * M_10hr + w_d100 * M_100hr) / total_load;

    return M_f;
}

// Compute moisture damping coefficient eta_M (Rothermel 1972 Eq. 29).
// eta_M = 1 - 2.59*(M_f/M_x) + 5.11*(M_f/M_x)^2 - 3.52*(M_f/M_x)^3
// Returns 0 when M_f >= M_x (fuel above extinction moisture).
// Clamped to [0, 1].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real compute_moisture_damping(amrex::Real M_f, amrex::Real M_x) noexcept
{
    if (M_x <= 0.0 || M_f >= M_x) {
        return 0.0;  // Fuel is above extinction moisture
    }

    const amrex::Real ratio = M_f / M_x;
    const amrex::Real ratio2 = ratio * ratio;
    const amrex::Real ratio3 = ratio2 * ratio;

    // Rothermel (1972) moisture damping coefficient
    amrex::Real eta_M = 1.0 - 2.59 * ratio + 5.11 * ratio2 - 3.52 * ratio3;

    // Clamp to [0, 1]
    eta_M = amrex::max(0.0, amrex::min(eta_M, 1.0));

    return eta_M;
}

#endif  // ERF_MOISTURE_EXTINCTION_H
