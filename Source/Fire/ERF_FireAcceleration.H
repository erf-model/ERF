// ============================================================================
// Fire Acceleration Models (Phase 12)
// ============================================================================
// Small fires do not immediately achieve their quasi-steady-state (QSS)
// rate of spread. This module implements two optional acceleration models:
//
// 1. SIZE-BASED MODEL (Rothermel 1983 / Catchpole et al. 1992)
//    - Simple exponential approach based on current fire size
//    - Applies global scaling factor to ROS field
//    - α_acc = 1 − exp(−r_fire / L_acc)
//
// 2. FARSITE TEMPORAL MODEL (McAlpine & Wakimoto 1991 / VanWagner)
//    - Per-cell temporal tracking of acceleration state
//    - Based on elapsed time since ignition or wind change
//    - R(t) = R_E × (1 − exp(−A × t))
//    - Includes optional wind-onset time-lag capability
//
// Both models are disabled by default (enable = false) with zero runtime cost
// when disabled.
//
// References:
//   Rothermel, R.C. (1983). "How to predict the spread and intensity of forest
//       and range fires." USDA Forest Service GTR INT-143.
//   Catchpole, E.A., de Mestre, N.J. & Gill, A.M. (1992). "Intensity of fire
//       at its perimeter." Australian Journal of Ecology, 17(1), 1–4.
//   McAlpine, R.S. & Wakimoto, R.H. (1991). "The acceleration of fire from
//       point source to equilibrium spread." Forest Science, 37(5), 1314–1337.
//   Alexander, M.E., Stocks, B.J. & Lawson, B.D. (1992). "Fire behavior in
//       Black Spruce-lichen woodland." Info. Rep. NOR-X-310. USDA/CFS.
//   Finney, M.A. (1998/2004). "FARSITE: Fire Area Simulator." USDA Forest
//       Service RMRS-RP-4.
// ============================================================================

#ifndef ERF_FIRE_ACCELERATION_H
#define ERF_FIRE_ACCELERATION_H

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_ParallelDescriptor.H>
#include <cmath>

// ============================================================================
// AccelerationParams - Fire acceleration model parameters
// ============================================================================
struct AccelerationParams {
    bool          enable         = false;   ///< Master on/off switch (default: disabled)
    bool          use_temporal   = false;   ///< false = size-based, true = temporal per-cell
    amrex::Real   L_acc          = 50.0;    ///< Size-based: acceleration length scale [m]
    amrex::Real   A_point        = 0.115;   ///< Temporal: acceleration constant, point fire [1/min]
    amrex::Real   A_line         = 0.886;   ///< Temporal: acceleration constant, line fire [1/min]
    amrex::Real   perim_limit    = 500.0;   ///< Temporal: perimeter length [m] switching point vs line
    bool          enable_wind_lag = false;  ///< Temporal: apply exponential lag on wind increase
    amrex::Real   tau_wind       = 60.0;    ///< Temporal: wind-lag time constant [s]
};

// ============================================================================
// apply_fire_acceleration_size_based
//
// SIZE-BASED MODEL: Computes fire area from phi and scales R_mf by
// α_acc = 1 − exp(−r_fire / L_acc).
//
// When the fire is large (r_fire >> L_acc), α_acc → 1 and R_mf is unchanged.
// When the fire is small (r_fire << L_acc), α_acc → 0 and R_mf → 0.
//
// Parameters:
//   R_mf      – rate-of-spread field [m/s] (modified in place)
//   phi       – level-set field (phi < 0 = burned)
//   geom      – AMReX geometry (provides cell sizes)
//   L_acc     – acceleration length scale [m]
//   fire_debug – whether to print debug diagnostics
// ============================================================================
inline void apply_fire_acceleration_size_based(amrex::MultiFab&       R_mf,
                                               const amrex::MultiFab& phi,
                                               const amrex::Geometry& geom,
                                               const amrex::Real      L_acc,
                                               bool                   fire_debug)
{
    const auto dx = geom.CellSize();

    // Count burned cells (CPU loop; same pattern as perimeter counting)
    amrex::Long n_burned = 0;
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const p         = phi.const_array(mfi);
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (p(i, j, k) < 0.0_rt) ++n_burned;
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_burned);

    if (n_burned == 0) return;  // no fire yet — skip

    // Burned area and effective fire radius
    const amrex::Real A_fire     = static_cast<amrex::Real>(n_burned) * dx[0] * dx[1];
    const amrex::Real r_fire     = std::sqrt(A_fire / amrex::Math::pi<amrex::Real>());
    const amrex::Real L_acc_safe = amrex::max(L_acc, 1.0_rt);  // floor prevents divide-by-zero
    const amrex::Real alpha      = 1.0_rt - std::exp(-r_fire / L_acc_safe);

    // Nothing to do when essentially at full speed (saves a kernel launch)
    if (alpha >= 1.0_rt - 1.0e-4_rt) return;

    if (fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Fire acceleration (size-based): r_fire=" << r_fire
                       << " m  L_acc=" << L_acc << " m  alpha=" << alpha << "\n";
    }

    // Scale R_mf by alpha
    for (amrex::MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto R               = R_mf.array(mfi);
        const amrex::Real a  = alpha;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            R(i, j, k) *= a;
        });
    }
}

// ============================================================================
// apply_fire_acceleration_temporal
//
// FARSITE TEMPORAL MODEL: Per-cell acceleration using VanWagner's equation
// R(t) = R_E × (1 − exp(−A × t)).
//
// Each cell tracks:
//   - ROS_current: current instantaneous ROS [m/s]
//   - elapsed_time: time since entering current state [s]
//   - ROS_equilibrium: target equilibrium ROS [m/s]
//
// When R_mf (quasi-steady Rothermel ROS) changes (e.g., due to wind change),
// the acceleration resets and ramps up exponentially.
//
// Parameters:
//   R_mf            – equilibrium ROS from Rothermel model [m/s] (read; output written to R_mf)
//   phi             – level-set field (phi < 0 = burned)
//   geom            – AMReX geometry
//   params          – acceleration parameters
//   dt              – timestep [s]
//   accel_state_mf  – 3-component state (ROS_current, elapsed_time, ROS_equilibrium)
//   fire_debug      – whether to print debug diagnostics
// ============================================================================
inline void apply_fire_acceleration_temporal(amrex::MultiFab&       R_mf,
                                             const amrex::MultiFab& phi,
                                             const amrex::Geometry& geom,
                                             const AccelerationParams& params,
                                             const amrex::Real      dt,
                                             amrex::MultiFab*       accel_state_mf,
                                             bool                   fire_debug)
{
    if (accel_state_mf == nullptr) {
        amrex::Abort("apply_fire_acceleration_temporal: accel_state_mf is null");
    }

    const auto dx = geom.CellSize();

    // Compute fire perimeter to determine acceleration constant
    // (point ignition: A_point; line ignition: A_line)
    amrex::Long n_perimeter = 0;
    for (amrex::MFIter mfi(phi); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto const p         = phi.const_array(mfi);
        const auto lo        = amrex::lbound(bx);
        const auto hi        = amrex::ubound(bx);
        
        amrex::LoopOnCpu(bx, [&](int i, int j, int k) {
            if (p(i, j, k) < 0.0_rt) {
                // Check if this is a perimeter cell (has unburned neighbor)
                // Guard all stencil accesses with validbox bounds
                bool is_perimeter = false;
                if (i+1 <= hi.x && p(i+1, j, k) >= 0.0_rt) is_perimeter = true;
                if (i-1 >= lo.x && p(i-1, j, k) >= 0.0_rt) is_perimeter = true;
                if (j+1 <= hi.y && p(i, j+1, k) >= 0.0_rt) is_perimeter = true;
                if (j-1 >= lo.y && p(i, j-1, k) >= 0.0_rt) is_perimeter = true;
                if (is_perimeter) ++n_perimeter;
            }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(n_perimeter);

    // Estimate perimeter length in meters
    const amrex::Real perim_length = static_cast<amrex::Real>(n_perimeter) * 
                                     std::sqrt(dx[0] * dx[0] + dx[1] * dx[1]);

    // Select acceleration constant based on perimeter size
    const amrex::Real A = (perim_length > params.perim_limit) ? params.A_line : params.A_point;
    
    // Convert A from [1/min] to [1/s]
    const amrex::Real A_per_sec = A / 60.0_rt;

    if (fire_debug) {
        amrex::Print() << "[FIRE DEBUG] Fire acceleration (temporal): perim=" << perim_length 
                       << " m  A=" << A << " 1/min  (using " 
                       << (perim_length > params.perim_limit ? "line" : "point") 
                       << " constant)\n";
    }

    // Update acceleration state for each cell
    for (amrex::MFIter mfi(R_mf); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        auto R_eq            = R_mf.const_array(mfi);          // equilibrium ROS from Rothermel
        auto state           = accel_state_mf->array(mfi);     // [ROS_current, elapsed_time, ROS_equilibrium]
        auto p               = phi.const_array(mfi);           // level set
        auto R_out           = R_mf.array(mfi);                // output (will be modified)

        const amrex::Real dt_loc        = dt;
        const amrex::Real A_loc         = A_per_sec;
        const int  enable_wind_lag_loc  = params.enable_wind_lag ? 1 : 0;
        const amrex::Real tau_wind_loc  = params.tau_wind;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            // Only accelerate burned cells
            if (p(i, j, k) >= 0.0_rt) {
                state(i, j, k, 0) = 0.0_rt;  // ROS_current
                state(i, j, k, 1) = 0.0_rt;  // elapsed_time
                state(i, j, k, 2) = 0.0_rt;  // ROS_equilibrium
                return;
            }

            const amrex::Real R_E = R_eq(i, j, k);  // equilibrium ROS (target)
            amrex::Real R_curr    = state(i, j, k, 0);
            amrex::Real t_elapsed = state(i, j, k, 1);
            amrex::Real R_prev_eq = state(i, j, k, 2);

            // Detect equilibrium ROS change (wind change or fuel change)
            const amrex::Real EQUILIBRIUM_CHANGE_THRESHOLD = 1.0e-6_rt;
            const amrex::Real dR_eq = amrex::Math::abs(R_E - R_prev_eq);
            const bool eq_changed = (dR_eq > EQUILIBRIUM_CHANGE_THRESHOLD * 
                                     amrex::max(R_E, R_prev_eq, 1.0_rt));

            if (eq_changed) {
                // Equilibrium changed — apply wind-lag or reset
                if (enable_wind_lag_loc == 1 && R_E > R_prev_eq) {
                    // Wind increase: apply exponential lag
                    // New target approaches R_E with time constant tau_wind
                    // R_target(t) = R_prev + (R_E - R_prev) × (1 - exp(-dt/tau_wind))
                    const amrex::Real lag_factor = 1.0_rt - std::exp(-dt_loc / tau_wind_loc);
                    const amrex::Real R_target = R_prev_eq + (R_E - R_prev_eq) * lag_factor;
                    
                    // Start new acceleration towards lagged target
                    R_prev_eq = R_target;
                    t_elapsed = 0.0_rt;
                } else {
                    // Wind decrease or lag disabled: reset immediately to new equilibrium
                    R_prev_eq = R_E;
                    t_elapsed = 0.0_rt;
                }
            }

            // VanWagner acceleration equation: R(t) = R_E × (1 - exp(-A × t))
            t_elapsed += dt_loc;
            R_curr = R_prev_eq * (1.0_rt - std::exp(-A_loc * t_elapsed));

            // Clamp to equilibrium (avoid overshoot due to numerical error)
            if (R_curr > R_prev_eq) R_curr = R_prev_eq;

            // Update state
            state(i, j, k, 0) = R_curr;
            state(i, j, k, 1) = t_elapsed;
            state(i, j, k, 2) = R_prev_eq;

            // Write accelerated ROS to output
            R_out(i, j, k) = R_curr;
        });
    }
}

// ============================================================================
// apply_fire_acceleration
//
// Main dispatch function: selects size-based or temporal model based on params.
// Returns immediately when enable == false (zero cost when disabled).
//
// Parameters:
//   R_mf            – rate-of-spread field [m/s] (modified in place)
//   phi             – level-set field (phi < 0 = burned)
//   geom            – AMReX geometry
//   params          – AccelerationParams from FireParams
//   dt              – timestep [s] (only used for temporal model)
//   accel_state_mf  – per-cell state for temporal model (nullptr for size-based)
//   fire_debug      – whether to print debug diagnostics
// ============================================================================
inline void apply_fire_acceleration(amrex::MultiFab&       R_mf,
                                    const amrex::MultiFab& phi,
                                    const amrex::Geometry& geom,
                                    const AccelerationParams& params,
                                    const amrex::Real      dt,
                                    amrex::MultiFab*       accel_state_mf,
                                    bool                   fire_debug)
{
    if (!params.enable) return;

    if (params.use_temporal) {
        // FARSITE temporal model
        apply_fire_acceleration_temporal(R_mf, phi, geom, params, dt, accel_state_mf, fire_debug);
    } else {
        // Simple size-based model
        apply_fire_acceleration_size_based(R_mf, phi, geom, params.L_acc, fire_debug);
    }
}

#endif // ERF_FIRE_ACCELERATION_H
