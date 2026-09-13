#ifndef ERF_CLOUD_CHAMBER_WALL_DT_GUARD_H_
#define ERF_CLOUD_CHAMBER_WALL_DT_GUARD_H_

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <ostream>

namespace erf_cloud_chamber_wall_dt_guard {

inline bool
fixed_dt_exceeds_limit (double fixed_dt, double wall_dt) noexcept
{
    return fixed_dt > 0.0 && fixed_dt > wall_dt;
}

inline void
enforce_fixed_dt_limit (int level, double fixed_dt, double wall_dt,
                        amrex::Real max_wall_rate)
{
    if (!fixed_dt_exceeds_limit(fixed_dt, wall_dt)) { return; }

    amrex::Print() << "Cloud Chamber wall-transfer timestep violation at level "
                   << level
                   << ": fixed_dt=" << fixed_dt
                   << ", wall_dt=" << wall_dt
                   << ", max_wall_rate=" << max_wall_rate << std::endl;
    amrex::Abort(
        "Cloud Chamber wall-transfer timestep exceeds dt_wall = 0.5 / max_wall_rate");
}

} // namespace erf_cloud_chamber_wall_dt_guard

#endif
