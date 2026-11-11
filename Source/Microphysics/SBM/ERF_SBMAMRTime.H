#ifndef ERF_SBM_AMR_TIME_H_
#define ERF_SBM_AMR_TIME_H_

#include <algorithm>
#include <cmath>
#include <limits>

namespace erf_sbm {

/** Compare values carrying one semantic time identity without admitting a
 * neighboring physical substep.  Callers must provide the smallest physical
 * step represented by the comparison.
 */
[[nodiscard]] inline bool semantic_time_matches (const double lhs,
                                                  const double rhs,
                                                  const double dt_min) noexcept
{
    if (!std::isfinite(lhs) || !std::isfinite(rhs) ||
        !std::isfinite(dt_min) || !(dt_min > 0.0)) {
        return false;
    }
    const double scale = std::max({1.0, std::abs(lhs), std::abs(rhs)});
    const double roundoff =
        16.0 * std::numeric_limits<double>::epsilon() * scale;
    const double tolerance = std::min(roundoff, dt_min / 32.0);
    return std::abs(lhs - rhs) <= tolerance;
}

/** Transfer tuples without a physical-step scale must share the exact finite
 * timestamp supplied by their owner.
 */
[[nodiscard]] inline bool semantic_time_matches_exact (const double lhs,
                                                        const double rhs) noexcept
{
    return std::isfinite(lhs) && std::isfinite(rhs) && lhs == rhs;
}

/** Whether a level's accepted state is at the host level time: to roundoff of
 * the level's step, or exactly before any step has been taken.
 */
[[nodiscard]] inline bool level_time_matches (const double state_time,
                                              const double level_time,
                                              const double level_dt) noexcept
{
    return level_dt > 0.0
        ? semantic_time_matches(state_time, level_time, level_dt)
        : semantic_time_matches_exact(state_time, level_time);
}

} // namespace erf_sbm

#endif
