#ifndef ERF_SBM_OWNERSHIP_H_
#define ERF_SBM_OWNERSHIP_H_

#include "ERF_IndexDefines.H"

#include <string>

namespace erf_sbm {

enum class HostWritePath {
    Advection,
    Diffusion,
    Source,
    Positivity,
    Microphysics,
    Wall,
    Boundary
};

[[nodiscard]] const char* host_write_path_name(HostWritePath path) noexcept;

/** The native host retains qv; SBM exclusively owns the compact liquid lanes. */
[[nodiscard]] bool host_write_allowed(bool sbm_active, int component) noexcept;

/** Check every component in a native writer's complete, non-empty write range. */
[[nodiscard]] bool host_write_range_allowed(bool sbm_active, int start_component,
                                             int num_components) noexcept;

/** Fail closed when a native host writer reaches a component owned by SBM. */
void require_host_write_allowed(bool sbm_active, int component, HostWritePath path);

/** Fail closed when any component in a native writer's range is owned by SBM. */
void require_host_write_range_allowed(bool sbm_active, int start_component,
                                      int num_components, HostWritePath path);

} // namespace erf_sbm

#endif
