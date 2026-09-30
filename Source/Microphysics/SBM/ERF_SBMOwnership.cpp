#include "ERF_SBMOwnership.H"

#include <AMReX.H>

#include <limits>

namespace erf_sbm {

const char* host_write_path_name(const HostWritePath path) noexcept
{
    switch (path) {
    case HostWritePath::Advection: return "advection";
    case HostWritePath::Diffusion: return "diffusion";
    case HostWritePath::Source: return "source";
    case HostWritePath::Positivity: return "positivity";
    case HostWritePath::Microphysics: return "microphysics";
    case HostWritePath::Wall: return "wall";
    case HostWritePath::Boundary: return "boundary";
    }
    return "unknown";
}

bool host_write_allowed(const bool sbm_active, const int component) noexcept
{
    return host_write_range_allowed(sbm_active, component, 1);
}

bool host_write_range_allowed(const bool sbm_active, const int start_component,
                              const int num_components) noexcept
{
    if (start_component < 0 || num_components <= 0 ||
        num_components > std::numeric_limits<int>::max() - start_component) {
        return false;
    }
    if (!sbm_active) return true;
    const int end_component = start_component + num_components;
    return !(start_component <= RhoQ2_comp && RhoQ2_comp < end_component) &&
           !(start_component <= RhoQ3_comp && RhoQ3_comp < end_component);
}

void require_host_write_allowed(const bool sbm_active, const int component,
                               const HostWritePath path)
{
    require_host_write_range_allowed(sbm_active, component, 1, path);
}

void require_host_write_range_allowed(const bool sbm_active, const int start_component,
                                      const int num_components, const HostWritePath path)
{
    if (!host_write_range_allowed(sbm_active, start_component, num_components)) {
        const std::string range = num_components > 0 && start_component >= 0 &&
            num_components <= std::numeric_limits<int>::max() - start_component
            ? " [" + std::to_string(start_component) + "," +
              std::to_string(start_component + num_components) + ")"
            : " (invalid range)";
        amrex::Abort(std::string("SBM owns compact qc/qr; native ") +
                     host_write_path_name(path) + " write denied for component range" + range);
    }
}

} // namespace erf_sbm
