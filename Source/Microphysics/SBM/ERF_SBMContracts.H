#ifndef ERF_SBM_CONTRACTS_H_
#define ERF_SBM_CONTRACTS_H_

#include "ERF_SBMLayout.H"

#include <string>
#include <vector>

namespace erf_sbm {

struct CapabilityInput {
    int max_level{0};
    bool diffusion{false};
    bool implicit_moisture_diffusion{false};
    bool shoc_or_macrophysics{false};
    bool static_cartesian{true};
    bool moving_terrain{false};
    bool embedded_boundary{false};
    bool high_order_or_fct{false};
    bool sedimentation{false};
    bool condensation{false};
    bool activation{false};
    bool collision{false};
    bool dynamic_grid{false};
    bool restart_schema_conversion{false};
    bool two_moment_transport{false};
    bool custom_moisture_forcing{false};
    bool large_scale_forcing{false};
    bool sounding_nudging{false};
    bool sponge_or_wall_modification{false};
    bool periodic_cartesian{true};
    bool double_precision{true};
    bool p2_requested{false};
    bool explicit_sbm_diffusion{false};
    int chunk_size{256};
};

struct CapabilityReport {
    bool supported{false};
    // `flags` is retained as the implementation inventory for source and
    // restart compatibility.  Qualification is intentionally reported
    // separately: an implemented path must not be presented as qualified
    // until its decisive production evidence exists.
    std::vector<std::string> flags;
    std::vector<std::string> qualified_flags;
    std::vector<std::string> qualification_limitations;
    std::vector<std::string> invariant_ids;
    std::vector<std::string> rejected_reasons;
    std::string qualification_status{"unqualified"};
    [[nodiscard]] std::string stable_description() const;
};

CapabilityReport evaluate_p1_capabilities(const CapabilityInput& input);
CapabilityReport evaluate_p2_capabilities(const CapabilityInput& input);

/** Validate a runtime bin count before it can reach a vector resize. */
[[nodiscard]] std::string validate_runtime_bin_count(int nbins);

/**
 * Host-only, stable conformance description assembled from production layout
 * and capability objects.  The optional identities are provenance strings;
 * no parser or JSON dependency is introduced in the transport path.
 */
[[nodiscard]] std::string stable_inspection(
    const SBMLayout& layout,
    const CapabilityInput& input,
    const std::string& baseline_identity = {},
    const std::string& amrex_identity = {},
    const std::string& design_identity = {});

enum class NativeWritePath {
    Advection,
    Diffusion,
    Source,
    PositivityClip,
    Microphysics,
    Wall
};

/** Host-side ownership registry used to gate all native compact-field writers. */
class OwnershipRegistry {
public:
    explicit OwnershipRegistry(bool provider_active = false) : m_provider_active(provider_active) {}
    [[nodiscard]] bool provider_active() const noexcept { return m_provider_active; }
    [[nodiscard]] bool owns_cloud_or_rain(const int component) const noexcept;
    [[nodiscard]] bool owns(const int component, NativeWritePath path) const noexcept;

private:
    bool m_provider_active{false};
};

} // namespace erf_sbm

#endif
