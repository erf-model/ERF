#include "ERF_SBMContracts.H"

#include "ERF_IndexDefines.H"

#include <algorithm>
#include <sstream>

namespace erf_sbm {

CapabilityReport evaluate_p1_capabilities(const CapabilityInput& input)
{
    CapabilityReport report;
    report.flags = {"single_level", "static_cartesian", "periodic_manufactured",
                    "runtime_bins", "first_order_donor", "qv_qc_qr", "double"};
    report.invariant_ids = {"SBM-AUX-NONNEGATIVE", "SBM-BULK-PROJECTION",
                            "SBM-ACCEPTED-TRANSFER", "SBM-OLD-BASELINE",
                            "SBM-ONE-VAPOR"};
    auto reject = [&report](const bool condition, const char* reason) {
        if (condition) report.rejected_reasons.emplace_back(reason);
    };
    reject(input.max_level > 0, "AMR levels greater than zero are unsupported in P1");
    reject(input.diffusion, "auxiliary/projection diffusion is unsupported in P1");
    reject(input.implicit_moisture_diffusion, "implicit moisture diffusion is unsupported in P1");
    reject(input.shoc_or_macrophysics, "SHOC/macrophysics is unsupported in P1");
    reject(!input.static_cartesian, "non-Cartesian geometry is unsupported in P1");
    reject(input.moving_terrain, "moving terrain is unsupported in P1");
    reject(input.embedded_boundary, "embedded boundaries are unsupported in P1");
    reject(input.high_order_or_fct, "high-order/FCT transport is unsupported in P1");
    reject(input.sedimentation, "sedimentation is unsupported in P1");
    reject(input.condensation, "condensation/evaporation is unsupported in P1");
    reject(input.activation, "activation/regeneration is unsupported in P1");
    reject(input.collision, "collision/coalescence is unsupported in P1");
    reject(input.dynamic_grid, "dynamic spectral grids are unsupported in P1");
    reject(input.restart_schema_conversion, "restart schema conversion is unsupported in P1");
    reject(input.two_moment_transport, "two-moment transport is a P0 contract and is unsupported in P1");
    reject(input.custom_moisture_forcing, "custom moisture forcing is unsupported in P1");
    reject(input.large_scale_forcing, "large-scale forcing is unsupported in P1");
    reject(input.sounding_nudging, "sounding nudging is unsupported in P1");
    reject(input.sponge_or_wall_modification, "sponge and wall modification are unsupported in P1");
    reject(!input.periodic_cartesian, "only static Cartesian periodic manufactured cases are supported in P1");
    reject(!input.double_precision, "P1 manufactured transport is currently double precision only");
    report.supported = report.rejected_reasons.empty();
    return report;
}

CapabilityReport evaluate_p2_capabilities(const CapabilityInput& input)
{
    CapabilityReport report;
    report.flags = {"multi_level", "static_cartesian", "runtime_bins", "complete_groups",
                    "two_moment_endpoints", "explicit_density_weighted_diffusion",
                    "weno_z3_fct", "conservative_amr", "strict_restart_schema", "double"};
    report.invariant_ids = {"SBM-P2-GROUP-COMPLETE", "SBM-P2-ENDPOINT-REALIZABLE",
                            "SBM-P2-PHYSICAL-TRANSFER", "SBM-P2-BOUNDARY-BUDGET",
                            "SBM-P2-AMR-CONSERVATIVE", "SBM-P2-RESTART-STRICT"};
    auto reject = [&report](const bool condition, const char* reason) {
        if (condition) report.rejected_reasons.emplace_back(reason);
    };
    reject(!input.p2_requested, "P2 capability evaluation requires a P2 transport request");
    reject(input.shoc_or_macrophysics, "SHOC/macrophysics is unsupported by the P2 SBM hook");
    reject(!input.static_cartesian, "P2 auxiliary transfer currently requires static Cartesian geometry");
    reject(input.moving_terrain, "moving terrain is unsupported by P2 auxiliary transfer");
    reject(input.embedded_boundary, "embedded boundaries are unsupported by P2 auxiliary transfer");
    reject(input.sedimentation, "sedimentation is reserved for P3+");
    reject(input.condensation, "condensation/evaporation is reserved for P3+");
    reject(input.activation, "activation/regeneration is reserved for P3+");
    reject(input.collision, "collision/coalescence is reserved for P3+");
    reject(input.dynamic_grid, "dynamic spectral grids are unsupported by P2");
    reject(input.restart_schema_conversion, "restart schema conversion is unsupported by P2");
    reject(input.custom_moisture_forcing, "custom moisture forcing is unsupported by P2");
    reject(input.large_scale_forcing, "large-scale forcing is unsupported by P2");
    reject(input.sounding_nudging, "sounding nudging is unsupported by P2");
    reject(input.sponge_or_wall_modification, "sponge/wall modification is unsupported by P2");
    reject(input.diffusion && !input.explicit_sbm_diffusion,
           "native moisture diffusion must not write provider-owned SBM components");
    reject(input.implicit_moisture_diffusion, "implicit moisture diffusion is unsupported by P2");
    reject(!input.periodic_cartesian, "P2 production transport currently requires periodic Cartesian boundaries");
    reject(!input.double_precision, "P2 manufactured transport is currently double precision only");
    reject(input.chunk_size <= 0, "P2 scratch chunk size must be positive");
    report.supported = report.rejected_reasons.empty();
    return report;
}

std::string validate_runtime_bin_count(const int nbins)
{
    if (nbins < 2 || nbins > 1000000) {
        return "SBM sbm_nbins must be in [2, 1000000]";
    }
    return {};
}

std::string CapabilityReport::stable_description() const
{
    std::ostringstream out;
    out << "supported=" << (supported ? 1 : 0) << "\nflags=";
    for (const auto& flag : flags) out << flag << ',';
    out << "\ninvariants=";
    for (const auto& invariant : invariant_ids) out << invariant << ',';
    out << "\nrejected=";
    for (const auto& reason : rejected_reasons) out << reason << '|';
    return out.str();
}

std::string stable_inspection(const SBMLayout& layout,
                              const CapabilityInput& input,
                              const std::string& baseline_identity,
                              const std::string& amrex_identity,
                              const std::string& design_identity)
{
    const auto report = input.p2_requested ?
        evaluate_p2_capabilities(input) : evaluate_p1_capabilities(input);
    std::ostringstream out;
    out << "format=erf-sbm-inspection-v1\n";
    out << layout.inspection();
    out << "capabilities\n" << report.stable_description() << "\n";
    out << "identity.baseline=" << baseline_identity << "\n";
    out << "identity.amrex=" << amrex_identity << "\n";
    out << "identity.design=" << design_identity << "\n";
    return out.str();
}

bool OwnershipRegistry::owns_cloud_or_rain(const int component) const noexcept
{
    return m_provider_active && (component == RhoQ2_comp || component == RhoQ3_comp);
}

bool OwnershipRegistry::owns(const int component, const NativeWritePath path) const noexcept
{
    // The registry intentionally has one provider-owned result for all native
    // write classes.  qv (RhoQ1_comp) is never owned by SBM.
    (void) path;
    return owns_cloud_or_rain(component);
}

} // namespace erf_sbm
