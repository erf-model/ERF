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
    reject(!input.periodic_cartesian, "only static Cartesian periodic manufactured cases are supported in P1");
    reject(!input.double_precision, "P1 manufactured transport is currently double precision only");
    report.supported = report.rejected_reasons.empty();
    return report;
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
    const auto report = evaluate_p1_capabilities(input);
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
