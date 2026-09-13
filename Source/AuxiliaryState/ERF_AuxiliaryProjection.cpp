#include "ERF_AuxiliaryProjection.H"

#include <set>

namespace erf_auxiliary {

ProjectionValidation AuxiliaryProjection::validate(const int source_components) const
{
    if (source_components < 0) return {false, "negative auxiliary component count"};
    std::set<std::string> targets;
    std::vector<int> seen(static_cast<std::size_t>(source_components), 0);
    for (const auto& rule : m_rules) {
        if (rule.target.empty()) return {false, "projection target is empty"};
        if (rule.target_kind != ProjectionTargetKind::BulkCoupling) {
            return {false, "non-bulk target cannot be used by this coupling projection"};
        }
        if (rule.source_kind != ProjectionSourceKind::AuxiliaryMass) {
            return {false, "only auxiliary mass components may enter this coupling projection"};
        }
        if (!targets.insert(rule.target).second) return {false, "duplicate projection target"};
        if (rule.source_begin < 0 || rule.source_count <= 0 ||
            rule.source_begin + rule.source_count > source_components) {
            return {false, "projection source range is outside the auxiliary state"};
        }
        for (int c = rule.source_begin; c < rule.source_begin + rule.source_count; ++c) {
            if (seen[static_cast<std::size_t>(c)] != 0) return {false, "overlapping projection source ranges"};
            seen[static_cast<std::size_t>(c)] = 1;
        }
    }
    return {true, {}};
}

} // namespace erf_auxiliary
