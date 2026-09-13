#ifndef ERF_AUXILIARY_PROJECTION_H_
#define ERF_AUXILIARY_PROJECTION_H_

#include <string>
#include <utility>
#include <vector>

namespace erf_auxiliary {

enum class ProjectionTargetKind { BulkCoupling, Vapor, Number, AttachedMaterial };
enum class ProjectionSourceKind { AuxiliaryMass, AuxiliaryNumber, AttachedMaterial };

/** A fixed linear range projection from auxiliary components to a coupling lane. */
struct ProjectionRule {
    std::string target;
    int source_begin{0};
    int source_count{0};
    ProjectionTargetKind target_kind{ProjectionTargetKind::BulkCoupling};
    ProjectionSourceKind source_kind{ProjectionSourceKind::AuxiliaryMass};
};

struct ProjectionValidation {
    bool valid{false};
    std::string message;
};

/** Generic, fixed-range projection metadata. */
class AuxiliaryProjection {
public:
    AuxiliaryProjection() = default;
    explicit AuxiliaryProjection(std::vector<ProjectionRule> rules)
        : m_rules(std::move(rules)) {}

    [[nodiscard]] const std::vector<ProjectionRule>& rules() const noexcept { return m_rules; }
    [[nodiscard]] ProjectionValidation validate(const int source_components) const;

private:
    std::vector<ProjectionRule> m_rules;
};

} // namespace erf_auxiliary

#endif
