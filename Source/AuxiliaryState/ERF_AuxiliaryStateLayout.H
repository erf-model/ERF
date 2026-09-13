#ifndef ERF_AUXILIARY_STATE_LAYOUT_H_
#define ERF_AUXILIARY_STATE_LAYOUT_H_

#include <string>
#include <utility>
#include <vector>

namespace erf_auxiliary {

/** A fixed component in an ERF-owned auxiliary prognostic state. */
struct ComponentDescriptor {
    std::string name;
    std::string semantic_id;
    std::string units;
};

/**
 * Generic metadata for an auxiliary state.  The generic layer deliberately
 * knows nothing about microphysics, bins, or water species.
 */
class AuxiliaryStateLayout {
public:
    AuxiliaryStateLayout() = default;
    AuxiliaryStateLayout(std::string schema_id,
                         std::vector<ComponentDescriptor> components)
        : m_schema_id(std::move(schema_id)), m_components(std::move(components)) {}

    [[nodiscard]] const std::string& schema_id() const noexcept { return m_schema_id; }
    [[nodiscard]] int ncomp() const noexcept { return static_cast<int>(m_components.size()); }
    [[nodiscard]] const std::vector<ComponentDescriptor>& components() const noexcept {
        return m_components;
    }

private:
    std::string m_schema_id;
    std::vector<ComponentDescriptor> m_components;
};

} // namespace erf_auxiliary

#endif
