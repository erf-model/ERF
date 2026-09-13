#ifndef ERF_SBM_LAYOUT_H_
#define ERF_SBM_LAYOUT_H_

#include "ERF_SpectralGrid.H"
#include "ERF_AuxiliaryStateLayout.H"
#include "ERF_AuxiliaryProjection.H"

#include <string>
#include <vector>

namespace erf_sbm {

enum class PropertyKind { NumberCarried, MassBoundedSubset };
enum class SupportRequirement { None, PositiveMass, PositiveNumber };

struct AttachedPropertyDescriptor {
    std::string name;
    std::string semantic_id;
    std::string units;
    int carrier_population{0};
    PropertyKind kind{PropertyKind::NumberCarried};
    SupportRequirement support{SupportRequirement::None};
    bool remap_with_mass{true};
    bool transported{true};
    bool may_overlap_mass{false};
};

struct PopulationLayout {
    int population_id{0};
    SpectralGrid grid;
    MomentMode moment_mode{MomentMode::OneMoment};
    int liquid_mass_offset{0};
    int number_offset{-1};
    int component_count{0};
};

struct SBMLayoutSpec {
    std::vector<SpectralGridSpec> populations;
    std::vector<MomentMode> moment_modes;
    std::vector<AttachedPropertyDescriptor> attached_properties;
};

struct LayoutValidation {
    bool valid{false};
    std::string message;
};

/** Immutable runtime SBM component and projection layout. */
class SBMLayout {
public:
    static LayoutValidation validate(const SBMLayoutSpec& spec);
    explicit SBMLayout(SBMLayoutSpec spec);

    [[nodiscard]] const std::vector<PopulationLayout>& populations() const noexcept { return m_populations; }
    [[nodiscard]] const std::vector<AttachedPropertyDescriptor>& attached_properties() const noexcept { return m_properties; }
    [[nodiscard]] const ::erf_auxiliary::AuxiliaryProjection& bulk_projection() const noexcept { return m_projection; }
    [[nodiscard]] int ncomp() const noexcept { return m_ncomp; }
    [[nodiscard]] int liquid_mass_offset(const int population) const;
    [[nodiscard]] int property_offset(const int property) const;
    [[nodiscard]] std::string schema_identity() const noexcept { return m_schema_identity; }
    [[nodiscard]] std::string inspection() const;
    [[nodiscard]] ::erf_auxiliary::AuxiliaryStateLayout auxiliary_layout() const;

private:
    std::vector<PopulationLayout> m_populations;
    std::vector<AttachedPropertyDescriptor> m_properties;
    std::vector<int> m_property_offsets;
    ::erf_auxiliary::AuxiliaryProjection m_projection;
    int m_ncomp{0};
    std::string m_schema_identity;
};

} // namespace erf_sbm

#endif
