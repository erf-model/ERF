#ifndef ERF_SBM_LAYOUT_H_
#define ERF_SBM_LAYOUT_H_

#include "ERF_SpectralGrid.H"
#include "ERF_AuxiliaryStateLayout.H"
#include "ERF_AuxiliaryProjection.H"

#include <string>
#include <limits>
#include <vector>

namespace erf_sbm {

/**
 * Semantic class of an attached carrier-bin extensive quantity.
 *
 * ExtensiveMass is nonnegative additional mass and is not part of the water
 * projection. NumberCarried is a nonnegative count-like quantity. A
 * MassBoundedSubset is nonnegative and may be constrained by its carrier mass
 * when that future physics is implemented.
 */
enum class PropertyKind { ExtensiveMass, NumberCarried, MassBoundedSubset };
/** P1/P2 remap contract; independent remapping is reserved for future work. */
enum class PropertyRemapPolicy { CarrierBinConservative, Independent };
enum class SupportRequirement { None, PositiveMass, PositiveNumber };
enum class PopulationPhase { Liquid, Ice, Aerosol, Generic };

struct SpectralPopulationSpec {
    int population_id{0};
    std::string semantic_id{"liquid_mass"};
    PopulationPhase phase{PopulationPhase::Liquid};
    SpectralGridSpec grid;
    MomentMode moment_mode{MomentMode::OneMoment};
    // These are units of the transported prognostic components, not of the
    // particle-coordinate grid.
    std::string mass_state_units{"kg m^-3"};
    std::string number_state_units{"m^-3"};
};

struct LiquidProjectionSpec {
    int population_id{0};
    int cloud_rain_split{-1};
};

struct AttachedPropertyDescriptor {
    std::string name;
    std::string semantic_id;
    std::string units;
    int carrier_population{0};
    PropertyKind kind{PropertyKind::NumberCarried};
    SupportRequirement support{SupportRequirement::None};
    // All P1/P2 property kinds are aligned to carrier bins and conservatively
    // remapped with the carrier. The kind, not this policy, distinguishes
    // additional mass from count-like or mass-bounded semantics.
    PropertyRemapPolicy remap_policy{PropertyRemapPolicy::CarrierBinConservative};
    bool transported{true};
    bool may_overlap_mass{false};
    // Optional finite per-carrier support bounds.  NaN means that the
    // descriptor only requests nonnegativity; no arbitrary carrier floor is
    // manufactured by the generic transport layer.
    amrex::Real support_min{0.0};
    amrex::Real support_max{std::numeric_limits<amrex::Real>::quiet_NaN()};
};

struct PopulationLayout {
    int population_id{0};
    std::string semantic_id;
    PopulationPhase phase{PopulationPhase::Generic};
    SpectralGrid grid;
    MomentMode moment_mode{MomentMode::OneMoment};
    std::string mass_state_units;
    std::string number_state_units;
    int mass_offset{0};
    int number_offset{-1};
    int component_count{0};
};

struct SBMLayoutSpec {
    std::vector<SpectralPopulationSpec> populations;
    LiquidProjectionSpec liquid_projection;
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
    [[nodiscard]] const LiquidProjectionSpec& liquid_projection() const noexcept { return m_liquid_projection; }
    [[nodiscard]] int ncomp() const noexcept { return m_ncomp; }
    [[nodiscard]] int mass_offset(const int population) const;
    [[nodiscard]] int property_offset(const int property) const;
    [[nodiscard]] std::string schema_identity() const noexcept { return m_schema_identity; }
    [[nodiscard]] std::string inspection() const;
    [[nodiscard]] ::erf_auxiliary::AuxiliaryStateLayout auxiliary_layout() const;

private:
    std::vector<PopulationLayout> m_populations;
    std::vector<AttachedPropertyDescriptor> m_properties;
    std::vector<int> m_property_offsets;
    LiquidProjectionSpec m_liquid_projection;
    ::erf_auxiliary::AuxiliaryProjection m_projection;
    int m_ncomp{0};
    std::string m_schema_identity;
};

} // namespace erf_sbm

#endif
