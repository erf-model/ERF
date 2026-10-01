#include "ERF_SBMLayout.H"
#include "ERF_SBMCanonicalIdentity.H"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

LayoutValidation SBMLayout::validate(const SBMLayoutSpec& spec)
{
    if (spec.populations.empty()) return {false, "at least one SBM population is required"};
    std::vector<int> ids;
    for (std::size_t i = 0; i < spec.populations.size(); ++i) {
        const auto& population = spec.populations[i];
        const auto grid_result = SpectralGrid::validate(population.grid);
        if (!grid_result.valid) return {false, grid_result.message};
        if (population.population_id < 0 || population.semantic_id.empty() ||
            population.mass_state_units.empty() || population.number_state_units.empty()) {
            return {false, "population id, semantic id, and state units are required"};
        }
        if (std::find(ids.begin(), ids.end(), population.population_id) != ids.end()) {
            return {false, "population ids must be unique"};
        }
        ids.push_back(population.population_id);
    }
    const auto liquid = std::find_if(spec.populations.begin(), spec.populations.end(),
        [&](const SpectralPopulationSpec& population) {
            return population.population_id == spec.liquid_projection.population_id;
        });
    if (liquid == spec.populations.end()) return {false, "liquid projection refers to an unknown population"};
    if (liquid->phase != PopulationPhase::Liquid) {
        return {false, "bulk cloud/rain projection must refer to a liquid population"};
    }
    if (spec.liquid_projection.cloud_rain_split <= 0 ||
        spec.liquid_projection.cloud_rain_split >= static_cast<int>(liquid->grid.edges.size()) - 1) {
        return {false, "cloud/rain split must be an interior liquid-population bin index"};
    }
    for (const auto& property : spec.attached_properties) {
        if (property.name.empty() || property.semantic_id.empty() || property.units.empty()) {
            return {false, "attached properties need name, semantic id, and units"};
        }
        if (std::find(ids.begin(), ids.end(), property.carrier_population) == ids.end()) {
            return {false, "attached property refers to an unknown population"};
        }
        if (property.may_overlap_mass && property.kind == PropertyKind::MassBoundedSubset) {
            return {false, "mass-bounded subset cannot overlap liquid mass"};
        }
        if (property.remap_policy != PropertyRemapPolicy::CarrierBinConservative) {
            return {false, "P1 attached properties require carrier-bin conservative remapping"};
        }
        if (!property.transported) {
            return {false, "P1 attached properties must be transported when registered"};
        }
        if (!std::isfinite(property.support_min) || property.support_min < amrex::Real(0.0)) {
            return {false, "attached-property support_min must be finite and nonnegative"};
        }
        if (!std::isnan(property.support_max) &&
            (!std::isfinite(property.support_max) || property.support_max < property.support_min)) {
            return {false, "attached-property support_max must be NaN or finite and >= support_min"};
        }
    }
    return {true, {}};
}

SBMLayout::SBMLayout(SBMLayoutSpec spec)
{
    const auto result = validate(spec);
    if (!result.valid) throw std::invalid_argument("invalid SBM layout: " + result.message);

    // Input order is part of the deterministic schema.  Offsets are assigned
    // once and never changed after construction.
    int offset = 0;
    for (std::size_t i = 0; i < spec.populations.size(); ++i) {
        const auto& input = spec.populations[i];
        PopulationLayout population{input.population_id, input.semantic_id, input.phase,
                                    SpectralGrid(input.grid), input.moment_mode,
                                    input.mass_state_units, input.number_state_units,
                                    offset, -1, 0};
        population.component_count = population.grid.nbins();
        offset += population.component_count;
        if (population.moment_mode == MomentMode::TwoMoment) {
            population.number_offset = offset;
            population.component_count += population.grid.nbins();
            offset += population.grid.nbins();
        }
        m_populations.push_back(std::move(population));
    }
    for (const auto& property : spec.attached_properties) {
        m_properties.push_back(property);
        const auto& carrier = std::find_if(m_populations.begin(), m_populations.end(),
            [&](const PopulationLayout& candidate) {
                return candidate.population_id == m_properties.back().carrier_population;
            });
        if (carrier == m_populations.end()) {
            throw std::invalid_argument("attached property refers to an unknown population");
        }
        // An attached extensive property is represented per carrier bin.  A
        // number-carried property is not folded into the number block and a
        // mass-bounded subset is not folded into liquid mass; both remain
        // typed, independently addressable transport components.
        m_property_offsets.push_back(offset);
        offset += carrier->grid.nbins();
    }
    m_ncomp = offset;

    std::vector<::erf_auxiliary::ProjectionRule> rules;
    m_liquid_projection = spec.liquid_projection;
    const auto liquid = std::find_if(m_populations.begin(), m_populations.end(),
        [&](const PopulationLayout& population) {
            return population.population_id == m_liquid_projection.population_id;
        });
    const int split = m_liquid_projection.cloud_rain_split;
    rules.push_back({"qc", liquid->mass_offset, split});
    rules.push_back({"qr", liquid->mass_offset + split, liquid->grid.nbins() - split});
    m_projection = ::erf_auxiliary::AuxiliaryProjection(std::move(rules));
    if (!m_projection.validate(m_ncomp).valid) {
        throw std::invalid_argument("SBM projection does not match layout");
    }

    std::ostringstream schema;
    schema << "sbm-layout-m1-v1|ncomp=" << m_ncomp
           << "|constraint_policy=nonnegative-bin-mass-and-moments-v1"
           << "|representation=bin-mass-density-v1"
           << "|projection=liquid-mass-sum-to-qc-qr-v1"
           << "|projection_population=" << m_liquid_projection.population_id
           << "|projection_cloud_rain_split=" << m_liquid_projection.cloud_rain_split << '|';
    for (const auto& p : m_populations) {
        schema << "population=" << p.population_id << ':' << p.semantic_id
               << ":phase=" << static_cast<int>(p.phase) << ':' << p.grid.identity()
               << ":moment=" << static_cast<int>(p.moment_mode)
               << ":mass_units=" << p.mass_state_units
               << ":number_units=" << p.number_state_units
               << ":mass=" << p.mass_offset << ":number=" << p.number_offset << '|';
    }
    for (std::size_t i = 0; i < m_properties.size(); ++i) {
        const auto& p = m_properties[i];
        schema << "property=" << p.name << ':' << p.semantic_id << ':' << p.units
               << ':' << p.carrier_population << ':' << static_cast<int>(p.kind)
               << ":remap=" << static_cast<int>(p.remap_policy)
               << ":support=" << static_cast<int>(p.support)
               << ":transported=" << (p.transported ? 1 : 0)
               << ":may_overlap_mass=" << (p.may_overlap_mass ? 1 : 0)
               << ":support_min=" << canonical_real(p.support_min)
               << ":support_max=" << canonical_real(p.support_max, true)
               << ":offset=" << m_property_offsets[i] << '|';
    }
    m_schema_identity = schema.str();
}

int SBMLayout::property_offset(const int property) const
{
    if (property < 0 || property >= static_cast<int>(m_property_offsets.size())) {
        throw std::out_of_range("unknown SBM attached property");
    }
    return m_property_offsets[static_cast<std::size_t>(property)];
}

int SBMLayout::mass_offset(const int population) const
{
    for (const auto& p : m_populations) if (p.population_id == population) return p.mass_offset;
    throw std::out_of_range("unknown SBM population");
}

std::string SBMLayout::inspection() const
{
    std::ostringstream out;
    out << "schema=" << m_schema_identity << "\ncomponents=" << m_ncomp << "\n";
    for (const auto& p : m_populations) {
        out << "population " << p.population_id << " bins=" << p.grid.nbins()
            << " semantic_id=" << p.semantic_id
            << " moment_mode=" << static_cast<int>(p.moment_mode)
            << " coordinate_units=" << p.grid.coordinate_units()
            << " mass_state_units=" << p.mass_state_units
            << " number_state_units=" << p.number_state_units
            << " mass_offset=" << p.mass_offset
            << " number_offset=" << p.number_offset << "\n";
    }
    const auto liquid = std::find_if(m_populations.begin(), m_populations.end(),
        [&](const PopulationLayout& population) {
            return population.population_id == m_liquid_projection.population_id;
        });
    out << "projection qc=population" << liquid->population_id << "[0:" << m_liquid_projection.cloud_rain_split
        << "] qr=population" << liquid->population_id << "[" << m_liquid_projection.cloud_rain_split << ':'
        << liquid->grid.nbins() << "]\n";
    for (std::size_t i = 0; i < m_properties.size(); ++i) {
        out << "property " << m_properties[i].name
            << " kind=" << static_cast<int>(m_properties[i].kind)
            << " remap=" << static_cast<int>(m_properties[i].remap_policy)
            << " carrier=" << m_properties[i].carrier_population
            << " offset=" << m_property_offsets[i] << '\n';
    }
    return out.str();
}

::erf_auxiliary::AuxiliaryStateLayout SBMLayout::auxiliary_layout() const
{
    std::vector<::erf_auxiliary::ComponentDescriptor> components;
    components.reserve(static_cast<std::size_t>(m_ncomp));
    for (const auto& p : m_populations) {
        for (int b = 0; b < p.grid.nbins(); ++b) {
            components.push_back({"population" + std::to_string(p.population_id) + ".mass." + std::to_string(b),
                                  p.semantic_id + ".mass", p.mass_state_units});
        }
        if (p.number_offset >= 0) {
            for (int b = 0; b < p.grid.nbins(); ++b) {
                components.push_back({"population" + std::to_string(p.population_id) + ".number." + std::to_string(b),
                                      p.semantic_id + ".number", p.number_state_units});
            }
        }
    }
    for (std::size_t i = 0; i < m_properties.size(); ++i) {
        const auto& property = m_properties[i];
        const auto& carrier = std::find_if(m_populations.begin(), m_populations.end(),
            [&](const PopulationLayout& candidate) {
                return candidate.population_id == property.carrier_population;
            });
        for (int b = 0; b < carrier->grid.nbins(); ++b) {
            components.push_back({property.name + "." + std::to_string(b),
                                  property.semantic_id, property.units});
        }
    }
    return {m_schema_identity, std::move(components)};
}

} // namespace erf_sbm
