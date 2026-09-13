#include "ERF_SBMLayout.H"

#include <algorithm>
#include <sstream>
#include <stdexcept>

namespace erf_sbm {

LayoutValidation SBMLayout::validate(const SBMLayoutSpec& spec)
{
    if (spec.populations.empty()) return {false, "at least one SBM population is required"};
    if (spec.populations.size() != spec.moment_modes.size()) return {false, "every population needs a moment mode"};
    std::vector<int> ids;
    for (std::size_t i = 0; i < spec.populations.size(); ++i) {
        const auto grid_result = SpectralGrid::validate(spec.populations[i]);
        if (!grid_result.valid) return {false, grid_result.message};
        if (std::find(ids.begin(), ids.end(), spec.populations[i].population_id) != ids.end()) {
            return {false, "population ids must be unique"};
        }
        ids.push_back(spec.populations[i].population_id);
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
        if (!property.transported) {
            return {false, "P1 attached properties must be transported when registered"};
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
        PopulationLayout population{spec.populations[i].population_id,
                                    SpectralGrid(spec.populations[i]), spec.moment_modes[i], offset,
                                    -1, 0};
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
        auto p = property;
        p.remap_with_mass = (p.kind == PropertyKind::MassBoundedSubset);
        m_properties.push_back(std::move(p));
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
    const auto& liquid = m_populations.front();
    rules.push_back({"qc", liquid.liquid_mass_offset, liquid.grid.cloud_bin_count()});
    rules.push_back({"qr", liquid.liquid_mass_offset + liquid.grid.cloud_bin_count(),
                     liquid.grid.nbins() - liquid.grid.cloud_bin_count()});
    m_projection = ::erf_auxiliary::AuxiliaryProjection(std::move(rules));
    if (!m_projection.validate(m_ncomp).valid) {
        throw std::invalid_argument("SBM projection does not match layout");
    }

    std::ostringstream schema;
    schema << "sbm-layout-v1|ncomp=" << m_ncomp << '|';
    for (const auto& p : m_populations) {
        schema << "population=" << p.population_id << ":" << p.grid.identity()
               << ":moment=" << static_cast<int>(p.moment_mode)
               << ":mass=" << p.liquid_mass_offset << ":number=" << p.number_offset << '|';
    }
    for (std::size_t i = 0; i < m_properties.size(); ++i) {
        const auto& p = m_properties[i];
        schema << "property=" << p.name << ':' << p.semantic_id << ':' << p.units
               << ':' << p.carrier_population << ':' << static_cast<int>(p.kind)
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

int SBMLayout::liquid_mass_offset(const int population) const
{
    for (const auto& p : m_populations) if (p.population_id == population) return p.liquid_mass_offset;
    throw std::out_of_range("unknown SBM population");
}

std::string SBMLayout::inspection() const
{
    std::ostringstream out;
    out << "schema=" << m_schema_identity << "\ncomponents=" << m_ncomp << "\n";
    for (const auto& p : m_populations) {
        out << "population " << p.population_id << " bins=" << p.grid.nbins()
            << " moment_mode=" << static_cast<int>(p.moment_mode)
            << " mass_offset=" << p.liquid_mass_offset
            << " number_offset=" << p.number_offset << "\n";
    }
    out << "projection qc=population0[0:" << m_populations.front().grid.cloud_bin_count()
        << "] qr=population0[" << m_populations.front().grid.cloud_bin_count() << ':'
        << m_populations.front().grid.nbins() << "]\n";
    for (std::size_t i = 0; i < m_properties.size(); ++i) {
        out << "property " << m_properties[i].name
            << " kind=" << static_cast<int>(m_properties[i].kind)
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
            components.push_back({"population" + std::to_string(p.population_id) + ".liquid_mass." + std::to_string(b),
                                  "liquid_mass", p.grid.units()});
        }
        if (p.number_offset >= 0) {
            for (int b = 0; b < p.grid.nbins(); ++b) {
                components.push_back({"population" + std::to_string(p.population_id) + ".number." + std::to_string(b),
                                      "number", "m^-3"});
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
