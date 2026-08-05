/**
 * \file ERF_Plotfile2DSampledField.cpp
 */
#include "ERF_Plotfile2DSampledField.H"

#include <unordered_set>

namespace plotfile2d
{

namespace
{

const amrex::Vector<SampledFieldDescriptor>& field_catalog_storage ()
{
    static const amrex::Vector<SampledFieldDescriptor> catalog{
        {SampledFieldID::Rho,       "rho",       "Density sampled at the target level",                  "kg/m^3"},
        {SampledFieldID::Theta,     "theta",     "Potential temperature sampled at the target level",    "K"},
        {SampledFieldID::Temp,      "temp",      "Temperature sampled at the target level",              "K"},
        {SampledFieldID::Pressure,  "pressure",  "Pressure sampled at the target level",                 "Pa"},
        {SampledFieldID::HeightMSL, "height_msl","Height above mean sea level sampled at the target level","m"},
        {SampledFieldID::HeightAGL, "height_agl","Height above ground level sampled at the target level", "m"},
        {SampledFieldID::Qv,        "qv",        "Water vapor mixing ratio sampled at the target level", "kg/kg"},
        {SampledFieldID::Qc,        "qc",        "Cloud liquid water mixing ratio sampled at the target level", "kg/kg"},
        {SampledFieldID::Qi,        "qi",        "Cloud ice mixing ratio sampled at the target level",    "kg/kg"},
        {SampledFieldID::Qr,        "qr",        "Rain mixing ratio sampled at the target level",         "kg/kg"},
        {SampledFieldID::Qs,        "qs",        "Snow mixing ratio sampled at the target level",         "kg/kg"},
        {SampledFieldID::Qg,        "qg",        "Graupel mixing ratio sampled at the target level",      "kg/kg"},
        {SampledFieldID::UEast,     "u_east",    "Eastward wind sampled at the target level",             "m/s"},
        {SampledFieldID::VNorth,    "v_north",   "Northward wind sampled at the target level",            "m/s"},
        {SampledFieldID::W,         "w",         "Vertical wind sampled at the target level",             "m/s"},
        {SampledFieldID::WindSpeed,  "wind_speed","Horizontal wind speed sampled at the target level",     "m/s"},
        {SampledFieldID::WindDir,    "wind_dir",  "Meteorological wind direction sampled at the target level", "degrees"},
    };

    return catalog;
}

bool field_is_available (SampledFieldID field_id, const SolverChoice& solver_choice) noexcept
{
    const auto& mi = solver_choice.moisture_indices;
    switch (field_id) {
    case SampledFieldID::Rho:
    case SampledFieldID::Theta:
    case SampledFieldID::Temp:
    case SampledFieldID::Pressure:
    case SampledFieldID::HeightMSL:
    case SampledFieldID::HeightAGL:
        return true;
    case SampledFieldID::Qv:
        return mi.qv >= 0;
    case SampledFieldID::Qc:
        return mi.qc >= 0;
    case SampledFieldID::Qi:
        return mi.qi >= 0;
    case SampledFieldID::Qr:
        return mi.qr >= 0;
    case SampledFieldID::Qs:
        return mi.qs >= 0;
    case SampledFieldID::Qg:
        return mi.qg >= 0;
    case SampledFieldID::UEast:
    case SampledFieldID::VNorth:
    case SampledFieldID::W:
    case SampledFieldID::WindSpeed:
    case SampledFieldID::WindDir:
        return true;
    }

    return false;
}

} // namespace

const amrex::Vector<SampledFieldDescriptor>&
sampled_field_catalog ()
{
    return field_catalog_storage();
}

const SampledFieldDescriptor*
find_sampled_field (const std::string& name)
{
    for (const auto& descriptor : sampled_field_catalog()) {
        if (name == descriptor.name) {
            return &descriptor;
        }
    }

    return nullptr;
}

amrex::Vector<std::string>
available_sampled_field_names (const SolverChoice& solver_choice)
{
    amrex::Vector<std::string> names;
    names.reserve(sampled_field_catalog().size());

    for (const auto& descriptor : sampled_field_catalog()) {
        if (field_is_available(descriptor.id, solver_choice)) {
            names.push_back(descriptor.name);
        }
    }

    return names;
}

SampledFieldSelection
select_requested_sampled_fields (const amrex::Vector<std::string>& requested,
                                 const SolverChoice& solver_choice)
{
    SampledFieldSelection selection;

    std::unordered_set<std::string> seen;

    for (const auto& name : requested) {
        const auto* descriptor = find_sampled_field(name);
        if (descriptor == nullptr) {
            selection.unavailable.push_back(name);
            continue;
        }

        if (!field_is_available(descriptor->id, solver_choice)) {
            selection.unavailable.push_back(name);
            continue;
        }

        if (seen.insert(name).second) {
            selection.accepted.push_back(*descriptor);
        }
    }

    return selection;
}

std::string
sampled_field_id_to_string (SampledFieldID field_id)
{
    return sampled_field_name(field_id);
}

} // namespace plotfile2d
