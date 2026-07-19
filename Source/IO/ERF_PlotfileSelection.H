#ifndef ERF_PLOTFILE_SELECTION_H_
#define ERF_PLOTFILE_SELECTION_H_

#include <algorithm>
#include <string>

#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"

namespace erf_plotfile
{

struct Plot3DSelectionCapabilities
{
    int conserved_state_size = 0;
    int moist_state_size = 0;
    int moist_numconc_size = 0;
    int qmoist_size = 0;
    MoistureType moisture_type = MoistureType::None;
};

inline int
plot3d_conserved_component_index (const std::string& name)
{
    if (name == "density")  return 0;
    if (name == "rhotheta") return 1;
    if (name == "rhoKE")    return 2;
    if (name == "rhoadv_0") return 3;

    if (name.size() > 4 && name.compare(0, 4, "rhoQ") == 0) {
        const std::string suffix = name.substr(4);
        if (!suffix.empty() && std::all_of(suffix.cbegin(), suffix.cend(),
                                           [] (const char c) { return c >= '0' && c <= '9'; })) {
            const int q = std::stoi(suffix);
            if (q >= 1 && q <= 11) return 3 + q;
        }
    }

    return -1;
}

inline bool
plot3d_has_moist_state (const Plot3DSelectionCapabilities& caps)
{
    return caps.moisture_type != MoistureType::None && caps.moist_state_size > 0;
}

inline bool
plot3d_fixed_variable_available (const std::string& name,
                                 const Plot3DSelectionCapabilities& caps)
{
    const int cons_comp = plot3d_conserved_component_index(name);
    if (cons_comp >= 0) return cons_comp < caps.conserved_state_size;

    const bool has_moist_state = plot3d_has_moist_state(caps);
    const int nq = caps.moist_state_size;

    if (name == "moist_density") return has_moist_state && nq >= 2;
    if (name == "qv")           return has_moist_state && nq >= 1;
    if (name == "qc")           return has_moist_state && nq >= 2;
    if (name == "qrain")        return has_moist_state && nq >= 3;
    if (name == "qi")           return has_moist_state && nq >= 4;
    if (name == "qsnow")        return has_moist_state && nq >= 5;
    if (name == "qgraup")       return has_moist_state && nq >= 6;
    if (name == "nc")           return has_moist_state && nq >= 7;
    if (name == "ni")           return has_moist_state && nq >= 8;
    if (name == "nr")           return has_moist_state && nq >= 9;
    if (name == "ns")           return has_moist_state && nq >= 10;
    if (name == "ng")           return has_moist_state && nq >= 11;
    if (name == "qt")           return has_moist_state && nq >= 1;
    if (name == "qn")           return has_moist_state && nq >= 2;
    if (name == "qp")           return has_moist_state && nq >= 3;
    if (name == "qsat")         return has_moist_state && nq >= 1;
    if (name == "precipitable") return has_moist_state;

    if (name == "rain_accum") {
        const bool supported_model =
            caps.moisture_type == MoistureType::Kessler ||
            caps.moisture_type == MoistureType::Morrison_NoIce ||
            caps.moisture_type == MoistureType::SAM_NoIce ||
            caps.moisture_type == MoistureType::SAM ||
            caps.moisture_type == MoistureType::Morrison ||
            caps.moisture_type == MoistureType::WSM6 ||
            caps.moisture_type == MoistureType::SuperDroplets;
        return supported_model &&
            ((caps.moisture_type == MoistureType::SuperDroplets && caps.qmoist_size >= 7) ||
             (caps.moisture_type != MoistureType::SuperDroplets && caps.qmoist_size >= 1));
    }

    if (name == "snow_accum" || name == "graup_accum") {
        const bool supported_model =
            caps.moisture_type == MoistureType::SAM ||
            caps.moisture_type == MoistureType::Morrison ||
            caps.moisture_type == MoistureType::WSM6;
        const int required_qmoist = name == "snow_accum" ? 2 : 3;
        return supported_model && caps.qmoist_size >= required_qmoist;
    }

    if (name == "rel_humidity") {
        return caps.moisture_type == MoistureType::SuperDroplets && caps.qmoist_size >= 6;
    }
    if (name == "condensation_rate") {
        return caps.moisture_type == MoistureType::SuperDroplets && caps.qmoist_size >= 4;
    }

    // The writer deliberately supplies a zero fallback for unsupported
    // reflectivity schemes, so these names remain valid fixed diagnostics.
    return true;
}

inline bool
plot3d_needs_pressure (const amrex::Vector<std::string>& names)
{
    for (const auto& name : names) {
        if (name == "pressure" || name == "pert_pres" ||
            name == "dpdx" || name == "dpdy" || name == "dpdz" ||
            name == "eq_pot_temp" || name == "VPD" || name == "qsat") {
            return true;
        }
#ifdef ERF_COMPUTE_ERROR
        if (name == "pp_err") return true;
#endif
    }
    return false;
}

inline amrex::Vector<std::string>
plot3d_selected_particle_count_names (const amrex::Vector<std::string>& requested,
                                      const amrex::Vector<std::string>& configured)
{
    amrex::Vector<std::string> selected;
    for (const auto& count_name : requested) {
        constexpr std::size_t suffix_size = 6; // "_count"
        if (count_name.size() > suffix_size &&
            count_name.compare(count_name.size() - suffix_size, suffix_size, "_count") == 0) {
            const std::string container = count_name.substr(0, count_name.size() - suffix_size);
            if (std::find(configured.cbegin(), configured.cend(), container) == configured.cend() ||
                std::find(selected.cbegin(), selected.cend(), container) != selected.cend()) {
                continue;
            }
            selected.push_back(container);
        }
    }
    return selected;
}

} // namespace erf_plotfile

#endif
