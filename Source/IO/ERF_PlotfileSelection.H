#ifndef ERF_PLOTFILE_SELECTION_H_
#define ERF_PLOTFILE_SELECTION_H_

#include <algorithm>
#include <string>

#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"

namespace erf_plotfile
{

struct Plot3DMoistureCapabilities
{
    bool vapor = false;
    bool cloud_liquid = false;
    bool cloud_ice = false;
    bool rain = false;
    bool snow = false;
    bool graupel = false;

    bool cloud_number = false;
    bool ice_number = false;
    bool rain_number = false;
    bool snow_number = false;
    bool graupel_number = false;

    bool rain_accumulation = false;
    bool snow_accumulation = false;
    bool graupel_accumulation = false;
    bool relative_humidity_diagnostic = false;
    bool condensation_rate_diagnostic = false;
};

struct Plot3DSelectionCapabilities
{
    int conserved_state_size = 0;
    int moist_state_size = 0;
    int moist_numconc_size = 0;
    int qmoist_size = 0;
    MoistureType moisture_type = MoistureType::None;

    Plot3DMoistureCapabilities moisture;

    bool time_average_storage = false;
    bool radiation_heating_storage = false;
    bool eddy_diffusivity_storage = false;
    bool dissipation_storage = false;
    bool wall_distance_storage = false;
};

inline Plot3DMoistureCapabilities
plot3d_moisture_capabilities (const MoistureType moisture_type)
{
    Plot3DMoistureCapabilities caps;

    switch (moisture_type) {
    case MoistureType::None:
        break;
    case MoistureType::MoistNoCondensation:
        caps.vapor = true;
        break;
    case MoistureType::SatAdj:
        caps.vapor = true;
        caps.cloud_liquid = true;
        break;
    case MoistureType::Kessler_NoRain:
        caps.vapor = true;
        caps.cloud_liquid = true;
        break;
    case MoistureType::Kessler:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.rain = true;
        caps.rain_accumulation = true;
        break;
    case MoistureType::SAM_NoPrecip_NoIce:
        caps.vapor = true;
        caps.cloud_liquid = true;
        break;
    case MoistureType::SAM_NoIce:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.rain = true;
        caps.rain_accumulation = true;
        break;
    case MoistureType::SAM:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.cloud_ice = true;
        caps.rain = true;
        caps.snow = true;
        caps.graupel = true;
        caps.rain_accumulation = true;
        caps.snow_accumulation = true;
        caps.graupel_accumulation = true;
        break;
    case MoistureType::Morrison_NoIce:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.rain = true;
        caps.cloud_number = true;
        caps.rain_number = true;
        caps.rain_accumulation = true;
        break;
    case MoistureType::Morrison:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.cloud_ice = true;
        caps.rain = true;
        caps.snow = true;
        caps.graupel = true;
        caps.cloud_number = true;
        caps.ice_number = true;
        caps.rain_number = true;
        caps.snow_number = true;
        caps.graupel_number = true;
        caps.rain_accumulation = true;
        caps.snow_accumulation = true;
        caps.graupel_accumulation = true;
        break;
    case MoistureType::WSM6:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.cloud_ice = true;
        caps.rain = true;
        caps.snow = true;
        caps.graupel = true;
        caps.rain_accumulation = true;
        caps.snow_accumulation = true;
        caps.graupel_accumulation = true;
        break;
    case MoistureType::SuperDroplets:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.rain_accumulation = true;
        caps.relative_humidity_diagnostic = true;
        caps.condensation_rate_diagnostic = true;
        break;
    default:
        amrex::Abort("Unknown MoistureType in plot3d_moisture_capabilities");
    }

    return caps;
}

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
plot3d_source_component_available (const Plot3DSelectionCapabilities& caps,
                                    const int q)
{
    return plot3d_has_moist_state(caps) && q >= 1 && q <= caps.moist_state_size &&
        3 + q < caps.conserved_state_size;
}

inline bool
plot3d_fixed_unconditional_diagnostic (const std::string& name)
{
    static constexpr const char* names[] = {
        "temp", "theta", "KE", "scalar", "soundspeed", "reflectivity",
        "max_reflectivity", "mucape", "vorticity_x", "vorticity_y", "vorticity_z",
        "helicity", "local_helicity", "magvel", "divU", "pres_hse", "dens_hse",
        "theta_hse", "qv_hse", "pressure", "pert_pres", "pert_dens", "buoyancy",
        "eq_pot_temp", "VPD", "dpdx", "dpdy", "dpdz", "pres_hse_x", "pres_hse_y",
        "pres_hse_z", "z_phys", "detJ", "mapfac", "lat_m", "lon_m", "pblh",
        "shoc_cldfrac", "shoc_ql", "shoc_ql2", "shoc_cond", "wqls_sec", "wthv_sec",
        "w_sec", "thl_sec", "qw_sec", "qwthl_sec", "wthl_sec", "wqw_sec", "w3",
        "brunt", "isotropy", "shear_prod", "buoy_prod", "diss_tke",
        "terrain_IB_mask", "volfrac"
    };

    for (const char* candidate : names) {
        if (name == candidate) return true;
    }

#ifdef ERF_COMPUTE_ERROR
    if (name == "xvel_err" || name == "yvel_err" || name == "zvel_err" ||
        name == "pp_err") return true;
#endif

    return false;
}

inline bool
plot3d_fixed_variable_available (const std::string& name,
                                 const Plot3DSelectionCapabilities& caps)
{
    const int cons_comp = plot3d_conserved_component_index(name);
    if (cons_comp >= 0) return cons_comp < caps.conserved_state_size;

    if (name == "u_t_avg" || name == "v_t_avg" || name == "w_t_avg" ||
        name == "umag_t_avg") {
        return caps.time_average_storage;
    }
    if (name == "qsrc_sw" || name == "qsrc_lw") {
        return caps.radiation_heating_storage;
    }
    if (name == "nut" || name == "Kmv" || name == "Kmh" ||
        name == "Khv" || name == "Khh" || name == "Lturb") {
        return caps.eddy_diffusivity_storage;
    }
    if (name == "diss") return caps.dissipation_storage;
    if (name == "walldist") return caps.wall_distance_storage;

    const auto& moisture = caps.moisture;
    if (name == "moist_density") {
        return moisture.vapor && moisture.cloud_liquid &&
            plot3d_source_component_available(caps, 2);
    }
    if (name == "qv")    return moisture.vapor && plot3d_source_component_available(caps, 1);
    if (name == "qc")    return moisture.cloud_liquid && plot3d_source_component_available(caps, 2);
    if (name == "qi")    return moisture.cloud_ice && plot3d_source_component_available(caps, 4);
    if (name == "qrain") return moisture.rain && plot3d_source_component_available(caps, caps.moist_state_size > 3 ? 4 : 3);
    if (name == "qsnow") return moisture.snow && plot3d_source_component_available(caps, 5);
    if (name == "qgraup") return moisture.graupel && plot3d_source_component_available(caps, 6);
    if (name == "nc") return moisture.cloud_number && plot3d_source_component_available(caps, 7);
    if (name == "ni") return moisture.ice_number && plot3d_source_component_available(caps, 8);
    if (name == "nr") return moisture.rain_number && plot3d_source_component_available(caps, 9);
    if (name == "ns") return moisture.snow_number && plot3d_source_component_available(caps, 10);
    if (name == "ng") return moisture.graupel_number && plot3d_source_component_available(caps, 11);
    if (name == "qt") return moisture.vapor && plot3d_source_component_available(caps, 1);
    if (name == "qn") {
        return moisture.cloud_liquid && plot3d_source_component_available(caps, 2);
    }
    if (name == "qp") {
        const int precip_start_q = caps.moist_state_size > 3 ? 4 : 3;
        const int precip_start_comp = 3 + precip_start_q;
        const int precip_end_comp = caps.conserved_state_size - 1 - caps.moist_numconc_size;
        return moisture.rain &&
            plot3d_source_component_available(caps, precip_start_q) &&
            caps.moist_numconc_size >= 0 && precip_start_comp <= precip_end_comp &&
            precip_end_comp < caps.conserved_state_size;
    }
    if (name == "qsat") return moisture.vapor && plot3d_source_component_available(caps, 1);
    if (name == "precipitable") return moisture.vapor && plot3d_source_component_available(caps, 1);

    if (name == "rain_accum") return moisture.rain_accumulation &&
        caps.qmoist_size >= (caps.moisture_type == MoistureType::SuperDroplets ? 7 : 1);
    if (name == "snow_accum") return moisture.snow_accumulation && caps.qmoist_size >= 2;
    if (name == "graup_accum") return moisture.graupel_accumulation && caps.qmoist_size >= 3;
    if (name == "rel_humidity") return moisture.relative_humidity_diagnostic && caps.qmoist_size >= 6;
    if (name == "condensation_rate") return moisture.condensation_rate_diagnostic && caps.qmoist_size >= 4;

    return plot3d_fixed_unconditional_diagnostic(name);
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
