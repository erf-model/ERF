/**
 * \file ERF_PlotfileSelection.H
 */
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

struct Plot3DQRange
{
    int first_q = 0;
    int last_q = -1;

    [[nodiscard]] bool valid () const noexcept
    {
        return first_q >= 1 && last_q >= first_q;
    }
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
        // The second moist conserved component is the canonical cloud-water
        // state even when this model supplies no condensation source.
        caps.cloud_liquid = true;
        break;
    case MoistureType::SatAdj:
        caps.vapor = true;
        caps.cloud_liquid = true;
        caps.relative_humidity_diagnostic = true;
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
        caps.cloud_ice = true;
        caps.rain = true;
        caps.snow = true;
        caps.graupel = true;
        caps.rain_accumulation = true;
        // graupel accumulation is not advertised because nothing writes it
        caps.snow_accumulation = true;
        caps.relative_humidity_diagnostic = true;
        caps.condensation_rate_diagnostic = true;
        break;
    default:
        amrex::Abort("Unknown MoistureType in plot3d_moisture_capabilities");
    }

    return caps;
}

inline int
plot3d_q_conserved_component_index (const int q)
{
    if (q < 1 || q > 11) return -1;
    return 3 + q; // RhoQ1_comp + q - 1 in the fixed conserved layout.
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
            return plot3d_q_conserved_component_index(q);
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
    const int comp = plot3d_q_conserved_component_index(q);
    return plot3d_has_moist_state(caps) && comp >= 0 && comp < caps.conserved_state_size;
}

inline bool
plot3d_moist_q_component_available (const Plot3DSelectionCapabilities& caps,
                                     const int q)
{
    return q >= 1 && q <= caps.moist_state_size &&
        plot3d_source_component_available(caps, q);
}

inline Plot3DQRange
plot3d_total_mass_q_range (const int moist_state_size,
                           const int moist_numconc_size)
{
    if (moist_state_size <= 0 || moist_numconc_size < 0 ||
        moist_numconc_size > moist_state_size) {
        return {};
    }

    const int mass_count = moist_state_size - moist_numconc_size;
    if (mass_count <= 0) return {};
    return {1, mass_count};
}

inline Plot3DQRange
plot3d_nonprecipitating_q_range (const int moist_state_size)
{
    if (moist_state_size < 2) return {};
    return {1, moist_state_size > 3 ? 3 : 2};
}

inline Plot3DQRange
plot3d_precipitating_q_range (const int moist_state_size,
                              const int moist_numconc_size)
{
    const Plot3DQRange total =
        plot3d_total_mass_q_range(moist_state_size, moist_numconc_size);
    if (!total.valid()) return {};

    const int first_precip_q = moist_state_size > 3 ? 4 : 3;
    if (total.last_q < first_precip_q) return {};
    return {first_precip_q, total.last_q};
}

inline bool
plot3d_q_range_available (const Plot3DSelectionCapabilities& caps,
                          const Plot3DQRange range)
{
    return range.valid() &&
        plot3d_source_component_available(caps, range.first_q) &&
        plot3d_source_component_available(caps, range.last_q);
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
            plot3d_q_range_available(
                caps, plot3d_nonprecipitating_q_range(caps.moist_state_size));
    }
    if (name == "qv")    return moisture.vapor && plot3d_moist_q_component_available(caps, 1);
    if (name == "qc")    return moisture.cloud_liquid && plot3d_moist_q_component_available(caps, 2);
    if (name == "qi")    return moisture.cloud_ice && plot3d_moist_q_component_available(caps, 3);
    if (name == "qrain") return moisture.rain && plot3d_moist_q_component_available(
        caps, caps.moist_state_size > 3 ? 4 : 3);
    if (name == "qsnow") return moisture.snow && plot3d_moist_q_component_available(caps, 5);
    if (name == "qgraup") return moisture.graupel && plot3d_moist_q_component_available(caps, 6);
    if (name == "nc") return moisture.cloud_number && plot3d_moist_q_component_available(caps, 7);
    if (name == "ni") return moisture.ice_number && plot3d_moist_q_component_available(caps, 8);
    if (name == "nr") return moisture.rain_number && plot3d_moist_q_component_available(caps, 9);
    if (name == "ns") return moisture.snow_number && plot3d_moist_q_component_available(caps, 10);
    if (name == "ng") return moisture.graupel_number && plot3d_moist_q_component_available(caps, 11);
    if (name == "qt") {
        return moisture.vapor &&
            plot3d_q_range_available(
                caps, plot3d_total_mass_q_range(
                    caps.moist_state_size, caps.moist_numconc_size));
    }
    if (name == "qn") {
        return moisture.cloud_liquid &&
            plot3d_q_range_available(
                caps, plot3d_nonprecipitating_q_range(caps.moist_state_size));
    }
    if (name == "qp") {
        return moisture.rain &&
            plot3d_q_range_available(
                caps, plot3d_precipitating_q_range(
                    caps.moist_state_size, caps.moist_numconc_size));
    }
    if (name == "qsat") return moisture.vapor && plot3d_moist_q_component_available(caps, 1);
    if (name == "precipitable") return moisture.vapor && plot3d_moist_q_component_available(caps, 1);

    if (name == "rain_accum") return moisture.rain_accumulation &&
        caps.qmoist_size >= (caps.moisture_type == MoistureType::SuperDroplets ? 7 : 1);
    if (name == "snow_accum") return moisture.snow_accumulation && caps.qmoist_size >= 2;
    if (name == "graup_accum") return moisture.graupel_accumulation && caps.qmoist_size >= 3;
    if (name == "rel_humidity") return moisture.relative_humidity_diagnostic &&
        (caps.moisture_type == MoistureType::SatAdj || caps.qmoist_size >= 6);
    if (name == "condensation_rate") return moisture.condensation_rate_diagnostic && caps.qmoist_size >= 4;

    return plot3d_fixed_unconditional_diagnostic(name);
}

inline bool
plot3d_needs_pressure (const amrex::Vector<std::string>& names)
{
    for (const auto& name : names) {
        if (name == "pressure" || name == "pert_pres" ||
            name == "dpdx" || name == "dpdy" || name == "dpdz" ||
            name == "eq_pot_temp" || name == "VPD" || name == "qsat" ||
            name == "rel_humidity") {
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
