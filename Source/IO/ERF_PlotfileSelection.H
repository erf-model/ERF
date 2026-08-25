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

/**
 * @brief What the active run can supply to a 3D plotfile.
 *
 * Availability of the *moist* components is not described here: it comes from
 * the active scheme's MoistureComponentIndices, which is the single source of
 * truth for which moisture variables exist (see the class comment on that
 * struct).  What remains in this struct is the conserved-state layout -- needed
 * to know where the moist window ends and the non-water species begin -- and the
 * optional non-moisture storage, which a plot variable can only be selected
 * against if it was actually allocated.
 */
struct Plot3DSelectionCapabilities
{
    //! Width of the allocated conserved state, i.e. NDRY + NSCALARS plus every
    //! component the microphysics interface adds.  This bounds the dry
    //! components and the non-water species some schemes carry above the moist
    //! window; it is deliberately *not* what gates the moist components.
    int conserved_state_size = 0;

    //! Number of conserved components in the moist window that starts at
    //! RhoQ1_comp.  Inside that window a component is gated by moisture_indices
    //! rather than by allocation, because several schemes allocate moist slots
    //! they never integrate.  Above the window sit the non-water species a
    //! scheme may add (SuperDroplets), which are allocated *and* integrated and
    //! so are gated by conserved_state_size like the dry components.
    //!
    //! The default is the full RhoQ1..RhoQ11 window, i.e. "let the moisture map
    //! decide about every rhoQn slot", which is the conservative answer for a
    //! caller that does not know the scheme's width.
    int moist_state_size = NMOIST_max;

    //! The active moisture scheme's component/diagnostic map.
    MoistureComponentIndices moisture_indices;

    bool time_average_storage = false;
    bool radiation_heating_storage = false;
    bool eddy_diffusivity_storage = false;
    bool dissipation_storage = false;
    bool wall_distance_storage = false;
};

/**
 * @brief Record the conserved-state layout of a run on its capabilities.
 *
 * Every output path that selects raw state names must describe the state the
 * same way, so the three fields are filled here rather than at each call site.
 *
 * @param caps              Capabilities to fill in.
 * @param mi                The active scheme's component/diagnostic map.
 * @param qstate_moist_size Conserved components in the moist window
 *                          (Microphysics::Get_Qstate_Moist_Size()).
 * @param qstate_size       All conserved components the microphysics adds,
 *                          moist window plus non-water species
 *                          (Microphysics::Get_Qstate_Size()).
 */
inline void
plot3d_set_state_capabilities (Plot3DSelectionCapabilities& caps,
                               const MoistureComponentIndices& mi,
                               const int qstate_moist_size,
                               const int qstate_size)
{
    caps.moisture_indices     = mi;
    caps.moist_state_size     = qstate_moist_size;
    caps.conserved_state_size = NDRY + NSCALARS + qstate_size;
}

/**
 * @brief Conserved-state component holding the q-th moist variable, or -1 if
 *        q is outside the fixed RhoQ1..RhoQ11 layout.
 */
inline int
plot3d_q_conserved_component_index (const int q)
{
    if (q < 1 || q > 11) return -1;
    return RhoQ1_comp + q - 1;
}

/**
 * @brief Conserved-state component behind a raw state variable name
 *        ("density", "rhotheta", "rhoQ3", ...), or -1 if the name is not one.
 */
inline int
plot3d_conserved_component_index (const std::string& name)
{
    if (name == "density")  return Rho_comp;
    if (name == "rhotheta") return RhoTheta_comp;
    if (name == "rhoKE")    return RhoKE_comp;
    if (name == "rhoadv_0") return RhoScalar_comp;

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

/**
 * @brief Diagnostics that every run can write, moist or dry.
 */
inline bool
plot3d_fixed_unconditional_diagnostic (const std::string& name)
{
    static constexpr const char* names[] = {
        "temp", "theta", "KE", "scalar", "soundspeed", "reflectivity",
        "max_reflectivity", "mucape", "vorticity_x", "vorticity_y", "vorticity_z",
        "helicity", "local_helicity", "magvel", "divU", "pres_hse", "dens_hse",
        "theta_hse", "pi_hse", "qv_hse", "pressure", "pert_pres", "pert_dens", "buoyancy",
        "eq_pot_temp", "VPD", "dpdx", "dpdy", "dpdz", "pres_hse_x", "pres_hse_y",
        "pres_hse_z", "z_phys", "detJ", "h_xi", "h_eta", "h_zeta",
        "mapfac", "lat_m", "lon_m", "pblh",
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

/**
 * @brief Whether a named 3D plot variable can be written by this run.
 *
 * The order of the tests below is the order of authority: raw state components
 * answer from the state layout, moisture variables answer from the active
 * scheme's index map, optional diagnostics answer from their storage flag, and
 * anything left is an always-available diagnostic.
 */
inline bool
plot3d_fixed_variable_available (const std::string& name,
                                 const Plot3DSelectionCapabilities& caps)
{
    const MoistureComponentIndices& mi = caps.moisture_indices;

    // Raw conserved-state components.  Inside the moist window the moisture map
    // decides, since allocated-but-never-integrated slots must not be published
    // as data.  Below the window the state is always allocated dry-or-moist, and
    // above it are the non-water species a scheme may add -- both of those are
    // real data whenever they were allocated, so a bounds check is the right
    // test for them.
    const int cons_comp = plot3d_conserved_component_index(name);
    if (cons_comp >= 0) {
        if ( (cons_comp >= RhoQ1_comp) &&
             (cons_comp <  RhoQ1_comp + caps.moist_state_size) ) {
            return mi.has_comp(cons_comp);
        }
        return (cons_comp < caps.conserved_state_size);
    }

    // Every moisture variable -- species, aggregates and moist diagnostics --
    // is decided by the scheme's index map.
    const auto moist_status = mi.query_var(name);
    if (moist_status.governed) return moist_status.available;

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

    return plot3d_fixed_unconditional_diagnostic(name);
}

/**
 * @brief Whether any requested variable needs the pressure field computed.
 */
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
