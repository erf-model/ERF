#ifndef ERF_FIRE_LAYER_H
#define ERF_FIRE_LAYER_H

#include <AMReX_MultiFab.h>
#include <AMReX_Geometry.H>
#include <memory>
#include <ERF_FireParams.H>
#include <ERF_FireGrid.H>
#include <ERF_Rothermel.H>

// Forward declarations
class SurfaceLayer;
class ERF;

/**
 * @file ERF_FireLayer.H
 * @brief Main fire layer container and orchestrator
 *
 * The FireLayer class is the primary interface between the ERF atmospheric
 * solver and the wildfire spread module. It owns all fire-specific data
 * (level-set, wind, ROS, slopes, etc.) and orchestrates the computation
 * pipeline each timestep.
 *
 * Owned by ERF as std::unique_ptr<FireLayer> m_fire_layer, created in
 * ERF::InitData_post() if erf.fire.enable = true.
 */

/**
 * @class FireLayer
 * @brief Manages fire simulation on a 2D level-set grid
 *
 * Maintains all fire-related MultiFabs on the refined fire grid and calls
 * the wind extraction, WAF, terrain, and Rothermel computations each timestep.
 *
 * Multifabs (all on FireGrid, size C*nx × C*ny × 1):
 *
 * | Name | ncomp | Description |
 * |------|-------|---|
 * | fire_phi | 1 | Level-set: <0 burned, >0 unburned, ≈0 front |
 * | fire_wind_ref | 2 | u,v at z_ref (e.g., 6.1 m) from MOST |
 * | fire_wind_eff | 2 | After WAF + terrain corrections |
 * | fire_slopes | 2 | dz/dx, dz/dy (static, computed at init) |
 * | fire_curvature | 1 | Terrain curvature (static, at init) |
 * | fire_ros | 1 | Rate of spread [m/s] |
 * | fire_fuel_load | 1 | Fuel load [kg/m²] |
 * | fire_fuel_mc | 3 | M_1hr, M_10hr, M_100hr [fraction] |
 * | fire_heat_flux | 1 | Heat flux [W/m²] (Phase 6, stub) |
 * | fire_spread_vec | 2 | Scratch for FARSITE spread direction (Phase 3) |
 */
class FireLayer
{
public:

    /**
     * @brief Constructor
     */
    FireLayer() = default;

    /**
     * @brief Destructor
     */
    ~FireLayer() = default;

    /**
     * @brief Initialize fire module after SurfaceLayer is ready
     *
     * Called from ERF::InitData_post(). Performs:
     * 1. Prerequisite verification (verify_fire_prerequisites)
     * 2. Create fire grid (create_fire_grid)
     * 3. Allocate all MultiFabs
     * 4. Compute static fields: terrain slopes, curvature
     * 5. Initialize fire level-set (ignition circle)
     * 6. Set initial fuel properties (moisture, fuel load)
     *
     * @param[in] erf             ERF instance (for grid, geometry, SurfaceLayer access)
     * @param[in] fire_params     Fire configuration
     *
     * @throw amrex::Error if prerequisites not met
     */
    void initialize(const ERF& erf, const FireParams& fire_params);

    /**
     * @brief Advance fire simulation by one timestep
     *
     * Called from ERF::Advance() each atmospheric timestep.
     * Executes the fire computation pipeline:
     * 1. Extract wind from MOST layer at reference height
     * 2. Apply Wind Adjustment Factor for canopy
     * 3. Apply terrain wind corrections (FARSITE)
     * 4. Compute effective midflame wind
     * 5. Compute rate-of-spread using Rothermel model
     * 6. (Phase 3) Advance level-set using ROS and wind
     * 7. (Phase 6) Compute heat flux to atmosphere
     *
     * @param[in] dt              Timestep [s]
     * @param[in] surface_layer   SurfaceLayer with MOST diagnostics
     */
    void advance(amrex::Real dt, const SurfaceLayer& surface_layer);

    /**
     * @brief Get heat flux for atmospheric coupling
     *
     * @return Reference to fire_heat_flux MultiFab
     *
     * Returns nullptr if heat flux not yet computed.
     * Used in Phase 6 for adding fire heat source to atmosphere.
     */
    const amrex::MultiFab* get_heat_flux() const
    {
        return fire_heat_flux.get();
    }

    /**
     * @brief Get ROS field for diagnostics
     *
     * @return Reference to fire_ros MultiFab [m/s]
     */
    const amrex::MultiFab* get_ros() const
    {
        return fire_ros.get();
    }

    /**
     * @brief Get effective wind for diagnostics
     *
     * @return Reference to fire_wind_eff MultiFab [m/s]
     */
    const amrex::MultiFab* get_wind_eff() const
    {
        return fire_wind_eff.get();
    }

    /**
     * @brief Get level-set for diagnostics
     *
     * @return Reference to fire_phi MultiFab
     */
    const amrex::MultiFab* get_levelset() const
    {
        return fire_phi.get();
    }

private:

    // Fire grid and geometry
    FireGrid                                 m_fg;

    // Fire state MultiFabs (all on m_fg grid, C*nx × C*ny × 1)
    std::unique_ptr<amrex::MultiFab>         fire_phi;          // Level-set
    std::unique_ptr<amrex::MultiFab>         fire_wind_ref;     // Wind at z_ref
    std::unique_ptr<amrex::MultiFab>         fire_wind_eff;     // Effective wind
    std::unique_ptr<amrex::MultiFab>         fire_slopes;       // Terrain slopes
    std::unique_ptr<amrex::MultiFab>         fire_curvature;    // Terrain curvature
    std::unique_ptr<amrex::MultiFab>         fire_ros;          // Rate of spread
    std::unique_ptr<amrex::MultiFab>         fire_fuel_load;    // Fuel load
    std::unique_ptr<amrex::MultiFab>         fire_fuel_mc;      // Fuel moisture (3 comp)
    std::unique_ptr<amrex::MultiFab>         fire_heat_flux;    // Heat flux (stub)
    std::unique_ptr<amrex::MultiFab>         fire_spread_vec;   // FARSITE scratch

    // Precomputed coefficients
    RothermelComputed                        m_rc;

    // Configuration
    FireParams                               m_params;

    /**
     * @brief Initialize fire level-set as a circle
     *
     * @param[in] center_x Ignition center x [m]
     * @param[in] center_y Ignition center y [m]
     * @param[in] radius   Ignition radius [m]
     *
     * Sets fire_phi to a circle:
     * - fire_phi < 0 inside (burned area)
     * - fire_phi > 0 outside (unburned fuel)
     * - fire_phi ≈ 0 on the front
     */
    void init_ignition(amrex::Real center_x, amrex::Real center_y,
                       amrex::Real radius);

    /**
     * @brief Apply Wind Adjustment Factor to wind field
     *
     * Multiplies wind speed by WAF per cell. Used in advance().
     */
    void apply_waf_to_wind();

};

#endif
