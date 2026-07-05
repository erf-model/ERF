#ifndef ERF_FIRE_LAYER_H
#define ERF_FIRE_LAYER_H

#include <AMReX_MultiFab.h>
#include <AMReX_Geometry.H>
#include <memory>
#include <ERF_FireParams.H>
#include <ERF_FireGrid.H>
#include <ERF_Rothermel.H>
#include <ERF_FarsiteEllipse.H>
#include <ERF_FireIgnition.H>
#include <ERF_FuelMoisture.H>
#include <ERF_MoistureExtinction.H>

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
 * | Name           | ncomp | Description                                  |
 * |----------------|-------|----------------------------------------------|
 * | fire_phi       | 1     | Level-set: <0 burned, >0 unburned, ≈0 front  |
 * | fire_wind_ref  | 2     | u,v at z_ref (e.g., 6.1 m) from MOST        |
 * | fire_wind_eff  | 2     | After WAF + terrain corrections              |
 * | fire_slopes    | 2     | dz/dx, dz/dy (static, computed at init)      |
 * | fire_curvature | 1     | Terrain curvature (static, at init)          |
 * | fire_ros       | 1     | Rate of spread [m/s]                         |
 * | fire_fuel_load | 1     | Fuel load [kg/m²]                            |
 * | fire_fuel_mc   | 3     | M_1hr, M_10hr, M_100hr [fraction]           |
 * | fire_mext      | 1     | Moisture of extinction [fraction] (Phase 4)  |
 * | fire_heat_flux | 1     | Heat flux [W/m²] (Phase 6, stub)             |
 * | fire_spread_vec| 2     | Scratch for FARSITE spread direction (Phase 3)|
 */
class FireLayer
{
public:

    /** @brief Constructor */
    FireLayer() = default;

    /** @brief Destructor */
    ~FireLayer() = default;

    /**
     * @brief Initialize fire module after SurfaceLayer is ready.
     *
     * Called from ERF::InitData_post(). Performs:
     * 1. Prerequisite verification (verify_fire_prerequisites)
     * 2. Create fire grid (create_fire_grid)
     * 3. Allocate all MultiFabs
     * 4. Compute static fields: terrain slopes, curvature
     * 5. Initialize fire level-set (ignition circle)
     * 6. Set initial fuel properties (moisture, fuel load)
     *
     * The caller (ERF, a friend of itself) passes the private members as
     * explicit arguments so that FireLayer never needs to access ERF's
     * private/protected data directly.
     *
     * @param[in] erf                ERF instance — used only for public
     *                               AmrCore accessors (boxArray, DistributionMap,
     *                               Geom).
     * @param[in] surface_layer_ptr  Pointer to the SurfaceLayer owned by ERF
     *                               (pass m_SurfaceLayer.get() at the call site).
     *                               May be nullptr when called before the surface
     *                               layer is constructed (prereq check only).
     * @param[in] z_phys_nd_atm      Atmospheric nodal terrain-height MultiFab
     *                               (pass *z_phys_nd[0] at the call site).
     * @param[in] fire_params        Fire configuration struct.
     *
     * @throw amrex::Error if prerequisites not met.
     */
    void initialize(const ERF&            erf,
                    const SurfaceLayer*   surface_layer_ptr,
                    const amrex::MultiFab& z_phys_nd_atm,
                    const FireParams&     fire_params);

    /**
     * @brief Advance fire simulation by one atmospheric timestep.
     *
     * Called from ERF::Advance() each timestep.
     * Executes the fire computation pipeline:
     * 1. Extract wind from atmosphere at reference height
     * 2. Apply Wind Adjustment Factor for canopy
     * 3. Apply terrain wind corrections (FARSITE)
     * 4. Update fuel moisture from atmospheric state (Phase 2)
     * 5. Compute effective midflame wind
     * 6. Compute rate-of-spread using Rothermel model
     * 7. (Phase 3) Advance level-set using ROS and wind
     * 8. (Phase 6) Compute heat flux to atmosphere
     *
     * @param[in]     dt             Timestep [s]
     * @param[in,out] surface_layer  SurfaceLayer with MOST diagnostics.
     *                               Non-const because get_u_star / get_z0 /
     *                               get_olen are not marked const in
     *                               ERF_SurfaceLayer.H.
     * @param[in]     xvel           x-velocity on atmospheric grid [m/s]
     * @param[in]     yvel           y-velocity on atmospheric grid [m/s]
     * @param[in]     z_phys_cc      Cell-center terrain height on atmospheric grid [m]
     * @param[in]     T_atm_k0       Atmospheric potential temperature at k=0 [K]
     *                               On atmospheric grid (coarser than fire grid).
     * @param[in]     RH_atm_k0      Atmospheric relative humidity at k=0 [fraction 0-1]
     *                               On atmospheric grid (coarser than fire grid).
     */
    void advance(amrex::Real time,
                 amrex::Real dt, 
                 SurfaceLayer& surface_layer,
                 const amrex::MultiFab& xvel,
                 const amrex::MultiFab& yvel,
                 const amrex::MultiFab& z_phys_cc,
                 const amrex::MultiFab& T_atm_k0,
                 const amrex::MultiFab& RH_atm_k0);

    /**
     * @brief Get heat flux for atmospheric coupling.
     * @return Pointer to fire_heat_flux MultiFab (nullptr until Phase 6).
     */
    const amrex::MultiFab* get_heat_flux() const
    {
        return fire_heat_flux.get();
    }

    /**
     * @brief Get rate-of-spread field for diagnostics.
     * @return Pointer to fire_ros MultiFab [m/s].
     */
    const amrex::MultiFab* get_ros() const
    {
        return fire_ros.get();
    }

    /**
     * @brief Get effective wind field for diagnostics.
     * @return Pointer to fire_wind_eff MultiFab [m/s].
     */
    const amrex::MultiFab* get_wind_eff() const
    {
        return fire_wind_eff.get();
    }

    /**
     * @brief Get level-set field for diagnostics.
     * @return Pointer to fire_phi MultiFab.
     */
    const amrex::MultiFab* get_levelset() const
    {
        return fire_phi.get();
    }

    /**
     * @brief Get mutable level-set field for checkpoint restart.
     * @return Mutable pointer to fire_phi MultiFab.
     */
    amrex::MultiFab* get_levelset_mut()
    {
        return fire_phi.get();
    }

    /**
     * @brief Get fire geometry for plotfile output.
     * @return Const reference to the fire grid geometry.
     */
    const amrex::Geometry& get_fire_geom() const
    {
        return m_fg.geom;
    }

    /**
     * @brief Get reference wind field for plotfile output.
     * @return Pointer to fire_wind_ref MultiFab [m/s].
     */
    const amrex::MultiFab* get_wind_ref() const
    {
        return fire_wind_ref.get();
    }

    /**
     * @brief Get terrain slopes for plotfile output.
     * @return Pointer to fire_slopes MultiFab [dz/dx, dz/dy].
     */
    const amrex::MultiFab* get_slopes() const
    {
        return fire_slopes.get();
    }

    /**
     * @brief Get fuel moisture content for plotfile output.
     * @return Pointer to fire_fuel_mc MultiFab [M_1hr, M_10hr, M_100hr].
     */
    const amrex::MultiFab* get_fuel_mc() const
    {
        return fire_fuel_mc.get();
    }

    /**
     * @brief Get fire grid reference for plotfile output.
     * @return Const reference to the FireGrid.
     */
    const FireGrid& get_fire_grid() const
    {
        return m_fg;
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
    std::unique_ptr<amrex::MultiFab>         fire_mext;         // Moisture of extinction (Phase 4)
    std::unique_ptr<amrex::MultiFab>         fire_heat_flux;    // Heat flux (stub)
    std::unique_ptr<amrex::MultiFab>         fire_spread_vec;   // FARSITE scratch
    std::unique_ptr<amrex::MultiFab>         fire_arrival_time; // Fire arrival time [s]

    // Precomputed Rothermel coefficients
    RothermelComputed                        m_rc;

    // FARSITE parameters for level-set propagation
    FarsiteParams                            m_fp;

    // Time tracking for fire advancement
    amrex::Real                              m_current_time = 0.0;
    amrex::Real                              m_dt_atm       = 0.0;

    // Vertical levels
    int                                      m_nz = 0;

    // Fuel bed depth [ft] (fuel-model dependent)
    amrex::Real                              m_fuel_bed_depth_ft = 6.5;

    // Ignition parameters
    amrex::Real                              m_ignition_x = 0.0;
    amrex::Real                              m_ignition_y = 0.0;
    amrex::Real                              m_ignition_r = 20.0;

    // Configuration
    FireParams                               m_params;

    /**
     * @brief Initialize fire level-set as a signed-distance circle.
     *
     * @param[in] center_x  Ignition centre x [m]
     * @param[in] center_y  Ignition centre y [m]
     * @param[in] radius    Ignition radius [m]
     *
     * Resulting fire_phi convention:
     * - fire_phi < 0  →  burned (inside circle)
     * - fire_phi > 0  →  unburned fuel (outside circle)
     * - fire_phi ≈ 0  →  fire front
     */
    void init_ignition(amrex::Real center_x, amrex::Real center_y,
                       amrex::Real radius);

    /**
     * @brief Apply Wind Adjustment Factor to the effective wind field.
     *
     * Multiplies fire_wind_eff by a scalar WAF. Used inside advance().
     */
    void apply_waf_to_wind();

    /**
     * @brief Update fire_fuel_mc and fire_mext for one atmospheric timestep.
     *
     * Called from advance() before ComputeRothermellSpreadRate. Updates
     * fuel moisture using time-lag ODEs and recomputes moisture of extinction.
     *
     * @param[in] dt_s       Atmospheric timestep [seconds]
     * @param[in] T_atm_k0   Atmospheric potential temperature at k=0 [K]
     * @param[in] RH_atm_k0  Atmospheric relative humidity at k=0 [fraction 0-1]
     */
    void advance_fuel_moisture(amrex::Real dt_s,
                               const amrex::MultiFab& T_atm_k0,
                               const amrex::MultiFab& RH_atm_k0);

};

#endif