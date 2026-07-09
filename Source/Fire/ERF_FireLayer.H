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
#include <ERF_FireStatsOutput.H>
#include <ERF_FireHeatFlux.H>
#include <ERF_FireDiagnostics.H>
#include <ERF_FireAtmCoupling.H>
#include <ERF_AlbiniSpotting.H>
#include <ERF_ScottSpottingTable.H>
#include <ERF_CrownFire.H>
#include <ERF_FuelMap.H>
#include <ERF_FireBreak.H>
#include <ERF_FuelBlending.H>
#include <ERF_IgnitionSchedule.H>
#include <ERF_PolygonIgnition.H>
#include <ERF_FireAcceleration.H>
// Phase 13A: Per-fuel wind height
#include <ERF_FuelWindHeight.H>
// Phase 13B: ROS models
#include <ERF_BalbiModel.H>
#include <ERF_CheneyGouldModel.H>
#include <ERF_MacArthurModel.H>
// Phase 15: BEHAVE multi-class Rothermel
#include <ERF_BehaveModel.H>
// Level-set propagation headers
#include <ERF_NumericalSchemes.H>
#include <ERF_Reinitialize.H>
#include <ERF_LevelSetAdvection.H>

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
 * | fire_fuel_mc   | 5     | M_1hr, M_10hr, M_100hr, M_lh, M_lw [fraction] |
 * | fire_mext      | 1     | Moisture of extinction [fraction] (Phase 4)  |
 * | fire_heat_flux | 1     | Heat flux [W/m²] (Phase 6, stub)             |
 * | fire_spread_vec| 2     | Scratch for FARSITE spread direction (Phase 3)|
 * | fire_crown_active  | 1 | Crown fire status: 1=active crown, 0=surface fire (Phase 9) |
 * | fire_crown_load    | 1 | Remaining crown fuel load [kg/m²] (Phase 9)                 |
 * | fire_crown_fraction_burned | 1 | Crown fraction burned [0-1] (Phase 9)            |
 * | fire_flame_tilt    | 1 | Flame tilt angle [deg from vertical] (Phase 9)             |
 * | fire_flame_temp    | 1 | Flame temperature [K] (Phase 9)                             |
 *
 * **fire_fuel_mc component layout (Phase 15, 5 components):**
 *   - comp 0: M_1hr      [fraction] — 1-hour dead fuel moisture (single-class path)
 *   - comp 1: M_10hr     [fraction] — 10-hour dead fuel moisture (single-class path)
 *   - comp 2: M_100hr    [fraction] — 100-hour dead fuel moisture (single-class path)
 *   - comp 3: M_lh       [fraction] — Live herbaceous moisture (BEHAVE multi-class path only)
 *   - comp 4: M_lw       [fraction] — Live woody moisture (BEHAVE multi-class path only)
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
     * @brief Get fire arrival time field for diagnostics.
     * @return Pointer to fire_arrival_time MultiFab.
     */
    const amrex::MultiFab* get_arrival_time() const
    {
        return fire_arrival_time.get();
    }

    /**
     * @brief Get mutable fire arrival time field for checkpoint restart.
     * @return Mutable pointer to fire_arrival_time MultiFab.
     */
    amrex::MultiFab* get_arrival_time_mut()
    {
        return fire_arrival_time.get();
    }

    /**
     * @brief Get mutable fuel load field for checkpoint restart.
     * @return Mutable pointer to fire_fuel_load MultiFab.
     */
    amrex::MultiFab* get_fuel_load_mut()
    {
        return fire_fuel_load.get();
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
     * @brief Get wind extraction height for plotfile output.
     * @return Pointer to fire_wind_extract_z MultiFab [m].
     */
    const amrex::MultiFab* get_wind_extract_z() const
    {
        return fire_wind_extract_z.get();
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
     * @brief Get fuel load for plotfile output.
     * @return Pointer to fire_fuel_load MultiFab [kg/m²].
     */
    const amrex::MultiFab* get_fuel_load() const
    {
        return fire_fuel_load.get();
    }

    /**
     * @brief Get surface temperature for plotfile output.
     * @return Pointer to fire_surface_temp MultiFab [K].
     */
    const amrex::MultiFab* get_surface_temp() const
    {
        return fire_surface_temp.get();
    }

    /**
     * @brief Get surface relative humidity for plotfile output.
     * @return Pointer to fire_surface_rh MultiFab [0-1].
     */
    const amrex::MultiFab* get_surface_rh() const
    {
        return fire_surface_rh.get();
    }

    /**
     * @brief Get fireline intensity for diagnostics and plotfile output.
     * @return Pointer to fire_fireline_intensity MultiFab [kW/m].
     */
    const amrex::MultiFab* get_fireline_intensity() const
    {
        return fire_fireline_intensity.get();
    }

    /**
     * @brief Get flame length for diagnostics and plotfile output.
     * @return Pointer to fire_flame_length MultiFab [m].
     */
    const amrex::MultiFab* get_flame_length() const
    {
        return fire_flame_length.get();
    }


    /**
     * @brief Get flame temperature for diagnostics and plotfile output.
     * @return Pointer to fire_flame_temp MultiFab [K].
     */
    const amrex::MultiFab* get_flame_temp() const
    {
        return fire_flame_temp.get();
    }

    /**
     * @brief Get flame tilt angle for diagnostics and plotfile output.
     * @return Pointer to fire_flame_tilt MultiFab [deg from vertical], or nullptr when disabled.
     */
    const amrex::MultiFab* get_flame_tilt() const
    {
        return fire_flame_tilt.get();
    }

    /**
     * @brief Get crown fire activation field for diagnostics and plotfile output.
     * @return Pointer to fire_crown_active MultiFab, or nullptr when crown fire is disabled.
     */
    const amrex::MultiFab* get_crown_active() const
    {
        return fire_crown_active.get();
    }

    /**
     * @brief Get remaining crown fuel load for diagnostics and plotfile output.
     * @return Pointer to fire_crown_load MultiFab [kg/m²], or nullptr when crown fire is disabled.
     */
    const amrex::MultiFab* get_crown_load() const
    {
        return fire_crown_load.get();
    }

    /**
     * @brief Get crown fraction burned diagnostic for plotfile output.
     * @return Pointer to fire_crown_fraction_burned MultiFab [0-1].
     */
    const amrex::MultiFab* get_crown_fraction_burned() const
    {
        return fire_crown_fraction_burned.get();
    }

    /**
     * @brief Get Albini spotting diagnostics for plotfile output.
     * @return Pointer to fire_albini_data MultiFab (nullptr when spotting disabled).
     */
    const amrex::MultiFab* get_albini_data() const
    {
        return fire_albini_data.get();
    }

    /**
     * @brief Get fire grid reference for plotfile output.
     * @return Const reference to the FireGrid.
     */
    const FireGrid& get_fire_grid() const
    {
        return m_fg;
    }

    /**
     * @brief Get coarsened sensible heat flux on atmospheric grid from previous timestep.
     * @return Pointer to m_Q_atm_prev MultiFab (nullptr until first advance).
     */
    const amrex::MultiFab* get_Q_atm_prev() const
    {
        return m_Q_atm_prev.get();
    }

    /**
     * @brief Get coarsened latent heat flux on atmospheric grid from previous timestep.
     * @return Pointer to m_Q_lat_atm_prev MultiFab (nullptr until first advance).
     */
    const amrex::MultiFab* get_Q_lat_atm_prev() const
    {
        return m_Q_lat_atm_prev.get();
    }

    /**
     * @brief Get fire latent heat flux on fire grid.
     * @return Pointer to fire_latent_flux MultiFab [W/m²].
     */
    const amrex::MultiFab* get_latent_flux() const
    {
        return fire_latent_flux.get();
    }

    /**
     * @brief Store the current fire heat flux into the previous-step buffer.
     *
     * Called from ERF::Advance() at the END of the fire advance, after
     * compute_heat_flux_and_diagnostics() has filled fire_heat_flux.
     * Coarsens fire_heat_flux from fire grid to atmospheric grid and
     * computes latent flux.
     *
     * @param[in] geom_atm  Atmospheric geometry (used for coarsening)
     */
    void update_atm_flux_buffer(const amrex::Geometry& geom_atm);

    /**
     * @brief Apply the previous-step fire fluxes to cc_source using WRF-Fire profile.
     *
     * Injects the one-step-lagged fire flux from previous timestep into cc_source
     * using the exponential distribution profile. Called from ERF::Advance()
     * BEFORE advance_dycore().
     *
     * @param[in,out] cc_source    3D source array for conserved variables
     * @param[in]     S_old        atmospheric conserved state at time n (for density)
     * @param[in]     z_phys_cc    cell-centre terrain height [m]
     * @param[in]     geom_atm     atmospheric geometry
     * @param[in]     has_moisture whether moisture is active
     */
    void apply_fire_coupling_to_cc_source(
        amrex::MultiFab& cc_source,
        const amrex::MultiFab& S_old,
        const amrex::MultiFab& z_phys_cc,
        const amrex::Geometry& geom_atm,
        bool has_moisture);

    /**
     * @brief Get the fire parameters for access to coupling type and other settings.
     * @return Const reference to the FireParams.
     */
    const FireParams& get_params() const { return m_params; }

    /**
     * @brief Get the per-cell fuel model MultiFab (Phase 10).
     * @return Const pointer to fire_fuel_model MultiFab (nullptr if uniform fuel).
     */
    const amrex::MultiFab* get_fuel_model() const { return fire_fuel_model.get(); }

    /**
     * @brief Check if spatial fuel map is active (Phase 10).
     * @return true if per-cell fuel model is loaded and active.
     */
    bool has_spatial_fuel() const { return m_has_spatial_fuel; }

private:

    // Fire grid and geometry
    FireGrid                                 m_fg;

    // Fire state MultiFabs (all on m_fg grid, C*nx × C*ny × 1)
    std::unique_ptr<amrex::MultiFab>         fire_phi;          // Level-set
    std::unique_ptr<amrex::MultiFab>         fire_wind_ref;     // Wind at z_ref
    std::unique_ptr<amrex::MultiFab>         fire_wind_eff;     // Effective wind
    std::unique_ptr<amrex::MultiFab>         fire_wind_extract_z; // Wind extraction height [m] (z_surf + z_ref)
    std::unique_ptr<amrex::MultiFab>         fire_slopes;       // Terrain slopes
    std::unique_ptr<amrex::MultiFab>         fire_curvature;    // Terrain curvature
    std::unique_ptr<amrex::MultiFab>         fire_ros;          // Rate of spread
    std::unique_ptr<amrex::MultiFab>         fire_fuel_load;    // Fuel load
    std::unique_ptr<amrex::MultiFab>         fire_fuel_mc;      // Fuel moisture (3 comp)
    std::unique_ptr<amrex::MultiFab>         fire_mext;         // Moisture of extinction (Phase 4)
    std::unique_ptr<amrex::MultiFab>         fire_heat_flux;    // Heat flux (stub)
    std::unique_ptr<amrex::MultiFab>         fire_spread_vec;   // FARSITE scratch
    std::unique_ptr<amrex::MultiFab>         fire_arrival_time; // Fire arrival time [s]
    std::unique_ptr<amrex::MultiFab>         fire_disp_accum;   // Persistent displacement accumulator (2 comp)
    std::unique_ptr<amrex::MultiFab>         fire_surface_temp; // Surface temperature [K] (Phase 2)
    std::unique_ptr<amrex::MultiFab>         fire_surface_rh;   // Surface relative humidity [0-1] (Phase 2)
    std::unique_ptr<amrex::MultiFab>         fire_fireline_intensity; // Byram fireline intensity [kW/m]
    std::unique_ptr<amrex::MultiFab>         fire_flame_length;       // Thomas flame length [m]
    std::unique_ptr<amrex::MultiFab>         fire_flame_tilt;         // Flame tilt angle [deg from vertical] (Phase 9)
    std::unique_ptr<amrex::MultiFab>         fire_flame_temp;         // Flame temperature [K] (Phase 9)

    // Phase 9: Crown fire diagnostics (2 comp, on fire grid)
    // comp 0: crown fire status (1.0 = active crown, 0.0 = surface),
    // comp 1: remaining crown fuel load [kg/m²]
    std::unique_ptr<amrex::MultiFab>         fire_crown_active;       // Crown fire status: 1=active, 0=surface (Phase 9)
    std::unique_ptr<amrex::MultiFab>         fire_crown_load;         // Remaining crown fuel load [kg/m²] (Phase 9)
    std::unique_ptr<amrex::MultiFab>         fire_crown_fraction_burned; // Crown fraction burned [0-1] (Phase 9)

    // Phase 8: Albini spotting diagnostics (4 comp, on fire grid)
    // comp 0: lofting height H_z [m], comp 1: brand count,
    // comp 2: max landing dist [m], comp 3: active spot flag
    std::unique_ptr<amrex::MultiFab>         fire_albini_data;
     
    // Phase 12: Per-cell acceleration state for temporal model (3 components).
    // comp 0: ROS_current [m/s], comp 1: elapsed_time [s], comp 2: ROS_equilibrium [m/s].
    // Only allocated when accel.enable = true AND accel.use_temporal = true.
    // nullptr when acceleration is disabled or size-based model is used.
    std::unique_ptr<amrex::MultiFab>  fire_accel_state;
     
    // Phase 10: Per-cell fuel model (non-null when fuel_map.file is non-empty)
    std::unique_ptr<amrex::MultiFab>         fire_fuel_model;  // integer code stored as Real
    amrex::Gpu::DeviceVector<int>            m_d_fuel_codes;   // device copy of integer codes
    bool                                     m_has_spatial_fuel = false;
     
    // Phase 11: Multi-ignition schedule runtime state
    IgnitionSchedule  m_ignition_schedule;   // loaded from ignition_schedule_file
    bool              m_has_schedule = false;
     
    // Phase 6: fire-atmosphere coupling MultiFabs
    // Coarsened fire fluxes on the atmospheric level-0 grid (nx × ny × 1).
    // These hold fire fluxes from the previous timestep for injection into
    // cc_source at the start of the next timestep (one-step explicit lag).
    std::unique_ptr<amrex::MultiFab>         m_Q_atm_prev;            // sensible [W/m²]
    std::unique_ptr<amrex::MultiFab>         m_Q_lat_atm_prev;        // latent   [W/m²]
    std::unique_ptr<amrex::MultiFab>         fire_latent_flux;        // latent heat flux on fire grid [W/m²]

    // Precomputed Rothermel coefficients
    RothermelComputed                        m_rc;

    // FARSITE parameters for level-set propagation
    FarsiteParams                            m_fp;

    // Time tracking for fire advancement
    amrex::Real                              m_current_time = 0.0;
    amrex::Real                              m_dt_atm       = 0.0;
    int                                      m_step         = 0;  ///< Simulation step counter

    // Level-set propagation tracking
    int                                      m_levelset_subcycle_count = 0;  ///< Counts fire subcycles for reinit scheduling

    // Vertical levels
    int                                      m_nz = 0;

    // Fuel bed depth [ft] (fuel-model dependent)
    amrex::Real                              m_fuel_bed_depth_ft = 6.5;

    // Initial fuel load [kg/m²] for diagnostics (Phase 5)
    amrex::Real                              m_fuel_load_initial_kg_m2 = 0.0;

    // Phase 13A: Per-fuel wind height device tables
    amrex::Gpu::DeviceVector<amrex::Real>    m_d_fcwh;          // fcwh[0..13], device copy
    amrex::Gpu::DeviceVector<amrex::Real>    m_d_fcz0;          // fcz0[0..13], device copy
    bool                                     m_use_per_fuel_wind_ht = false;

    // Phase 13B: Pre-computed ROS model coefficients (set in initialize, used in advance)
    BalbiComputed                            m_bc_default;      // Balbi global default (uniform fuel)
    CheneyGouldComputed                      m_cgc;             // Cheney-Gould computed params
    amrex::Gpu::DeviceVector<BalbiComputed>  m_d_balbi_table;   // per-fuel Balbi table (spatial fuel)

    // Phase 15: BEHAVE pre-computed state (set in initialize from domain-avg moisture,
    // updated in advance when moisture_dynamic = false)
    BehaveState                              m_bs_default;

    // Ignition parameters
    amrex::Real                              m_ignition_x = 0.0;
    amrex::Real                              m_ignition_y = 0.0;
    amrex::Real                              m_ignition_r = 20.0;

    // Configuration
    FireParams                               m_params;


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

    /**
     * @brief Compute heat flux and deplete fuel load for one atmospheric timestep (Phase 5).
     *
     * Called from advance() after fire_subcycle. Computes:
     * - Sensible heat flux [W/m²] in fire_heat_flux
     * - Fuel load depletion via exponential burnout model
     * - Byram fireline intensity [kW/m] in fire_fireline_intensity
     * - Thomas flame length [m] in fire_flame_length
     *
     * @param[in] dt_fire_s  Last fire subcycle timestep [s]
     */
    void compute_heat_flux_and_diagnostics(amrex::Real dt_fire_s);

};

#endif