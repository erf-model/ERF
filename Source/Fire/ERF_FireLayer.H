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

// Forward declarations
class SurfaceLayer;
class ERF;

class FireLayer
{
public:

    FireLayer() = default;
    ~FireLayer() = default;

    void initialize(const ERF&             erf,
                    const SurfaceLayer*    surface_layer_ptr,
                    const amrex::MultiFab& z_phys_nd_atm,
                    const FireParams&      fire_params);

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

    FireGrid                                 m_fg;

    std::unique_ptr<amrex::MultiFab>         fire_phi;           // Level-set
    std::unique_ptr<amrex::MultiFab>         fire_wind_ref;      // Wind at z_ref
    std::unique_ptr<amrex::MultiFab>         fire_wind_eff;      // Effective wind
    std::unique_ptr<amrex::MultiFab>         fire_slopes;        // Terrain slopes
    std::unique_ptr<amrex::MultiFab>         fire_curvature;     // Terrain curvature
    std::unique_ptr<amrex::MultiFab>         fire_ros;           // Rate of spread
    std::unique_ptr<amrex::MultiFab>         fire_fuel_load;     // Fuel load
    std::unique_ptr<amrex::MultiFab>         fire_fuel_mc;       // Fuel moisture (3 comp)
    std::unique_ptr<amrex::MultiFab>         fire_mext;          // Moisture of extinction
    std::unique_ptr<amrex::MultiFab>         fire_heat_flux;     // Heat flux (stub)
    std::unique_ptr<amrex::MultiFab>         fire_spread_vec;    // FARSITE scratch
    std::unique_ptr<amrex::MultiFab>         fire_arrival_time;  // Fire arrival time [s]
    std::unique_ptr<amrex::MultiFab>         fire_disp_accum;    // Persistent displacement accumulator (2 comp)

    RothermelComputed                        m_rc;
    FarsiteParams                            m_fp;

    amrex::Real                              m_current_time      = 0.0;
    amrex::Real                              m_dt_atm            = 0.0;
    int                                      m_nz                = 0;
    amrex::Real                              m_fuel_bed_depth_ft = 6.5;
    amrex::Real                              m_ignition_x        = 0.0;
    amrex::Real                              m_ignition_y        = 0.0;
    amrex::Real                              m_ignition_r        = 20.0;
    FireParams                               m_params;


    /**
     * @brief Apply Wind Adjustment Factor to the effective wind field.
     *
     * Multiplies fire_wind_eff by a scalar WAF. Used inside advance().
     */
    void apply_waf_to_wind();
    void advance_fuel_moisture(amrex::Real dt_s,
                               const amrex::MultiFab& T_atm_k0,
                               const amrex::MultiFab& RH_atm_k0);
};

#endif