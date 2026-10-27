#ifndef ERF_FIRE_ENGINE_H
#define ERF_FIRE_ENGINE_H

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <string>
#include <ERF_FireGrid.H>

class ERF;
class SurfaceLayer;

/**
 * @file ERF_FireEngine.H
 * @brief The interface between the fire coupling and a fire engine
 *
 * A fire engine advances a fire front on the fire grid and reports the heat
 * it releases. The shared coupling layer (FireCoupler) owns everything that
 * touches the atmosphere: the fire grid's placement, the atmospheric fields
 * handed to the engine, the area average of the engine's fluxes onto the
 * atmospheric columns, their one-step lag, their injection into the
 * atmosphere's sources and the smoke tracer. An engine sees only this
 * interface, the two structs below and the shared helpers it chooses to use
 * (the wind sampler of ERF_FireWindExtract.H, the column grounds of
 * ERF_FireSurfaceHeight.H).
 *
 * The engines built into this ERF are listed in ERF_FireEngineFactory.cpp.
 * ERF_FireEngineTemplate.H declares every method below for new engines: a
 * method added to FireEngine is declared there too (a pure virtual one that is
 * missing fails the build; an optional one does not).
 */

/**
 * @brief The atmospheric state a fire engine reads at the start of its step
 *
 * Every field lives on the level the fire grid refines (erf.fire.anchor_level),
 * in that level's index space and with its ghost cells; the engine samples
 * them onto its own grid (FireGrid::atm_lo gives the corner of the fire
 * region). The k = 0 slabs hold one cell in z per box of the level. The
 * slabs and the precipitation are temporaries of FireCoupler::advance(): an
 * engine reads them during its advance() and keeps no pointer to them.
 */
struct FireForcing
{
    amrex::Real            time;          ///< Start of the step [s]
    amrex::Real            dt;            ///< Step of the fire grid's level [s]
    const amrex::MultiFab& xvel;          ///< x-velocity on x-faces [m/s]: before the dycore (passive, lagged) or after it (synchronous)
    const amrex::MultiFab& yvel;          ///< y-velocity on y-faces [m/s], the same time level as xvel
    const amrex::MultiFab& z_phys_cc;     ///< Cell-centre height [m] of the terrain-following mesh
    const amrex::MultiFab& T_k0;          ///< Air temperature at k = 0 [K] (theta times Exner), before the dycore
    const amrex::MultiFab& RH_k0;         ///< Relative humidity at k = 0 [0-1], before the dycore
    const amrex::MultiFab* precip_accum_k0; ///< Surface precipitation accumulated since t = 0 [kg/m2] after this step's microphysics, or nullptr unless the engine asks for it (FireEngine::wants_surface_precip)
};

/**
 * @brief The fluxes a fire engine hands to the atmosphere after its step
 *
 * Both fields live on the fire grid, per unit ground area, as the mean over
 * the engine's step; the coupling averages them onto the atmospheric columns
 * and injects them at the next step. The pointers stay valid until the
 * engine's next advance. Sensible and latent pointers are independently
 * optional: nullptr means no flux for that channel on this step, without
 * suppressing a valid flux from the other channel.
 */
struct FireFluxes
{
    const amrex::MultiFab* sensible = nullptr;  ///< Sensible heat flux into the air [W/m2]; nullptr = none this step
    const amrex::MultiFab* latent   = nullptr;  ///< Latent heat flux into the air [W/m2]; nullptr = no water released
    /// With erf.fire.fire_debug: text the coupling appends to its debug line of
    /// the sensible flux (for example how the engine split its heat)
    std::string            debug_note;
};

/**
 * @brief A fire engine behind the shared coupling layer
 *
 * The coupling calls, in this order:
 * - verify_prerequisites() once at start-up, after its own checks of the
 *   ERF configuration and before it builds the fire grid;
 * - initialize() once, with the fire grid built, then start_run();
 * - read_checkpoint() on a restart, right after start_run();
 * - advance() then fill_fluxes() once per step of the fire grid's level,
 *   after that level's dycore;
 * - smoke_heat_per_kg(), heat_open_fraction_atm() and roof_height_atm() when
 *   it injects the lagged fluxes, at every Runge-Kutta stage of the next step;
 * - write_plotfile() on the fire plotfile schedule (erf.fire_plot_int,
 *   erf.fire_plot_per), and write_checkpoint() from FireCoupler::write_checkpoint
 *   when ERF writes a checkpoint.
 *
 * The getters at the end (level set, reference wind, heat release) are what
 * other surface models may read from the fire (a dust or smoke source model,
 * for example); an engine that has none returns nullptr and those couplings
 * do nothing.
 */
class FireEngine
{
public:
    virtual ~FireEngine() = default;

    /// Name of the engine, as erf.fire.engine selects it.
    virtual const char* name () const = 0;

    /**
     * @brief Check the engine's own inputs against the ERF configuration
     *
     * Called once at start-up, after the coupling's checks of the grid
     * (verify_fire_prerequisites) and before the fire grid is built. Aborts,
     * naming the input, on a configuration the engine cannot run.
     *
     * @param[in] erf  ERF instance (public AmrCore accessors only)
     * @param[in] lev  The level the fire grid refines
     */
    virtual void verify_prerequisites (const ERF& erf, int lev) = 0;

    /**
     * @brief Set the engine up on the fire grid
     *
     * @param[in] erf            ERF instance (public AmrCore accessors only)
     * @param[in] lev            The level the fire grid refines
     * @param[in] fg             The fire grid the coupling built
     * @param[in] surface_layer  The bottom surface layer (never null)
     * @param[in] z_phys_nd_atm  Nodal height of level lev [m] (ERF always passes it)
     */
    virtual void initialize (const ERF& erf,
                             int lev,
                             const FireGrid& fg,
                             const SurfaceLayer* surface_layer,
                             const amrex::MultiFab* z_phys_nd_atm) = 0;

    /**
     * @brief Start-up after initialize(), once the run's end is known
     *
     * @param[in] stop_time              ERF's stop_time [s] (numeric_limits max when unset)
     * @param[in] max_step               ERF's max_step, counting level-0 steps (negative when unset)
     * @param[in] steps_per_coarse_step  Steps of the fire grid's level per level-0 step
     */
    virtual void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step) = 0;

    /// Advance the fire over one step of its level.
    virtual void advance (const FireForcing& forcing) = 0;

    /**
     * @brief The fluxes of the step just taken, on the fire grid
     *
     * Called after advance(), only when the coupling injects heat or smoke.
     */
    virtual void fill_fluxes (FireFluxes& fluxes) = 0;

    /**
     * @brief Heat released per kilogram of fuel burned [J/kg], matching the
     *        sensible flux last handed over, to turn it into the fuel burnt
     *        that the smoke emission factor applies to
     */
    virtual amrex::Real smoke_heat_per_kg () const = 0;

    /// The engine needs the surface precipitation accumulation in FireForcing.
    virtual bool wants_surface_precip () const { return false; }

    /**
     * @brief Optional heat placement per atmospheric column, for
     *        erf.fire.heat_open_fraction: the fraction of the column's
     *        footprint open to the air [0-1] and the mean roof height of the
     *        rest [m], on the k = 0 slabs of the fire grid's level with one
     *        ghost cell: one slab per box of the level, in the level's box
     *        order and on its DistributionMapping (the coupling reads them
     *        with the level's MFIter and stops otherwise). nullptr when the
     *        engine has no such blocking.
     */
    virtual const amrex::MultiFab* heat_open_fraction_atm () const { return nullptr; }
    virtual const amrex::MultiFab* roof_height_atm () const { return nullptr; }

    /**
     * @brief Write the engine's own state to the checkpoint directory
     *
     * Everything the engine carries from one step to the next, under
     * checkpointname/Level_<lev>/ of the fire grid's level. Called by
     * FireCoupler::write_checkpoint, which then writes the lagged fluxes.
     */
    virtual void write_checkpoint (const std::string& checkpointname) const = 0;

    /**
     * @brief Restore the engine's state from a checkpoint on a restart
     *
     * @param[in] restart_chkfile  The checkpoint directory
     * @param[in] max_level        ERF's max_level (the levels to look for the fire state on)
     * @param[in] step             The step count of the fire grid's level at the checkpoint
     */
    virtual void read_checkpoint (const std::string& restart_chkfile, int max_level, int step) = 0;

    /**
     * @brief Write the engine's plotfile of the fire grid
     *
     * @param[in] plotfile_prefix  erf.fire_plot_file (the step is appended)
     * @param[in] time             Simulation time [s]
     * @param[in] step             Level-0 step number
     */
    virtual void write_plotfile (const std::string& plotfile_prefix, amrex::Real time, int step) const = 0;

    /// Level set on the fire grid, negative where burned; nullptr if the engine has none.
    virtual const amrex::MultiFab* level_set () const { return nullptr; }
    /// Wind (u, v) [m/s] the engine sampled at reference_wind_height() on the fire grid; nullptr if none.
    virtual const amrex::MultiFab* reference_wind () const { return nullptr; }
    /// Height above ground [m] of reference_wind().
    virtual amrex::Real reference_wind_height () const { return amrex::Real(0.0); }
    /// Heat released by the fire on the fire grid [W/m2], before any partition into
    /// sensible and latent flux; nullptr if the engine has none.
    virtual const amrex::MultiFab* heat_release () const { return nullptr; }
};

#endif
