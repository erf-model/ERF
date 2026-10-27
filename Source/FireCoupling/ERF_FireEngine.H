#ifndef ERF_FIRE_ENGINE_H
#define ERF_FIRE_ENGINE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <string>
#include <vector>
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
 * ERF itself contains no engine. An engine lives in its own directory outside
 * ERF (FIRE_HOME in GNU make, ERF_FIRE_ENGINE_DIRS in CMake) and registers
 * itself with ERF_REGISTER_FIRE_ENGINE (ERF_FireEngineRegistry.H).
 * ERF_FireEngineTemplate.H declares every method below for new engines: a
 * method added to FireEngine is declared there too (a pure virtual one that is
 * missing fails the build; an optional one does not).
 *
 * Memory contract. Every MultiFab the coupling hands an engine (FireForcing,
 * the FireGrid layout) and every one the engine hands back (FireFluxes, the
 * getters below) is an ordinary AMReX MultiFab in the default arena, so on a
 * GPU build its data is on the device. The coupling reads the engine's fields
 * on the device right after the call that returned them. An engine that runs
 * on the host:
 * - copies what it reads from FireForcing to the host inside advance()
 *   (fire_host_copy(), ERF_FireUtils.H, does the copy and the stream
 *   synchronisation), because those fields are temporaries of that call;
 * - owns the MultiFabs it returns, allocated on fg.ba and fg.dm (or the
 *   level's k = 0 slabs where stated), and keeps them alive until its next
 *   advance(); it fills them on the host and copies them back to device
 *   memory (fire_host_copy() the other way) before returning from
 *   fill_fluxes() or the getter;
 * - may hold managed or pinned memory internally, but synchronises the GPU
 *   stream (amrex::Gpu::streamSynchronize()) before it touches the host side
 *   of data the device may still be writing.
 * The coupling never frees, resizes or keeps across steps a pointer an engine
 * returns.
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
 * Both fields live on the fire grid (FireGrid::ba, FireGrid::dm), per unit
 * ground area, and are positive into the atmosphere. Each is the time mean
 * over exactly the interval of the last advance(), [forcing.time,
 * forcing.time + forcing.dt], whatever steps the engine takes inside it: an
 * engine that subcycles accumulates its fluxes over its substeps and divides
 * by forcing.dt. The coupling averages them onto the atmospheric columns and
 * injects them, unchanged, at every stage of the next ERF step of the fire
 * grid's level (the one-step lag of WRF-Fire; erf.fire.coupling_type =
 * synchronous only changes the wind the engine saw).
 *
 * The latent flux is the energy of the water vapour released, so the vapour
 * mass flux is latent / L_v [kg/m2/s] with L_v of ERF_Constants.H; an engine
 * that computes a mass flux multiplies by that same L_v. The pointers stay
 * valid until the engine's next advance(). Sensible and latent pointers are
 * independently optional: nullptr means no flux for that channel on this
 * step, without suppressing a valid flux from the other channel.
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
 * - write_plotfile() on the fire plotfile schedule (erf.fire.plot_int,
 *   erf.fire.plot_per), and write_checkpoint() from FireCoupler::write_checkpoint
 *   when ERF writes a checkpoint;
 * - plot_var_names() once, when the coupler is constructed (before
 *   verify_prerequisites), and plot_var() each time ERF writes its own
 *   plotfile, for the fire fields that file carries.
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
     * @param[in] surface_layer  The bottom surface layer; never null when
     *                           needs_surface_layer() is true, and nullptr when
     *                           it is false and the bottom boundary is not one
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
    virtual void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step)
    {
        amrex::ignore_unused(stop_time, max_step, steps_per_coarse_step);
    }

    /// Advance the fire over one step of its level.
    virtual void advance (const FireForcing& forcing) = 0;

    /**
     * @brief The fluxes of the step just taken, on the fire grid
     *
     * Means over exactly [forcing.time, forcing.time + forcing.dt] of the last
     * advance() (see FireFluxes). Called after advance(), only when the
     * coupling injects heat or smoke.
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
     * @brief The engine needs a surface layer at the bottom boundary
     *
     * True by default: the coupling then stops at start-up unless zlo.type =
     * surface_layer. An engine that reads nothing from the surface layer (it
     * takes the wind and heights from FireForcing alone) returns false and may
     * run over any bottom boundary.
     */
    virtual bool needs_surface_layer () const { return true; }

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
     * FireCoupler::write_checkpoint, which then writes the lagged fluxes. The
     * default writes nothing; such an engine cannot restart (read_checkpoint).
     */
    virtual void write_checkpoint (const std::string& checkpointname) const
    {
        amrex::ignore_unused(checkpointname);
    }

    /**
     * @brief Restore the engine's state from a checkpoint on a restart
     *
     * @param[in] restart_chkfile  The checkpoint directory
     * @param[in] max_level        ERF's max_level (the levels to look for the fire state on)
     * @param[in] step             The step count of the fire grid's level at the checkpoint
     *
     * The default stops the run: an engine that does not implement restart
     * says so instead of continuing from its initial state.
     */
    virtual void read_checkpoint (const std::string& restart_chkfile, int max_level, int step)
    {
        amrex::ignore_unused(restart_chkfile, max_level, step);
        amrex::Abort(std::string("[FIRE] fire engine \"") + name()
                     + "\" does not support restart from a checkpoint");
    }

    /**
     * @brief Write the engine's plotfile of the fire grid
     *
     * @param[in] plotfile_prefix  erf.fire.plot_file (the step is appended)
     * @param[in] time             Simulation time [s]
     * @param[in] step             Level-0 step number
     *
     * The default writes nothing (the engine's fields may still reach the ERF
     * plotfile through plot_var_names()).
     */
    virtual void write_plotfile (const std::string& plotfile_prefix, amrex::Real time, int step) const
    {
        amrex::ignore_unused(plotfile_prefix, time, step);
    }

    /**
     * @brief Fire fields the ERF plotfile can carry, by name
     *
     * The names the deck may list in erf.plot_vars_* (for example "fire_phi",
     * "fire_heat_flux"); they should start with "fire_" so they never collide
     * with an ERF variable. The coupler asks once, when it is constructed:
     * before initialize() and before ERF reads its plot variable lists, so the
     * answer depends on the engine's inputs only. The default lists none.
     */
    virtual std::vector<std::string> plot_var_names () const { return {}; }

    /**
     * @brief The field plot_var_names() lists as name, on the fire grid
     *
     * One component, on FireGrid::ba and FireGrid::dm, valid until the next
     * advance(). The coupler averages it over each atmospheric column (the
     * grid_ratio^2 fire cells under it) and writes it into the k = 0 layer of
     * the ERF plotfile, zero elsewhere. nullptr writes zeros.
     */
    virtual const amrex::MultiFab* plot_var (const std::string& name) const
    {
        amrex::ignore_unused(name);
        return nullptr;
    }

    /// Area burned so far [m2], for the run log; negative when the engine does not track it.
    virtual amrex::Real burned_area () const { return amrex::Real(-1.0); }

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
