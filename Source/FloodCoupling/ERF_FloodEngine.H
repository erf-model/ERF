#ifndef ERF_FLOOD_ENGINE_H
#define ERF_FLOOD_ENGINE_H

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <string>
#include <vector>
#include <ERF_FloodGrid.H>

class ERF;

/**
 * @file ERF_FloodEngine.H
 * @brief The interface between the flood coupling and a flood engine
 *
 * A flood engine is a surface-water model (shallow water, a runoff or routing
 * model) that runs on the flood grid and takes the atmosphere's surface
 * precipitation as its source. The shared coupling layer (FloodCoupler) owns
 * everything that touches the atmosphere: the flood grid's placement on one
 * level, the precipitation of each step taken from the microphysics
 * accumulators, the wind handed to the engine, the engine's fields in the ERF
 * plotfile and the check that a restart uses the same flood grid. An engine
 * sees only this interface, FloodForcing and the helpers of ERF_FloodGrid.H.
 *
 * The flood coupling is one-way in this version: the engine reads the
 * atmosphere and returns nothing to it.
 *
 * ERF itself contains no engine. An engine lives in its own directory outside
 * ERF (FLOOD_HOME in GNU make, ERF_FLOOD_ENGINE_DIRS in CMake) and registers
 * itself with ERF_REGISTER_FLOOD_ENGINE (ERF_FloodEngineRegistry.H).
 * ERF_FloodEngineTemplate.H declares every method below for new engines: a
 * method added to FloodEngine is declared there too (a pure virtual one that
 * is missing fails the build; an optional one does not).
 *
 * Memory contract. Every MultiFab the coupling hands an engine (FloodForcing,
 * the FloodGrid layout) and every one the engine hands back (plot_var) is an
 * ordinary AMReX MultiFab in the default arena, so on a GPU build its data is
 * on the device. An engine that runs on the host:
 * - copies what it reads from FloodForcing to the host inside advance(),
 *   because those fields are temporaries of that call, and synchronizes the
 *   GPU stream (amrex::Gpu::streamSynchronize()) before it reads them there;
 * - owns the MultiFabs it returns, allocated on fg.ba and fg.dm, keeps them
 *   alive until its next advance(), and copies them back to device memory
 *   before returning them.
 * The coupling never frees, resizes or keeps across steps a pointer an engine
 * returns.
 */

/**
 * @brief The atmospheric state a flood engine reads for one step
 *
 * Every field lives on the level the flood grid refines
 * (erf.flood.anchor_level), in that level's index space; the engine maps them
 * onto its own grid (FloodGrid::atm_lo gives the corner of the flood region,
 * fill_flood_from_atm_k0 copies a column field onto the flood grid). The k = 0
 * slabs hold one cell in z per box of the level. The slabs are temporaries of
 * FloodCoupler::advance(): an engine reads them during its advance() and keeps
 * no pointer to them.
 */
struct FloodForcing
{
    double                 time;       ///< Start of the step [s], in double as ERF's own clock
    amrex::Real            dt;         ///< Step of the flood grid's level [s]
    const amrex::MultiFab& xvel;       ///< x-velocity on x-faces [m/s], after this step's dycore
    const amrex::MultiFab& yvel;       ///< y-velocity on y-faces [m/s], after this step's dycore
    const amrex::MultiFab& z_phys_cc;  ///< Cell-center height [m] of the terrain-following mesh
    /**
     * Surface precipitation of this step on the k = 0 slabs [kg/m2, equal to
     * mm of liquid water]: component 0 the liquid part, component 1 the frozen
     * part (snow, graupel and hail as liquid-water equivalent; Morrison and
     * WSM6 count their cloud-water fallout, and Morrison its cloud-ice
     * fallout, as liquid). It is the change of the microphysics accumulators
     * over the step, so the sum over the steps of a run is the change of ERF's
     * own accumulation over the run. Not bounded at zero: the liquid part is
     * the total less the frozen part for a scheme with a total accumulator and
     * can be slightly negative from round-off; an engine that bounds it should
     * carry the negative part to the next step. Where a finer level covers the
     * flood grid's level, ERF averages the finer accumulators down onto it
     * after the flood's step (and, with super-droplets, at plotfiles and
     * regrids); the next step carries that correction, of either sign, in the
     * columns under the finer level. In single precision the per-step value is
     * resolved to the round-off of the accumulator. nullptr unless the engine
     * asks for it (FloodEngine::wants_surface_precip).
     */
    const amrex::MultiFab* precip_k0;
};

/**
 * @brief A flood engine behind the shared coupling layer
 *
 * The coupling calls, in this order:
 * - plot_var_names() once, when the coupler is constructed (before
 *   verify_prerequisites and before ERF reads its plot variable lists);
 * - verify_prerequisites() once at start-up, after its own checks of the
 *   ERF configuration and before it builds the flood grid;
 * - initialize() once, with the flood grid built, then start_run();
 * - read_checkpoint() on a restart from a checkpoint that holds flood state,
 *   right after start_run();
 * - advance() once per step of the flood grid's level, after that level's
 *   dycore and microphysics;
 * - write_plotfile() on the flood plotfile schedule (erf.flood.plot_int,
 *   erf.flood.plot_per), write_checkpoint() when ERF writes a checkpoint, and
 *   plot_var() each time ERF writes its own plotfile, for the flood fields
 *   that file carries.
 */
class FloodEngine
{
public:
    virtual ~FloodEngine() = default;

    /// Name of the engine, as erf.flood.engine selects it.
    virtual const char* name () const = 0;

    /**
     * @brief Check the engine's own inputs against the ERF configuration
     *
     * Called once at start-up, after the coupling's checks of the grid
     * (verify_flood_prerequisites) and before the flood grid is built. Aborts,
     * naming the input, on a configuration the engine cannot run.
     *
     * @param[in] erf  ERF instance (public AmrCore accessors only)
     * @param[in] lev  The level the flood grid refines
     */
    virtual void verify_prerequisites (const ERF& erf, int lev) = 0;

    /**
     * @brief Set the engine up on the flood grid
     *
     * @param[in] erf            ERF instance (public AmrCore accessors only)
     * @param[in] lev            The level the flood grid refines
     * @param[in] fg             The flood grid the coupling built
     * @param[in] z_phys_nd_atm  Nodal height of level lev [m]: the terrain
     *                           the atmosphere sees, at k = 0
     */
    virtual void initialize (const ERF& erf,
                             int lev,
                             const FloodGrid& fg,
                             const amrex::MultiFab* z_phys_nd_atm) = 0;

    /**
     * @brief Start-up after initialize(), once the run's end is known
     *
     * @param[in] stop_time              Elapsed time at which the run stops [s], on the clock the
     *                                   engine's step times use (ERF's stop_time less start_time);
     *                                   numeric_limits<amrex::Real>::max() when no stop time is set
     * @param[in] max_step               ERF's max_step, counting level-0 steps;
     *                                   numeric_limits<int>::max() when unset (never negative)
     * @param[in] steps_per_coarse_step  Steps of the flood grid's level per level-0 step
     */
    virtual void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step)
    {
        amrex::ignore_unused(stop_time, max_step, steps_per_coarse_step);
    }

    /// Advance the flood over one step of its level.
    virtual void advance (const FloodForcing& forcing) = 0;

    /**
     * @brief The engine takes the surface precipitation in FloodForcing
     *
     * True by default: the coupling then stops at start-up unless the
     * microphysics scheme accumulates surface precipitation. An engine driven
     * by other sources only (an inflow hydrograph) returns false.
     */
    virtual bool wants_surface_precip () const { return true; }

    /**
     * @brief Write the engine's own state to the checkpoint directory
     *
     * Everything the engine carries from one step to the next, under
     * checkpointname/Level_<lev>/ of the flood grid's level. The default
     * writes nothing; such an engine cannot restart (read_checkpoint).
     */
    virtual void write_checkpoint (const std::string& checkpointname) const
    {
        amrex::ignore_unused(checkpointname);
    }

    /**
     * @brief Restore the engine's state from a checkpoint on a restart
     *
     * Called only when the checkpoint holds flood state (it was written with
     * erf.flood.enable = true) and the coupling has checked that it belongs to
     * this run's flood grid. A restart from a checkpoint without flood state
     * (an atmosphere-only spin-up) skips this call, and the engine continues
     * from its initial state.
     *
     * @param[in] restart_chkfile  The checkpoint directory
     * @param[in] step             The step count of the flood grid's level at the checkpoint
     *
     * The default stops the run: an engine that does not implement restart
     * says so instead of continuing from its initial state.
     */
    virtual void read_checkpoint (const std::string& restart_chkfile, int step)
    {
        amrex::ignore_unused(restart_chkfile, step);
        amrex::Abort(std::string("[FLOOD] flood engine \"") + name()
                     + "\" does not support restart from a checkpoint");
    }

    /**
     * @brief Write the engine's plotfile of the flood grid
     *
     * @param[in] plotfile_prefix  erf.flood.plot_file (the step is appended)
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
     * @brief Flood fields the ERF plotfile can carry, by name
     *
     * The names the deck may list in erf.plot_vars_*. Each must start with
     * "flood_" (the coupler stops otherwise), so they never collide with an
     * ERF variable. The coupler asks once, when it is constructed, so the
     * answer depends on the engine's inputs only. The default lists none.
     */
    virtual std::vector<std::string> plot_var_names () const { return {}; }

    /**
     * @brief The field plot_var_names() lists as name, on the flood grid
     *
     * One component, on FloodGrid::ba and FloodGrid::dm, valid until the next
     * advance(). The coupler averages it over each atmospheric column (the
     * grid_ratio^2 flood cells under it) and writes it into the k = 0 layer of
     * the ERF plotfile, zero elsewhere. nullptr writes zeros.
     */
    virtual const amrex::MultiFab* plot_var (const std::string& name) const
    {
        amrex::ignore_unused(name);
        return nullptr;
    }
};

#endif
