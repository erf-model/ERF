#ifndef ERF_DUST_MODEL_H
#define ERF_DUST_MODEL_H

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <string>

class ERF;
class SurfaceLayer;
#ifdef ERF_USE_FIRE
class FireEngine;
struct FireGrid;
#endif

/**
 * @file ERF_DustModel.H
 * @brief The interface between ERF and a dust model
 *
 * A dust model emits dust from the ground into the dust scalars of the
 * conserved state, settles and deposits it, and writes its own output. ERF
 * owns the dust scalars: ERF_NUM_DUST_SCALARS (CMake; NUM_DUST_SCALARS in GNU
 * make) of them, starting at RhoDust_comp (ERF_IndexDefines.H), advected and
 * diffused like every passive scalar and zero at the start of a run. A model
 * adds its tendencies to them in add_sources() and keeps everything else (its
 * surface grid, its fields, its particles) itself.
 *
 * ERF itself contains no model. A model lives in its own directory outside ERF
 * (DUST_HOME in GNU make, ERF_DUST_MODEL_DIRS in CMake) and registers itself
 * with ERF_REGISTER_DUST_MODEL (ERF_DustModelRegistry.H). erf.dust.enable =
 * true turns dust on; erf.dust.model names the model, and may be left out
 * when one model is compiled in.
 *
 * ERF calls, in this order:
 * - initialize() once, at the end of start-up (after the fire coupling's
 *   initialization), which on a restart also restores the model's state;
 * - add_sources() at every Runge-Kutta stage of every level, right after
 *   ERF's own sources are built and before the fire coupling adds its fluxes;
 * - advance() once per step of level 0, after the step;
 * - write_output() at start-up (step 0, before the first plotfile), after
 *   each step and at the end of the run, each time before ERF writes a
 *   checkpoint, so a checkpoint records the model's output of that step as
 *   written;
 * - write_checkpoint() each time ERF writes a checkpoint.
 *
 * Memory contract, as for fire engines (ERF_FireEngine.H): every MultiFab in
 * DustInit, DustStep and the add_sources() arguments is an ordinary AMReX
 * MultiFab in the default arena (device memory on a GPU build) and a
 * temporary of the call; a model keeps no pointer to it.
 */

#ifdef ERF_USE_FIRE
/**
 * @brief The fire, read-only, for a dust model that couples to it
 *
 * The engine's getters (level_set(), reference_wind(),
 * reference_wind_height(), heat_release()) are what a dust model may read;
 * an engine without one returns nullptr. The fields live on the fire grid.
 */
struct DustFire
{
    const FireEngine& engine;   ///< The fire engine (erf.fire.engine)
    const FireGrid&   grid;     ///< The fire grid: its layout and refinement factor C (erf.fire.grid_ratio)
    double            time;     ///< The time the engine's fields belong to [s]
};
#endif

/**
 * @brief What a dust model reads when it is set up
 */
struct DustInit
{
    const ERF&             erf;              ///< ERF instance (public accessors only)
    SurfaceLayer*          surface_layer;    ///< The bottom surface layer, or nullptr when zlo is not one
    const amrex::MultiFab& z_phys_nd;        ///< Nodal height [m] of level 0
    const std::string&     restart_chkfile;  ///< The checkpoint the run restarts from; empty on a fresh start
    int                    step;             ///< Level 0's step: 0, or the step of the checkpoint
    double                 time;             ///< Level 0's time [s]: 0, or the time of the checkpoint
#ifdef ERF_USE_FIRE
    const DustFire*        fire;             ///< The fire, when the run has one; nullptr otherwise
#endif
};

/**
 * @brief The state a dust model reads after each step of level 0
 */
struct DustStep
{
    int                    step;             ///< The step just taken (0 for the first one)
    double                 time;             ///< Start of the step [s]
    amrex::Real            dt;               ///< Step of level 0 [s]
    const amrex::Geometry& geom;             ///< Geometry of level 0
    const amrex::MultiFab& cons;             ///< Conserved state of level 0 at the end of the step
    const amrex::MultiFab& xvel;             ///< x-velocity on x-faces [m/s], end of the step
    const amrex::MultiFab& yvel;             ///< y-velocity on y-faces [m/s], end of the step
    const amrex::MultiFab& zvel;             ///< z-velocity on z-faces [m/s], end of the step
    const amrex::MultiFab* z_phys_cc;        ///< Cell-centre height [m]; nullptr when ERF holds none
    const amrex::MultiFab* surface_q1_flux;  ///< Surface vapour flux of level 0 (the turbulence scheme's z flux of q1); nullptr when there is none
    SurfaceLayer*          surface_layer;    ///< The bottom surface layer, or nullptr
#ifdef ERF_USE_FIRE
    const DustFire*        fire;             ///< The fire, when the run has one; nullptr otherwise
#endif
};

/**
 * @brief A dust model behind ERF's dust coupling
 */
class DustModel
{
public:
    virtual ~DustModel () = default;

    /// Name of the model, as erf.dust.model selects it.
    virtual const char* name () const = 0;

    /// Number of dust scalars the model fills, from RhoDust_comp on; at most
    /// NSCALARS_DUST (ERF stops the run otherwise).
    virtual int num_scalars () const { return 1; }

    /**
     * @brief Set the model up, and on a restart restore its state
     *
     * Reads the model's inputs (erf.dust.*). On a restart
     * (init.restart_chkfile not empty) the model reads back what its
     * write_checkpoint() wrote there. After this call, on a fresh start, ERF
     * zeroes the dust scalars on every level.
     */
    virtual void initialize (const DustInit& init) = 0;

    /**
     * @brief Add the dust tendencies to the source of the conserved state
     *
     * @param[in]    lev      The level; a model that runs on level 0 only returns on the others
     * @param[inout] cc_src   Source of the conserved state [per second]; the model adds to
     *                        the dust components only
     * @param[in]    state    The conserved state of this Runge-Kutta stage
     * @param[in]    detJ_cc  Cell Jacobian (1 on a mesh without terrain): a cell's height is detJ dz
     * @param[in]    geom     Geometry of the level
     * @param[in]    dt       The step of the stage's slow update [s]
     */
    virtual void add_sources (int lev, amrex::MultiFab& cc_src, const amrex::MultiFab& state,
                              const amrex::MultiFab& detJ_cc, const amrex::Geometry& geom,
                              amrex::Real dt) = 0;

    /// Advance the model's own state over the step of level 0 just taken.
    virtual void advance (const DustStep& step) = 0;

    /// Write the model's output for step (0 at start-up); is_final at the end of the run.
    virtual void write_output (int /*step*/, double /*time*/, bool /*is_final*/) {}

    /// Write the model's state into the checkpoint directory, for initialize() to restore.
    virtual void write_checkpoint (const std::string& /*checkpointname*/) const {}
};

#endif
