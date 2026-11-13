#ifndef ERF_STRUCTURE_MODEL_H
#define ERF_STRUCTURE_MODEL_H

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_RealBox.H>

#include <functional>
#include <string>
#include <vector>

/**
 * @file ERF_StructureModel.H
 * @brief The interface between ERF and a model of structures in the flow
 *
 * A structure model represents bodies that the air loads and that push back on
 * the air: power-line conductors and their towers, for example. It samples the
 * wind at its own points, advances its own dynamics, and returns a momentum
 * source that ERF adds to the air at every Runge-Kutta stage. The velocity
 * sampler and the force spreading of the actuator core (Source/MovingBodies/Core,
 * ERF_ActuatorSampling.H and ERF_ActuatorSpreading.H) are the tools a model is
 * expected to use for both.
 *
 * ERF itself contains no structure model. A model lives in its own directory
 * outside ERF (STRUCTURE_HOME in GNU make, ERF_STRUCTURE_MODEL_DIRS in CMake)
 * and registers itself with ERF_REGISTER_STRUCTURE_MODEL
 * (ERF_StructureModelRegistry.H); erf.structures.model selects it by name.
 * Tests/StructureModelExample is a complete model directory.
 *
 * ERF calls, in this order:
 * - the registered creator once, in the ERF constructor (the model reads its
 *   own inputs there, under its own prefix);
 * - level() each time ERF needs the level the model couples to;
 * - initialize() once, after the initial data (or the checkpoint) exist and
 *   before the step-0 plotfile; on a restart the model restores its own state
 *   there, from StructureInit::restart_chkfile. It is not called again after a
 *   regrid: what the model needs of the mesh it reads from StructureStep each step;
 * - advance() once per step of level(), with the state at the start of the step,
 *   before that step's Runge-Kutta stages;
 * - support() right after each advance(): when a finer level than level() exists,
 *   ERF stops the run if that level's boxes overlap the support (the force would
 *   be lost there, see level());
 * - add_momentum_sources() at every Runge-Kutta stage of level(), and, with
 *   erf.coupling_type = TwoWay, once at the start of each step of every coarser
 *   level, into zeroed MultiFabs on the grids of level(), which ERF restricts
 *   onto that level (see add_momentum_sources());
 * - write_checkpoint() each time ERF writes a checkpoint;
 * - finalize() once, when the run ends normally, after the last plotfile and checkpoint;
 * - plot_var_names() each time ERF builds its plotfile variable list, and
 *   fill_plot_var() for each listed name on every level of each plotfile.
 *
 * Memory contract. Every MultiFab ERF hands a model is an ordinary AMReX
 * MultiFab in the default arena (device memory on a GPU build) and is valid
 * only during the call. A model keeps no pointer to it after the call returns.
 */

/**
 * @brief What a structure model reads once, at start-up
 *
 * All fields live on the level the model couples to (StructureModel::level()).
 */
struct StructureInit
{
    int                    lev;               ///< The level the model couples to
    double                 time;              ///< Simulation time [s]: 0, or the checkpoint's time on a restart
    const amrex::Geometry& geom;              ///< Geometry of level lev
    const amrex::MultiFab& xvel;              ///< x-velocity on x-faces [m/s]
    const amrex::MultiFab* z_phys_nd;         ///< Nodal height [m] of a stretched or terrain-following mesh; nullptr on a uniform-dz mesh (flat ground at prob_lo z)
    const amrex::MultiFab* detJ_cc;           ///< Cell Jacobian of a stretched or terrain-following mesh; nullptr on a uniform-dz mesh
    /// With erf.terrain_type = ImmersedForcing: the terrain height [m] at the nodes of level lev
    /// (a k = 0 slab over the domain grown by 3 nodes), from the surface the immersed boundary
    /// was built from, in host-readable memory and valid only during initialize(); nullptr
    /// otherwise, or when wants_immersed_terrain_surface() is false
    const amrex::FArrayBox* immersed_terrain_surface;
    bool                   unit_map_factors;  ///< Every map factor of level lev is 1 (the actuator core assumes so)
    bool                   keqn_rans;         ///< Level lev runs the k-equation RANS closure, so cons holds rho k
    amrex::Real            Cmu0;              ///< The closure's C_mu^0 on level lev
    std::string            restart_chkfile;   ///< The checkpoint ERF restarts from; empty on a fresh start
    /// The nominal step of level lev [s]: erf.fixed_dt divided by the sub-cycling ratios down to lev.
    /// Not positive when the run has no fixed step (ERF then chooses each step from erf.cfl).
    /// StructureStep::dt is the step taken: it differs on the first steps with erf.init_shrink
    /// below 1 or after a restart (growth limited by erf.change_max), and on a last step that
    /// ERF fits to the stop time
    double                 fixed_dt;
    /// The time [s], on the clock of time, at which the run stops: the input stop_time, or
    /// stop_datetime less start_datetime; the largest double when neither is given. It is ERF's
    /// own stop, which can differ from the input by about 1e-7 s when start_datetime is given,
    /// so compare with a tolerance. max_step may end the run earlier
    double                 stop_time;
    /// The run stops after this many level-0 steps, counted from the start of the first run (a
    /// restart continues the count); the largest int when max_step is not given or negative
    int                    max_step;
    /// A level of the run, not only lev, is FullyCompressible (erf.compressibility_type); false
    /// when every level is Anelastic or PseudoIncompressible. It names the solver, not the time
    /// scheme: a FullyCompressible level with erf.timestepping_type = RK3_NoSubstepping counts
    bool                   fully_compressible;
};

/**
 * @brief What a structure model may read from the fire, when the run has one
 *
 * The fire engine's fields (the getters of ERF_FireEngine.H) on the fire grid: a
 * horizontal grid, one cell in z, that refines one atmospheric level
 * (erf.fire.anchor_level). That level need not be the one the model couples to.
 * Fire cell (i, j) has its centre at geom.ProbLo() + ((i + 1/2) dx, (j + 1/2) dy),
 * in the coordinates of the model's own points.
 *
 * The fields are the engine's state at time: the time ERF assigns to the fire's
 * level at the end of the fire's last step, or the time ERF started or restarted
 * from before that step. The fire steps after the dycore of its own level. A
 * model on the fire's level or a coarser one therefore sees the fire at its own
 * step's time, exactly; a model
 * on a finer level sees it up to one step of the fire's level ahead of its own
 * step. On a restart, until the fire's first step, they are what the engine
 * restored from the checkpoint (an engine that does not restore a field holds
 * its initial value there). Each field is spread over the ranks as the fire grid is: a model that
 * needs the fire at its own points reduces over the ranks itself (every rank
 * calls advance()). Like the flow fields, they are temporaries of the call.
 */
struct StructureFire
{
    const amrex::Geometry& geom;          ///< Geometry of the fire grid [m] (one cell in z)
    int                    lev;           ///< The atmospheric level the fire grid refines
    double                 time;          ///< The time the fields belong to [s] (see above)
    const amrex::MultiFab* heat_release;  ///< Heat released [W/m2] on the fire grid (FireEngine::heat_release); nullptr if the engine has none
    const amrex::MultiFab* level_set;     ///< Level set on the fire grid, negative where burned (FireEngine::level_set); nullptr if the engine has none
    /// Area burned so far [m2] (FireEngine::burned_area), negative when the engine does not
    /// track it. It may reduce over the ranks: every rank calls it, or none.
    std::function<amrex::Real()> burned_area;
};

/**
 * @brief The flow a structure model reads at the start of each step of its level
 *
 * The state at the start of the step, before its Runge-Kutta stages, on level lev
 * with the level's ghost cells filled. Temporaries of the call: the model samples
 * them during advance() and keeps no pointer to them.
 */
struct StructureStep
{
    int                    lev;        ///< The level the model couples to
    int                    finest_level; ///< ERF's finest level now
    double                 time;       ///< Start of the step [s]
    double                 dt;         ///< Step of level lev [s], in double as ERF keeps it, also in a single-precision build
    const amrex::Geometry& geom;       ///< Geometry of level lev
    const amrex::MultiFab& xvel;       ///< x-velocity on x-faces [m/s]
    const amrex::MultiFab& yvel;       ///< y-velocity on y-faces [m/s]
    const amrex::MultiFab& zvel;       ///< z-velocity on z-faces [m/s]
    const amrex::MultiFab& cons;       ///< The conserved state (rho, rho theta, the scalars; rho k with the k-equation RANS)
    const amrex::MultiFab* z_phys_nd;  ///< Nodal height [m]; nullptr on a uniform-dz mesh, as in StructureInit
    const amrex::MultiFab* detJ_cc;    ///< Cell Jacobian; nullptr on a uniform-dz mesh, as in StructureInit
    bool                   unit_map_factors; ///< Every map factor of level lev is 1 now (checked each step, so also after a regrid)
    const StructureFire*   fire;       ///< The fire (see StructureFire) when the run has one; nullptr otherwise
};

/**
 * @brief A model of structures in the flow, behind ERF's structure coupling
 */
class StructureModel
{
public:
    virtual ~StructureModel () = default;

    /// Name of the model, as erf.structures.model selects it.
    virtual const char* name () const = 0;

    /**
     * @brief The level the model couples to
     *
     * ERF calls advance() on this level only, and add_momentum_sources() for
     * this level only (see add_momentum_sources() for the coarser levels).
     * Aborts, naming the model's input, when the model cannot couple to a level
     * of this run (for example a level above finest_level).
     *
     * The level must be the finest level that covers the model's support. On a
     * refined level the model's points and their kernels must lie inside the
     * level's boxes, also after a regrid (erf_actuator::points_covered_by checks
     * it). The force reaches the coarser levels through ERF's averaging down and,
     * with TwoWay coupling, as a source restricted from this level.
     * It never reaches a finer level: averaging down would overwrite the forced
     * cells with that level's unforced data, so ERF stops the run when a finer
     * level's boxes overlap support().
     *
     * @param[in] finest_level  ERF's finest level now
     */
    virtual int level (int finest_level) const = 0;

    /**
     * @brief Where the model's force can be nonzero this step
     *
     * Boxes in physical coordinates [m] (z as the height of z_phys_nd on a
     * terrain-following mesh) that hold every point together with its kernel's
     * reach, as of the last advance(). The same on every rank. ERF checks them
     * against the next finer level's boxes after each advance() (see level()).
     * The boxes may extend past the domain; a periodic direction wraps.
     * Every box needs its low end at or below its high end in each direction and no
     * NaN (an infinite end is allowed); ERF stops the run otherwise, naming the model.
     */
    virtual std::vector<amrex::RealBox> support () const = 0;

    /// The model wants StructureInit::immersed_terrain_surface (ERF builds it only on request).
    virtual bool wants_immersed_terrain_surface () const { return false; }

    /**
     * @brief Set the model up once the mesh and the initial state exist
     *
     * On a restart (StructureInit::restart_chkfile not empty) the model restores
     * what write_checkpoint() wrote and continues from it. Aborts, naming the
     * input, on a configuration the model cannot run.
     */
    virtual void initialize (const StructureInit& init) = 0;

    /**
     * @brief Advance the structures over one step of level()
     *
     * Samples the flow, advances the structures from step.time to step.time +
     * step.dt, and sets the momentum source add_momentum_sources() returns at the
     * stages of this step.
     */
    virtual void advance (const StructureStep& step) = 0;

    /**
     * @brief Add the structures' force on the air to the momentum sources
     *
     * Face-centred sources on level lev [N/m^3, i.e. kg/(m^2 s^2)], added (not set),
     * held fixed over the step. Called at every Runge-Kutta stage, after ERF's
     * own sources. ERF applies no mask to them: the model keeps its kernels clear
     * of embedded boundaries, immersed terrain and buildings.
     *
     * With erf.coupling_type = TwoWay, ERF also calls it at the start of each step
     * of every coarser level, with lev = level() and zeroed MultiFabs on the grids
     * of level() as of the last advance() (they may have ghost faces, which the
     * model leaves alone). It adds the source of its last advance(), as at a
     * stage. ERF restricts that source onto the coarser level, keeping the force
     * it puts on the air (its integral over the faces' control volumes, with the
     * mesh Jacobian on a stretched or terrain-following mesh), and holds it over
     * that level's step: the air under the finer level then feels the force while
     * the coarser level steps, and the coarse/fine faces, which the finer level
     * takes from the coarser one, get the force once.
     * Before the first advance() and after a regrid of level(), until the next
     * advance(), ERF does not call it for the coarser levels. ERF writes the
     * source to its checkpoint, so after a restart, until the first advance(), it
     * restricts the checkpointed source instead: the model need not hold one.
     * ERF also calls it when it writes a checkpoint; being const, it leaves the
     * model as it was.
     */
    virtual void add_momentum_sources (int lev,
                                       amrex::MultiFab& xmom_src,
                                       amrex::MultiFab& ymom_src,
                                       amrex::MultiFab& zmom_src) const = 0;

    /**
     * @brief Write the model's state under the checkpoint directory
     *
     * Everything the model carries from one step to the next, under
     * checkpointname/<name()>/ (the model creates that directory on the I/O rank).
     * The default writes nothing, and supports_restart() must then stay false.
     */
    virtual void write_checkpoint (const std::string& checkpointname) const
    {
        amrex::ignore_unused(checkpointname);
    }

    /**
     * @brief The model restores its state in initialize() on a restart
     *
     * False by default: ERF then stops a restart instead of letting the model
     * continue from its initial state.
     */
    virtual bool supports_restart () const { return false; }

    /**
     * @brief The run has ended normally
     *
     * Called once, on every rank, after the last step, plotfile and checkpoint: the place to
     * flush and close the model's own files or to shut down a solver it drives. Not called
     * when the run aborts. The default does nothing.
     */
    virtual void finalize () {}

    /**
     * @brief Fields the ERF plotfile can carry, by name
     *
     * The names the deck may list in erf.plot_vars_*; they start with the model's
     * own prefix (for example "conductor_") so they never collide with an ERF
     * variable. ERF stops the run on an empty name, a name listed twice, or a
     * name ERF or another model already uses. The answer depends on the model's
     * inputs only. The default lists none.
     */
    virtual std::vector<std::string> plot_var_names () const { return {}; }

    /**
     * @brief Fill the field plot_var_names() lists as name on level lev
     *
     * One cell-centred component written into dst at component comp, on dst's
     * boxes, with no ghost cells; zero where the model has nothing (for example
     * on levels other than level()).
     */
    virtual void fill_plot_var (const std::string& name, int lev,
                                amrex::MultiFab& dst, int comp) const
    {
        amrex::ignore_unused(name, lev);
        dst.setVal(amrex::Real(0.0), comp, 1, 0);
    }
};

#endif
