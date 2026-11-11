#ifndef ERF_FIRE_ENGINE_TEMPLATE_H
#define ERF_FIRE_ENGINE_TEMPLATE_H

#include <ERF_FireEngine.H>

/**
 * @file ERF_FireEngineTemplate.H
 * @brief A starting point for a new fire engine: declarations only
 *
 * This class declares every FireEngine method a new fire model provides and
 * says what ERF hands it and expects back. Nothing here has a body, no input
 * selects it and nothing creates one, so it adds no code to ERF. The stand-in
 * engine of the FireCoupler unit tests (Tests/Unit/Fire/ERF_GTestFireCoupler.cpp)
 * is the smallest engine that works.
 *
 * To start a new engine:
 *  1. Copy this file into the new model's engine directory (outside ERF),
 *     rename the class and the include guard, and write the bodies in a .cpp
 *     file there. Delete the
 *     declarations of the optional methods the engine does not provide (their
 *     defaults in FireEngine then apply); a declaration left without a body
 *     fails to link.
 *  2. Read the model's inputs under its own prefix (fire.<engine>.* is
 *     recommended), so they never mix with the coupling's keys
 *     (FireCouplingParams, erf.fire.*) or another engine's.
 *  3. Register the engine in that .cpp, at namespace scope:
 *         ERF_REGISTER_FIRE_ENGINE(MyEngine, "my_fire")
 *     (ERF_FireEngineRegistry.H), and give the engine directory a
 *     Make.package (GNU make, FIRE_HOME) and an ERFFireEngine.cmake defining
 *     erf_fire_engine_add(target) (CMake, ERF_FIRE_ENGINE_DIRS).
 *     Tests/FireEngineExample is a complete engine directory to copy.
 *  4. Run the FireCoupler unit tests with the new engine in place of their
 *     stand-in, and a restart of a coupled run against a straight run: the
 *     atmosphere and fire plotfiles must agree at zero tolerance.
 *
 * What the coupling layer already does for every engine: it builds the fire
 * grid over erf.fire.anchor_level, hands the engine the atmosphere's state
 * each step, averages the engine's fluxes onto the atmospheric columns, holds
 * them one step, injects heat, vapour and smoke, and saves those lagged fluxes
 * in the checkpoint. An engine owns only its fire physics and its own state.
 *
 * The call order is the one documented on FireEngine (ERF_FireEngine.H).
 * ERF_FireCoupler.cpp checks at compile time that this class declares every
 * method FireEngine requires (the pure virtual ones), and override checks the
 * signature of every method declared here. The optional methods are listed by
 * hand: a new optional method in FireEngine is added here too.
 */
class FireEngineTemplate : public FireEngine
{
public:
    /// The name erf.fire.engine selects this engine by, for example "my_fire".
    const char* name () const override;

    /**
     * Check the engine's own inputs against the ERF configuration, once at
     * start-up, before the fire grid is built. Abort naming the input on
     * anything the engine cannot run (a missing fuel file, a time step its
     * front tracking cannot take, a grid ratio it does not support).
     */
    void verify_prerequisites (const ERF& erf, int lev) override;

    /**
     * Allocate the engine's fields on the fire grid fg (fg.ba, fg.dm, fg.geom;
     * fg.atm_lo is the fire region's corner in the level's index space) and
     * set the initial fire: fuels, terrain slope, ignition. surface_layer is
     * never null; ERF passes z_phys_nd_atm too.
     */
    void initialize (const ERF& erf,
                     int lev,
                     const FireGrid& fg,
                     const SurfaceLayer* surface_layer,
                     const amrex::MultiFab* z_phys_nd_atm) override;

    /**
     * Called after initialize(), once the run's end is known: anything that
     * depends on stop_time or max_step (output schedules, a check that the
     * ignition happens before the run ends). On a restart, read_checkpoint()
     * comes after this.
     */
    void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step) override;   // optional

    /**
     * Advance the fire over one step of the fire grid's level (forcing.dt from
     * forcing.time). forcing holds the face velocities, the cell-centre
     * heights, the lowest cell's temperature and humidity, and the surface
     * rain if wants_surface_precip() is true; sample them onto the fire grid
     * here and keep no pointer to them after this call.
     */
    void advance (const FireForcing& forcing) override;

    /**
     * Point fluxes.sensible (and fluxes.latent, if the fuel releases water) at
     * the heat flux into the air on the fire grid [W/m2], positive upward,
     * averaged over exactly [forcing.time, forcing.time + forcing.dt] of the
     * last advance() (sum over any substeps, divided by forcing.dt); leave a
     * pointer nullptr for none. A vapour mass flux becomes a latent flux by
     * multiplying by L_v (ERF_Constants.H). The fields must stay valid until the next
     * advance(). Called only when the coupling injects heat or smoke.
     */
    void fill_fluxes (FireFluxes& fluxes) override;

    /**
     * Heat released per kilogram of fuel burned [J/kg], matching the sensible
     * flux last handed over; the coupling divides by it to get the fuel burnt
     * that the smoke emission factor applies to.
     */
    amrex::Real smoke_heat_per_kg () const override;

    /// Optional, false unless declared: true if advance() needs forcing.precip_accum_k0.
    bool wants_surface_precip () const override;

    /// Optional, true unless declared: false if the engine reads nothing from the
    /// surface layer, so it may run over any bottom boundary.
    bool needs_surface_layer () const override;

    /**
     * Optional, nullptr unless declared: for erf.fire.heat_open_fraction, the
     * fraction of each atmospheric column open to the air [0-1] and the mean
     * roof height of the rest [m], on the level's k = 0 slabs with one ghost cell.
     */
    const amrex::MultiFab* heat_open_fraction_atm () const override;
    const amrex::MultiFab* roof_height_atm () const override;

    /**
     * Optional (the default writes nothing and the engine cannot restart).
     * Write everything the engine carries from one step to the next under
     * checkpointname/Level_<lev>/ (lev the fire grid's level). The coupling
     * writes its lagged fluxes after this; use names that do not start with
     * FireQ.
     */
    void write_checkpoint (const std::string& checkpointname) const override;

    /**
     * Optional (the default stops the run, saying restart is not supported).
     * Read back what write_checkpoint() wrote, on a restart, after start_run()
     * and before the first step; step is the fire grid level's step count at
     * the checkpoint. Abort if the checkpoint's fire grid is not this run's.
     */
    void read_checkpoint (const std::string& restart_chkfile, int max_level, int step) override;

    /// Optional (the default writes nothing): the engine's plotfile of the fire grid,
    /// named plotfile_prefix plus the step.
    void write_plotfile (const std::string& plotfile_prefix, amrex::Real time, int step) const override;

    /**
     * Optional, none unless declared: the fire fields the ERF plotfile can
     * carry (names starting with "fire_", listed in erf.plot_vars_*), known
     * from the inputs alone (asked before initialize()), and each field on the
     * fire grid; the coupler writes its column means at k = 0.
     */
    std::vector<std::string> plot_var_names () const override;
    const amrex::MultiFab* plot_var (const std::string& name) const override;

    /// Optional, -1 unless declared: area burned so far [m2], printed with erf.fire.fire_debug.
    amrex::Real burned_area () const override;

    /**
     * Optional, nullptr (or 0) unless declared: what other surface models read
     * from the fire (a dust or smoke source model, for example). A level set negative where burned;
     * the wind (u, v) the engine sampled at reference_wind_height() above
     * ground; and the heat released before any split into sensible and latent.
     */
    const amrex::MultiFab* level_set () const override;
    const amrex::MultiFab* reference_wind () const override;
    amrex::Real reference_wind_height () const override;
    const amrex::MultiFab* heat_release () const override;
};

#endif
