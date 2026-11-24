#ifndef ERF_FLOOD_ENGINE_TEMPLATE_H
#define ERF_FLOOD_ENGINE_TEMPLATE_H

#include <ERF_FloodEngine.H>

/**
 * @file ERF_FloodEngineTemplate.H
 * @brief A starting point for a new flood engine: declarations only
 *
 * This class declares every FloodEngine method a new flood model provides.
 * Nothing here has a body, no input selects it and nothing creates one, so it
 * adds no code to ERF. Tests/FloodEngineExample is the smallest engine that
 * works, and a complete engine directory to copy.
 *
 * To start a new engine:
 *  1. Copy this file into the new model's engine directory (outside ERF),
 *     rename the class and the include guard, and write the bodies in a .cpp
 *     file there. Delete the declarations of the optional methods the engine
 *     does not provide (their defaults in FloodEngine then apply); a
 *     declaration left without a body fails to link.
 *  2. Read the model's inputs under its own prefix (flood.<engine>.* is
 *     recommended), so they never mix with the coupling's keys
 *     (FloodCouplingParams, erf.flood.*) or another engine's.
 *  3. Register the engine in that .cpp, at namespace scope:
 *         ERF_REGISTER_FLOOD_ENGINE(MyEngine, "my_flood")
 *     (ERF_FloodEngineRegistry.H), and give the engine directory a
 *     Make.package (GNU make, FLOOD_HOME) and an ERFFloodEngine.cmake defining
 *     erf_flood_engine_add(target) (CMake, ERF_FLOOD_ENGINE_DIRS).
 *  4. Run a restart of a coupled run against a straight run: the atmosphere
 *     and flood plotfiles must agree at zero tolerance.
 *
 * The call order is the one documented on FloodEngine (ERF_FloodEngine.H).
 * ERF_FloodCoupler.cpp checks at compile time that this class declares every
 * method FloodEngine requires (the pure virtual ones), and override checks the
 * signature of every method declared here. The optional methods are listed by
 * hand: a new optional method in FloodEngine is added here too.
 */
class FloodEngineTemplate : public FloodEngine
{
public:
    /// The name erf.flood.engine selects this engine by, for example "my_flood".
    const char* name () const override;

    /**
     * Check the engine's own inputs against the ERF configuration, once at
     * start-up, before the flood grid is built. Abort naming the input on
     * anything the engine cannot run (a missing terrain file, a grid ratio it
     * does not support).
     */
    void verify_prerequisites (const ERF& erf, int lev) override;

    /**
     * Allocate the engine's fields on fg.ba and fg.dm and set its initial
     * state. z_phys_nd_atm is the terrain the atmosphere sees.
     */
    void initialize (const ERF& erf, int lev, const FloodGrid& fg,
                     const amrex::MultiFab* z_phys_nd_atm) override;

    /// Optional: the run's end, once known (stop time, max_step, substeps).
    void start_run (amrex::Real stop_time, int max_step, int steps_per_coarse_step) override;

    /**
     * Advance over [forcing.time, forcing.time + forcing.dt]. forcing.precip_k0
     * holds this step's liquid and frozen precipitation per atmospheric column
     * [kg/m2]; fill_flood_from_atm_k0 (ERF_FloodGrid.H) copies a component onto
     * the flood grid.
     */
    void advance (const FloodForcing& forcing) override;

    /// Optional: false for an engine that reads no precipitation (default true).
    bool wants_surface_precip () const override;

    /// Optional: the engine's state under checkpointname/Level_<lev>/.
    void write_checkpoint (const std::string& checkpointname) const override;

    /// Optional: restore what write_checkpoint wrote (the default aborts).
    void read_checkpoint (const std::string& restart_chkfile, int step) override;

    /// Optional: the engine's own plotfile of the flood grid.
    void write_plotfile (const std::string& plotfile_prefix, amrex::Real time, int step) const override;

    /// Optional: fields for the ERF plotfile, each named "flood_...".
    std::vector<std::string> plot_var_names () const override;

    /// Optional: the field plot_var_names() lists as name, on the flood grid.
    const amrex::MultiFab* plot_var (const std::string& name) const override;
};

#endif
