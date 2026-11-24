#ifndef ERF_FLOOD_COUPLER_H
#define ERF_FLOOD_COUPLER_H

#include <memory>
#include <string>
#include <vector>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <ERF_FloodCouplingParams.H>
#include <ERF_FloodEngine.H>
#include <ERF_FloodGrid.H>
#include <ERF_SurfacePrecipitation.H>

class ERF;

/**
 * @file ERF_FloodCoupler.H
 * @brief The atmosphere side of the flood coupling, shared by every flood engine
 *
 * ERF holds one FloodCoupler (std::unique_ptr<FloodCoupler> m_flood_coupler),
 * which owns the flood engine selected by erf.flood.engine. The coupler:
 * - checks the ERF configuration a flood grid needs and builds that grid on
 *   the level it refines (erf.flood.anchor_level);
 * - keeps the microphysics surface precipitation accumulation it saw last, and
 *   each step of that level hands the engine the change since then (liquid and
 *   frozen), with the wind after the dycore (FloodForcing);
 * - writes the engine's fields into the ERF plotfile on request;
 * - records the flood grid in the checkpoint and refuses a restart onto
 *   another one.
 *
 * The coupling is one-way: nothing returns to the atmosphere.
 */
class FloodCoupler
{
public:
    /**
     * @param[in] params  The coupling inputs (FloodCouplingParams::read())
     * @param[in] engine  The flood engine erf.flood.engine selects
     *
     * Aborts when the engine lists a plot field whose name does not start with
     * "flood_" or lists a name twice.
     */
    FloodCoupler (const FloodCouplingParams& params, std::unique_ptr<FloodEngine> engine);

    /**
     * @brief Check the configuration, build the flood grid and set the engine up
     *
     * Called once from ERF::InitData_post(), on a fresh start and on a restart
     * alike. Runs the coupling's checks (verify_flood_prerequisites), then the
     * engine's, builds the flood grid on level lev and initializes the engine.
     *
     * @param[in] erf            ERF instance (public AmrCore accessors only)
     * @param[in] lev            The level the flood grid refines (flood_anchor_level)
     * @param[in] z_phys_nd_atm  Nodal height of level lev [m]
     */
    void initialize (const ERF& erf, int lev, const amrex::MultiFab* z_phys_nd_atm);

    /**
     * @brief Lay the coupling out on a flood grid: what initialize() does after
     *        its start-up checks and before it sets the engine up
     *
     * Takes the flood grid and the k = 0 slabs of ba_atm. Public so that a test
     * can drive the coupling with an engine of its own without an ERF object.
     *
     * @param[in] lev       The level the flood grid refines
     * @param[in] fg        The flood grid (create_flood_grid)
     * @param[in] ba_atm    The level's boxes, each spanning the level's height
     * @param[in] dm_atm    Their owners
     * @param[in] geom_atm  The level's geometry
     */
    void set_flood_grid (int lev, const FloodGrid& fg,
                         const amrex::BoxArray& ba_atm,
                         const amrex::DistributionMapping& dm_atm,
                         const amrex::Geometry& geom_atm);

    /**
     * @brief Take the current surface precipitation accumulation as the
     *        reference the next step's precipitation is measured from
     *
     * Called once at start-up, after initialize(), when the engine reads the
     * precipitation (FloodEngine::wants_surface_precip). On a fresh start the
     * accumulators are zero; on a restart they hold the checkpoint's values, so
     * the first restarted step hands the engine the same precipitation as the
     * straight run.
     *
     * @param[in] precip_sources  The microphysics surface precipitation accumulators of the level
     */
    void set_precip_reference (const SurfacePrecipAccumulationSources& precip_sources);

    /**
     * @brief Advance the flood one step of its level
     *
     * Called from ERF::Advance() once per step of the flood grid's level, after
     * that level's dycore and microphysics.
     *
     * @param[in] time            Start of the step [s]
     * @param[in] dt              Step of the level [s], ERF's double dt[lev]: fields_time()
     *                            becomes time + dt, the time ERF assigns to the level
     * @param[in] xvel, yvel      Face velocities after the dycore [m/s]
     * @param[in] z_phys_cc       Cell-center height of the level [m]
     * @param[in] precip_sources  The microphysics surface precipitation accumulators of the level
     */
    void advance (double time, double dt,
                  const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
                  const amrex::MultiFab& z_phys_cc,
                  const SurfacePrecipAccumulationSources& precip_sources);

    /**
     * @brief Write the flood to a checkpoint: the engine's state
     *        (FloodEngine::write_checkpoint), then the flood grid's level,
     *        ratio and region in checkpointname/Level_<lev>/FloodGrid, and the
     *        precipitation reference in Level_<lev>/FloodPrecipRef
     */
    void write_checkpoint (const std::string& checkpointname) const;

    /**
     * @brief Why a checkpoint's flood state does not belong to this run's flood
     *        grid, or an empty string when it does
     *
     * @param[in] chk_lev     The level the checkpoint's flood grid refines
     * @param[in] chk_ratio   Its grid ratio
     * @param[in] chk_region  The atmospheric columns it covers
     */
    std::string checkpoint_mismatch (int chk_lev, int chk_ratio, const amrex::Box& chk_region) const;

    /**
     * @brief Restore the engine's state on a restart
     *
     * A checkpoint written with a flood engine holds Level_<lev>/FloodGrid: the
     * coupler aborts when it names another level, ratio or region than this run,
     * and otherwise calls FloodEngine::read_checkpoint. A checkpoint without it
     * (an atmosphere-only spin-up) leaves the engine in its initial state, said
     * once.
     *
     * @param[in] restart_chkfile  The checkpoint directory
     * @param[in] step             The step count of the flood grid's level at the checkpoint
     */
    void read_checkpoint (const std::string& restart_chkfile, int step);

    /// Flood fields the ERF plotfile can carry (FloodEngine::plot_var_names(), read at construction).
    const std::vector<std::string>& plot_var_names () const { return m_plot_var_names; }

    /**
     * @brief Fill component dcomp of an ERF plotfile MultiFab of level lev
     *        with the flood field name
     *
     * On the flood grid's level, the k = 0 cell of each column holds the mean
     * of the engine's field over the flood cells under it; every other cell,
     * every other level and a field the engine returns as nullptr hold zero.
     */
    void fill_plot_var (const std::string& name, amrex::MultiFab& dst, int dcomp, int lev) const;

    /// The level the flood grid refines (erf.flood.anchor_level, resolved at initialize).
    int level () const { return m_lev; }

    /// The coupling inputs.
    const FloodCouplingParams& params () const { return m_params; }

    /// The flood engine.
    FloodEngine&       engine ()       { return *m_engine; }
    const FloodEngine& engine () const { return *m_engine; }

    /// The flood grid.
    const FloodGrid& flood_grid () const { return m_fg; }

    /**
     * @brief The time the engine's fields belong to [s], in ERF's double time
     *
     * ERF sets it at start-up and on a restart, advance() to the end of the
     * step, and ERF again to the coarser level's time once the finer levels'
     * substeps have reached it.
     */
    double fields_time () const { return m_fields_time; }
    /// Set fields_time() to t [s] (ERF only: start-up, restart, end of a coarser step).
    void set_fields_time (double t) { m_fields_time = t; }

private:
    FloodCouplingParams          m_params;
    std::unique_ptr<FloodEngine> m_engine;
    FloodGrid                    m_fg;
    int                          m_lev = 0;
    bool                         m_grid_set = false;
    amrex::Geometry              m_geom_atm;   ///< geometry of the flood grid's level
    amrex::BoxArray              m_ba_cols;    ///< the k = 0 slab of each box of the level
    amrex::DistributionMapping   m_dm_cols;    ///< their owners (the level's)
    double                       m_fields_time = 0.0;

    /// Liquid (0) and frozen (1) accumulation [kg/m2] at the end of the last
    /// step, on m_ba_cols; null until set_precip_reference()
    std::unique_ptr<amrex::MultiFab> m_precip_prev;

    std::vector<std::string> m_plot_var_names;
};

#endif
