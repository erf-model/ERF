#ifndef ERF_FIRE_COUPLER_H
#define ERF_FIRE_COUPLER_H

#include <memory>
#include <string>
#include <vector>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <ERF_FireCouplingParams.H>
#include <ERF_FireEngine.H>
#include <ERF_FireGrid.H>
#include <ERF_SurfacePrecipitation.H>

class ERF;
class SurfaceLayer;

/**
 * @file ERF_FireCoupler.H
 * @brief The atmosphere side of the fire coupling, shared by every fire engine
 *
 * ERF holds one FireCoupler (std::unique_ptr<FireCoupler> m_fire_coupler),
 * which owns the fire engine selected by erf.fire.engine. The coupler:
 * - checks the ERF configuration a fire grid needs and builds that grid on
 *   the level it refines (erf.fire.anchor_level);
 * - each step of that level, hands the engine the atmospheric state
 *   (FireForcing) and advances it; with lagged or synchronous coupling, or
 *   with the smoke tracer, it averages the fluxes the engine reports
 *   (FireFluxes) onto the atmospheric columns, where they wait for the next
 *   step (the one-step explicit lag of WRF-Fire);
 * - at every Runge-Kutta stage of that next step, injects the lagged heat and
 *   vapour into the cell sources with the WRF-Fire exponential profile
 *   (ERF_FireAtmCoupling.H; lagged and synchronous coupling only) and the
 *   smoke tracer (ERF_FireSmokeEmission.H; erf.fire.smoke_enable);
 * - holds the lagged fluxes for the boundary-layer scheme and the checkpoint.
 *
 * Coarser levels see the fire through average-down, finer levels not at all.
 */
class FireCoupler
{
public:
    /**
     * @param[in] params  The coupling inputs (FireCouplingParams::read())
     * @param[in] engine  The fire engine erf.fire.engine selects
     */
    FireCoupler (const FireCouplingParams& params, std::unique_ptr<FireEngine> engine);

    /**
     * @brief Check the configuration, build the fire grid and set the engine up
     *
     * Called once from ERF::InitData_post() after every surface layer exists.
     * Runs the coupling's checks (verify_fire_prerequisites), then the
     * engine's, builds the fire grid on level lev, allocates the lagged flux
     * buffers (zero) and initialises the engine.
     *
     * @param[in] erf            ERF instance (public AmrCore accessors only)
     * @param[in] lev            The level the fire grid refines (fire_anchor_level)
     * @param[in] surface_layer  The bottom surface layer, or nullptr (aborts when
     *                           the engine needs one, FireEngine::needs_surface_layer)
     * @param[in] z_phys_nd_atm  Nodal height of level lev [m] (ERF always passes it)
     */
    void initialize (const ERF& erf,
                     int lev,
                     const SurfaceLayer* surface_layer,
                     const amrex::MultiFab* z_phys_nd_atm);

    /**
     * @brief Lay the coupling out on a fire grid: what initialize() does after
     *        its start-up checks and before it sets the engine up
     *
     * Takes the fire grid and allocates the lagged flux buffers (zero) on the
     * k = 0 slabs of ba_atm, with 5 ghost columns. Public so that a test can
     * drive the coupling with an engine of its own without an ERF object
     * (Tests/Unit/Fire/ERF_GTestFireCoupler.cpp).
     *
     * @param[in] lev       The level the fire grid refines
     * @param[in] fg        The fire grid (create_fire_grid)
     * @param[in] ba_atm    The level's boxes
     * @param[in] dm_atm    Their owners
     * @param[in] geom_atm  The level's geometry
     */
    void set_fire_grid (int lev, const FireGrid& fg,
                        const amrex::BoxArray& ba_atm,
                        const amrex::DistributionMapping& dm_atm,
                        const amrex::Geometry& geom_atm);

    /**
     * @brief Advance the fire one step of its level and refill the lagged fluxes
     *
     * Called from ERF::Advance() once per step of the fire grid's level, after
     * that level's dycore. The air temperature and humidity at k = 0 come from
     * the state before the dycore; the wind from before the dycore (passive and
     * lagged coupling) or after it (synchronous).
     *
     * @param[in] time            Start of the step [s]
     * @param[in] dt              Step of the level [s]
     * @param[in] S_old           Conserved state before the dycore
     * @param[in] xvel_old, yvel_old  Face velocities before the dycore [m/s]
     * @param[in] xvel_new, yvel_new  Face velocities after the dycore [m/s]
     * @param[in] z_phys_cc       Cell-centre height of the level [m]
     * @param[in] ba2d            The level's boxes coarsened to one cell in z
     * @param[in] geom            Geometry of the level
     * @param[in] precip_sources  The microphysics surface precipitation accumulators of the level
     */
    void advance (amrex::Real time, amrex::Real dt,
                  const amrex::MultiFab& S_old,
                  const amrex::MultiFab& xvel_old, const amrex::MultiFab& yvel_old,
                  const amrex::MultiFab& xvel_new, const amrex::MultiFab& yvel_new,
                  const amrex::MultiFab& z_phys_cc,
                  const amrex::BoxArray& ba2d,
                  const amrex::Geometry& geom,
                  const SurfacePrecipAccumulationSources& precip_sources);

    /**
     * @brief Inject the lagged fire heat, vapour and smoke into the cell sources
     *
     * Called from the slow right-hand side at every Runge-Kutta stage on the
     * fire grid's level, after make_sources has rebuilt cc_source. The lagged
     * fluxes are fixed over the step, so every stage receives the same
     * tendency (erf.fire.source_mode decides whether it replaces or adds to the
     * theta and vapour sources). Heat and vapour only with lagged or
     * synchronous coupling; smoke with erf.fire.smoke_enable.
     *
     * @param[in,out] cc_source     Cell-centred sources of the level
     * @param[in]     S_old         Conserved state at the start of the step (density)
     * @param[in]     z_phys_cc     Cell-centre height [m]
     * @param[in]     z_phys_nd     Nodal height [m] (ERF always passes it): the cell
     *                              thickness and roof datum of the heat, and the
     *                              smoke cell's thickness
     * @param[in]     geom          Geometry of the level
     * @param[in]     has_moisture  The vapour source exists
     * @param[in]     step          Step of the level, for the debug output
     */
    void add_sources (amrex::MultiFab& cc_source,
                      const amrex::MultiFab& S_old,
                      const amrex::MultiFab& z_phys_cc,
                      const amrex::MultiFab* z_phys_nd,
                      const amrex::Geometry& geom,
                      bool has_moisture,
                      int step);

    /**
     * @brief Write the fire to a checkpoint: the engine's state
     *        (FireEngine::write_checkpoint), then the coupling's lagged fluxes
     *        (FireQAtmPrev and FireQLatAtmPrev) under checkpointname/Level_<lev>/
     */
    void write_checkpoint (const std::string& checkpointname) const;

    /**
     * @brief Why a checkpoint's lagged fluxes do not belong to this run's fire
     *        grid, or an empty string when they do
     *
     * @param[in] chk_lev     The level the checkpoint holds them on
     * @param[in] chk_region  The columns they cover (minimal box of their boxes)
     */
    std::string checkpoint_mismatch (int chk_lev, const amrex::Box& chk_region) const;

    /**
     * @brief Restore the engine's state, then the lagged fluxes, on a restart
     *
     * Aborts when the checkpoint holds the lagged fluxes on another level or
     * over another region than this run's fire grid. A checkpoint without them
     * (no fire state, or written before they were checkpointed) leaves them at
     * zero: the first restarted step injects no fire heat.
     *
     * @param[in] restart_chkfile  The checkpoint directory
     * @param[in] max_level        ERF's max_level
     * @param[in] step             The step count of the fire grid's level at the checkpoint
     */
    void read_checkpoint (const std::string& restart_chkfile, int max_level, int step);

    /// Fire fields the ERF plotfile can carry (FireEngine::plot_var_names(), read at construction,
    /// before ERF reads erf.plot_vars_*).
    const std::vector<std::string>& plot_var_names () const { return m_plot_var_names; }

    /**
     * @brief Fill component dcomp of an ERF plotfile MultiFab of level lev
     *        with the fire field name
     *
     * On the fire grid's level, the k = 0 cell of each column holds the mean
     * of the engine's field (FireEngine::plot_var) over the fire cells under
     * it; every other cell, every other level and a field the engine returns
     * as nullptr hold zero.
     *
     * @param[in]  name  One of plot_var_names()
     * @param[out] dst   Cell-centred MultiFab on the boxes of level lev
     * @param[in]  dcomp Component of dst to fill
     * @param[in]  lev   The level dst belongs to
     */
    void fill_plot_var (const std::string& name, amrex::MultiFab& dst, int dcomp, int lev) const;

    /// The level the fire grid refines (erf.fire.anchor_level, resolved at initialize).
    int level () const { return m_lev; }

    /// The coupling inputs.
    const FireCouplingParams& params () const { return m_params; }

    /// The fire engine.
    FireEngine&       engine ()       { return *m_engine; }
    const FireEngine& engine () const { return *m_engine; }

    /// The fire grid.
    const FireGrid& fire_grid () const { return m_fg; }

    /**
     * @brief Sensible heat flux [W/m2] of the previous fire step on the
     *        atmospheric columns of the fire grid's level, with 5 ghost columns
     *        (for a reader that samples a halo around each tile, such as a
     *        boundary-layer scheme taking the fire's heat into its mixing); zero
     *        outside the fire region. Allocated at initialize, so never null afterwards.
     */
    const amrex::MultiFab* get_Q_atm_prev () const { return m_Q_atm_prev.get(); }
    /// Latent heat flux [W/m2] of the previous fire step, laid out as get_Q_atm_prev().
    const amrex::MultiFab* get_Q_lat_atm_prev () const { return m_Q_lat_atm_prev.get(); }

private:
    /// Average the engine's fluxes onto the atmospheric columns for the next step.
    void update_flux_buffer (const amrex::Geometry& geom_atm);

    FireCouplingParams          m_params;
    std::unique_ptr<FireEngine> m_engine;
    FireGrid                    m_fg;
    int                         m_lev = 0;
    amrex::Geometry             m_geom_atm;   ///< geometry of the fire grid's level

    // The engine's fluxes averaged onto the atmospheric columns of the fire
    // grid's level (one slab per box of that level), held from the previous
    // step for injection at every stage of the next (one-step explicit lag).
    // Their ghost columns hold the neighbouring box's value for a reader that
    // samples a halo around every tile; outside the refined region of a finer
    // level they hold zero.
    std::unique_ptr<amrex::MultiFab> m_Q_atm_prev;      ///< sensible [W/m2]
    std::unique_ptr<amrex::MultiFab> m_Q_lat_atm_prev;  ///< latent   [W/m2]

    std::vector<std::string> m_plot_var_names;  ///< FireEngine::plot_var_names(), read at construction
};

#endif
