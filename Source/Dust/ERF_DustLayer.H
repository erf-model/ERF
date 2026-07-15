/**
 * @file ERF_DustLayer.H
 * @brief Container class for the ERF-Dust 2D surface emission layer.
 *
 * DustLayer owns all dust-state MultiFabs on the 2D dust grid and
 * provides the initialize() and advance() interface to the ERF solver.
 *
 * In Phase 2, initialize() allocates MultiFabs and sets them to their
 * default values. advance() is a stub that returns without computation.
 * Subsequent phases replace the stub advance() with physics.
 *
 * The dust grid is a 2D slab (k-extent = 1) at k=0 of the atmospheric
 * domain, optionally refined by DustParams::grid_ratio in x and y.
 *
 * MultiFab layout (all on DustGrid, grid_ratio*nx x grid_ratio*ny x 1):
 *
 * | Name                 | ncomp           | Description                             |
 * |----------------------|-----------------|-----------------------------------------|
 * | dust_ustar_t         | 1               | Threshold friction velocity [m/s]       |
 * | dust_soil_type       | 1               | Soil type index [integer stored as Real]|
 * | dust_silt_fraction   | 1               | Surface silt mass fraction [-]          |
 * | dust_crust_index     | 1               | Mineral crust strength index [0,1]      |
 * | dust_moisture_flag   | 1               | Surface moisture inhibition flag [0,1]  |
 * | dust_suppression     | 1               | Suppression agent coverage fraction [-] |
 * | dust_emission_flux   | n_size_bins     | Emission flux per size bin [kg/m^2/s]  |
 * | dust_efflor          | 1               | Efflorescence fraction [0,1]            |
 * | dust_ustar_base      | 1               | Base u*_t before chemistry [m/s]        |
 * | dust_ustar_in        | 1               | Friction velocity u* [m/s]              |
 *
 * dust_soil_type stores an integer code as amrex::Real. Codes follow the
 * STATSGO classification for natural soils (1-16) and reserve codes
 * >= 100 for mine-specific surface types (defined in Phase 3).
 *
 * References:
 *   Marticorena, B., & Bergametti, G. (1995). J. Geophys. Res., 100, 16415.
 *     https://doi.org/10.1029/95JD00690
 *   Shao, Y., & Lu, H. (2000). J. Geophys. Res., 105, 22437.
 *     https://doi.org/10.1029/2000JD900304
 *   Bagnold, R. A. (1941). The Physics of Blown Sand and Desert Dunes.
 *     Methuen, London.
 */

#ifndef ERF_DUST_LAYER_H
#define ERF_DUST_LAYER_H

#include <ERF_DustGrid.H>
#include <ERF_DustParams.H>
#include <memory>
#include <AMReX_MultiFab.H>

#ifdef ERF_USE_DUST
#include <ERF_DustBlastSchedule.H>
#include <ERF_DustRoadSchedule.H>
#include <ERF_DustSettling.H>
#include <ERF_DustDeposition.H>
#include <ERF_DustAtmReturn.H>
#include <ERF_DustPlotfile.H>
#include <ERF_DustStatsOutput.H>
#include <ERF_DustPM.H>
#include <ERF_DustNAAQSOutput.H>
#include <ERF_DustSiteRegistry.H>
#include <ERF_DustPHREEQCWriter.H>
#include <ERF_DustCriticalMaterials.H>
#endif

#if defined(ERF_USE_DUST) && defined(ERF_USE_PARTICLES)
#include <ERF_DustPC.H>
#endif


// Forward declarations
class SurfaceLayer;
class ERF;

class DustLayer
{
public:
    DustLayer()  = default;
    ~DustLayer() = default;

    /**
     * @brief Initialize the dust layer after SurfaceLayer is ready.
     *
     * Performs prerequisite verification, creates the 2D dust grid,
     * allocates all MultiFabs listed in the class-level table, and
     * fills them with default values from dust_params.
     *
     * Called once from ERF::InitData_post() when erf.dust.enable = true.
     *
     * @param[in] erf              ERF instance (public AmrCore accessors only).
     * @param[in] surface_layer    Pointer to the initialized SurfaceLayer.
     * @param[in] dust_params      Dust configuration parameters.
     */
    void initialize(const ERF&          erf,
                    const SurfaceLayer* surface_layer,
                    const DustParams&   dust_params);

    /**
     * @brief Advance the dust layer by one atmospheric timestep.
     *
     * Phases 5-13 replace the stub with emission, coupling, and deposition physics.
    * Phase 19 adds Lagrangian particle advection (requires xvel, yvel, zvel, geom_atm).
    *
    * @param[in] dt                 Atmospheric timestep [s].
    * @param[in] dust_params        Dust configuration parameters.
    * @param[in] surface_layer      Pointer to SurfaceLayer (nullptr = placeholder fallback).
    * @param[in] xvel_mf            3D u velocity field (nullptr = placeholder fallback).
    * @param[in] yvel_mf            3D v velocity field (nullptr = placeholder fallback).
    * @param[in] zvel_mf            3D w velocity field (nullptr = placeholder fallback). Phase 19.
    * @param[in] z_phys_cc_mf       3D terrain height field (nullptr = placeholder fallback).
    * @param[in] geom_atm           Atmospheric geometry. Phase 19.
    * @param[in] nz                 Number of vertical levels (0 = placeholder fallback).
    */
    void advance(amrex::Real            dt,
                 const DustParams&      dust_params,
                 SurfaceLayer*          surface_layer = nullptr,   // non-const: getters not const
                 const amrex::MultiFab* xvel_mf = nullptr,
                 const amrex::MultiFab* yvel_mf = nullptr,
                 const amrex::MultiFab* zvel_mf = nullptr,
                 const amrex::MultiFab* z_phys_cc_mf = nullptr,
                 const amrex::Geometry* geom_atm = nullptr,
                 int                    nz = 0);

    // Const accessors for plotfile and diagnostic output
    const amrex::MultiFab* get_ustar_t()         const { return dust_ustar_t.get(); }
    const amrex::MultiFab* get_soil_type()        const { return dust_soil_type.get(); }
    const amrex::MultiFab* get_silt_fraction()    const { return dust_silt_fraction.get(); }
    const amrex::MultiFab* get_crust_index()      const { return dust_crust_index.get(); }
    const amrex::MultiFab* get_moisture_flag()    const { return dust_moisture_flag.get(); }
    const amrex::MultiFab* get_suppression()      const { return dust_suppression.get(); }
    const amrex::MultiFab* get_efflor()           const { return dust_efflor.get(); }
    const amrex::MultiFab* get_emission_flux()    const { return dust_emission_flux.get(); }
    const amrex::MultiFab* get_surf_moist()       const { return dust_surf_moist.get(); }
    const amrex::MultiFab* get_ustar_in()         const { return dust_ustar_in.get(); }
    const amrex::MultiFab* get_retreat_flag()     const { return dust_retreat_flag.get(); }
    const amrex::MultiFab* get_wind_ref()         const { return dust_wind_ref.get(); }
    const amrex::MultiFab* get_pblh()             const { return dust_pblh.get(); }
    const amrex::MultiFab* get_tsfc()             const { return dust_tsfc.get(); }
    const amrex::MultiFab* get_deposition_rate()  const { return dust_deposition_rate.get(); }
    const amrex::MultiFab* get_conc_sfc()         const { return dust_conc_sfc.get(); }
    const amrex::MultiFab* get_pm25()             const { return dust_pm25.get(); }
    const amrex::MultiFab* get_pm10()             const { return dust_pm10.get(); }
    const amrex::MultiFab* get_pm25_24h()         const { return dust_pm25_24h.get(); }
    const amrex::MultiFab* get_pm10_24h()         const { return dust_pm10_24h.get(); }
    const amrex::MultiFab* get_pm25_exceed()      const { return dust_pm25_exceed.get(); }
    const amrex::MultiFab* get_pm10_exceed()      const { return dust_pm10_exceed.get(); }
    const amrex::MultiFab* get_msha_dose()        const { return dust_msha_dose.get(); }
    const amrex::MultiFab* get_msha_twa()         const { return dust_msha_twa.get(); }
    const amrex::MultiFab* get_msha_exceed()      const { return dust_msha_exceed.get(); }
    const amrex::MultiFab* get_msha_shift_twa()   const { return dust_msha_shift_twa.get(); }
    /**
     * @brief Get the dust source map (Phase 19).
     * @return Pointer to dust_source_map, nullptr if not initialized.
     */
    const amrex::MultiFab* get_source_map()        const { return dust_source_map.get(); }

    const amrex::MultiFab* get_site_id()           const { return dust_site_id.get(); }

    const amrex::MultiFab* get_cm_flux()           const { return dust_cm_flux.get(); }

    const DustGrid&        get_dust_grid()        const { return m_dg; }
    const amrex::Geometry& get_dust_geom()        const { return m_dg.geom; }
    const DustParams&      get_params()           const { return m_params; }

    /**
     * @brief Get the dust scalar component index.
     * @return Component index in cc_source for dust mass density.
     */
    int get_dust_scalar_comp() const { return m_dust_scalar_comp; }

    /**
     * @brief Inject dust emission flux into the ERF cc_source MultiFab.
     *
     * Called from ERF_Advance.cpp before advance_dycore(), following the same
     * call site as the fire module's apply_fire_tendency_to_cc_source.
     *
     * @param[in,out] cc_source    Cell-centred source MultiFab on atm grid.
     * @param[in]     z_phys_cc    Cell-centre terrain height [m].
     * @param[in]     geom_atm     Atmospheric geometry at level 0.
     */
    void apply_to_cc_source(amrex::MultiFab&       cc_source,
                            const amrex::MultiFab& z_phys_cc,
                            const amrex::Geometry& geom_atm);

    /**
     * @brief Apply Stokes settling tendency to cc_source for all active bins.
     *
     * Called from ERF_Advance.cpp at the same site as apply_to_cc_source.
     * When transport_bins_separately = false, uses bin 0 diameter for the
     * single scalar. When true, applies per-bin settling to each slot.
     *
     * @param[in,out] cc_source   Cell-centred source MultiFab.
     * @param[in]     S_old       Conserved state at time n.
     * @param[in]     z_phys_cc   Cell-centre terrain height [m].
     * @param[in]     geom_atm    Atmospheric geometry.
     */
    void apply_settling_to_cc_source(amrex::MultiFab&       cc_source,
                                      const amrex::MultiFab& S_old,
                                      const amrex::MultiFab& z_phys_cc,
                                      const amrex::Geometry& geom_atm);

    /**
     * @brief Apply dry deposition BC at k=klo and accumulate deposited mass.
     *
     * Called from ERF_TI_slow_rhs_pre.H after apply_settling_to_cc_source (Phase 11),
     * following the same call-site pattern as apply_to_cc_source (Phase 10).
     * Uses SurfaceLayer u* extracted in Phase 9 (dust_ustar_in).
     *
     * @param[in,out] cc_source  Cell-centred source MultiFab on atm grid.
     * @param[in]     S_old      Conserved state at time n.
     * @param[in]     z_phys_cc  Cell-centre terrain height [m].
     * @param[in]     geom_atm   Atmospheric geometry.
     * @param[in]     dt         Timestep [s].
     */
    void apply_deposition_bc(amrex::MultiFab&       cc_source,
                             const amrex::MultiFab& S_old,
                             const amrex::MultiFab& z_phys_cc,
                             const amrex::Geometry& geom_atm,
                             amrex::Real            dt);

    /**
     * @brief Extract return fields from the 3D solver to the 2D dust grid.
     *
     * Called from ERF_Advance.cpp once per timestep after erf_slow_rhs_post
     * and before the next DustLayer::advance() call. Extracts:
     *   - RhoDust at k=klo -> dust_conc_sfc (loading feedback for Phase 5)
     *   - Q1fx3 at z-face klo -> dust_surf_moist (moisture inhibition, Phase 5)
     *   - Prints dust_deposition_rate domain integral (Phase 12 confirmation)
     *
     * Pattern: follows SurfaceLayer extraction in Phase 9 (ERF_DustWindExtract.H),
     * which adapts fill_fire_wind_from_interpolation from Source/Fire/.
     *
     * @param[in] S_new_cons  Conserved state MultiFab after dycore advance (level 0).
     * @param[in] Q1fx3       Vertical moisture flux (face-staggered z). Null when
     *                        moisture_type == None. From SFS_q1fx3_lev[0].get().
     * @param[in] geom_atm    Atmospheric geometry at level 0.
     */
    void extract_atm_return_fields(const amrex::MultiFab& S_new_cons,
                                    const amrex::MultiFab* Q1fx3,
                                    const amrex::Geometry& geom_atm);

    /**
     * @brief Compute PM2.5/PM10 concentrations, running averages, and exceedance flags.
     *
     * Called from advance() each timestep after extract_atm_return_fields (Phase 13).
     * Uses dust_conc_sfc from Phase 13 as the concentration source.
     * Updates running 24-hour averages and exceedance flags.
     * Appends one row to dust_naaqs.csv.
     *
     * @param dt       Timestep [s].
     * @param cur_time Simulation time [s].
     * @param nstep    Current step number.
     */
    void compute_naaqs_diagnostics(amrex::Real dt, amrex::Real cur_time, int nstep);

    /**
     * @brief Compute MSHA worker exposure: dose accumulator, TWA, and exceedance flag.
     *
     * Called from advance() each timestep after compute_naaqs_diagnostics (Phase 18).
     * Updates dose and TWA accumulators, checks PEL exceedance, handles shift resets,
     * writes per-step and per-shift CSV outputs, and samples receptor points.
     *
     * @param dt       Timestep [s].
     * @param cur_time Simulation time [s].
     * @param nstep    Current step number.
     */
    void compute_msha_exposure(amrex::Real dt, amrex::Real cur_time, int nstep);

    /**
     * @brief Release and advance dust super-particles one timestep (Phase 19).
     *
     * Called from advance() when enable_particles=true.
     * Releases one particle per emission cell with flux > 0.
     * Advects with nearest-cell interpolated ERF face velocities + Stokes settling.
     * Deposits when z < z_lo + 0.5*dz_atm.
     *
     * @param xvel,yvel,zvel  ERF level-0 face velocities.
     * @param geom_atm        Atmospheric geometry.
     * @param dt              Timestep [s].
     * @param nstep           Current step number.
     */
    void advance_particles(const amrex::MultiFab& xvel,
                           const amrex::MultiFab& yvel,
                           const amrex::MultiFab& zvel,
                           const amrex::Geometry& geom_atm,
                           amrex::Real dt, int nstep);


    /**
     * @brief Write dust plotfile and CSV diagnostics.
     *
     * Called from ERF::WriteAtIntermediateTime and ERF::WriteAtFinalTime.
     * Plotfile uses VisMF::Write + WriteGenericPlotfileHeader on the native
     * dust grid (pattern: ERF_FirePlotfile.cpp).
     * CSV diagnostics use append pattern (pattern: ERF_FireStatsOutput.H).
     *
     * Final-step guard: is_final=true and nstep > m_last_dust_plot_step
     * ensures a write at simulation end regardless of dust_plot_int.
     *
     * @param[in] nstep     Current step number.
     * @param[in] cur_time  Current simulation time [s].
     * @param[in] is_final  True when called from WriteAtFinalTime.
     */
    void write_output(int nstep, double cur_time, bool is_final = false);

    /**
     * @brief Write PHREEQC deposition feedback files if interval has elapsed.
     *
     * Called from write_output() each timestep. Writes:
     *   (A) phreeqc_feedback_file: full grid deposition [kg/m^2] (overwritten).
     *   (B) phreeqc_site_summary_file: per-site totals [kg] (appended).
     *
     * Triggered when:
     *   phreeqc_feedback_interval_s > 0 AND
     *   (cur_time - m_last_phreeqc_write_time) >= phreeqc_feedback_interval_s
     *   OR is_final = true AND nstep > m_last_phreeqc_write_step.
     *
     * Pattern: WriteAtFinalTime guard in Source/ERF.cpp.
     * Fire reference: ERF_FirePrerequisites.H/cpp.
     *
     * @param nstep     Current step.
     * @param cur_time  Simulation time [s].
     * @param is_final  True when called from WriteAtFinalTime.
     */
    void write_phreeqc_feedback(int nstep, amrex::Real cur_time,
                                 bool is_final = false);

private:
    DustGrid                              m_dg;

    std::unique_ptr<amrex::MultiFab>      dust_ustar_t;
    std::unique_ptr<amrex::MultiFab>      dust_soil_type;
    std::unique_ptr<amrex::MultiFab>      dust_silt_fraction;
    std::unique_ptr<amrex::MultiFab>      dust_crust_index;
    std::unique_ptr<amrex::MultiFab>      dust_moisture_flag;
    std::unique_ptr<amrex::MultiFab>      dust_suppression;
    std::unique_ptr<amrex::MultiFab>      dust_emission_flux;
    std::unique_ptr<amrex::MultiFab>      dust_efflor;
    std::unique_ptr<amrex::MultiFab>      dust_ustar_base;
    std::unique_ptr<amrex::MultiFab>      dust_surf_moist;
    std::unique_ptr<amrex::MultiFab>      dust_ustar_in;
    std::unique_ptr<amrex::MultiFab>      dust_retreat_flag;

    // Phase 9: atmospheric fields extracted each timestep
    std::unique_ptr<amrex::MultiFab>      dust_wind_ref;   ///< u,v at zref [m/s] (ncomp=2)
    std::unique_ptr<amrex::MultiFab>      dust_pblh;       ///< MRF PBL height [m]
    std::unique_ptr<amrex::MultiFab>      dust_tsfc;       ///< Surface temperature [K]

    // Phase 10: coarsened emission flux on atmospheric grid for injection.
    std::unique_ptr<amrex::MultiFab>      dust_flux_atm;   ///< Total emission flux coarsened
                                                             ///< to atm level-0 grid [kg/m^2/s].

    // Phase 12: deposition accumulation and flux buffer
    std::unique_ptr<amrex::MultiFab>      dust_deposition_rate; ///< Accumulated deposition [kg/m^2].
    std::unique_ptr<amrex::MultiFab>      dep_flux_atm;          ///< Deposition flux buffer on atm 2D slab.

    // Phase 13: return fields from 3D solver
    std::unique_ptr<amrex::MultiFab>      dust_conc_sfc;   ///< Near-surface dust [kg/m^3] on dust grid.

    // Phase 17: PM size fraction MultiFabs on dust grid.
    std::unique_ptr<amrex::MultiFab>      dust_pm25;       ///< Instantaneous PM2.5 [µg/m^3].
    std::unique_ptr<amrex::MultiFab>      dust_pm10;       ///< Instantaneous PM10 [µg/m^3].
    std::unique_ptr<amrex::MultiFab>      dust_pm25_24h;   ///< 24-hour running average PM2.5.
    std::unique_ptr<amrex::MultiFab>      dust_pm10_24h;   ///< 24-hour running average PM10.
    std::unique_ptr<amrex::MultiFab>      dust_pm25_exceed; ///< Exceedance flag PM2.5.
    std::unique_ptr<amrex::MultiFab>      dust_pm10_exceed; ///< Exceedance flag PM10.

    // Phase 18: MSHA worker exposure MultiFabs on dust grid.
    std::unique_ptr<amrex::MultiFab>      dust_msha_dose;       ///< Cumulative dose [mg/m³/h].
    std::unique_ptr<amrex::MultiFab>      dust_msha_twa;        ///< 8h TWA [mg/m³].
    std::unique_ptr<amrex::MultiFab>      dust_msha_exceed;     ///< PEL exceedance [0/1].
    std::unique_ptr<amrex::MultiFab>      dust_msha_shift_twa;  ///< TWA at last shift end [mg/m³].

#if defined(ERF_USE_PARTICLES)
     // Phase 19: Lagrangian super-particles for source-receptor attribution
     std::unique_ptr<class ERFDustPC>      m_dust_pc;         ///< Particle container
     std::unique_ptr<amrex::MultiFab>     dust_source_map;   ///< Cumulative source deposition [kg/m^2]
     amrex::Geometry                       m_geom_atm;        ///< Atmospheric geometry (stored for particle use)
#endif

    // Phase 20: multi-mine site assignment
    std::unique_ptr<amrex::MultiFab>      dust_site_id;      ///< Per-cell site index [0=unassigned, 1..N].
    amrex::Vector<amrex::Real>            m_site_ustar_t_factors; ///< Per-site u*_t reduction factors

    // Phase 23: DOE Critical Materials Assessment
    std::unique_ptr<amrex::MultiFab>      dust_cm_flux;      ///< Critical material flux [kg_CM/m^2/s].

    int                                   m_dust_scalar_comp = -1; ///< cc_source component index for dust mass.
                                                                     ///< Set during initialize from ERF_IndexDefines.H.

    DustParams                            m_params;
    int                                   m_step = 0;
    amrex::Real                           m_time = 0.0;
    int                                   m_last_dust_plot_step = -1;
    amrex::Real                           m_last_phreeqc_write_time = 0.0; ///< Time of last PHREEQC feedback write.
    int                                   m_last_phreeqc_write_step = -1;      ///< Step of last PHREEQC feedback write.
    amrex::Real                           m_last_phreeqc_update = -1.0e30;
    int                                   m_msha_shift_count = 0; ///< Phase 18: shift counter for MSHA.
#ifdef ERF_USE_DUST
    BlastSchedule                         m_blast_schedule;
    bool                                  m_has_blast_schedule = false;
    RoadSchedule                          m_road_schedule;
#endif
};

#endif // ERF_DUST_LAYER_H
