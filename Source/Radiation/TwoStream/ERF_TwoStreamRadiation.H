#ifndef ERF_TWO_STREAM_RADIATION_H_
#define ERF_TWO_STREAM_RADIATION_H_

#include <memory>
#include <string>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

#include <ERF_RadStruct.H>
#include <ERF_LandSurface.H>
#include <ERF_RadiationDiagnostics.H>

/**
 * @brief The two-stream radiation model (erf.radiation_model = TwoStream).
 *
 * Owns everything the model keeps between steps: the per-level surface
 * property fields it falls back on when no land-surface model supplies
 * them, the surface-energy-balance (SEB) fields and force-restore state,
 * and the domain-mean flux diagnostics of the last column sweep. ERF holds
 * one instance and calls it only after testing
 * `solverChoice.rad_type == RadiationType::TwoStream`; every method
 * here checks the same predicate and returns without acting otherwise.
 *
 * Lifecycle:
 *  - resize(nlevs) once at construction (before the inputs are read);
 *  - define_level(lev, radChoice, rdOcp, ba2d, dm, ba, domain) whenever a level is
 *    (re)built;
 *  - advance(...) twice per step: "pre_dycore" runs the column sweep on the
 *    old state and writes the heating rates, "post_dycore" reuses that
 *    sweep's flux diagnostics and advances the SEB state by dt_step;
 *  - write_checkpoint / read_checkpoint for the prognostic surface state.
 *
 * The heating-rate MultiFab itself stays in ERF, since RRTMGP writes the
 * same 2-component (SW, LW) array and the RhoTheta source term and the
 * plotfile read it without caring which model filled it.
 */
class TwoStreamRadiation
{
public:
    TwoStreamRadiation () = default;

    /** @brief Size the per-level containers; no allocation happens here. */
    void resize (int nlevs_max);

    /** @brief True once a level has been defined with the model selected. */
    [[nodiscard]] bool active () const {
        return (m_rad != nullptr) && m_rad->enabled;
    }

    /**
     * @brief Allocate the 2D surface fields of one level on the given
     * horizontal BoxArray and fill them with the scalar defaults. A no-op
     * unless rad_choice.rad_type == TwoStream.
     *
     * Also enforces the model's one grid requirement: every box of @p ba must
     * span @p domain in z, because the column sweep applies the
     * top-of-atmosphere and surface boundary conditions at the ends of a box.
     * Aborting here means a bad refinement patch is rejected when the level is
     * built (including on regrid), not at the first radiation call.
     *
     * @param lev        level
     * @param rad_choice the model's inputs
     * @param ba2d       horizontal BoxArray the 2D surface fields live on
     * @param dm         distribution mapping shared with the level's state
     * @param ba         the level's 3D BoxArray, checked for the z span
     * @param domain     the level's domain Box
     */
    void define_level (int lev,
                       const RadChoice& rad_choice,
                       amrex::Real rdOcp,
                       const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm,
                       const amrex::BoxArray& ba,
                       const amrex::Box& domain);

    /**
     * @brief One radiation call at one level.
     *
     * @param lev        level
     * @param nstep      step number, used as the diagnostics row index
     * @param time       time logged with the diagnostics row [s]
     * @param dt_step    step size, used by the SEB Euler update [s]
     * @param call_site  "pre_dycore" (column sweep + diagnostics) or
     *                   "post_dycore" (cached diagnostics + SEB update)
     * @param cons_old   old-state conserved variables of this level
     * @param z_phys_nd  nodal heights (layer interfaces), or nullptr on a uniform grid
     * @param geom       geometry of this level
     * @param lsm        land-surface model (queried for per-column fields)
     * @param qheating   2-component (SW, LW) heating-rate MultiFab of this
     *                   level, written at the pre-dycore call; may be nullptr
     * @param rad_fluxes 4-component (SW up, SW down, LW up, LW down) interface
     *                   fluxes of this level in RRTMGP's level layout (index k
     *                   the lower interface of layer k, the top of the
     *                   atmosphere in the z-ghost cell above the top layer),
     *                   written at the pre-dycore call; may be nullptr
     * @param t_surf     the surface layer's surface temperature, or nullptr;
     *                   used as the longwave boundary condition where no
     *                   land-surface field exists, as RRTMGP does
     * @param lat_m      per-column latitude [degrees], or nullptr (then
     *                   erf.rad_cons_lat); only the WRF/metgrid inits fill it
     * @param lon_m      per-column longitude [degrees], likewise
     * @param epoch_time seconds since 1970-01-01 00:00 UTC of this call
     *                   (start_time + time), the RRTMGP time base; double, as
     *                   RRTMGP keeps it, since a float cannot resolve a time
     *                   of day at ~1.7e9 s
     * @param have_datetime true when start_datetime (or a WRF/metgrid start
     *                   date) fixed the calendar; without it a dynamic sun
     *                   cannot be placed and the call aborts
     */
    void advance (int lev,
                  int nstep,
                  amrex::Real time,
                  amrex::Real dt_step,
                  const std::string& call_site,
                  const amrex::MultiFab& cons_old,
                  const amrex::MultiFab* z_phys_nd,
                  const amrex::Geometry& geom,
                  LandSurface& lsm,
                  amrex::MultiFab* qheating,
                  amrex::MultiFab* rad_fluxes,
                  const amrex::MultiFab* t_surf,
                  const amrex::MultiFab* lat_m,
                  const amrex::MultiFab* lon_m,
                  double epoch_time,
                  bool have_datetime);

    /** @brief Write the prognostic surface state of one level (SEB on). */
    void write_checkpoint (int lev, const std::string& checkpointname) const;

    /** @brief Read it back when the checkpoint carries it; else keep defaults. */
    void read_checkpoint (int lev, const std::string& restart_chkfile);

    // Orbital elements of the current year (Berger 1978, or the erf.rad_orbital_*
    // overrides), computed once per year rather than every step, and the
    // one-time note about an unscaled irradiance. Public so the file-local
    // solar-state helper can fill it.
    struct OrbitalCache {
        int year = -9999;
        double eccen = 0.0;
        double obliqr = 0.0;
        double lambm0 = 0.0;
        double mvelpp = 0.0;
        bool noted_unscaled = false;
    };

private:
    // Domain-mean flux diagnostics of the last column sweep, per level, so the
    // post-dycore call can report them without sweeping again.
    struct FluxDiag {
        // False until a column sweep has filled this level's entry. The post-dycore call
        // reports the cached values rather than sweeping again, and a level whose
        // pre-dycore sweep was skipped -- the interp-from-coarse-pending path, which
        // returns before the sweep -- would otherwise report the zero-initialised struct
        // as though the fluxes really were zero.
        bool valid = false;
        amrex::Real SW_surface = 0.0;
        amrex::Real SW_TOA = 0.0;      // domain-mean incident SW at the top
        amrex::Real SW_up_TOA = 0.0;
        amrex::Real LW_net_surface = 0.0;
        amrex::Real LW_up_TOA = 0.0;
        amrex::Real heating_rate_max = 0.0;
    };

    const RadChoice* m_rad = nullptr;
    amrex::Real m_rdOcp = RdoCp;
    OrbitalCache m_orbit;

    // CSV / console diagnostics writer, one per level, created on that level's
    // first advance() so its header and duplicate-row state persist across
    // calls. Per level because every level of a hierarchy reports the same
    // step, call_site and time, and a shared duplicate guard would suppress
    // all but the first to report.
    amrex::Vector<std::unique_ptr<RadiationDiagnostics>> m_diag;

    // Fallback surface properties when no LSM supplies them (2D: i,j).
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_alb_sw;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_emiss_lw;
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_t_sfc;      // prognostic when the SEB evolves it

    // Surface-energy-balance fields (2D: i,j), allocated when seb_enable.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_sw_flux_sfc; // net absorbed shortwave [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_lw_flux_sfc; // net absorbed longwave [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_hfx_sfc;     // sensible heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_lh_sfc;      // latent heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_grdflx_sfc;  // ground heat flux [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_sfc;       // surface moisture [kg/kg], prognostic
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_t_deep;      // deep-soil temperature [K]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_deep;      // deep-soil moisture [kg/kg]

    amrex::Vector<FluxDiag> m_flux_diag;
};

#endif // ERF_TWO_STREAM_RADIATION_H_
