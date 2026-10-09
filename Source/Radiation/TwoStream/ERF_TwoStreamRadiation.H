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
 *  - define_level(lev, radChoice, rdOcp, ba2d, dm, ba, domain, supply_land_forcing)
 *    whenever a level is (re)built;
 *  - advance(...) twice per step: "pre_dycore" runs the column sweep on the
 *    old state and writes the heating rates, "post_dycore" reuses that
 *    sweep's flux diagnostics and advances the SEB state by dt_step;
 *  - write_land_forcing(...) after a sweep, when a land-surface model integrates
 *    on the radiation this model computes (see supply_land_forcing);
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

    /** @brief Read-only view of the per-level prognostic surface-temperature state. */
    [[nodiscard]] const amrex::MultiFab* prognostic_surface_temperature_state (int lev) const noexcept {
        return (lev >= 0 && lev < static_cast<int>(m_t_sfc.size()))
            ? m_t_sfc[lev].get() : nullptr;
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
     * @param supply_land_forcing allocate the downwelling SW, downwelling LW and
     *                   cos(zenith) fields the sweep fills for a land-surface model
     *                   (SolverChoice::rad_feeds_lsm()); off, the sweep does no
     *                   extra work and the fields do not exist
     */
    void define_level (int lev,
                       const RadChoice& rad_choice,
                       amrex::Real rdOcp,
                       const amrex::BoxArray& ba2d,
                       const amrex::DistributionMapping& dm,
                       const amrex::BoxArray& ba,
                       const amrex::Box& domain,
                       bool supply_land_forcing = false);

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
     * @param radiation_inputs canonical surface-property inputs, in fixed order
     * @param noahmp_active whether Noah-MP diagnostics may be queried
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
     * @param sfc_sens_flux the conservative heat flux (rho w'theta', z faces)
     *                   the zlo surface layer applies to this level, or
     *                   nullptr when it applies none; the surface energy
     *                   balance's H where no land-surface model supplies one
     *                   (see ERF_SEBTurbulentFlux.H)
     * @param sfc_laten_flux the moisture flux (rho w'qv', z faces) likewise,
     *                   for LE; nullptr without a moisture model
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
                  const amrex::Vector<const amrex::MultiFab*>& radiation_inputs,
                  bool noahmp_active,
                  amrex::MultiFab* qheating,
                  amrex::MultiFab* rad_fluxes,
                  const amrex::MultiFab* t_surf,
                  const amrex::MultiFab* sfc_sens_flux,
                  const amrex::MultiFab* sfc_laten_flux,
                  const amrex::MultiFab* lat_m,
                  const amrex::MultiFab* lon_m,
                  double epoch_time,
                  bool have_datetime);

    /**
     * @brief The prognostic surface state of one level, or nullptr when the SEB is
     * off or the level has not been built.
     *
     * Exposed so ERF can average t_sfc and q_sfc down after the finer levels have
     * advanced. Every level evolves its own copy, and without that transfer the two
     * estimates of one surface drift apart -- which would give the longwave boundary
     * condition two different values for the same ground. ERF owns refRatio and the
     * level geometries; this class does not, and should not grow them for one call.
     */
    [[nodiscard]] amrex::MultiFab* seb_t_sfc (int lev) {
        return (active() && m_rad->seb_enable && lev < m_t_sfc.size()) ? m_t_sfc[lev].get() : nullptr;
    }
    [[nodiscard]] amrex::MultiFab* seb_q_sfc (int lev) {
        return (active() && m_rad->seb_enable && lev < m_q_sfc.size()) ? m_q_sfc[lev].get() : nullptr;
    }

    /**
     * @brief The sensible and latent heat fluxes [W/m^2] the surface energy balance
     * of one level last used, or nullptr when the SEB is off or the level has not
     * been built. Read-only: for output.
     */
    [[nodiscard]] const amrex::MultiFab* seb_hfx (int lev) const {
        return (active() && m_rad->seb_enable && lev < m_hfx_sfc.size()) ? m_hfx_sfc[lev].get() : nullptr;
    }
    [[nodiscard]] const amrex::MultiFab* seb_lh (int lev) const {
        return (active() && m_rad->seb_enable && lev < m_lh_sfc.size()) ? m_lh_sfc[lev].get() : nullptr;
    }

    /**
     * @brief What level lev's surface hands the surface layer for its surface mixing ratio
     * (erf.radiation.seb_surface_layer_uses_moisture), per column, in the components of
     * erf_surface_moisture::Component: the soil-water factor (seb_moisture_availability
     * of the current q_s); with erf.radiation.seb_soil_type, the bare-soil resistance and
     * the relative humidity of its pore air at the skin temperature; and
     * with erf.radiation.seb_vegetation_type as well, the vegetated fraction and the canopy
     * resistance without its vapour-deficit factor (from the net shortwave the balance
     * holds divided by (1 - surface albedo) -- the sweep's with seb_use_radiation_fluxes,
     * seb_sw_flux_default otherwise -- and the skin temperature). Absent pieces are a
     * fraction of 0 and resistances of 1e6 s/m. nullptr
     * when the option is off or the level has not been built. Recomputed on every call,
     * on the valid cells; the halo is not filled.
     */
    [[nodiscard]] const amrex::MultiFab* seb_surface_moisture (int lev);

    /** @brief True when the force-restore state is being evolved and so needs transferring. */
    [[nodiscard]] bool seb_prognostic_active () const {
        return active() && m_rad->seb_enable && m_rad->seb_prognostic_enable;
    }

    /** @brief True when define_level allocated the land-model forcing fields on @p lev. */
    [[nodiscard]] bool supplies_land_forcing (int lev) const noexcept {
        return lev >= 0 && lev < static_cast<int>(m_sw_dn_sfc.size()) && m_sw_dn_sfc[lev] != nullptr;
    }

    /**
     * @brief The surface forcing of the last sweep on @p lev (2D: i,j), or nullptr when
     * the level does not supply it: total downwelling shortwave [W/m^2], downwelling
     * longwave [W/m^2] and the cosine of the solar zenith angle floored at zero [1].
     * lsm_undefined until a sweep has run on the level, so a copy made before one is
     * caught by the land model's missing-input check instead of reading as zero.
     */
    [[nodiscard]] const amrex::MultiFab* land_forcing_sw_dn (int lev) const noexcept {
        return supplies_land_forcing(lev) ? m_sw_dn_sfc[lev].get() : nullptr;
    }
    [[nodiscard]] const amrex::MultiFab* land_forcing_lw_dn (int lev) const noexcept {
        return supplies_land_forcing(lev) ? m_lw_dn_sfc[lev].get() : nullptr;
    }
    [[nodiscard]] const amrex::MultiFab* land_forcing_cos_zenith (int lev) const noexcept {
        return supplies_land_forcing(lev) ? m_cos_zenith[lev].get() : nullptr;
    }

    /**
     * @brief Copy the surface forcing of the last sweep on @p lev into a land-surface
     * model's radiation inputs (SurfaceModel's "sw_flux_dn", "lw_flux_dn" and
     * "cos_zenith_angle" outputs). A null destination is skipped. Each destination
     * must lie on this level's horizontal grids; see copy_surface_plane. Before the
     * level's first sweep this copies the lsm_undefined sentinel, not zero.
     */
    void write_land_forcing (int lev,
                             amrex::MultiFab* sw_flux_dn,
                             amrex::MultiFab* lw_flux_dn,
                             amrex::MultiFab* cos_zenith_angle) const;

    /** @brief Write the prognostic surface state of one level (SEB on). */
    void write_checkpoint (int lev, const std::string& checkpointname) const;

    /** @brief Read it back when the checkpoint carries it; else keep defaults. */
    void read_checkpoint (int lev, const std::string& restart_chkfile);

    // Orbital elements of the current year (Berger 1978, or the erf.rad_orbital_*
    // overrides), computed once per year rather than every step, and the
    // one-time note about an unscaled irradiance. Public so the file-local
    // solar-state helper and two_stream_sun_date() can fill it.
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
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_sfc;       // surface moisture, prognostic (see RadChoice::seb_q_sfc_default)
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_t_deep;      // deep-soil temperature [K]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_q_deep;      // deep-soil moisture q_s restores to
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_moisture;    // seb_surface_moisture (seb_surface_layer_uses_moisture)

    // What the sweep hands to a land-surface model (2D: i,j), allocated when
    // define_level is asked to supply it. These are the DOWNWELLING fluxes, as opposed to
    // the net fluxes above that the SEB consumes: a land model applies its own albedo
    // and emissivity, so they cannot be recovered from the net fields.
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_sw_dn_sfc;   // downwelling shortwave at the surface [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_lw_dn_sfc;   // downwelling longwave at the surface [W/m^2]
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> m_cos_zenith;  // cosine of the solar zenith angle, >= 0 [1]

    amrex::Vector<FluxDiag> m_flux_diag;
};

/**
 * @brief Copy a 2D field (k = 0 only) into the k = 0 plane of a surface-model field.
 *
 * Surface-model fields keep their surface values in the k = 0 plane, which may lie in
 * the field's z ghost region (SLM stores soil layers at k < 0), and carry x/y ghost cells
 * the source need not have. Only the i,j extent of each destination VALID box is
 * written, so ghost columns and every other k are left alone. The two fields must share
 * a horizontal layout: the destination's boxes flattened to k = 0 equal the source's,
 * on the same DistributionMapping -- which holds for ERF's ba2d and a land model's
 * data built from the same level grids, and is asserted rather than assumed.
 */
void copy_surface_plane (const amrex::MultiFab& src2d, amrex::MultiFab& dst);

/** @brief RRTMGP's reference total solar irradiance [W/m^2], scaled by the date's Earth-Sun distance factor. */
inline constexpr amrex::Real two_stream_tsi_reference = 1360.9;

/** @brief The calendar sun of the two-stream model at one time (see two_stream_sun_date()). */
struct TwoStreamSunDate {
    double calday = 1.0;   ///< Day of the year (1 = 1 January), UTC time of day as the fraction
    double declin = 0.0;   ///< Solar declination [rad]
    double eccf   = 1.0;   ///< Earth-Sun distance factor (d0/d)^2
};

/**
 * @brief The calendar sun of the two-stream model at a UTC epoch time (start_time + t).
 *
 * The date as the RRTMGP interface forms it, the orbit of the year (or the
 * erf.rad_orbital_* overrides) cached in @p orbit, and the declination and distance
 * factor of ERF_OrbCosZenith.H. The two-stream sweep places its sun from this, with the
 * column's position, and so does the immersed-boundary balance with
 * erf.ibseb.sun_mode = two_stream.
 */
TwoStreamSunDate two_stream_sun_date (const RadChoice& rc, TwoStreamRadiation::OrbitalCache& orbit,
                                      double epoch_time);

#endif // ERF_TWO_STREAM_RADIATION_H_
