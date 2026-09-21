#include <ERF_TwoStreamRadiation.H>
#include <ERF_RadStruct.H>
#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_Math.H>
#include <ERF_RadiationDiagnostics.H>
#include <ERF_TwoStreamColumn.H>
#include <ERF_PrognosticCloudFraction.H>
#include <ERF_AerosolOpticalDepth.H>
#include <ERF_SimplifiedSEB.H>
#include <ERF_OrbCosZenith.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Gpu.H>
#include <ERF_IndexDefines.H>
#include <ERF_EOS.H>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>

using namespace amrex;

namespace {
// RRTMGP's reference total solar irradiance [W/m^2]; the date's Earth-Sun
// distance factor scales it unless erf.fixed_total_solar_irradiance is set.
constexpr amrex::Real tsi_reference = 1360.9;

/**
 * Set the sun of this call from the inputs shared with RRTMGP. A fixed
 * cosine (erf.fixed_solar_zenith_angle > 0) with a fixed irradiance needs no
 * calendar; a fixed cosine alone takes the date-scaled irradiance when a
 * start date is known and the unscaled reference otherwise. A calendar sun
 * needs the date: it comes from start_datetime (epoch_time = start_time + t,
 * as RRTMGP forms it, in double so a time of day at ~1.7e9 s is resolved) and
 * goes through the same orbital code (orbital_params once per year,
 * orbital_decl per call, both of ERF_OrbCosZenith.H) for the declination and
 * the Earth-Sun distance factor; the column sweep then evaluates the position
 * over each column. Without a start date there is no sun to place, so the run
 * stops and says what to set.
 */
void set_solar_state (TwoStreamParams& p, const RadChoice& rc,
                      TwoStreamRadiation::OrbitalCache& orbit,
                      double epoch_time, bool have_datetime, int lev, int nstep)
{
    const bool fixed_sun = (rc.fixed_solar_zenith_angle > 0.0);
    const bool fixed_tsi = (rc.fixed_total_solar_irradiance >= 0.0);
    p.solar_dynamic    = !fixed_sun;
    p.cos_zenith_fixed = rc.fixed_solar_zenith_angle;
    p.lat_cons_rad     = rc.rad_cons_lat * PI / 180.0;
    p.lon_cons_rad     = rc.rad_cons_lon * PI / 180.0;
    p.S0               = fixed_tsi ? rc.fixed_total_solar_irradiance : tsi_reference;
    if (!rc.sw_enabled || (fixed_sun && fixed_tsi)) { return; }

    if (fixed_sun && !have_datetime) {
        // A fixed sun with the default irradiance and no calendar to scale it
        // by: the unscaled reference. Say so once.
        if (!orbit.noted_unscaled) {
            if (ParallelDescriptor::IOProcessor()) {
                Print() << "NOTE: erf.fixed_solar_zenith_angle is set, erf.fixed_total_solar_irradiance "
                           "is not, and no start_datetime is known: the two-stream model uses the "
                           "unscaled reference irradiance of " << tsi_reference << " W/m^2.\n";
            }
            orbit.noted_unscaled = true;
        }
        return;
    }
    if (!have_datetime) {
        amrex::Abort("TwoStreamRadiation (level " + std::to_string(lev) + ", step " +
                     std::to_string(nstep) + "): the sun follows the calendar because "
                     "erf.fixed_solar_zenith_angle is not set, and no start date is known. "
                     "Set start_datetime = \"YYYY-MM-DD HH:MM:SS\" (UTC), as for RRTMGP, or fix "
                     "the sun with erf.fixed_solar_zenith_angle (the cosine of the angle).");
    }

    // Calendar date of this call (UTC), as the RRTMGP interface forms it.
    time_t timestamp = time_t(epoch_time);
    struct tm timeinfo{};
#if defined(_WIN32)
    gmtime_s(&timeinfo, &timestamp);
#else
    gmtime_r(&timestamp, &timeinfo);
#endif
    const int year = (rc.rad_orbital_year >= 0) ? rc.rad_orbital_year : timeinfo.tm_year + 1900;
    const int mon = timeinfo.tm_mon + 1;
    const int day = timeinfo.tm_mday;
    const int sec = timeinfo.tm_hour*3600 + timeinfo.tm_min*60 + timeinfo.tm_sec;

    // Orbital parameters of the year (Berger 1978) unless overridden: a few
    // hundred series terms and six vectors, so once per year, not per step.
    if (orbit.year != year) {
        int  iyear = year;
        double eccen = rc.rad_orbital_eccentricity;
        double obliq = rc.rad_orbital_obliquity;
        double mvelp = rc.rad_orbital_mvelp;
        double obliqr = 0.0, lambm0 = 0.0, mvelpp = 0.0;
        orbital_params(iyear, eccen, obliq, mvelp, obliqr, lambm0, mvelpp);
        orbit.year = year;
        orbit.eccen = eccen; orbit.obliqr = obliqr; orbit.lambm0 = lambm0; orbit.mvelpp = mvelpp;
    }

    // Declination and Earth-Sun distance factor of the day.
    double calday = orbital_calday(year, mon, day, sec);
    double eccen = orbit.eccen, mvelpp = orbit.mvelpp, lambm0 = orbit.lambm0, obliqr = orbit.obliqr;
    double delta = 0.0, eccf = 1.0;
    orbital_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, eccf);

    p.calday = static_cast<amrex::Real>(calday);
    p.declin = static_cast<amrex::Real>(delta);
    if (!fixed_tsi) { p.S0 = tsi_reference * static_cast<amrex::Real>(eccf); }
}
} // namespace


namespace {
// Fill a 2D surface-energy-balance field from the LSM field of the given
// name, scaled by `scale` (Noah-MP's fira is positive upward, the SEB wants
// absorbed fluxes positive), plus an optional second field added on top
// (Noah-MP splits absorbed shortwave into sav and sag). Falls back to the
// scalar default when the LSM does not expose the field.
bool lsm_has_field(LandSurface& lsm, int lev, const char* field_name)
{
    std::string varname(field_name);
    const int lsm_idx = lsm.Get_DataIdx(lev, varname);
    return (lsm_idx >= 0) && (lsm.Get_Data_Ptr(lev, lsm_idx) != nullptr);
}

void fill_or_copy_seb_field(
    MultiFab* seb_mf,
    LandSurface& lsm,
    int lev,
    const char* field_name,
    amrex::Real fallback_value,
    amrex::Real scale = 1.0,
    const char* add_field_name = nullptr)
{
    if (seb_mf == nullptr) return;

    std::string varname(field_name);
    int lsm_idx = lsm.Get_DataIdx(lev, varname);
    if (lsm_idx >= 0) {
        if (MultiFab* lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx)) {
            MultiFab::Copy(*seb_mf, *lsm_ptr, 0, 0, 1, 0);
            if (scale != 1.0) seb_mf->mult(scale, 0, 1, 0);
            if (add_field_name != nullptr) {
                std::string addname(add_field_name);
                int add_idx = lsm.Get_DataIdx(lev, addname);
                if (add_idx >= 0) {
                    if (MultiFab* add_ptr = lsm.Get_Data_Ptr(lev, add_idx)) {
                        MultiFab::Add(*seb_mf, *add_ptr, 0, 0, 1, 0);
                    }
                }
            }
            return;
        }
    }
    seb_mf->setVal(fallback_value);
}
}

/**
 * @file ERF_TwoStreamRadiation.cpp
 * @brief TwoStreamRadiation: the two-stream radiation model and its surface state.
 *
 * Computes SW/LW fluxes and per-level heating rates using real per-column
 * vertical sweeps over the atmospheric grid. Reads temperature and density
 * from the state and properly accumulates optical depth through all
 * vertical levels.
 *
 * Capabilities include:
 * - Height-varying optical depth: optional "cloud_layer" tau_profile_type
 *   adds cloud_tau_per_layer on top of the clear-sky background within
 *   [cloud_base_height_m, cloud_top_height_m].
 * - Cloud fraction masking: blends clear-sky and cloudy-column fluxes via
 *   F = (1 - cloud_fraction) * F_clear + cloud_fraction) * F_cloudy.
 * - Diffuse (scattered) SW flux: two-stream reflectance/transmittance per
 *   layer combined with the surface albedo by the adding method, giving
 *   upward and downward diffuse streams (see ERF_TwoStreamSW.H).
 * - Per-level heating rate output: writes the SW/LW radiative tendencies of
 *   potential temperature, (dT/dt)/pi, to a 2-component MultiFab
 *   (component 0 = SW, component 1 = LW), mirroring the RRTMGP convention
 *   expected by the RhoTheta source term.
 *
 * Includes support for:
 * - Height-varying surface properties (albedo, emissivity, temperature)
 * - Dynamic cloud fraction from relative humidity and cloud water
 * - Prescribed bulk aerosol optical depth
 * - Dynamic solar geometry (time-varying solar position)
 * - Surface energy balance diagnostics and prognostic updates
 *
 * Vertical orientation follows ERF: k = kmin is the surface layer and
 * k = kmax the top layer. SW sweeps downward from kmax to kmin; LW sweeps
 * downward (TOA -> surface) and then upward (surface -> TOA). Layer
 * temperature is obtained from rho*theta through the Exner function.
 *
 * CSV diagnostics (domain means unless noted): SW_surface is the SW absorbed
 * by the surface, SW_TOA the incident SW at the top of the atmosphere,
 * SW_up_TOA the reflected SW leaving the top, LW_net_surface the net
 * (up - down) LW at the surface, LW_up_TOA the outgoing LW at the top, and
 * heating_rate_max the max(|Q_sw|+|Q_lw|) over the column evaluations.
 */

void
TwoStreamRadiation::resize (int nlevs_max)
{
    m_alb_sw.resize(nlevs_max);
    m_emiss_lw.resize(nlevs_max);
    m_t_sfc.resize(nlevs_max);
    m_sw_flux_sfc.resize(nlevs_max);
    m_lw_flux_sfc.resize(nlevs_max);
    m_hfx_sfc.resize(nlevs_max);
    m_lh_sfc.resize(nlevs_max);
    m_grdflx_sfc.resize(nlevs_max);
    m_q_sfc.resize(nlevs_max);
    m_t_deep.resize(nlevs_max);
    m_q_deep.resize(nlevs_max);
    m_flux_diag.resize(nlevs_max);
}

void
TwoStreamRadiation::define_level (int lev,
                                  const RadChoice& rad_choice,
                                  const amrex::Real rdOcp,
                                  const BoxArray& ba2d,
                                  const DistributionMapping& dm)
{
    if (!rad_choice.enabled) { return; }
    m_rad = &rad_choice;
    m_rdOcp = rdOcp;

    // 2D surface fields on the horizontal BoxArray, one ghost cell in x and y
    const IntVect ng_sfc{1,1,0};
    m_alb_sw[lev]   = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
    m_emiss_lw[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
    m_t_sfc[lev]    = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
    m_alb_sw[lev]->setVal(rad_choice.surface_albedo_sw);
    m_emiss_lw[lev]->setVal(rad_choice.surface_emissivity_lw);
    m_t_sfc[lev]->setVal(rad_choice.rad_t_sfc);

    if (rad_choice.seb_enable) {
        m_sw_flux_sfc[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_lw_flux_sfc[lev] = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_hfx_sfc[lev]     = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_lh_sfc[lev]      = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_grdflx_sfc[lev]  = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_q_sfc[lev]       = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_t_deep[lev]      = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_q_deep[lev]      = std::make_unique<MultiFab>(ba2d, dm, 1, ng_sfc);
        m_sw_flux_sfc[lev]->setVal(rad_choice.seb_sw_flux_default);
        m_lw_flux_sfc[lev]->setVal(rad_choice.seb_lw_flux_default);
        m_hfx_sfc[lev]->setVal(rad_choice.seb_hfx_default);
        m_lh_sfc[lev]->setVal(rad_choice.seb_lh_default);
        m_grdflx_sfc[lev]->setVal(rad_choice.seb_grdflx_default);
        m_q_sfc[lev]->setVal(rad_choice.seb_q_sfc_default);
        m_t_deep[lev]->setVal(rad_choice.seb_t_deep_default);
        m_q_deep[lev]->setVal(rad_choice.seb_q_deep_default);
    }
    m_flux_diag[lev] = FluxDiag{};
}

void
TwoStreamRadiation::write_checkpoint (int lev, const std::string& checkpointname) const
{
    // The force-restore state is the only part of this model a restart must
    // carry; without it T_s and q_s restart from the scalar defaults.
    if (!active() || !m_rad->seb_enable) { return; }
    if (m_t_sfc[lev]) {
        VisMF::Write(*m_t_sfc[lev],
                     MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TwoStream_TSfc"));
    }
    if (m_q_sfc[lev]) {
        VisMF::Write(*m_q_sfc[lev],
                     MultiFabFileFullPrefix(lev, checkpointname, "Level_", "TwoStream_QSfc"));
    }
}

void
TwoStreamRadiation::read_checkpoint (int lev, const std::string& restart_chkfile)
{
    // Older checkpoints do not carry the state; then the defaults set by
    // define_level stand.
    if (!active() || !m_rad->seb_enable) { return; }
    const std::string tsfc_name =
        MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "TwoStream_TSfc");
    if (m_t_sfc[lev] && amrex::FileExists(tsfc_name + "_H")) {
        VisMF::Read(*m_t_sfc[lev], tsfc_name);
    }
    const std::string qsfc_name =
        MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "TwoStream_QSfc");
    if (m_q_sfc[lev] && amrex::FileExists(qsfc_name + "_H")) {
        VisMF::Read(*m_q_sfc[lev], qsfc_name);
    }
}

void
TwoStreamRadiation::advance (int lev,
                             int nstep,
                             amrex::Real time,
                             amrex::Real dt_step,
                             const std::string& call_site,
                             const MultiFab& cons_old,
                             const MultiFab* z_phys_nd,
                             const Geometry& geom,
                             LandSurface& lsm,
                             MultiFab* qheating,
                            MultiFab* rad_fluxes,
                            const MultiFab* t_surf,
                            const MultiFab* lat_m,
                            const MultiFab* lon_m,
                            double epoch_time,
                            bool have_datetime)
{
    BL_PROFILE("TwoStreamRadiation::advance()");

    // Only proceed if TwoStream radiation is enabled
    if (!active()) { return; }
    const RadChoice& rad_choice = *m_rad;

    // ---- Contract checks. Each of these would otherwise surface as a wrong
    // heating rate or an out-of-bounds read several routines away.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(call_site == "pre_dycore" || call_site == "post_dycore",
        "TwoStreamRadiation::advance: call_site must be pre_dycore or post_dycore");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev < m_alb_sw.size() && m_alb_sw[lev] != nullptr,
        "TwoStreamRadiation::advance called on a level define_level has not built");
    // The sweep indexes the 2D surface fields with the MFIter of the 3D state,
    // so both must have the same number of boxes on the same ranks.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        m_alb_sw[lev]->boxArray().size() == cons_old.boxArray().size() &&
        m_alb_sw[lev]->DistributionMap() == cons_old.DistributionMap(),
        "TwoStreamRadiation: the 2D surface fields and the 3D state are laid out differently");
    if (qheating != nullptr) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(qheating->nComp() == 2 &&
                                         qheating->boxArray() == cons_old.boxArray() &&
                                         qheating->DistributionMap() == cons_old.DistributionMap(),
            "TwoStreamRadiation: qheating must be the 2-component (SW, LW) field on the state's grids");
    }
    if (rad_choice.seb_enable) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_sw_flux_sfc[lev] && m_lw_flux_sfc[lev] && m_hfx_sfc[lev] &&
                                         m_lh_sfc[lev] && m_grdflx_sfc[lev] && m_q_sfc[lev] &&
                                         m_t_deep[lev] && m_q_deep[lev],
            "TwoStreamRadiation: seb_enable is set but the SEB fields were not allocated");
    }
    if (call_site == "post_dycore" && rad_choice.seb_prognostic_enable) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dt_step > 0.0 && std::isfinite(dt_step),
            "TwoStreamRadiation: the force-restore update needs a positive, finite dt_step");
    }
    // The column kernel would substitute placeholders (rho = 1, rho*theta of
    // 288 K) for a non-finite or non-positive density or rho*theta and carry
    // on, hiding a corrupt state behind plausible heating rates. Refuse such
    // a state here instead.
    //
    // This runs every step, so it is one pass over the state and one
    // collective, not the three that contains_nan() plus two MultiFab::min()
    // calls would cost. A non-finite value is mapped to -infinity so that the
    // same minimum answers both questions: -infinity means non-finite, and any
    // other value <= 0 means non-positive. The reduction is collective, so
    // every rank takes the same branch.
    if (call_site == "pre_dycore") {
        ReduceOps<ReduceOpMin, ReduceOpMin> state_ops;
        ReduceData<amrex::Real, amrex::Real> state_data(state_ops);
        using StateTuple = typename decltype(state_data)::Type;
        for (MFIter mfi(cons_old, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto& arr = cons_old.const_array(mfi);
            state_ops.eval(mfi.tilebox(), state_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> StateTuple {
                    constexpr amrex::Real neg_inf = -std::numeric_limits<amrex::Real>::infinity();
                    const amrex::Real rho = arr(i,j,k,Rho_comp);
                    const amrex::Real rth = arr(i,j,k,RhoTheta_comp);
                    return {amrex::Math::isfinite(rho) ? rho : neg_inf,
                            amrex::Math::isfinite(rth) ? rth : neg_inf};
                });
        }
        auto state_tuple = state_data.value(state_ops);
        amrex::Real mins[2] = {amrex::get<0>(state_tuple), amrex::get<1>(state_tuple)};
        ParallelDescriptor::ReduceRealMin(mins, 2);
        const amrex::Real rho_min = mins[0];
        const amrex::Real rth_min = mins[1];
        if (!std::isfinite(rho_min) || !std::isfinite(rth_min)) {
            amrex::Abort("TwoStreamRadiation: the state handed to the column sweep at level " +
                         std::to_string(lev) + ", step " + std::to_string(nstep) +
                         " has a non-finite density or rho*theta");
        }
        if (!(rho_min > 0.0) || !(rth_min > 0.0)) {
            amrex::Abort("TwoStreamRadiation: the state handed to the column sweep at level " +
                         std::to_string(lev) + ", step " + std::to_string(nstep) +
                         " has a non-positive density (min " + std::to_string(rho_min) +
                         ") or rho*theta (min " + std::to_string(rth_min) + ")");
        }
    }

    // The column sweep runs once per step, at the pre-dycore call. The
    // post-dycore call would sweep the same old state (vars_old is not
    // swapped until the next step), so it reuses the cached flux diagnostics
    // and only advances the surface-energy-balance state.
    const bool do_sweep = (call_site != "post_dycore");

    // One diagnostics writer for the life of the run, so its header-written
    // flag and (step, call_site, time) duplicate guard actually carry over
    // between calls.
    if (!m_diag) {
        m_diag = std::make_unique<RadiationDiagnostics>(
            rad_choice.verbosity, rad_choice.diag_file,
            rad_choice.diag_enable, rad_choice.diag_stdout_enable,
            rad_choice.diag_tagged_enable, rad_choice.diag_regtest_line_enable,
            rad_choice.diag_csv_enable, rad_choice.diag_callsite_mode,
            rad_choice.diag_dedup_tol);
    }
    RadiationDiagnostics& rad_diag = *m_diag;

    // ========================================
    // GPU-Safe ParallelFor Implementation with Cloud Fraction
    // Blending, Diffuse (Scattering) SW Flux, and Per-Level Heating Rate
    // Output (qheating_rates MultiFab)
    // ========================================

    // Initialize global diagnostics
    amrex::Real SW_surface = 0.0;
    amrex::Real SW_TOA = 0.0;
    amrex::Real SW_up_TOA = 0.0;
    amrex::Real LW_net_surface = 0.0;
    amrex::Real LW_up_TOA = 0.0;
    amrex::Real heating_rate_max = 0.0;
    amrex::Real seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real seb_residual_max  = std::numeric_limits<amrex::Real>::quiet_NaN();

    //  Prognostic SEB surface temperature and moisture diagnostics
    amrex::Real t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real t_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
    amrex::Real q_s_max  = std::numeric_limits<amrex::Real>::quiet_NaN();

    // Get state at this level (conservative variables: density, RhoTheta, etc.)
    const MultiFab& state_cons = cons_old;

    // Only compute radiation if we have valid state data
    if (state_cons.nComp() > 0 ) {

    // Trivially copyable parameter set for the device lambdas below
    // (RadChoice itself holds std::string members and cannot be captured),
    // with the sun of this call from the inputs shared with RRTMGP. The
    // incident SW at the top (SW_TOA) is a domain mean formed by the sweep,
    // since with a calendar sun it varies across the columns.
    TwoStreamParams ts_params = make_two_stream_params(rad_choice, m_rdOcp);
    if (do_sweep) { set_solar_state(ts_params, rad_choice, m_orbit, epoch_time, have_datetime, lev, nstep); }

        // Host-side storage for reduction results (will be set by device-side reduction)
        amrex::Real max_heating_global = 0.0;
        amrex::Real sw_surface_sum = 0.0;
        amrex::Real sw_up_toa_sum = 0.0;
        amrex::Real lw_net_sum = 0.0;
        amrex::Real lw_up_toa_sum = 0.0;
        amrex::Real sw_toa_sum = 0.0;
        amrex::Long n_columns_total = 0;

        // SEB residual diagnostics
        amrex::Real seb_residual_sum = 0.0;
        amrex::Long n_seb_columns = 0;

        // cloud fraction used to blend clear-sky and cloudy-column results.
        // cloud_fraction == 0.0 (default) means only the clear-sky column is
        // ever evaluated, and the blend below reduces to F = F_clear exactly.
        amrex::Real cloud_fraction = rad_choice.cloud_fraction;

        // qheating is ERF's 2-component (SW, LW) heating-rate MultiFab of this
        // level. If it is not allocated yet, the sweep still runs for the
        // diagnostics and skips the per-level heating write.
        MultiFab* qheating_mf = qheating;

        // Surface radiative fluxes for the SEB. Precedence: an LSM field when
        // the LSM exposes one; otherwise, with seb_use_radiation_fluxes, the
        // fluxes the column sweep below computes at the surface (written per
        // column by the sweep and left in place for the post-dycore call);
        // otherwise the scalar defaults.
        const bool sw_flux_from_rad = rad_choice.seb_enable &&
                                      rad_choice.seb_use_radiation_fluxes &&
                                      !lsm_has_field(lsm, lev, "sav");
        const bool lw_flux_from_rad = rad_choice.seb_enable &&
                                      rad_choice.seb_use_radiation_fluxes &&
                                      !lsm_has_field(lsm, lev, "fira");

        if (rad_choice.seb_enable) {
            fill_or_copy_seb_field(m_alb_sw[lev].get(), lsm, lev, "sfc_alb_dir_vis", rad_choice.surface_albedo_sw);
            fill_or_copy_seb_field(m_emiss_lw[lev].get(), lsm, lev, "sfc_emis", rad_choice.surface_emissivity_lw);

            // Gate t_sfc fill on prognostic mode: when seb_prognostic_enable is true,
            // t_sfc is owned and evolved by the prognostic update, not reset by fill_or_copy.
            // This prevents silently overwriting the prognostic state before the update reads it.
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(m_t_sfc[lev].get(), lsm, lev, "t_sfc", rad_choice.rad_t_sfc);
            }

            // Net absorbed shortwave: Noah-MP splits it into the canopy (sav)
            // and ground (sag) parts. Net longwave: Noah-MP's fira is the
            // net flux to the atmosphere (positive up); the SEB wants the
            // absorbed flux, so the sign flips.
            if (!sw_flux_from_rad) {
                fill_or_copy_seb_field(m_sw_flux_sfc[lev].get(), lsm, lev, "sav",
                                       rad_choice.seb_sw_flux_default, 1.0, "sag");
            }
            if (!lw_flux_from_rad) {
                fill_or_copy_seb_field(m_lw_flux_sfc[lev].get(), lsm, lev, "fira",
                                       rad_choice.seb_lw_flux_default, -1.0);
            }
            // The LSM data lists carry no sensible or latent heat flux under
            // these names, so H and LE come from the scalar defaults unless a
            // model exposes them; G is Noah-MP's grdflx when present.
            fill_or_copy_seb_field(m_hfx_sfc[lev].get(), lsm, lev, "hfx", rad_choice.seb_hfx_default);
            fill_or_copy_seb_field(m_lh_sfc[lev].get(), lsm, lev, "lh", rad_choice.seb_lh_default);
            fill_or_copy_seb_field(m_grdflx_sfc[lev].get(), lsm, lev, "grdflx", rad_choice.seb_grdflx_default);

            // Gate q_sfc fill on prognostic mode: same reasoning as t_sfc.
            if (!rad_choice.seb_prognostic_enable) {
                fill_or_copy_seb_field(m_q_sfc[lev].get(), lsm, lev, "noahmp_water_vapor_mixing_ratio_2m_vegetated", rad_choice.seb_q_sfc_default);
            }
            // No LSM exposes a deep-soil temperature or moisture in kg/kg by
            // name (Noah-MP's smstav / smstot are soil-moisture availability
            // and total column water), so the reservoir values are the
            // scalar defaults.
            m_t_deep[lev]->setVal(rad_choice.seb_t_deep_default);
            m_q_deep[lev]->setVal(rad_choice.seb_q_deep_default);
        }

        // The column sweep integrates the whole atmospheric column in one
        // kernel, so it needs every k of a box in a single pass. MFIter tiling
        // would hand it partial columns (the default CPU tile size splits z),
        // so this loop is deliberately untiled and works on valid boxes; the
        // horizontal ParallelFor below still provides the parallelism.
        const Box& rad_domain = geom.Domain();
        if (do_sweep) {
        for (MFIter mfi(state_cons, false); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.validbox();

            // A box that does not span the domain vertically would give this
            // column solver a truncated atmosphere: no beam from above, no
            // cooling to space. ERF only decomposes in z when max_grid_size_z
            // is smaller than the domain, so refuse that configuration rather
            // than return heating rates that look plausible and are wrong.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                bx.smallEnd(2) == rad_domain.smallEnd(2) &&
                bx.bigEnd(2)   == rad_domain.bigEnd(2),
                "TwoStream radiation requires grids that span the domain in z; "
                "set amr.max_grid_size_z to at least amr.n_cell in z");
            const auto& state_arr = state_cons.const_array(mfi);
            // Geometry::CellSize() is host-only; read it here and hand the
            // value to the device sweep.
            const amrex::Real dz_uniform_lev = geom.CellSize(2);

            // Nodal heights for the layer thicknesses on a stretched or
            // terrain-following grid (nullptr on a uniform grid).
            Array4<const amrex::Real> z_phys_nd_arr;
            if (z_phys_nd != nullptr) {
                z_phys_nd_arr = z_phys_nd->const_array(mfi);
            }

            // Per-column scratch of the sweep (interfaces need nlev + 1 entries),
            // shared by the clear and cloudy evaluations of the same column.
            FArrayBox scratch_fab(two_stream_scratch_box(bx), TwoStreamScratch::NCOMP);
            Array4<amrex::Real> scratch_arr = scratch_fab.array();

            //  Wire LSM surface property fields or use standalone fallback MultiFabs
            // Priority:
            // 1. If LSM is active (lsm.Get_DataIdx() returns >=0), use real LSM fields
            // 2. Otherwise, use standalone fallback MultiFabs (allocated and constant-filled from RadChoice scalars)
            // The resolve_surface_*() helpers implement the full precedence chain with finite guards

            // SW albedo: Try LSM field "sfc_alb_dir_vis" (simplified: broadband approx from vis-direct only;
            // future work: full 4-band vis/nir dir/dif support is planned).
            bool has_hetero_alb_sw = false;
            Array4<const amrex::Real> hetero_alb_sw_arr;
            {
                std::string varname_alb = "sfc_alb_dir_vis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_alb);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_alb_sw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_alb_sw = true;
                    }
                } else if (m_alb_sw[lev]) {
                    hetero_alb_sw_arr = m_alb_sw[lev]->const_array(mfi);
                    has_hetero_alb_sw = true;
                }
            }

            // LW emissivity: Try LSM field "sfc_emis"
            bool has_hetero_emiss_lw = false;
            Array4<const amrex::Real> hetero_emiss_lw_arr;
            {
                std::string varname_emiss = "sfc_emis";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_emiss);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        hetero_emiss_lw_arr = lsm_ptr->const_array(mfi);
                        has_hetero_emiss_lw = true;
                    }
                } else if (m_emiss_lw[lev]) {
                    hetero_emiss_lw_arr = m_emiss_lw[lev]->const_array(mfi);
                    has_hetero_emiss_lw = true;
                }
            }

            // Keep the surface-temperature sources separate until the column
            // kernel resolves them per cell. An LSM field can exist globally
            // while carrying an undefined sentinel in an individual column;
            // that column must still fall through to SurfaceLayer theta.
            // The kernel's precedence is valid LSM absolute temperature,
            // valid prognostic SEB absolute temperature, SurfaceLayer theta,
            // then the scalar absolute-temperature fallback.
            bool has_lsm_t_sfc = false;
            Array4<const amrex::Real> lsm_t_sfc_arr;
            bool has_seb_t_sfc = false;
            Array4<const amrex::Real> seb_t_sfc_arr;
            bool has_surface_layer = false;
            Array4<const amrex::Real> surface_layer_theta_arr;
            {
                std::string varname_t_sfc = "t_sfc";
                int lsm_idx = lsm.Get_DataIdx(lev, varname_t_sfc);
                if (lsm_idx >= 0) {
                    auto lsm_ptr = lsm.Get_Data_Ptr(lev, lsm_idx);
                    if (lsm_ptr) {
                        lsm_t_sfc_arr = lsm_ptr->const_array(mfi);
                        has_lsm_t_sfc = true;
                    }
                }
                // Prognostic SEB is a level-wide alternative to Noah t_sfc.
                // When Noah exposes t_sfc on this level, m_t_sfc is not advanced,
                // so it must not be offered as a per-cell fallback.
                if (!has_lsm_t_sfc && rad_choice.seb_prognostic_enable &&
                    rad_choice.seb_enable && m_t_sfc[lev]) {
                    seb_t_sfc_arr = m_t_sfc[lev]->const_array(mfi);
                    has_seb_t_sfc = true;
                }
                if (t_surf != nullptr) {
                    surface_layer_theta_arr = t_surf->const_array(mfi);
                    has_surface_layer = true;
                }
            }

            // Per-column latitude and longitude for the calendar sun (only
            // the WRF and metgrid initialisations fill them; otherwise the
            // erf.rad_cons_lat/lon constants apply).
            const bool has_latlon = (lat_m != nullptr) && (lon_m != nullptr);
            Array4<const amrex::Real> lat_arr, lon_arr;
            if (has_latlon) {
                lat_arr = lat_m->const_array(mfi);
                lon_arr = lon_m->const_array(mfi);
            }

            // Interface fluxes for ERF's rad_fluxes (SW up, SW down, LW up,
            // LW down at each layer's lower interface, plus the top of the
            // atmosphere in the z-ghost cell above the column), the level
            // layout RRTMGP writes. The cloudy evaluation goes to a scratch
            // FArrayBox on the same grown box and is blended like the
            // heating rates.
            const bool write_fluxes = (rad_fluxes != nullptr);
            Array4<amrex::Real> rad_flux_clear_arr;
            FArrayBox rad_flux_cloudy_fab;
            Array4<amrex::Real> rad_flux_cloudy_arr;
            if (write_fluxes) {
                rad_flux_clear_arr = rad_fluxes->array(mfi);
                if (cloud_fraction > 0.0) {
                    rad_flux_cloudy_fab.resize(two_stream_scratch_box(bx), 4);
                    rad_flux_cloudy_arr = rad_flux_cloudy_fab.array();
                }
            }


            // Surface flux arrays the sweep fills for the SEB when asked to.
            Array4<amrex::Real> sw_sfc_out;
            Array4<amrex::Real> lw_sfc_out;
            if (sw_flux_from_rad) sw_sfc_out = m_sw_flux_sfc[lev]->array(mfi);
            if (lw_flux_from_rad) lw_sfc_out = m_lw_flux_sfc[lev]->array(mfi);

            // Create a 2D box for (i,j) iteration over the horizontal extent
            // One GPU thread per (i,j) column; k-loop is sequential within each thread
            const auto& lo = bx.loVect();
            const auto& hi = bx.hiVect();
            Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

            // Count columns in this box for later averaging
            amrex::Long n_cols = static_cast<amrex::Long>(bx.length(0)) *
                                 static_cast<amrex::Long>(bx.length(1));
            n_columns_total += n_cols;

            // Clear-sky-column heating rates are written directly
            // into the real qheating_rates MultiFab when available.
            // Fall back to a throwaway local FArrayBox otherwise (keeps the
            // kernel call GPU-safe even if qheating_rates isn't allocated).
            FArrayBox qheating_fallback_fab;
            Array4<amrex::Real> qheating_clear_arr;
            if (qheating_mf != nullptr) {
                qheating_clear_arr = qheating_mf->array(mfi);
            } else {
                qheating_fallback_fab.resize(bx, 2);
                qheating_clear_arr = qheating_fallback_fab.array();
            }

            // Cloudy-column heating rates always go into a scratch
            // FArrayBox; only used/blended when cloud_fraction > 0.
            FArrayBox qheating_cloudy_fab;
            Array4<amrex::Real> qheating_cloudy_arr;
            if (cloud_fraction > 0.0) {
                qheating_cloudy_fab.resize(bx, 2);
                qheating_cloudy_arr = qheating_cloudy_fab.array();
            }

            // GPU-safe reduction using ReduceOps (per-column results aggregated on device)
            amrex::Real max_heating_box = 0.0;
            amrex::Real sw_sum_box = 0.0;
            amrex::Real sw_up_sum_box = 0.0;
            amrex::Real lw_sum_box = 0.0;
            amrex::Real lw_up_sum_box = 0.0;
            amrex::Real sw_toa_sum_box = 0.0;

            // Device-side reduction: max heating, sums of the surface and
            // top-of-atmosphere fluxes
            ReduceOps<ReduceOpMax, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_ops;
            ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_ops);

            using ReduceTuple = typename decltype(reduce_data)::Type;

         // Launch parallel kernel over (i,j) columns
            reduce_ops.eval(xy_box, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ReduceTuple
                {
                    // Clear-sky column (always evaluated; this is the sole
                    // contributor when cloud_fraction == 0.0, matching earlier behavior)
                    amrex::Real max_heating_clear = 0.0;
                    amrex::Real sw_flux_clear = 0.0;
                    amrex::Real sw_up_clear = 0.0;
                    amrex::Real lw_net_clear = 0.0;
                    amrex::Real lw_up_clear = 0.0;
                    amrex::Real sw_toa_clear = 0.0;
                    vertical_two_stream_sweep(
                        i, j, bx, dz_uniform_lev, state_arr, ts_params, /*cloudy=*/false,
                        qheating_clear_arr,
                        max_heating_clear, sw_flux_clear, sw_up_clear, lw_net_clear, lw_up_clear,
                        sw_toa_clear,
                        z_phys_nd_arr, scratch_arr,
                        has_hetero_alb_sw, &hetero_alb_sw_arr,
                        has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                        has_lsm_t_sfc, &lsm_t_sfc_arr,
                        has_seb_t_sfc, &seb_t_sfc_arr,
                        has_surface_layer, &surface_layer_theta_arr,
                        has_latlon, &lat_arr, &lon_arr,
                        write_fluxes ? &rad_flux_clear_arr : nullptr);

                    amrex::Real max_heating_col = max_heating_clear;
                    amrex::Real sw_flux_col = sw_flux_clear;
                    amrex::Real sw_up_col = sw_up_clear;
                    amrex::Real lw_net_col = lw_net_clear;
                    amrex::Real lw_up_col = lw_up_clear;

                    // Cloudy column only needs to be evaluated when there is a
                    // nonzero cloud fraction; this keeps the cloud_fraction==0
                    // path numerically and computationally identical.
                    if (cloud_fraction > 0.0) {
                         amrex::Real max_heating_cloudy = 0.0;
                         amrex::Real sw_flux_cloudy = 0.0;
                         amrex::Real sw_up_cloudy = 0.0;
                         amrex::Real lw_net_cloudy = 0.0;
                         amrex::Real lw_up_cloudy = 0.0;
                         amrex::Real sw_toa_cloudy = 0.0;
                         vertical_two_stream_sweep(
                            i, j, bx, dz_uniform_lev, state_arr, ts_params, /*cloudy=*/true,
                            qheating_cloudy_arr,
                            max_heating_cloudy, sw_flux_cloudy, sw_up_cloudy, lw_net_cloudy, lw_up_cloudy,
                            sw_toa_cloudy,
                            z_phys_nd_arr, scratch_arr,
                            has_hetero_alb_sw, &hetero_alb_sw_arr,
                            has_hetero_emiss_lw, &hetero_emiss_lw_arr,
                            has_lsm_t_sfc, &lsm_t_sfc_arr,
                            has_seb_t_sfc, &seb_t_sfc_arr,
                            has_surface_layer, &surface_layer_theta_arr,
                            has_latlon, &lat_arr, &lon_arr,
                            write_fluxes ? &rad_flux_cloudy_arr : nullptr);

                        // Blend clear-sky and cloudy-column results
                        sw_flux_col = (1.0 - cloud_fraction) * sw_flux_clear +
                                      cloud_fraction * sw_flux_cloudy;
                        sw_up_col = (1.0 - cloud_fraction) * sw_up_clear +
                                    cloud_fraction * sw_up_cloudy;
                        lw_net_col = (1.0 - cloud_fraction) * lw_net_clear +
                                     cloud_fraction * lw_net_cloudy;
                        lw_up_col = (1.0 - cloud_fraction) * lw_up_clear +
                                    cloud_fraction * lw_up_cloudy;
                        max_heating_col = std::max(max_heating_clear, max_heating_cloudy);

                        // Blend per-level heating rates in place
                        // into qheating_clear_arr (which is the real output
                        // MultiFab when qheating_mf != nullptr).
                        int kmin = bx.smallEnd(2);
                        int kmax = bx.bigEnd(2);
                        for (int k = kmin; k <= kmax; ++k) {
                            for (int comp = 0; comp < 2; ++comp) {
                                amrex::Real q_clear_val = qheating_clear_arr(i, j, k, comp);
                                amrex::Real q_cloudy_val = qheating_cloudy_arr(i, j, k, comp);
                                qheating_clear_arr(i, j, k, comp) =
                                    (1.0 - cloud_fraction) * q_clear_val +
                                    cloud_fraction * q_cloudy_val;
                            }
                        }
                        if (write_fluxes) {
                            // nlev + 1 levels: the top interface sits at kmax + 1
                            for (int k = kmin; k <= kmax + 1; ++k) {
                                for (int comp = 0; comp < 4; ++comp) {
                                    rad_flux_clear_arr(i, j, k, comp) =
                                        (1.0 - cloud_fraction) * rad_flux_clear_arr(i, j, k, comp) +
                                        cloud_fraction * rad_flux_cloudy_arr(i, j, k, comp);
                                }
                            }
                        }
                    }

                    // Surface fluxes for the SEB: absorbed shortwave, and the
                    // absorbed longwave, which is minus the net (up - down).
                    if (sw_flux_from_rad) sw_sfc_out(i, j, 0) = sw_flux_col;
                    if (lw_flux_from_rad) lw_sfc_out(i, j, 0) = -lw_net_col;

                    // The incident TOA flux is the same for both evaluations.
                    // Return tuple for reduction
                    return {max_heating_col, sw_flux_col, sw_up_col, lw_net_col, lw_up_col, sw_toa_clear};
                }
            );

            // Copy results from device to host
            amrex::Gpu::synchronize();
            auto reduce_tuple = reduce_data.value(reduce_ops);
            max_heating_box = amrex::get<0>(reduce_tuple);
            sw_sum_box = amrex::get<1>(reduce_tuple);
            sw_up_sum_box = amrex::get<2>(reduce_tuple);
            lw_sum_box = amrex::get<3>(reduce_tuple);
            lw_up_sum_box = amrex::get<4>(reduce_tuple);
            sw_toa_sum_box = amrex::get<5>(reduce_tuple);

            // Accumulate box results into global results
            max_heating_global = std::max(max_heating_global, max_heating_box);
            sw_surface_sum += sw_sum_box;
            sw_up_toa_sum += sw_up_sum_box;
            lw_net_sum += lw_sum_box;
            lw_up_toa_sum += lw_up_sum_box;
            sw_toa_sum += sw_toa_sum_box;
        }
        // Every accumulator above is rank-local. Reduce before forming means
        // and maxima, so the diagnostics describe the whole domain and do
        // not change with the decomposition. These are collective calls; the
        // conditions around them are input-driven and identical on all ranks.
        if (do_sweep) {
            amrex::Real sums[5] = {sw_surface_sum, sw_up_toa_sum, lw_net_sum, lw_up_toa_sum, sw_toa_sum};
            ParallelDescriptor::ReduceRealSum(sums, 5);
            sw_surface_sum = sums[0]; sw_up_toa_sum = sums[1];
            lw_net_sum = sums[2];     lw_up_toa_sum = sums[3];
            sw_toa_sum = sums[4];
            ParallelDescriptor::ReduceLongSum(n_columns_total);
            ParallelDescriptor::ReduceRealMax(max_heating_global);
        }
        } // do_sweep

        // The heating rates go straight into the RhoTheta source term, so a
        // non-finite value here becomes a non-finite state one step later.
        if (do_sweep && qheating_mf != nullptr &&
            (qheating_mf->contains_nan(0, 2, 0) || qheating_mf->contains_inf(0, 2, 0))) {
            amrex::Abort("TwoStreamRadiation: non-finite heating rate after the column sweep at level " +
                         std::to_string(lev) + ", step " + std::to_string(nstep));
        }

         // Warn if diagnostic is requested but SEB infrastructure isn't enabled
        if (rad_choice.seb_diagnostic_enable && !rad_choice.seb_enable) {
            static bool warned_seb_misconfig = false;
            if (!warned_seb_misconfig && ParallelDescriptor::IOProcessor()) {
                Print() << "WARNING: erf.radiation.seb_diagnostic_enable=true but "
                           "seb_enable=false; SEB residual diagnostics will report NaN. "
                           "Set erf.radiation.seb_enable=true to enable SEB field "
                           "population.\n";
                warned_seb_misconfig = true;
            }
        }
        // Compute SEB residual diagnostics if enabled
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable) {
            seb_residual_max = 0.0;
            // Second loop over boxes to compute SEB residual from populated SEB
            // MultiFabs. Untiled: the work below is per surface column, and a
            // tiled iteration would visit each (i,j) once per z tile.
            for (MFIter mfi(state_cons, false); mfi.isValid(); ++mfi) {
                const Box& bx = mfi.validbox();
                const auto& lo = bx.loVect();
                const auto& hi = bx.hiVect();
                Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

                // Get SEB field arrays
                Array4<const amrex::Real> sw_flux_arr = m_sw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lw_flux_arr = m_lw_flux_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> hfx_arr = m_hfx_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> lh_arr = m_lh_sfc[lev]->const_array(mfi);
                Array4<const amrex::Real> grdflx_arr = m_grdflx_sfc[lev]->const_array(mfi);

                // Count columns and compute residuals
                amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                        static_cast<amrex::Long>(bx.length(1));
                amrex::Real residual_sum_box = 0.0;
                amrex::Real residual_max_box = 0.0;

                // GPU-safe reduction for SEB residuals
                ReduceOps<ReduceOpSum, ReduceOpMax> seb_reduce_ops;
                ReduceData<amrex::Real, amrex::Real> seb_reduce_data(seb_reduce_ops);

                using SEBReduceTuple = typename decltype(seb_reduce_data)::Type;

                seb_reduce_ops.eval(xy_box, seb_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> SEBReduceTuple {
                        amrex::Real sw_net = sw_flux_arr(i, j, 0);
                        amrex::Real lw_net = lw_flux_arr(i, j, 0);
                        amrex::Real hfx = hfx_arr(i, j, 0);
                        amrex::Real lh = lh_arr(i, j, 0);
                        amrex::Real grdflx = grdflx_arr(i, j, 0);

                        // Compute residual using helper function
                        amrex::Real residual = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);

                        // Return sum and abs(max) of residual
                        return {residual, std::abs(residual)};
                    }
                );

                // Copy results from device to host
                amrex::Gpu::synchronize();
                auto seb_reduce_tuple = seb_reduce_data.value(seb_reduce_ops);
                residual_sum_box = amrex::get<0>(seb_reduce_tuple);
                residual_max_box = amrex::get<1>(seb_reduce_tuple);

                // Accumulate into global results
                seb_residual_sum += residual_sum_box;
                seb_residual_max = std::max(seb_residual_max, residual_max_box);
                n_seb_columns += n_cols_box;
            }
            ParallelDescriptor::ReduceRealSum(seb_residual_sum);
            ParallelDescriptor::ReduceRealMax(seb_residual_max);
            ParallelDescriptor::ReduceLongSum(n_seb_columns);
        }

        //  Prognostic SEB surface temperature and moisture evolution
        // Only run if prognostic mode is enabled and Noah-MP is NOT driving LSM at this level
        if (rad_choice.seb_prognostic_enable && rad_choice.seb_enable &&
            call_site == "post_dycore") {
            // Check if Noah-MP is active at this level by attempting to get the LSM t_sfc field
            std::string varname_t_sfc_prog = "t_sfc";
            int lsm_idx_t_sfc = lsm.Get_DataIdx(lev, varname_t_sfc_prog);
            bool noahmp_active = (lsm_idx_t_sfc >= 0);

            if (!noahmp_active) {
                // Noah-MP is NOT active; proceed with prognostic update

                // Initialize diagnostics for T_s and q_s
                amrex::Real t_s_sum = 0.0;
                amrex::Real t_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Real q_s_sum = 0.0;
                amrex::Real q_s_max_val = -std::numeric_limits<amrex::Real>::max();
                amrex::Long n_prog_columns = 0;

                // Third loop over boxes for prognostic SEB update. Untiled for
                // the same reason as above, and here it also matters for
                // correctness: the force-restore update is applied in place, so
                // visiting a column twice would advance it twice in one step.
                for (MFIter mfi(state_cons, false); mfi.isValid(); ++mfi) {
                    const Box& bx = mfi.validbox();
                    const auto& lo = bx.loVect();
                    const auto& hi = bx.hiVect();
                    Box xy_box(IntVect(lo[0], lo[1], 0), IntVect(hi[0], hi[1], 0));

                    // Get SEB field arrays (read-only)
                    Array4<const amrex::Real> sw_flux_arr = m_sw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lw_flux_arr = m_lw_flux_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> hfx_arr = m_hfx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> lh_arr = m_lh_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> grdflx_arr = m_grdflx_sfc[lev]->const_array(mfi);
                    Array4<const amrex::Real> t_deep_arr = m_t_deep[lev]->const_array(mfi);
                    Array4<const amrex::Real> q_deep_arr = m_q_deep[lev]->const_array(mfi);

                    // Get SEB state arrays (read-write for prognostic update)
                    Array4<amrex::Real> t_s_arr = m_t_sfc[lev]->array(mfi);
                    Array4<amrex::Real> q_s_arr = m_q_sfc[lev]->array(mfi);
                    // Count columns and prepare for reductions
                    amrex::Long n_cols_box = static_cast<amrex::Long>(bx.length(0)) *
                                             static_cast<amrex::Long>(bx.length(1));
                    amrex::Real t_s_sum_box = 0.0;
                    amrex::Real t_s_max_box = -std::numeric_limits<amrex::Real>::max();
                    amrex::Real q_s_sum_box = 0.0;
                    amrex::Real q_s_max_box = -std::numeric_limits<amrex::Real>::max();

                    // GPU-safe update for prognostic T_s and q_s with reductions
                    ReduceOps<ReduceOpSum, ReduceOpMax, ReduceOpSum, ReduceOpMax> prog_reduce_ops;
                    ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real> prog_reduce_data(prog_reduce_ops);

                    using ProgReduceTuple = typename decltype(prog_reduce_data)::Type;

                    prog_reduce_ops.eval(xy_box, prog_reduce_data,
                            [=,  C_s=rad_choice.seb_surface_heat_capacity,
                            tau=rad_choice.seb_restore_timescale_s,
                            d_s=rad_choice.seb_moisture_layer_depth_m,
                            tau_q=rad_choice.seb_moisture_restore_timescale_s,
                            t_min=rad_choice.seb_prognostic_t_min_k,
                            t_max=rad_choice.seb_prognostic_t_max_k,
                            q_min=rad_choice.seb_prognostic_q_min,
                            q_max=rad_choice.seb_prognostic_q_max]
                            AMREX_GPU_DEVICE (int i, int j, int /*k_unused*/) -> ProgReduceTuple {
                            amrex::Real t_s_old = t_s_arr(i, j, 0);
                            amrex::Real q_s_old = q_s_arr(i, j, 0);

                            // Read forcing data
                            amrex::Real sw_net = sw_flux_arr(i, j, 0);
                            amrex::Real lw_net = lw_flux_arr(i, j, 0);
                            amrex::Real hfx = hfx_arr(i, j, 0);
                            amrex::Real lh = lh_arr(i, j, 0);
                            amrex::Real grdflx = grdflx_arr(i, j, 0);
                            amrex::Real t_deep_val = t_deep_arr(i, j, 0);
                            amrex::Real q_deep_val = q_deep_arr(i, j, 0);

                            // Compute SEB residual
                            amrex::Real seb_res = diagnose_seb_residual(sw_net, lw_net, hfx, lh, grdflx);

                            // Compute tendencies
                            amrex::Real dT_s_dt = prognostic_dTs_dt(seb_res, t_s_old, t_deep_val,
                                                                     C_s, tau);
                            amrex::Real dq_s_dt = prognostic_dqs_dt(lh, q_s_old, q_deep_val,
                                                                     d_s, tau_q);

                            // Explicit Euler update over the step size
                            amrex::Real t_s_new = t_s_old + dt_step * dT_s_dt;
                            amrex::Real q_s_new = q_s_old + dt_step * dq_s_dt;

                            // Clamp to valid ranges
                            t_s_new = amrex::max(t_min, amrex::min(t_max, t_s_new));
                            q_s_new = amrex::max(q_min, amrex::min(q_max, q_s_new));

                            // Write back updated values (this modifies the device array)
                            t_s_arr(i, j, 0) = t_s_new;
                            q_s_arr(i, j, 0) = q_s_new;

                            // Return for reduction: sum T_s, max T_s, sum q_s, max q_s
                            return {t_s_new, std::abs(t_s_new), q_s_new, std::abs(q_s_new)};
                        });

                    // Copy results from device to host
                    amrex::Gpu::synchronize();
                    auto prog_reduce_tuple = prog_reduce_data.value(prog_reduce_ops);
                    t_s_sum_box = amrex::get<0>(prog_reduce_tuple);
                    t_s_max_box = amrex::get<1>(prog_reduce_tuple);
                    q_s_sum_box = amrex::get<2>(prog_reduce_tuple);
                    q_s_max_box = amrex::get<3>(prog_reduce_tuple);

                    // Accumulate into global results
                    t_s_sum += t_s_sum_box;
                    t_s_max_val = std::max(t_s_max_val, t_s_max_box);
                    q_s_sum += q_s_sum_box;
                    q_s_max_val = std::max(q_s_max_val, q_s_max_box);
                    n_prog_columns += n_cols_box;
                }

                // The force-restore kernels return a zero tendency on any
                // non-finite input, so a NaN here means the state itself was
                // corrupted (a bad restart file, or an overwrite elsewhere).
                if (m_t_sfc[lev]->contains_nan(0, 1, 0) || m_q_sfc[lev]->contains_nan(0, 1, 0)) {
                    amrex::Abort("TwoStreamRadiation: non-finite surface temperature or moisture "
                                 "after the force-restore update at level " + std::to_string(lev) +
                                 ", step " + std::to_string(nstep));
                }

                {
                    amrex::Real sums[2] = {t_s_sum, q_s_sum};
                    ParallelDescriptor::ReduceRealSum(sums, 2);
                    t_s_sum = sums[0]; q_s_sum = sums[1];
                    ParallelDescriptor::ReduceRealMax(t_s_max_val);
                    ParallelDescriptor::ReduceRealMax(q_s_max_val);
                    ParallelDescriptor::ReduceLongSum(n_prog_columns);
                }

                // Compute mean values from sums
                if (n_prog_columns > 0) {
                    t_s_mean = t_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    t_s_max = t_s_max_val;
                    q_s_mean = q_s_sum / static_cast<amrex::Real>(n_prog_columns);
                    q_s_max = q_s_max_val;
                } else {
                    t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                    q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                }
            } else {
                // Noah-MP is active; skip prognostic update for this level
                // Leave t_s and q_s as populated by LSM passthrough
                t_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                t_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
                q_s_max = std::numeric_limits<amrex::Real>::quiet_NaN();
            }
        }
        // equivalent to a single-column value for spatially UNIFORM atmospheres
        // (as in the current SW_ClearSky_Analytical / LW_Isothermal RegTests).
        // Cloud layer and scattering tests are ALSO spatially uniform (identical
        // cloud/tau/scattering parameters applied to every column), so the
        // domain-averaged value still equals the true single-column flux there.
        // True horizontal heterogeneity (e.g., patchy clouds varying by column)
        // remains deferred to future work.
        if (do_sweep) {
            if (n_columns_total > 0) {
                const amrex::Real inv_n = 1.0 / static_cast<amrex::Real>(n_columns_total);
                SW_surface     = sw_surface_sum * inv_n;
                SW_TOA         = sw_toa_sum * inv_n;
                SW_up_TOA      = sw_up_toa_sum * inv_n;
                LW_net_surface = lw_net_sum * inv_n;
                LW_up_TOA      = lw_up_toa_sum * inv_n;
            }
            heating_rate_max = max_heating_global;
            m_flux_diag[lev] = FluxDiag{SW_surface, SW_TOA, SW_up_TOA,
                                                         LW_net_surface, LW_up_TOA,
                                                         heating_rate_max};
        } else {
            const FluxDiag& cached = m_flux_diag[lev];
            SW_surface       = cached.SW_surface;
            SW_TOA           = cached.SW_TOA;
            SW_up_TOA        = cached.SW_up_TOA;
            LW_net_surface   = cached.LW_net_surface;
            LW_up_TOA        = cached.LW_up_TOA;
            heating_rate_max = cached.heating_rate_max;
        }

        // Compute SEB residual mean from sum
        if (rad_choice.seb_diagnostic_enable && rad_choice.seb_enable && n_seb_columns > 0) {
            seb_residual_mean = seb_residual_sum / static_cast<amrex::Real>(n_seb_columns);
        } else {
            // When feature is disabled, use NaN for backward compatibility
            seb_residual_mean = std::numeric_limits<amrex::Real>::quiet_NaN();
            seb_residual_max = std::numeric_limits<amrex::Real>::quiet_NaN();
        }

    }

    // Logging output
    if (rad_choice.verbosity >= 1 && ParallelDescriptor::IOProcessor()) {
        Print() << "Radiation diagnostics at step " << nstep << ":\n"
                << "  SW TOA = " << SW_TOA << " W/m^2\n"
                << "  SW surface = " << SW_surface << " W/m^2\n"
                << "  SW up (TOA) = " << SW_up_TOA << " W/m^2\n"
                << "  LW net (surface) = " << LW_net_surface << " W/m^2\n"
                << "  LW up (TOA) = " << LW_up_TOA << " W/m^2\n"
                << "  Max heating rate = " << heating_rate_max << " K/s\n";
        if (rad_choice.seb_diagnostic_enable && std::isfinite(seb_residual_mean)) {
            Print() << "  SEB residual (mean) = " << seb_residual_mean << " W/m^2\n"
                    << "  SEB residual (max) = " << seb_residual_max << " W/m^2\n";
        }
        if (rad_choice.seb_prognostic_enable && std::isfinite(t_s_mean)) {
            Print() << "  Surface temperature (mean) = " << t_s_mean << " K\n"
                    << "  Surface temperature (max) = " << t_s_max << " K\n"
                    << "  Surface moisture (mean) = " << q_s_mean << " kg/kg\n"
                    << "  Surface moisture (max) = " << q_s_max << " kg/kg\n";
        }
    }

    rad_diag.append(nstep, time, call_site, SW_surface, SW_TOA,
                    SW_up_TOA, LW_net_surface, LW_up_TOA, heating_rate_max,
                    seb_residual_mean, seb_residual_max,
                    t_s_mean, t_s_max, q_s_mean, q_s_max);
}
