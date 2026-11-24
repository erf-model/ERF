#ifndef ERF_FORCE_RESTORE_HYDROLOGY_H_
#define ERF_FORCE_RESTORE_HYDROLOGY_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <cmath>

#include <ERF_MicrophysicsConstants.H>   // rhor, the density of water
#include <ERF_NumericalConstants.H>     // zero, one

/**
 * @file ERF_ForceRestoreHydrology.H
 * @brief Water of the two-stream force-restore surface: rain in, runoff, a prognostic deep
 *        reservoir with gravitational drainage, and a prognostic deep temperature.
 *
 * The force-restore balance (ERF_SimplifiedSEB.H) carries q_s, the volumetric water content
 * [m^3/m^3] of a surface layer d_s deep (erf.radiation.seb_moisture_layer_depth_m), restored
 * toward q_deep over tau_q. These functions add the terms of Noilhan and Planton (1989,
 * "NP89") and Mahfouf and Noilhan (1996, "MN96") that it leaves out, each behind its own
 * input:
 *
 * - Rain (erf.radiation.seb_precipitation_source = microphysics). P [kg/m^2/s] enters the
 *   surface layer as in NP89 Eq. (11) and Deardorff (1978) Eq. (12), with depth d_1 = d_s:
 *     dq_s/dt += P / (rho_w d_s).
 *   NP89 Eq. (11) multiplies (P - E) by one coefficient C1. The existing balance takes E out
 *   of q_s with C1 = 1 (prognostic_dqs_dt), so rain uses C1 = 1 too. NP89 Eq. (19) and Table 2
 *   give C1 = 0.08-0.34 at saturation, rising as the soil dries, and D78 Eq. (13) gives
 *   C1 = 0.5-14. With C1 = 1 the 10 cm layer fills faster than with NP89's C1 near saturation,
 *   so heavy rain runs off sooner than in NP89.
 * - Runoff (with rain). Water above the porosity w_sat runs off (NP89 Sec. 4a: "Runoff occurs
 *   when w_g or w_2 exceeds the saturation value w_sat"). The surface layer's runoff that
 *   leaves the column is at most this step's rain: water above w_sat that was already in q_s
 *   (a restore overshoot, or a restart from a run without hydrology) is clamped off q_s and
 *   does not come out of the column.
 * - Deep reservoir (erf.radiation.seb_deep_prognostic_enable). q_deep becomes w_2, the mean
 *   water content of the soil column of depth d_2 (erf.radiation.seb_deep_layer_depth_m),
 *   which contains the surface layer (NP89 Eq. 12, MN96 Eq. 4):
 *     dq_deep/dt = (P - E) / (rho_w d_2).
 *   Since w_2 is the mean of the whole column, the restore term of q_s moves water within the
 *   column and does not appear here; the column's water is rho_w d_2 q_deep.
 * - Deep temperature (with the deep reservoir). T_2 relaxes toward T_s over tau
 *   (erf.radiation.seb_restore_timescale_s), NP89 Eq. (6):
 *     dT_2/dt = (T_s - T_2) / tau.
 * - Gravitational drainage (erf.radiation.seb_drainage_enable). Water above field capacity
 *   w_fc leaves the bottom of the column, MN96 Eq. (4):
 *     dq_deep/dt += -(C_3 / tau_1d) max(0, q_deep - w_fc),
 *   with C_3 = (C_3 d_2) / d_2 (MN96 Eq. 16 and Table 1, which tabulates C_3 for d_2 = 1 m).
 *   MN96 Eq. (15) makes C_3 proportional to tau, and Table 1 is computed with tau = 1 day
 *   (MN96 appendix), so the drainage uses tau_1d = 86400 s, not seb_restore_timescale_s. The
 *   drainage rate then does not depend on that input.
 *
 * The relaxation terms (drainage, deep temperature) are linear in the state, so each is
 * integrated exactly over the step, x_new - x_eq = (x_old - x_eq) exp(-k dt); the update
 * cannot overshoot whatever the step. The source terms (rain, evaporation) are forward Euler,
 * as the existing balance is.
 *
 * References:
 * - Deardorff, J. W., 1978: Efficient prediction of ground surface temperature and moisture,
 *   with inclusion of a layer of vegetation. J. Geophys. Res., 83, 1889-1903.
 * - Noilhan, J., and S. Planton, 1989: A simple parameterization of land surface processes
 *   for meteorological models. Mon. Wea. Rev., 117, 536-549.
 * - Mahfouf, J.-F., and J. Noilhan, 1996: Inclusion of gravitational drainage in a land
 *   surface scheme based on the force-restore method. J. Appl. Meteor., 35, 987-992.
 */

/**
 * @brief Parameters of one force-restore hydrology step (from erf.radiation.seb_*).
 */
/// The restore period [s] with which MN96 Table 1 and Eq. (16) give C_3 (MN96 appendix).
constexpr amrex::Real force_restore_c3_tau = 86400.0;

struct ForceRestoreHydrologyParams
{
    amrex::Real d_s   = 0.1;      ///< Surface layer depth d_1 [m] (seb_moisture_layer_depth_m)
    amrex::Real d_2   = 1.0;      ///< Column depth d_2 [m] (seb_deep_layer_depth_m)
    amrex::Real tau   = 86400.0;  ///< Restore period of T_2 [s] (seb_restore_timescale_s)
    amrex::Real w_fc  = 0.0;      ///< Field capacity [m^3/m^3] (seb_soil_moisture_fc)
    amrex::Real w_sat = 1.0;      ///< Porosity [m^3/m^3] (seb_soil_porosity)
    amrex::Real c3    = 0.0;      ///< Drainage coefficient C_3 at depth d_2 [-] (seb_drainage_c3_m / d_2)
    amrex::Real q_min = 0.0;      ///< Lower bound of q_s (seb_prognostic_q_min)
    amrex::Real q_max = 1.0;      ///< Upper bound of q_s (seb_prognostic_q_max)
    bool rain  = false;           ///< seb_precipitation_source != none
    bool deep  = false;           ///< seb_deep_prognostic_enable
    bool drain = false;           ///< seb_drainage_enable (needs deep)
};

/**
 * @brief Water that left a column during one step [kg/m^2].
 */
struct ForceRestoreHydrologyFluxes
{
    amrex::Real runoff   = 0.0;  ///< Above-porosity water that left the column: rain that did not infiltrate, and column excess
    amrex::Real drainage = 0.0;  ///< Gravitational drainage out of the column bottom
    amrex::Real deficit  = 0.0;  ///< Water the column lacked when evaporation would take q_deep below 0
};

/**
 * @brief Rain rate [kg/m^2/s] from the change of the surface precipitation accumulation over
 *        one step.
 *
 * @param[in] accum_end   accumulation at the end of the step [kg/m^2]
 * @param[in] accum_start accumulation at the start of the step [kg/m^2]
 * @param[in] dt          the step [s]
 * @return (end - start) / dt; zero for a non-positive step, a non-positive change (round-off
 *         or a reset accumulator) and any non-finite input.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real force_restore_rain_rate (amrex::Real accum_end, amrex::Real accum_start,
                                     amrex::Real dt) noexcept
{
    if (!(dt > amrex::Real(0.0))) { return amrex::Real(0.0); }
    const amrex::Real d = accum_end - accum_start;
    if (!(d > amrex::Real(0.0)) || !amrex::Math::isfinite(d)) { return amrex::Real(0.0); }
    return d / dt;
}

/**
 * @brief One step of the force-restore hydrology for one column, after the balance has
 *        updated T_s and q_s (prognostic_dTs_dt, prognostic_dqs_dt).
 *
 * In order: rain into the surface layer and the surface runoff above the porosity (at most
 * the step's rain: water that never infiltrates); the infiltrated rain minus the evaporation
 * into the column, then its drainage; the deep temperature toward the updated T_s; runoff of
 * the column above the porosity; the bounds of q_s and q_deep. The surface runoff is taken
 * before the drainage, so water that never entered the soil does not drain.
 *
 * With the deep reservoir, the column then satisfies exactly
 *   rho_w d_2 (q_deep_new - q_deep_old) = (P - E) dt - drainage - runoff + deficit.
 *
 * @param[in,out] q_s    surface-layer water content [m^3/m^3]
 * @param[in,out] q_deep column water content [m^3/m^3]; changed only with prm.deep
 * @param[in,out] t_deep deep temperature [K]; changed only with prm.deep
 * @param[in]     t_s    surface temperature after this step's update [K]
 * @param[in]     p      rain reaching the soil [kg/m^2/s]; used only with prm.rain
 * @param[in]     e      evaporation LE / L_v [kg/m^2/s], positive upward
 * @param[in]     dt     the step [s]
 * @param[in]     prm    parameters
 * @return the water that left the column, and any deficit
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ForceRestoreHydrologyFluxes
force_restore_hydrology_step (amrex::Real& q_s, amrex::Real& q_deep, amrex::Real& t_deep,
                              amrex::Real t_s, amrex::Real p, amrex::Real e, amrex::Real dt,
                              const ForceRestoreHydrologyParams& prm) noexcept
{
    ForceRestoreHydrologyFluxes out;
    if (!prm.rain) { p = zero; }
    if (!amrex::Math::isfinite(p)) { p = zero; }
    if (!amrex::Math::isfinite(e)) { e = zero; }

    // Rain into the surface layer: NP89 Eq. (11), first term, with C1 = 1 and d_1 = d_s.
    // (The balance already removed E from q_s in prognostic_dqs_dt.)
    // Runoff above saturation (NP89 Sec. 4a) from the surface layer: at most this step's rain
    // runs off there, since that water never infiltrated. Water above w_sat that q_s held
    // before the rain is clamped off q_s and stays in the column.
    amrex::Real surface_runoff = zero;
    if (prm.rain) {
        q_s += dt * p / (rhor * prm.d_s);
        if (q_s > prm.w_sat) {
            surface_runoff = amrex::min(rhor * prm.d_s * (q_s - prm.w_sat), p * dt);
            q_s = prm.w_sat;
            out.runoff += surface_runoff;
        }
    }

    if (prm.deep) {
        // The column's water: NP89 Eq. (12), MN96 Eq. (4), first term, with the rain that
        // infiltrated.
        q_deep += (dt * (p - e) - surface_runoff) / (rhor * prm.d_2);

        // Gravitational drainage toward field capacity: MN96 Eq. (4), second term, exact
        // over the step for the linear relaxation.
        if (prm.drain && q_deep > prm.w_fc) {
            const amrex::Real removed = (q_deep - prm.w_fc) *
                                        (one - std::exp(-prm.c3 * dt / force_restore_c3_tau));
            q_deep -= removed;
            out.drainage = rhor * prm.d_2 * removed;
        }

        // Deep temperature: NP89 Eq. (6), exact over the step.
        t_deep += (t_s - t_deep) * (one - std::exp(-dt / prm.tau));

        // Runoff of the column above saturation.
        if (q_deep > prm.w_sat) {
            out.runoff += rhor * prm.d_2 * (q_deep - prm.w_sat);
            q_deep = prm.w_sat;
        }
    }

    // NP89 Eq. (12): 0 < w_2. An evaporation the column cannot supply is recorded, not hidden:
    // the latent heat flux has already taken it, so a deficit is water the coupled system did
    // not conserve (LE is not limited by q_deep).
    if (prm.deep && q_deep < zero) {
        out.deficit = -rhor * prm.d_2 * q_deep;
        q_deep = zero;
    }

    q_s = amrex::max(prm.q_min, amrex::min(prm.q_max, q_s));
    return out;
}

#endif // ERF_FORCE_RESTORE_HYDROLOGY_H_
