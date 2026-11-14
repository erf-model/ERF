#ifndef ERF_FORCE_RESTORE_SNOW_H_
#define ERF_FORCE_RESTORE_SNOW_H_

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>

#include <ERF_MicrophysicsConstants.H>   // lat_ice, the latent heat of fusion
#include <ERF_NumericalConstants.H>      // zero

/**
 * @file ERF_ForceRestoreSnow.H
 * @brief One bulk snow layer on the two-stream force-restore surface
 *        (erf.radiation.seb_snow_enable): the snow routine of the ISBA scheme as Douville,
 *        Royer and Mahfouf (1995, "DRM95") describe it in their Sec. 2.1, for bare ground
 *        (vegetation fraction 0).
 *
 * The state is the snow water equivalent W_n [kg/m^2] (1 kg/m^2 of water is 1 mm).
 *
 * - Mass (DRM95 Eq. 6):
 *     dW_n/dt = P_n - E_n - F_n,
 *   with P_n the frozen precipitation, E_n the sublimation and F_n the melt.
 * - Cover fraction (DRM95 Eq. 7):
 *     p_ng = W_n / (W_n + W_crn),   W_crn = 10 kg/m^2 (erf.radiation.seb_snow_critical_swe_kg_m2).
 * - Albedo (DRM95 Eq. 8 with veg = 0):
 *     alpha_t = (1 - p_ng) alpha_g + p_ng alpha_n,
 *   alpha_n = 0.70 (erf.radiation.seb_snow_albedo), "the uniform and constant value of the
 *   snow albedo".
 * - Sublimation. DRM95 states that p_ng partitions the surface evaporation between the soil
 *   and the snow, without an equation. Here E_n = p_ng E for an upward E = LE / L_v, at most
 *   the snow there is; the soil keeps the rest. A downward flux (dew, frost) stays with the
 *   soil, as in the balance without snow. The latent heat flux itself is not changed: it
 *   uses L_v over the snow too.
 * - Melt (DRM95 Eq. 12). The balance gives a first estimate T_s* of the surface temperature.
 *   Where snow lies and T_s* > T_00 = 273.16 K,
 *     F_n dt = (T_s* - T_00) / (C_t L_f),
 *   at most the snow there is. C_t [K m^2/J] is the inverse of the surface heat capacity
 *   C_s [J/m^2/K] (erf.radiation.seb_surface_heat_capacity). The melt takes L_f F_n dt from
 *   the surface, so T_s = T_s* - L_f F_n dt / C_s: T_00 while snow remains ("the surface
 *   temperature can only exceed the melting point when all the snow has melted on the
 *   ground", DRM95 Sec. 2.1).
 * - Melt water and the sublimated water the balance took from the soil go to the soil with
 *   the rain (ERF_ForceRestoreHydrology.H).
 *
 * The update is forward Euler, in the order snowfall, sublimation, melt, and satisfies
 * to round-off W_n_new - W_n_old = snowfall - sublimation - melt.
 *
 * Not here (DRM95 Sec. 2.2, their new scheme): prognostic snow density and albedo, the snow
 * heat capacity, the vegetation masking, and the snow roughness length (DRM95 Eqs. 9-10).
 *
 * Reference:
 * - Douville, H., J.-F. Royer, and J.-F. Mahfouf, 1995: A new snow parameterization for the
 *   Meteo-France climate model. Part I: validation in stand-alone experiments. Climate
 *   Dynamics, 12, 21-35.
 */

/// The melting point T_00 [K] of DRM95 Eq. (12).
constexpr amrex::Real force_restore_snow_t_melt = amrex::Real(273.16);

/**
 * @brief Parameters of one snow step (from erf.radiation.seb_*).
 */
struct ForceRestoreSnowParams
{
    amrex::Real w_crn  = amrex::Real(10.0);   ///< Critical snow water equivalent W_crn [kg/m^2] (seb_snow_critical_swe_kg_m2)
    amrex::Real albedo = amrex::Real(0.70);   ///< Snow albedo alpha_n [-] (seb_snow_albedo)
    amrex::Real c_s    = amrex::Real(2.0e4);  ///< Surface heat capacity C_s = 1 / C_t [J/m^2/K] (seb_surface_heat_capacity)
};

/**
 * @brief What one snow step moved in one column [kg/m^2].
 */
struct ForceRestoreSnowFluxes
{
    amrex::Real snowfall    = amrex::Real(0.0);  ///< Frozen precipitation taken in
    amrex::Real sublimation = amrex::Real(0.0);  ///< Snow sublimated
    amrex::Real melt        = amrex::Real(0.0);  ///< Snow melted
};

/**
 * @brief Snow cover fraction p_ng of the ground, DRM95 Eq. (7).
 *
 * @param[in] w_n   snow water equivalent [kg/m^2]
 * @param[in] w_crn critical snow water equivalent [kg/m^2], > 0
 * @return W_n / (W_n + W_crn) in [0, 1); zero for no snow or a non-finite W_n.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real force_restore_snow_fraction (amrex::Real w_n, amrex::Real w_crn) noexcept
{
    if (!(w_n > amrex::Real(0.0)) || !amrex::Math::isfinite(w_n)) { return amrex::Real(0.0); }
    return w_n / (w_n + w_crn);
}

/**
 * @brief Surface albedo with snow, DRM95 Eq. (8) with veg = 0.
 *
 * @param[in] alb_ground snow-free albedo alpha_g [-]
 * @param[in] w_n        snow water equivalent [kg/m^2]
 * @param[in] prm        parameters
 * @return (1 - p_ng) alpha_g + p_ng alpha_n
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real force_restore_snow_albedo (amrex::Real alb_ground, amrex::Real w_n,
                                       const ForceRestoreSnowParams& prm) noexcept
{
    const amrex::Real p_ng = force_restore_snow_fraction(w_n, prm.w_crn);
    return (amrex::Real(1.0) - p_ng) * alb_ground + p_ng * prm.albedo;
}

/**
 * @brief One step of the snow layer of one column, after the balance has updated T_s.
 *
 * @param[in,out] w_n snow water equivalent [kg/m^2]
 * @param[in,out] t_s surface temperature [K]: the balance's estimate in, the value after
 *                    the melt out
 * @param[in]     p_n frozen precipitation [kg/m^2/s]
 * @param[in]     e   evaporation LE / L_v [kg/m^2/s], positive upward
 * @param[in]     dt  the step [s]
 * @param[in]     prm parameters
 * @return the snow taken in, sublimated and melted
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ForceRestoreSnowFluxes
force_restore_snow_step (amrex::Real& w_n, amrex::Real& t_s, amrex::Real p_n, amrex::Real e,
                         amrex::Real dt, const ForceRestoreSnowParams& prm) noexcept
{
    ForceRestoreSnowFluxes out;
    if (!(p_n > zero) || !amrex::Math::isfinite(p_n)) { p_n = zero; }
    if (!(e > zero) || !amrex::Math::isfinite(e)) { e = zero; }

    // The cover fraction of the snow at the start of the step: DRM95 Eq. (7).
    const amrex::Real p_ng = force_restore_snow_fraction(w_n, prm.w_crn);

    // Each term below reports what W_n actually changed by, not the amount asked for: in
    // single precision an increment below the resolution of a deep pack leaves W_n as it
    // was, and the soil must then not be handed water the snow did not lose.

    // Snowfall: DRM95 Eq. (6), first term.
    amrex::Real w_new = w_n + p_n * dt;
    out.snowfall = w_new - w_n;
    w_n = w_new;

    // Sublimation from the covered fraction: DRM95 Eq. (6), second term. A negative W_n (a
    // corrupt state, which the caller refuses after the step) gives none.
    // amrex::max takes references: a namespace-scope constant passed to it is odr-used, which
    // CUDA refuses in device code, so the constants are passed as temporaries here and below.
    w_new = w_n - amrex::min(p_ng * e * dt, amrex::max(w_n, amrex::Real(0.0)));
    out.sublimation = w_n - w_new;
    w_n = w_new;

    // Melt: DRM95 Eq. (12), with C_t = 1 / C_s.
    if (w_n > zero && t_s > force_restore_snow_t_melt) {
        const amrex::Real potential = prm.c_s * (t_s - force_restore_snow_t_melt) / lat_ice;
        if (potential < w_n) {
            w_new = w_n - potential;
            out.melt = w_n - w_new;
            w_n = w_new;
            // The heat of the part W_n did not take stays in T_s for the next step; it is
            // zero to round-off in double precision.
            t_s = amrex::max(amrex::Real(force_restore_snow_t_melt),
                             force_restore_snow_t_melt + lat_ice * (potential - out.melt) / prm.c_s);
        } else {
            out.melt = w_n;
            w_n = zero;
            t_s -= lat_ice * out.melt / prm.c_s;
        }
    }
    return out;
}

#endif // ERF_FORCE_RESTORE_SNOW_H_
