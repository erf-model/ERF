#ifndef ERF_SURFACE_MOISTURE_H_
#define ERF_SURFACE_MOISTURE_H_

#include <AMReX_Math.H>
#include <AMReX_Algorithm.H>

#include <ERF_MOSTUtils.H>

#include <cmath>

/**
 * Pieces of the surface layer's land surface mixing ratio when the two-stream balance
 * supplies the surface's moisture (erf.radiation.seb_surface_layer_uses_moisture): the
 * surface evaporates at beta times the potential rate, and its mixing ratio is
 * beta * q_sat + (1 - beta) * q_air. With a vegetation type, beta comes from the
 * resistances of a canopy and of bare soil in series with the aerodynamic one.
 */
namespace erf_surface_moisture {

/**
 * Components of the per-column field the balance hands the surface layer
 * (TwoStreamRadiation::seb_surface_moisture): the soil-water factor beta in [0, 1], and,
 * with a vegetation type, the vegetated fraction, the canopy resistance without its
 * vapour-deficit factor [s/m] and the bare-soil resistance [s/m].
 */
enum Component : int {
    SoilFactor = 0,
    VegetationFraction,
    CanopyResistanceWithoutVPD,
    SoilResistance,
    NumComponents
};

/**
 * Surface mixing ratio of a surface with moisture availability beta: the saturation
 * mixing ratio at the surface temperature where beta = 1, the air's own where beta = 0
 * (no moisture flux), linear in between. A surface-layer moisture flux proportional to
 * (q_surface - q_air) is then beta times the potential one.
 *
 * @param[in] beta  moisture availability, clamped to [0, 1]
 * @param[in] q_sat saturation mixing ratio at the surface temperature and pressure [kg/kg]
 * @param[in] q_air mixing ratio of the air at the surface layer's reference height [kg/kg]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
surface_mixing_ratio (const amrex::Real beta,
                      const amrex::Real q_sat,
                      const amrex::Real q_air)
{
    const amrex::Real b = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), beta));
    return b * q_sat + (amrex::Real(1.0) - b) * q_air;
}

/**
 * The vapour-deficit factor of the Jarvis canopy resistance,
 * 1 / (1 + HS * max(0, q_sat - q_air)), floored at 0.01 as in Noah-MP
 * (ResistanceCanopyStomataJarvisMod).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
vapour_deficit_factor (const amrex::Real hs,
                       const amrex::Real q_sat,
                       const amrex::Real q_air)
{
    const amrex::Real deficit = amrex::max(amrex::Real(0.0), q_sat - q_air);
    return amrex::max(amrex::Real(0.01), amrex::Real(1.0) / (amrex::Real(1.0) + hs * deficit));
}

/**
 * Aerodynamic resistance to the moisture flux between the surface and the reference
 * height [s/m]: max(ln(z_ref / z0) - psi_h, 1) / (kappa u*).
 *
 * This is the limited CmPsih of the surface layer's surface_temp kernel
 * (ERF_MOSTStress.H), which uses z0 for heat and moisture too. Its q* is
 * kappa (q_air - q_surf) / CmPsih, so the flux is u* q* = (q_air - q_surf) / r_a.
 * u* is floored at 1e-3 m/s.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
aerodynamic_resistance (const amrex::Real z_ref,
                        const amrex::Real z0,
                        const amrex::Real psi_h,
                        const amrex::Real kappa,
                        const amrex::Real u_star)
{
    const amrex::Real cm_psi_h = amrex::max(std::log(z_ref / z0) - psi_h, amrex::Real(1.0));
    return cm_psi_h / (kappa * amrex::max(u_star, amrex::Real(1.0e-3)));
}

/**
 * The aerodynamic resistance [s/m] the surface layer's moisture flux sees over a land
 * column whose skin temperature the two-stream balance gives it.
 *
 * Skin coupling needs erf.most.surf_temp, so that flux comes from the surface_temp
 * kernel. Its stability function is Jimenez's psi_h2 at z_ref / L (calc_psi_h2), not
 * Businger-Dyer's psi_h, and it is limited as in aerodynamic_resistance.
 *
 * The resistance uses the last u* and Obukhov length L. Before the first flux, u* is
 * @p unset or non-finite, and a neutral first guess is used:
 * u* = kappa max(U, U_min) / max(ln(z_ref / z0), 1), with psi_h2 = 0.
 * An unset, non-finite or zero L also gives psi_h2 = 0.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
surface_layer_aerodynamic_resistance (const similarity_funs& sfuns,
                                      const amrex::Real z_ref,
                                      const amrex::Real z0,
                                      const amrex::Real u_star,
                                      const amrex::Real olen,
                                      const amrex::Real wind,
                                      const amrex::Real wind_min,
                                      const amrex::Real kappa,
                                      const amrex::Real unset)
{
    amrex::Real ustar = u_star;
    amrex::Real psi_h = amrex::Real(0.0);
    if (ustar == unset || !amrex::Math::isfinite(ustar)) {
        ustar = kappa * amrex::max(wind, wind_min) /
                amrex::max(std::log(z_ref / z0), amrex::Real(1.0));
    } else if (olen != unset && amrex::Math::isfinite(olen) && olen != amrex::Real(0.0)) {
        psi_h = sfuns.calc_psi_h2(z_ref / olen);
    }
    return aerodynamic_resistance(z_ref, z0, psi_h, kappa, ustar);
}

/**
 * Moisture availability of a surface whose vegetated fraction f_veg transpires through a
 * canopy resistance r_c and whose bare fraction evaporates through a soil resistance
 * r_soil, each in series with the aerodynamic resistance r_a:
 * f_veg r_a / (r_a + r_c) + (1 - f_veg) r_a / (r_a + r_soil).
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
two_source_availability (const amrex::Real r_a,
                         const amrex::Real r_c,
                         const amrex::Real r_soil,
                         const amrex::Real f_veg)
{
    const amrex::Real f = amrex::max(amrex::Real(0.0), amrex::min(amrex::Real(1.0), f_veg));
    return f * r_a / (r_a + r_c) + (amrex::Real(1.0) - f) * r_a / (r_a + r_soil);
}

} // namespace erf_surface_moisture

#endif
