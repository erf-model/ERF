#ifndef ERF_SURFACE_MOISTURE_H_
#define ERF_SURFACE_MOISTURE_H_

#include <AMReX_Math.H>
#include <AMReX_Algorithm.H>

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
 * height [s/m], (ln(z_ref / z0) - psi_h) / (kappa u*): the form of the surface layer's
 * moisture flux, which uses z0 for heat and moisture too. u* is floored at 1e-3 m/s and
 * the resistance at 1 s/m.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
aerodynamic_resistance (const amrex::Real z_ref,
                        const amrex::Real z0,
                        const amrex::Real psi_h,
                        const amrex::Real kappa,
                        const amrex::Real u_star)
{
    const amrex::Real r = (std::log(z_ref / z0) - psi_h) / (kappa * amrex::max(u_star, amrex::Real(1.0e-3)));
    return amrex::max(amrex::Real(1.0), r);
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
