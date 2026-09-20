#ifndef ERF_SURFACE_TEMPERATURE_H_
#define ERF_SURFACE_TEMPERATURE_H_

#include <AMReX_Math.H>

#include <ERF_Constants.H>
#include <ERF_EOS.H>

namespace erf_surface_temperature {

/**
 * Pressure at a boundary face from the adjacent atmospheric cell state.
 *
 * The state pressure is diagnosed at the cell centre. Adding rho*g*delta_z
 * gives the pressure at the boundary face, where delta_z is the signed
 * distance from the face to the cell centre.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
pressure_at_boundary_from_cell (const amrex::Real rho,
                                const amrex::Real rho_theta,
                                const amrex::Real qv,
                                const amrex::Real delta_z)
{
    return getPgivenRTh(rho_theta, qv) + rho * CONST_GRAV * delta_z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
pressure_at_surface (const amrex::Real rho,
                     const amrex::Real rho_theta,
                     const amrex::Real qv,
                     const amrex::Real delta_z)
{
    return pressure_at_boundary_from_cell(rho, rho_theta, qv, delta_z);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_pressure (const amrex::Real pressure)
{
    return amrex::Math::isfinite(pressure) && pressure > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_temperature (const amrex::Real temperature)
{
    return amrex::Math::isfinite(temperature) && temperature > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
theta_to_temperature (const amrex::Real theta,
                      const amrex::Real pressure,
                      const amrex::Real rdOcp,
                      amrex::Real& temperature)
{
    if (!valid_pressure(pressure) || !valid_temperature(theta) ||
        !amrex::Math::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
        return false;
    }
    const amrex::Real result = getTgivenPandTh(pressure, theta, rdOcp);
    if (!valid_temperature(result)) {
        return false;
    }
    temperature = result;
    return true;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
temperature_to_theta (const amrex::Real temperature,
                      const amrex::Real pressure,
                      const amrex::Real rdOcp,
                      amrex::Real& theta)
{
    if (!valid_pressure(pressure) || !valid_temperature(temperature) ||
        !amrex::Math::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
        return false;
    }
    const amrex::Real result = getThgivenTandP(temperature, pressure, rdOcp);
    if (!valid_temperature(result)) {
        return false;
    }
    theta = result;
    return true;
}

} // namespace erf_surface_temperature

#endif
