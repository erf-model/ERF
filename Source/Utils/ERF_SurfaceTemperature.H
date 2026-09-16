#ifndef ERF_SURFACE_TEMPERATURE_H_
#define ERF_SURFACE_TEMPERATURE_H_

#include <ERF_Constants.H>
#include <ERF_EOS.H>

namespace erf_surface_temperature {

/**
 * Pressure at the physical surface from the lowest atmospheric cell state.
 *
 * The state pressure is diagnosed at the cell centre. Adding rho*g*delta_z
 * gives the pressure at the physical surface, where delta_z is the
 * terrain-relative distance from the surface to that cell centre.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
pressure_at_surface (const amrex::Real rho,
                     const amrex::Real rho_theta,
                     const amrex::Real qv,
                     const amrex::Real delta_z)
{
    return getPgivenRTh(rho_theta, qv) + rho * CONST_GRAV * delta_z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_pressure (const amrex::Real pressure)
{
    return std::isfinite(pressure) && pressure > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
valid_temperature (const amrex::Real temperature)
{
    return std::isfinite(temperature) && temperature > amrex::Real(0.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool
theta_to_temperature (const amrex::Real theta,
                      const amrex::Real pressure,
                      const amrex::Real rdOcp,
                      amrex::Real& temperature)
{
    if (!valid_pressure(pressure) || !valid_temperature(theta) ||
        !std::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
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
        !std::isfinite(rdOcp) || rdOcp <= amrex::Real(0.0)) {
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
