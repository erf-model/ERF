#ifndef ERF_MORRISON_VAPOR_PRESSURE_H_
#define ERF_MORRISON_VAPOR_PRESSURE_H_

#include <cmath>

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

#include "ERF_NumericalConstants.H"

/**
 * Helper function to calculate saturation vapor pressure for water or ice.
 * This corresponds to the POLYSVP function in the Fortran code (line ~5580).
 *
 * @param[in] T Temperature in Kelvin
 * @param[in] type 0 for liquid water, 1 for ice
 * @return Saturation vapor pressure in Pascals
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
calc_saturation_vapor_pressure (const amrex::Real T, const int type)
{
    amrex::Real polysvp = amrex::Real(0);
    amrex::Real del_T = T - amrex::Real(273.15);  // Convert to Celsius

    if (type == 1) {  // Ice (lines ~5631-5644)
        if (T >= amrex::Real(195.8)) {
            // Flatau et al. formula for ice
            const amrex::Real a0i = amrex::Real(6.11147274);
            const amrex::Real a1i = amrex::Real(0.503160820);
            const amrex::Real a2i = amrex::Real(0.188439774e-1);
            const amrex::Real a3i = amrex::Real(0.420895665e-3);
            const amrex::Real a4i = amrex::Real(0.615021634e-5);
            const amrex::Real a5i = amrex::Real(0.602588177e-7);
            const amrex::Real a6i = amrex::Real(0.385852041e-9);
            const amrex::Real a7i = amrex::Real(0.146898966e-11);
            const amrex::Real a8i = amrex::Real(0.252751365e-14);

            polysvp = a0i + del_T*(a1i + del_T*(a2i + del_T*(a3i + del_T*(a4i + del_T*(a5i + del_T*(a6i + del_T*(a7i + a8i*del_T)))))));
            polysvp *= amrex::Real(100.0);  // Convert from hPa to Pa
        } else {
            // Goff-Gratch formula for ice at cold temperatures
            polysvp = std::pow(amrex::Real(10.0), (-amrex::Real(9.09718)*(amrex::Real(273.16)/T-one) - amrex::Real(3.56654)*std::log10(amrex::Real(273.16)/T) +
                                             amrex::Real(0.876793)*(one-T/amrex::Real(273.16)) + std::log10(amrex::Real(6.1071)))) * amrex::Real(100.0);
        } // T
    } else {  // Water (lines ~5648-5665)
      if (T >= amrex::Real(202.0)) {
        // Flatau et al. formula for liquid water
        const amrex::Real a0 = amrex::Real(6.11239921);
        const amrex::Real a1 = amrex::Real(0.443987641);
        const amrex::Real a2 = amrex::Real(0.142986287e-1);
        const amrex::Real a3 = amrex::Real(0.264847430e-3);
        const amrex::Real a4 = amrex::Real(0.302950461e-5);
        const amrex::Real a5 = amrex::Real(0.206739458e-7);
        const amrex::Real a6 = amrex::Real(0.640689451e-10);
        const amrex::Real a7 = -amrex::Real(0.952447341e-13);
        const amrex::Real a8 = -amrex::Real(0.976195544e-15);

        polysvp = a0 + del_T*(a1 + del_T*(a2 + del_T*(a3 + del_T*(a4 + del_T*(a5 + del_T*(a6 + del_T*(a7 + a8*del_T)))))));
        polysvp *= amrex::Real(100.0);  // Convert from hPa to Pa
      } else {
        // Goff-Gratch formula for water at cold temperatures
        polysvp = std::pow(amrex::Real(10.0), (-amrex::Real(7.90298)*(amrex::Real(373.16)/T-one) + amrex::Real(5.02808)*std::log10(amrex::Real(373.16)/T) -
                                  amrex::Real(1.3816e-7)*(std::pow(amrex::Real(10.0), (amrex::Real(11.344)*(one-T/amrex::Real(373.16))))-one) +
                                  amrex::Real(8.1328e-3)*(std::pow(amrex::Real(10.0), (-amrex::Real(3.49149)*(amrex::Real(373.16)/T-one)))-one) +
                                  std::log10(amrex::Real(1013.246)))) * amrex::Real(100.0);
      }
    }

    return polysvp;
}

#endif
