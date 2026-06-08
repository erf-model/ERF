#include <string>
#include <vector>
#include <memory>
#include <complex>
#include <cmath>

#include <AMReX_Math.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_TableData.H>
#include <AMReX_MultiFabUtil.H>

#include "ERF.H"
#include "ERF_Constants.H"
#include "ERF_MicrophysicsUtils.H"
#include "ERF_IndexDefines.H"
#include "ERF_DataStruct.H"
#include "ERF_NullMoist.H"
#include "ERF_Morrison.H"
#ifdef ERF_USE_MORR_FORT
#include <ERF_Morrison_Fortran_Interface.H>
#endif

using namespace amrex;

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
    amrex::Real polysvp = Real(0);
    amrex::Real del_T = T - Real(273.15);  // Convert to Celsius

    if (type == 1) {  // Ice (lines ~5631-5644)
        if (T >= Real(195.8)) {
            // Flatau et al. formula for ice
            const amrex::Real a0i = Real(6.11147274);
            const amrex::Real a1i = Real(0.503160820);
            const amrex::Real a2i = Real(0.188439774e-1);
            const amrex::Real a3i = Real(0.420895665e-3);
            const amrex::Real a4i = Real(0.615021634e-5);
            const amrex::Real a5i = Real(0.602588177e-7);
            const amrex::Real a6i = Real(0.385852041e-9);
            const amrex::Real a7i = Real(0.146898966e-11);
            const amrex::Real a8i = Real(0.252751365e-14);

            polysvp = a0i + del_T*(a1i + del_T*(a2i + del_T*(a3i + del_T*(a4i + del_T*(a5i + del_T*(a6i + del_T*(a7i + a8i*del_T)))))));
            polysvp *= Real(100.0);  // Convert from hPa to Pa
        } else {
            // Goff-Gratch formula for ice at cold temperatures
            polysvp = std::pow(Real(10.0), (-Real(9.09718)*(Real(273.16)/T-one) - Real(3.56654)*std::log10(Real(273.16)/T) +
                                             Real(0.876793)*(one-T/Real(273.16)) + std::log10(Real(6.1071)))) * Real(100.0);
        } // T
    } else {  // Water (lines ~5648-5665)
      if (T >= Real(202.0)) {
        // Flatau et al. formula for liquid water
        const amrex::Real a0 = Real(6.11239921);
        const amrex::Real a1 = Real(0.443987641);
        const amrex::Real a2 = Real(0.142986287e-1);
        const amrex::Real a3 = Real(0.264847430e-3);
        const amrex::Real a4 = Real(0.302950461e-5);
        const amrex::Real a5 = Real(0.206739458e-7);
        const amrex::Real a6 = Real(0.640689451e-10);
        const amrex::Real a7 = -Real(0.952447341e-13);
        const amrex::Real a8 = -Real(0.976195544e-15);

        polysvp = a0 + del_T*(a1 + del_T*(a2 + del_T*(a3 + del_T*(a4 + del_T*(a5 + del_T*(a6 + del_T*(a7 + a8*del_T)))))));
        polysvp *= Real(100.0);  // Convert from hPa to Pa
      } else {
        // Goff-Gratch formula for water at cold temperatures
        polysvp = std::pow(Real(10.0), (-Real(7.90298)*(Real(373.16)/T-one) + Real(5.02808)*std::log10(Real(373.16)/T) -
                                  Real(1.3816e-7)*(std::pow(Real(10.0), (Real(11.344)*(one-T/Real(373.16))))-one) +
                                  Real(8.1328e-3)*(std::pow(Real(10.0), (-Real(3.49149)*(Real(373.16)/T-one)))-one) +
                                  std::log10(Real(1013.246)))) * Real(100.0);
      }
    }

    return polysvp;
}
