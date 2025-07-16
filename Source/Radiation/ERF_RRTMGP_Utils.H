#ifndef ERF_RRTMGP_UTILS_H
#define ERF_RRTMGP_UTILS_H

#include "rrtmgp_const.h"
//#include "rrtmgp_conversion.h"

#include "YAKL.h"
#include "YAKL_Bounds_fortran.h"

#include "ERF_Constants.H"

namespace rrtmgp {

// Things we need from YAKL
using yakl::intrinsics::maxval;
using yakl::intrinsics::minval;
using yakl::intrinsics::count;
using yakl::intrinsics::sum;


template<class T, int myMem, int myStyle>
void
mixing_ratio_to_cloud_mass (yakl::Array<T,2,myMem,myStyle> const &mixing_ratio  ,
                            yakl::Array<T,2,myMem,myStyle> const &cloud_fraction,
                            yakl::Array<T,2,myMem,myStyle> const &rho           ,
                            yakl::Array<T,2,myMem,myStyle> const &dz            ,
                            yakl::Array<T,2,myMem,myStyle>       &cloud_mass)
{
    int ncol = mixing_ratio.dimension[0];
    int nlay = mixing_ratio.dimension[1];
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay, ncol), YAKL_LAMBDA(int ilay, int icol)
    {
        // Compute in-cloud mixing ratio (mixing ratio of the cloudy part of the layer)
        // NOTE: these thresholds (from E3SM) seem arbitrary, but included here for consistency
        // This limits in-cloud mixing ratio to 0.005 kg/kg. According to note in cloud_diagnostics
        // in EAM, this is consistent with limits in MG2. Is this true for P3?
        if (cloud_fraction(icol,ilay) > 0) {
            // Compute layer-integrated cloud mass (per unit area)
            auto incloud_mixing_ratio = std::min(mixing_ratio(icol,ilay) / std::max(0.0001, cloud_fraction(icol,ilay)), 0.005);
            cloud_mass(icol,ilay) = incloud_mixing_ratio * rho(icol,ilay) * dz(icol,ilay);
        } else {
            cloud_mass(icol,ilay) = 0;
        }
    });
}


template<class S, class T>
void
limit_to_bounds (S const &arr_in,
                 T const lower  ,
                 T const upper  ,
                 S &arr_out)
{
    yakl::c::parallel_for(arr_in.totElems(), YAKL_LAMBDA(int i)
    {
        arr_out.data()[i] = std::min(std::max(arr_in.data()[i], lower), upper);
    });
}


// Provide a routine to compute heating due to radiative fluxes. This is
// computed as net flux into a layer, converted to a heating rate. It is
// the responsibility of the user to ensure fields are passed with the
// proper units. I.e., pressure at level interfaces should be in Pa,
// fluxes in W m-2, Cpair in J kg-1 K-1, gravit in m s-2. This will give
// heating in units of K s-1.
// TODO: we should probably update this to use the pseudo-density pdel instead
// of approximating pdel by differencing the level interface pressures.
// We are leaving this for the time being for consistency with SCREAMv0,
// from which this code was directly ported.
template <class T, int myMem, int myStyle>
void
compute_heating_rate (yakl::Array<T,2,myMem,myStyle> const &flux_up,
                      yakl::Array<T,2,myMem,myStyle> const &flux_dn,
                      yakl::Array<T,2,myMem,myStyle> const &rho    ,
                      yakl::Array<T,2,myMem,myStyle> const &dz     ,
                      yakl::Array<T,2,myMem,myStyle> &heating_rate)
{
    auto ncol = flux_up.dimension[0];
    auto nlay = flux_up.dimension[1]-1;
    yakl::fortran::parallel_for(yakl::fortran::SimpleBounds<2>(nlay,ncol), YAKL_LAMBDA(int ilay, int icol)
    {
        heating_rate(icol,ilay) = ( flux_up(icol,ilay+1) - flux_up(icol,ilay)
                                  - flux_dn(icol,ilay+1) + flux_dn(icol,ilay) )
                                  / ( -rho(icol,ilay) * dz(icol,ilay) * Cp_d );
    });
}


inline
bool
radiation_do (const int irad,
              const int nstep)
{
    // If irad == 0, then never do radiation;
    // Otherwise, we always call radiation at the first step,
    // and afterwards we do radiation if the timestep is divisible
    // by irad
    if (irad == 0) {
        return false;
    } else {
        return ( (nstep == 0) || (nstep % irad == 0) );
    }
}

} // namespace rrtmgp
#endif
