#ifndef ERF_RRTMGP_UTILS_H
#define ERF_RRTMGP_UTILS_H

#include "rrtmgp_const.h"

#include "ERF_Constants.H"

namespace rrtmgp {

template<class View1, class View2, class View3, class View4, class View5>
void
mixing_ratio_to_cloud_mass (View1 const& mixing_ratio  ,
                            View2 const& cloud_fraction,
                            View3 const& rho           ,
                            View4 const& dz            ,
                            View5 const& cloud_mass)
{
    const int ncol = mixing_ratio.extent(0);
    const int nlay = mixing_ratio.extent(1);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                         KOKKOS_LAMBDA (int icol, int ilay)
    {
        // Compute in-cloud mixing ratio (mixing ratio of the cloudy part of the layer)
        // NOTE: these thresholds (from E3SM) seem arbitrary, but included here for consistency
        // This limits in-cloud mixing ratio to 0.005 kg/kg. According to note in cloud_diagnostics
        // in EAM, this is consistent with limits in MG2. Is this true for P3?
        if (cloud_fraction(icol,ilay) > 0) {
            // Compute layer-integrated cloud mass (per unit area)
            auto incloud_mixing_ratio = std::min(mixing_ratio(icol,ilay) /
                                                 std::max(0.0001, cloud_fraction(icol,ilay)), 0.005);
            cloud_mass(icol,ilay) = incloud_mixing_ratio * rho(icol,ilay) * dz(icol,ilay);
        } else {
            cloud_mass(icol,ilay) = 0;
        }
    });
}

// 1D limiting
template<typename InT, typename OutT, typename T>
void
limit_to_bounds_1d (InT const& arr_in,
                    T const lower    ,
                    T const upper    ,
                    OutT& arr_out)
{
    Kokkos::parallel_for(arr_out.size(), KOKKOS_LAMBDA(int i)
    {
        arr_out(i) = std::min(std::max(arr_in(i), lower), upper);
    });
}

// 2D limiting
template<typename InT, typename OutT, typename T>
void
limit_to_bounds_2d (InT const& arr_in,
                    T const lower    ,
                    T const upper    ,
                    OutT& arr_out)
{
    const int ex0 = (int) arr_out.extent(0);
    const int ex1 = (int) arr_out.extent(1);
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ex0, ex1}),
                         KOKKOS_LAMBDA (int i, int j)
    {
        arr_out(i, j) = std::min(std::max(arr_in(i, j), lower), upper);
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
template<class View1, class View2, class View3, class View4, class View5>
void
compute_heating_rate (View1 const& flux_up,
                      View2 const& flux_dn,
                      View3 const& rho    ,
                      View4 const& dz     ,
                      View5& heating_rate)
{
    const int ncol = (int)flux_up.extent(0);
    const int nlay = (int)flux_up.extent(1)-1;
    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {ncol, nlay}),
                         KOKKOS_LAMBDA (int icol, int ilay)
    {
        // NOTE: This calculation defines Fnet = Up - Dn
        //       H = dT/dt = g/Cp dF/dP
        //       Using dF/dz * dz/dP = dF/dz * (1/(-rho*g)) = dF/dP
        //       We have: dT/dt = dF/dz * (1/(-rho*Cp))
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


}// namespace rrtmgp
#endif
