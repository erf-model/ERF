#ifndef INTERPOLATE_WENO_H_
#define INTERPOLATE_WENO_H_

#include "DataStruct.H"

/**
 * Interpolation operators used in construction of advective fluxes using WENO schemes
 */

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInX_WENO (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                     int qty_index, int spatial_order)
{
    amrex::Real eps=1.0e-6;
    amrex::Real beta1=0., beta2=0., beta3=0.;
    amrex::Real w1=0.   , w2=0.   , w3=0.   , sum_wl=0.;
    amrex::Real phi1=0. , phi2=0. , phi3=0. ;

    // WENO3
    if (spatial_order == 3) {
        beta1  = ( qty(i-1,j,k,qty_index) - qty(i-2,j,k,qty_index) ) *
                 ( qty(i-1,j,k,qty_index) - qty(i-2,j,k,qty_index) );
        beta2  = ( qty(i  ,j,k,qty_index) - qty(i-1,j,k,qty_index) ) *
                 ( qty(i  ,j,k,qty_index) - qty(i-1,j,k,qty_index) );

        w1     = ( (1./3.) / ((eps+beta1)*(eps+beta1)) );
        w2     = ( (2./3.) / ((eps+beta2)*(eps+beta2)) );

        phi1   = 0.5 * ( -qty(i-2,j,k,qty_index) + 3.*qty(i-1,j,k,qty_index) );
        phi2   = 0.5 * (  qty(i-1,j,k,qty_index) +    qty(i  ,j,k,qty_index) );
    // WENO5
    } else {
        beta1  = (13./12.) * (    qty(i-3,j,k,qty_index) - 2.*qty(i-2,j,k,qty_index) +    qty(i-1,j,k,qty_index) ) *
                             (    qty(i-3,j,k,qty_index) - 2.*qty(i-2,j,k,qty_index) +    qty(i-1,j,k,qty_index) );
        beta1 += (1./4.)   * (    qty(i-3,j,k,qty_index) - 4.*qty(i-2,j,k,qty_index) + 3.*qty(i-1,j,k,qty_index) ) *
                             (    qty(i-3,j,k,qty_index) - 4.*qty(i-2,j,k,qty_index) + 3.*qty(i-1,j,k,qty_index) );

        beta2  = (13./12.) * (    qty(i-2,j,k,qty_index) - 2.*qty(i-1,j,k,qty_index) +    qty(i  ,j,k,qty_index) ) *
                             (    qty(i-2,j,k,qty_index) - 2.*qty(i-1,j,k,qty_index) +    qty(i  ,j,k,qty_index) );
        beta2 += (1./4.)   * (    qty(i-2,j,k,qty_index) -    qty(i  ,j,k,qty_index) ) *
                             (    qty(i-2,j,k,qty_index) -    qty(i  ,j,k,qty_index) );

        beta3  = (13./12.) * (    qty(i-1,j,k,qty_index) - 2.*qty(i  ,j,k,qty_index) +    qty(i+1,j,k,qty_index) ) *
                             (    qty(i-1,j,k,qty_index) - 2.*qty(i  ,j,k,qty_index) +    qty(i+1,j,k,qty_index) );
        beta3 += (1./4.)   * ( 3.*qty(i-1,j,k,qty_index) - 4.*qty(i  ,j,k,qty_index) +    qty(i+1,j,k,qty_index) ) *
                             ( 3.*qty(i-1,j,k,qty_index) - 4.*qty(i  ,j,k,qty_index) +    qty(i+1,j,k,qty_index) );

        w1   = ( (1./10.) / ((eps+beta1)*(eps+beta1)) );
        w2   = ( (3./5. ) / ((eps+beta2)*(eps+beta2)) );
        w3   = ( (3./10.) / ((eps+beta3)*(eps+beta3)) );

        phi1 =  (1./3.)*qty(i-3,j,k,qty_index) - (7./6.)*qty(i-2,j,k,qty_index) + (11./6.)*qty(i-1,j,k,qty_index);
        phi2 = -(1./6.)*qty(i-2,j,k,qty_index) + (5./6.)*qty(i-1,j,k,qty_index) + ( 1./3.)*qty(i  ,j,k,qty_index);
        phi3 =  (1./3.)*qty(i-1,j,k,qty_index) + (5./6.)*qty(i  ,j,k,qty_index) - ( 1./6.)*qty(i+1,j,k,qty_index);
    }

    sum_wl = w1 + w2 + w3;
    w1    /= sum_wl;
    w2    /= sum_wl;
    w3    /= sum_wl;
    return w1*phi1 + w2*phi2 + w3*phi3;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInY_WENO (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                     int qty_index, int spatial_order)
{
    amrex::Real eps=1.0e-6;
    amrex::Real beta1=0., beta2=0., beta3=0.;
    amrex::Real w1=0.   , w2=0.   , w3=0.   , sum_wl=0.;
    amrex::Real phi1=0. , phi2=0. , phi3=0. ;

    // WENO3
    if (spatial_order == 3) {
        beta1  = ( qty(i,j-1,k,qty_index) - qty(i,j-2,k,qty_index) ) *
                 ( qty(i,j-1,k,qty_index) - qty(i,j-2,k,qty_index) );
        beta2  = ( qty(i,j  ,k,qty_index) - qty(i,j-1,k,qty_index) ) *
                 ( qty(i,j  ,k,qty_index) - qty(i,j-1,k,qty_index) );

        w1     = ( (1./3.) / ((eps+beta1)*(eps+beta1)) );
        w2     = ( (2./3.) / ((eps+beta2)*(eps+beta2)) );

        phi1   = 0.5 * ( -qty(i,j-2,k,qty_index) + 3.*qty(i,j-1,k,qty_index) );
        phi2   = 0.5 * (  qty(i,j-1,k,qty_index) +    qty(i,j  ,k,qty_index) );
    // WENO5
    } else {
        beta1  = (13./12.) * (    qty(i,j-3,k,qty_index) - 2.*qty(i,j-2,k,qty_index) +    qty(i,j-1,k,qty_index) ) *
                             (    qty(i,j-3,k,qty_index) - 2.*qty(i,j-2,k,qty_index) +    qty(i,j-1,k,qty_index) );
        beta1 += (1./4.)   * (    qty(i,j-3,k,qty_index) - 4.*qty(i,j-2,k,qty_index) + 3.*qty(i,j-1,k,qty_index) ) *
                             (    qty(i,j-3,k,qty_index) - 4.*qty(i,j-2,k,qty_index) + 3.*qty(i,j-1,k,qty_index) );

        beta2  = (13./12.) * (    qty(i,j-2,k,qty_index) - 2.*qty(i,j-1,k,qty_index) +    qty(i,j  ,k,qty_index) ) *
                             (    qty(i,j-2,k,qty_index) - 2.*qty(i,j-1,k,qty_index) +    qty(i,j  ,k,qty_index) );
        beta2 += (1./4.)   * (    qty(i,j-2,k,qty_index) -    qty(i,j  ,k,qty_index) ) *
                             (    qty(i,j-2,k,qty_index) -    qty(i,j  ,k,qty_index) );

        beta3  = (13./12.) * (    qty(i,j-1,k,qty_index) - 2.*qty(i,j  ,k,qty_index) +    qty(i,j+1,k,qty_index) ) *
                             (    qty(i,j-1,k,qty_index) - 2.*qty(i,j  ,k,qty_index) +    qty(i,j+1,k,qty_index) );
        beta3 += (1./4.)   * ( 3.*qty(i,j-1,k,qty_index) - 4.*qty(i,j  ,k,qty_index) +    qty(i,j+1,k,qty_index) ) *
                             ( 3.*qty(i,j-1,k,qty_index) - 4.*qty(i,j  ,k,qty_index) +    qty(i,j+1,k,qty_index) );

        w1     = ( (1./10.) / ((eps+beta1)*(eps+beta1)) );
        w2     = ( (3./5. ) / ((eps+beta2)*(eps+beta2)) );
        w3     = ( (3./10.) / ((eps+beta3)*(eps+beta3)) );

        phi1 =  (1./3.)*qty(i,j-3,k,qty_index) - (7./6.)*qty(i,j-2,k,qty_index) + (11./6.)*qty(i,j-1,k,qty_index);
        phi2 = -(1./6.)*qty(i,j-2,k,qty_index) + (5./6.)*qty(i,j-1,k,qty_index) + ( 1./3.)*qty(i,j  ,k,qty_index);
        phi3 =  (1./3.)*qty(i,j-1,k,qty_index) + (5./6.)*qty(i,j  ,k,qty_index) - ( 1./6.)*qty(i,j+1,k,qty_index);
    }

    sum_wl = w1 + w2 + w3;
    w1    /= sum_wl;
    w2    /= sum_wl;
    w3    /= sum_wl;
    return w1*phi1 + w2*phi2 + w3*phi3;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
InterpolateInZ_WENO (int i, int j, int k, const amrex::Array4<const amrex::Real>& qty,
                     int qty_index, int spatial_order)
{
    amrex::Real eps=1.0e-6;
    amrex::Real beta1=0., beta2=0., beta3=0.;
    amrex::Real w1=0.   , w2=0.   , w3=0.   , sum_wl=0.;
    amrex::Real phi1=0. , phi2=0. , phi3=0. ;

    // WENO3
    if (spatial_order == 3) {
        beta1  = ( qty(i,j,k-1,qty_index) - qty(i,j,k-2,qty_index) ) *
                 ( qty(i,j,k-1,qty_index) - qty(i,j,k-2,qty_index) );
        beta2  = ( qty(i,j,k  ,qty_index) - qty(i,j,k-1,qty_index) ) *
                 ( qty(i,j,k  ,qty_index) - qty(i,j,k-1,qty_index) );

        w1     = ( (1./3.) / ((eps+beta1)*(eps+beta1)) );
        w2     = ( (2./3.) / ((eps+beta2)*(eps+beta2)) );

        phi1   = 0.5 * ( -qty(i,j,k-2,qty_index) + 3.*qty(i,j,k-1,qty_index) );
        phi2   = 0.5 * (  qty(i,j,k-1,qty_index) +    qty(i,j,k  ,qty_index) );
    // WENO5
    } else {
        beta1  = (13./12.) * (    qty(i,j,k-3,qty_index) - 2.*qty(i,j,k-2,qty_index) +    qty(i,j,k-1,qty_index) ) *
                             (    qty(i,j,k-3,qty_index) - 2.*qty(i,j,k-2,qty_index) +    qty(i,j,k-1,qty_index) );
        beta1 += (1./4.)   * (    qty(i,j,k-3,qty_index) - 4.*qty(i,j,k-2,qty_index) + 3.*qty(i,j,k-1,qty_index) ) *
                             (    qty(i,j,k-3,qty_index) - 4.*qty(i,j,k-2,qty_index) + 3.*qty(i,j,k-1,qty_index) );

        beta2  = (13./12.) * (    qty(i,j,k-2,qty_index) - 2.*qty(i,j,k-1,qty_index) +    qty(i,j,k  ,qty_index) ) *
                             (    qty(i,j,k-2,qty_index) - 2.*qty(i,j,k-1,qty_index) +    qty(i,j,k  ,qty_index) );
        beta2 += (1./4.)   * (    qty(i,j,k-2,qty_index) -    qty(i,j,k  ,qty_index) ) *
                             (    qty(i,j,k-2,qty_index) -    qty(i,j,k  ,qty_index) );

        beta3  = (13./12.) * (    qty(i,j,k-1,qty_index) - 2.*qty(i,j,k  ,qty_index) +    qty(i,j,k+1,qty_index) ) *
                             (    qty(i,j,k-1,qty_index) - 2.*qty(i,j,k  ,qty_index) +    qty(i,j,k+1,qty_index) );
        beta3 += (1./4.)   * ( 3.*qty(i,j,k-1,qty_index) - 4.*qty(i,j,k  ,qty_index) +    qty(i,j,k+1,qty_index) ) *
                             ( 3.*qty(i,j,k-1,qty_index) - 4.*qty(i,j,k  ,qty_index) +    qty(i,j,k+1,qty_index) );

        w1   = ( (1./10.) / ((eps+beta1)*(eps+beta1)) );
        w2   = ( (3./5. ) / ((eps+beta2)*(eps+beta2)) );
        w3   = ( (3./10.) / ((eps+beta3)*(eps+beta3)) );

        phi1 =  (1./3.)*qty(i,j,k-3,qty_index) - (7./6.)*qty(i,j,k-2,qty_index) + (11./6.)*qty(i,j,k-1,qty_index);
        phi2 = -(1./6.)*qty(i,j,k-2,qty_index) + (5./6.)*qty(i,j,k-1,qty_index) + ( 1./3.)*qty(i,j,k  ,qty_index);
        phi3 =  (1./3.)*qty(i,j,k-1,qty_index) + (5./6.)*qty(i,j,k  ,qty_index) - ( 1./6.)*qty(i,j,k+1,qty_index);
    }

    sum_wl = w1 + w2 + w3;
    w1    /= sum_wl;
    w2    /= sum_wl;
    w3    /= sum_wl;
    return w1*phi1 + w2*phi2 + w3*phi3;
}
#endif
