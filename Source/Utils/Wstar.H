#ifndef _UTILS_WSTAR_H_
#define _UTILS_WSTAR_H_

#include <ERF_Constants.H>

/*
 * Calculate the convective velocity scale
 */

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_wstar (const amrex::Real& ust,
            const amrex::Real& tst,
            const amrex::Real& qst,
            const amrex::Real& pblh,
            const amrex::Real& th,
            const amrex::Real& thv,
            const amrex::Real& qv = 0.0)
{
    amrex::Real hfx = -ust*tst; // <theta'w'>
    if (hfx <= 0) {
        return 0.0;
    }
    if (qv > 0) {
        // have moisture -- estimate buoyancy flux, <thetav'w'> (Stull 4.4.5d)
        hfx = hfx*(1 + 0.61*qv) + 0.61*th*(-ust*qst);
    }
    return std::cbrt(CONST_GRAV / thv * pblh * hfx);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_wstar (const amrex::Real& hfx,
            const amrex::Real& pblh,
            const amrex::Real& thv)
{
    if (hfx <= 0) {
        return 0.0;
    } else {
        return std::cbrt(CONST_GRAV / thv * pblh * hfx);
    }
}
#endif
