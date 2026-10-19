#ifndef ERF_UTILS_WSTAR_H_
#define ERF_UTILS_WSTAR_H_

#include <ERF_Constants.H>

/*
 * Whether a PBL height has actually been diagnosed.  The SurfaceLayer's pblh holds
 * bogus_large_value until a scheme writes one: on the first call to update_fluxes, on a
 * level that was just made or regridded, and in any cell no scheme reaches.  A w* built
 * from that is ~1e50 and destroys the surface fluxes, so we take w* = 0 there instead.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
pblh_is_diagnosed (const amrex::Real& pblh)
{
    return pblh < amrex::Real(0.5) * bogus_large_value;
}

/*
 * Calculate the convective velocity scale
 */

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_wstar (const amrex::Real& ust,
            const amrex::Real& tst,
            const amrex::Real& qst,
            const amrex::Real& pblh,
            const amrex::Real& th,
            const amrex::Real& thv,
            const amrex::Real& qv = amrex::Real(0))
{
    amrex::Real hfx = -ust*tst; // <theta'w'>
    if (hfx <= 0 || !pblh_is_diagnosed(pblh)) {
        return amrex::Real(0);
    }
    if (qv > 0) {
        // have moisture -- estimate buoyancy flux, <thetav'w'> (Stull 4.4.5d)
        hfx = hfx*(one + epsv*qv) + epsv*th*(-ust*qst);
    }
    return std::cbrt(CONST_GRAV / thv * pblh * hfx);
}

/**
 * Calculate the convective velocity scale from buoyancy flux.
 *
 * @param[in] hfx Buoyancy flux.
 * @param[in] pblh PBL height.
 * @param[in] thv Virtual potential temperature.
 * @return Convective velocity scale.
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_wstar (const amrex::Real& hfx,
            const amrex::Real& pblh,
            const amrex::Real& thv)
{
    if (hfx <= 0 || !pblh_is_diagnosed(pblh)) {
        return amrex::Real(0);
    } else {
        return std::cbrt(CONST_GRAV / thv * pblh * hfx);
    }
}
#endif
