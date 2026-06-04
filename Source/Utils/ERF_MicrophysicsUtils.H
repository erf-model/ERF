/*
 * utility tools for microphysics
 *
 */
#ifndef ERF_Microphysics_Utils_H
#define ERF_Microphysics_Utils_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <ERF_Constants.H>

// Positive-argument gamma wrapper. std::lgamma returns log(abs(Gamma(x))) for
// negative non-integers, so ERF keeps a positive-only contract here.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_gammafff (amrex::Real x){
    AMREX_ALWAYS_ASSERT(x > amrex::Real(0));
    return std::exp(std::lgamma(x));
}

// Saturation vapor pressure of ice [mbar / hPa].
//
// Flatau et al. (1992), Polynomial Fits to Saturation Vapor Pressure,
// J. Appl. Meteorol. Coefficients come from Table 4 and the value fit is valid
// over [-90, 0] C, i.e. down to about 183.16 K. ERF intentionally clamps the
// above-freezing branch to 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esati (amrex::Real t) {
    amrex::Real const a0 = amrex::Real(6.11147274);
    amrex::Real const a1 = amrex::Real(0.503160820);
    amrex::Real const a2 = amrex::Real(0.188439774e-1);
    amrex::Real const a3 = amrex::Real(0.420895665e-3);
    amrex::Real const a4 = amrex::Real(0.615021634e-5);
    amrex::Real const a5 = amrex::Real(0.602588177e-7);
    amrex::Real const a6 = amrex::Real(0.385852041e-9);
    amrex::Real const a7 = amrex::Real(0.146898966e-11);
    amrex::Real const a8 = amrex::Real(0.252751365e-14);

    amrex::Real dtt = t-amrex::Real(273.16);
    const amrex::Real tol = amrex::Real(16.0) * std::numeric_limits<amrex::Real>::epsilon();
    AMREX_ALWAYS_ASSERT(dtt >= -amrex::Real(90.0) - tol);

    amrex::Real esati;
    if (dtt > amrex::Real(0)) {
        esati = amrex::Real(0);
    } else {
        esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return esati;
}

// Magnus-style exponential saturation-pressure approximation over water
// [mbar / hPa]. This closed-form formula is commonly motivated by
// simplified integrations of the Clausius-Clapeyron relation, but ERF
// evaluates the empirical exponential expression below rather than
// performing a direct Clausius-Clapeyron integration. The `_cc` suffix is
// retained for API compatibility.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw_cc (amrex::Real t) {
    constexpr amrex::Real svp1  = amrex::Real(0.6112);
    constexpr amrex::Real svp2  = amrex::Real(17.67);
    constexpr amrex::Real svp3  = amrex::Real(29.65);
    constexpr amrex::Real svpt0 = amrex::Real(273.15);
    amrex::Real esatw = amrex::Real(10.0) * svp1 * std::exp(svp2 * (t - svpt0) / (t - svp3));
    return esatw;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw_flatau_poly (amrex::Real dtt)
{
    amrex::Real const a0 = amrex::Real(6.11239921);
    amrex::Real const a1 = amrex::Real(0.443987641);
    amrex::Real const a2 = amrex::Real(0.142986287e-1);
    amrex::Real const a3 = amrex::Real(0.264847430e-3);
    amrex::Real const a4 = amrex::Real(0.302950461e-5);
    amrex::Real const a5 = amrex::Real(0.206739458e-7);
    amrex::Real const a6 = amrex::Real(0.640689451e-10);
    amrex::Real const a7 = -amrex::Real(0.952447341e-13);
    amrex::Real const a8 = -amrex::Real(0.976195544e-15);

    return a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesatw_flatau_poly (amrex::Real dtt)
{
    amrex::Real const a0 = amrex::Real(0.443956472);
    amrex::Real const a1 = amrex::Real(0.285976452e-1);
    amrex::Real const a2 = amrex::Real(0.794747212e-3);
    amrex::Real const a3 = amrex::Real(0.121167162e-4);
    amrex::Real const a4 = amrex::Real(0.103167413e-6);
    amrex::Real const a5 = amrex::Real(0.385208005e-9);
    amrex::Real const a6 = -amrex::Real(0.604119582e-12);
    amrex::Real const a7 = -amrex::Real(0.792933209e-14);
    amrex::Real const a8 = -amrex::Real(0.599634321e-17);

    return a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool erf_use_positive_esatw_poly (amrex::Real t)
{
    amrex::Real const dtt = t - amrex::Real(273.16);
    amrex::Real const poly_min_dtt = -amrex::Real(70.0);
    return dtt > poly_min_dtt && dtt < amrex::Real(70.0) &&
           erf_esatw_flatau_poly(dtt) > amrex::Real(0.0);
}

// Saturation vapor pressure of water vapor [mbar / hPa].
//
// Flatau et al. (1992), Polynomial Fits to Saturation Vapor Pressure,
// J. Appl. Meteorol. The default branch uses the Flatau polynomial only on the
// current positive in-range interval selected by erf_use_positive_esatw_poly,
// currently about [-70, 70] C. Outside that interval, or whenever the
// polynomial would be non-positive, it falls back to the Magnus-style
// exponential approximation. The optional empirical fit is originally in Pa;
// ERF converts it to mbar / hPa before returning so both warm-water branches
// use the same units.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw (amrex::Real t, bool use_empirical=false)
{
    if (use_empirical) {
        return amrex::Real(0.01) *
               std::exp(amrex::Real(34.494) - amrex::Real(4924.99)/(t - amrex::Real(273.15) + amrex::Real(237.1))) /
               std::pow(t - amrex::Real(273.15) + amrex::Real(105.0), amrex::Real(1.57));
    }

    if (erf_use_positive_esatw_poly(t)) {
        return erf_esatw_flatau_poly(t - amrex::Real(273.16));
    }

    return erf_esatw_cc(t);
}

// Flatau et al. (1992), Polynomial Fits to Saturation Vapor Pressure,
// J. Appl. Meteorol. The derivative fit is valid over the narrower
// [-85, 0] C interval, i.e. down to about 188.16 K. ERF intentionally returns
// 0 above freezing to match erf_esati.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesati (amrex::Real t) {
    amrex::Real const a0 = amrex::Real(0.503223089);
    amrex::Real const a1 = amrex::Real(0.377174432e-1);
    amrex::Real const a2 = amrex::Real(0.126710138e-2);
    amrex::Real const a3 = amrex::Real(0.249065913e-4);
    amrex::Real const a4 = amrex::Real(0.312668753e-6);
    amrex::Real const a5 = amrex::Real(0.255653718e-8);
    amrex::Real const a6 = amrex::Real(0.132073448e-10);
    amrex::Real const a7 = amrex::Real(0.390204672e-13);
    amrex::Real const a8 = amrex::Real(0.497275778e-16);

    amrex::Real dtt = t-amrex::Real(273.16);
    const amrex::Real tol = amrex::Real(16.0) * std::numeric_limits<amrex::Real>::epsilon();
    AMREX_ALWAYS_ASSERT(dtt >= -amrex::Real(85.0) - tol);
    amrex::Real dtesati;
    if (dtt > amrex::Real(0)) {
        dtesati = amrex::Real(0);
    } else {
        dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return dtesati;
}

// Temperature derivative of the Magnus-style exponential water
// saturation-pressure approximation [mbar / hPa / K]. This differentiates
// the closed-form formula in erf_esatw_cc.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesatw_cc (amrex::Real t) {
    constexpr amrex::Real svp1  = amrex::Real(0.6112);
    constexpr amrex::Real svp2  = amrex::Real(17.67);
    constexpr amrex::Real svp3  = amrex::Real(29.65);
    constexpr amrex::Real svpt0 = amrex::Real(273.15);
    amrex::Real dtesatw = amrex::Real(10.0) * svp1 * svp2 * std::exp(svp2 * (t - svpt0) / (t - svp3))
                        * (svpt0 - svp3) / ((t - svp3) * (t - svp3));
    return dtesatw;
}

// Flatau et al. (1992), Polynomial Fits to Saturation Vapor Pressure,
// J. Appl. Meteorol. Use the same branch selection as erf_esatw so the value
// and derivative stay on the same function. Small derivative jumps can remain
// at the switch because the Flatau polynomial and Magnus-style exponential
// approximation are independent fits rather than a matched composite model.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesatw (amrex::Real t) {
    if (erf_use_positive_esatw_poly(t)) {
        return erf_dtesatw_flatau_poly(t - amrex::Real(273.16));
    }

    return erf_dtesatw_cc(t);
}

// Saturation mixing-ratio helper for a precomputed saturation pressure. qsat is
// capped at Rd_on_Rv whenever esat > p/2, i.e. when max(esat, p - esat) = esat.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_qsat_from_esat (amrex::Real esat, amrex::Real p)
{
    return Rd_on_Rv * esat / std::max(esat, p - esat);
}

// Temperature derivative of the capped saturation-mixing-ratio helper. The
// derivative is zero on the capped branch where qsat == Rd_on_Rv.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtqsat_from_esat (amrex::Real esat, amrex::Real dtesat, amrex::Real p)
{
    if ((p - esat) >= esat) {
        amrex::Real const denom = p - esat;
        return Rd_on_Rv * dtesat * p / (denom * denom);
    }

    return amrex::Real(0);
}

// Saturated ice vapor mixing ratio [-]. Input temperature is in K and input
// pressure is in mbar. The result is capped at Rd_on_Rv when esat > p/2.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsati (amrex::Real t, amrex::Real p, amrex::Real &qsati) {
    qsati = erf_qsat_from_esat(erf_esati(t), p);
}

/* saturated water vapor mixing ratio
 * t: temperature [K]
 * p: pressure [mbar]
 * return value is capped at Rd_on_Rv when esat > p/2.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsatw (amrex::Real t, amrex::Real p, amrex::Real &qsatw) {
    qsatw = erf_qsat_from_esat(erf_esatw(t), p);
}

// Temperature derivative of the capped erf_qsati relation. This is zero on the
// capped branch where qsat == Rd_on_Rv.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsati (amrex::Real t, amrex::Real p, amrex::Real &dtqsati) {
    amrex::Real esati   = erf_esati(t);
    amrex::Real dtesati = erf_dtesati(t);
    dtqsati = erf_dtqsat_from_esat(esati, dtesati, p);
}

// Temperature derivative of the capped erf_qsatw relation. This is zero on the
// capped branch where qsat == Rd_on_Rv.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsatw (amrex::Real t, amrex::Real p, amrex::Real &dtqsatw) {
    amrex::Real esatw   = erf_esatw(t);
    amrex::Real dtesatw = erf_dtesatw(t);
    dtqsatw = erf_dtqsat_from_esat(esatw, dtesatw, p);
}
#endif
