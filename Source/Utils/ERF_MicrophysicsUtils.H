/*
 * utility tools for microphysics
 *
 */
#ifndef ERF_Microphysics_Utils_H
#define ERF_Microphysics_Utils_H

#include <cmath>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <ERF_Constants.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_gammafff (amrex::Real x){
    return std::exp(std::lgamma(x));
}

// Saturation vapor pressure of ice [mbar]
//
// From Flatau et al. (1992):
// https://doi.org/amrex::Real(10.1175)/1520-0450(1992)031<1507:PFTSVP>two.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-90  0] C. Return 0 if above this temp range.
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
    AMREX_ALWAYS_ASSERT(dtt>-90);

    amrex::Real esati;
    if (dtt > amrex::Real(0)) {
        esati = amrex::Real(0);
    } else {
        esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return esati;
}

// From Clausius-Clapeyron
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw_cc (amrex::Real t) {
    constexpr amrex::Real svp1  = amrex::Real(0.6112);
    constexpr amrex::Real svp2  = amrex::Real(17.67);
    constexpr amrex::Real svp3  = amrex::Real(29.65);
    constexpr amrex::Real svpt0 = amrex::Real(273.15);
    // NOTE: units of
    amrex::Real esatw = amrex::Real(10.0) * svp1 * std::exp(svp2 * (t - svpt0) / (t - svp3));
    return esatw;
}

// Saturation vapor pressure of water vapor [mbar]
//
// From Flatau et al. (1992):
// https://doi.org/amrex::Real(10.1175)/1520-0450(1992)031<1507:PFTSVP>two.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-85  70] C. Assert we are in this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw (amrex::Real t, bool use_empirical=false)
{
    amrex::Real esatw;
    if (use_empirical) {
        esatw = std::exp(amrex::Real(34.494) - amrex::Real(4924.99)/(t - amrex::Real(273.15) + amrex::Real(237.1))) / std::pow(t - amrex::Real(273.15) + amrex::Real(105.0),amrex::Real(1.57));
    } else {
        amrex::Real const a0 = amrex::Real(6.11239921);
        amrex::Real const a1 = amrex::Real(0.443987641);
        amrex::Real const a2 = amrex::Real(0.142986287e-1);
        amrex::Real const a3 = amrex::Real(0.264847430e-3);
        amrex::Real const a4 = amrex::Real(0.302950461e-5);
        amrex::Real const a5 = amrex::Real(0.206739458e-7);
        amrex::Real const a6 = amrex::Real(0.640689451e-10);
        amrex::Real const a7 = -amrex::Real(0.952447341e-13);
        amrex::Real const a8 = -amrex::Real(0.976195544e-15);

        amrex::Real dtt = t-amrex::Real(273.16);
        if (dtt>-85 && dtt<amrex::Real(70.0)) {
            esatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
        } else {
            esatw = erf_esatw_cc(t);
        }
    }
    return esatw;
}

// From Flatau et al. (1992):
// https://doi.org/amrex::Real(10.1175)/1520-0450(1992)031<1507:PFTSVP>two.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-90  0] C. Return 0 if above this temp range.
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
    AMREX_ALWAYS_ASSERT(dtt>-85);
    amrex::Real dtesati;
    if (dtt > amrex::Real(0)) {
        dtesati = amrex::Real(0);
    } else {
        dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return dtesati;
}

// From Clausius-Clapeyron
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

// From Flatau et al. (1992):
// https://doi.org/amrex::Real(10.1175)/1520-0450(1992)031<1507:PFTSVP>two.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-85  70] C. Assert we are in this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesatw (amrex::Real t) {
    amrex::Real const a0 = amrex::Real(0.443956472);
    amrex::Real const a1 = amrex::Real(0.285976452e-1);
    amrex::Real const a2 = amrex::Real(0.794747212e-3);
    amrex::Real const a3 = amrex::Real(0.121167162e-4);
    amrex::Real const a4 = amrex::Real(0.103167413e-6);
    amrex::Real const a5 = amrex::Real(0.385208005e-9);
    amrex::Real const a6 = -amrex::Real(0.604119582e-12);
    amrex::Real const a7 = -amrex::Real(0.792933209e-14);
    amrex::Real const a8 = -amrex::Real(0.599634321e-17);

    amrex::Real dtt = t-amrex::Real(273.16); // Goff-Gratch
    amrex::Real dtesatw;
    if (dtt>-amrex::Real(85.0) && dtt<amrex::Real(70.)) {
        dtesatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    } else {
        dtesatw = erf_dtesatw_cc(t);
    }
    return dtesatw;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsati (amrex::Real t, amrex::Real p, amrex::Real &qsati) {
    amrex::Real esati;
    esati = erf_esati(t);
    qsati = Rd_on_Rv*esati/std::max(esati,p-esati);
}

/* saturated water vapor mixing ratio
 * t: temperature [K]
 * p: pressure [mbar]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsatw (amrex::Real t, amrex::Real p, amrex::Real &qsatw) {
    amrex::Real esatw;
    esatw = erf_esatw(t);
    qsatw = Rd_on_Rv*esatw/std::max(esatw,p-esatw);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsati (amrex::Real t, amrex::Real p, amrex::Real &dtqsati) {
    amrex::Real esati   = erf_esati(t);
    amrex::Real dtesati = erf_dtesati(t);
    amrex::Real denom   = p - esati;
    dtqsati = Rd_on_Rv * dtesati * p / (denom * denom);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsatw (amrex::Real t, amrex::Real p, amrex::Real &dtqsatw) {
    amrex::Real esatw   = erf_esatw(t);
    amrex::Real dtesatw = erf_dtesatw(t);
    amrex::Real denom   = p - esatw;
    dtqsatw = Rd_on_Rv * dtesatw * p / (denom * denom);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void z0_est (amrex::Real z, amrex::Real bflx, amrex::Real wnd, amrex::Real ustar, amrex::Real &z0) {
    amrex::Real vonk = amrex::Real(0.4);
    amrex::Real eps  = amrex::Real(1.0e-10);
    amrex::Real am   = amrex::Real(4.8);
    amrex::Real bm   = amrex::Real(19.3);
    amrex::Real c1   = amrex::Real(3.14159)/two - three*std::log(two);
    amrex::Real rlmo = -bflx*vonk/(ustar*ustar*ustar+eps);
    amrex::Real zeta = std::min(one,z*rlmo);
    amrex::Real x;
    amrex::Real psi1;
    if(zeta >= zero) {
        psi1 = -am*zeta;
    }
    else {
        x = std::sqrt(std::sqrt(one-bm*zeta));
        psi1 = two*std::log(one+x) + std::log(one+x*x) - two*std::atan(x) + c1;
    }
    amrex::Real lnz = std::max(amrex::Real(0), vonk*wnd/(ustar+eps) +psi1);
    z0 = z*std::exp(-lnz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real term_vel_qp (amrex::Real qploc,
                         amrex::Real vrain,
                         amrex::Real vsnow,
                         amrex::Real vgrau,
                         amrex::Real rho,
                         amrex::Real tabs)
{
    amrex::Real term_vel = amrex::Real(0);
    if(qploc > qp_threshold) {
        amrex::Real omp = std::max(amrex::Real(0),std::min(one,(tabs-tprmin)*a_pr));
        amrex::Real omg = std::max(amrex::Real(0),std::min(one,(tabs-tgrmin)*a_gr));
        amrex::Real qrr = omp*qploc;
        amrex::Real qss = (one-omp)*(one-omg)*qploc;
        amrex::Real qgg = (one-omp)*(omg)*qploc;
        term_vel = omp*vrain*std::pow(rho*qrr,crain)
                 + (one-omp)*( (one-omg)*vsnow*std::pow(rho*qss,csnow)
                             +      omg *vgrau*std::pow(rho*qgg,cgrau) );
    }
    return term_vel;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pp (amrex::Real y) {
    return std::max(amrex::Real(0),y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pn (amrex::Real y) {
    return -std::min(amrex::Real(0),y);
}
#endif
