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
    return std::exp(lgamma(x));
}

// From Flatau et al. (1992):
// https://doi.org/10.1175/1520-0450(1992)031<1507:PFTSVP>2.0.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-90  0] C. Return 0 if above this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esati (amrex::Real t) {
    amrex::Real const a0 = 6.11147274;
    amrex::Real const a1 = 0.503160820;
    amrex::Real const a2 = 0.188439774e-1;
    amrex::Real const a3 = 0.420895665e-3;
    amrex::Real const a4 = 0.615021634e-5;
    amrex::Real const a5 = 0.602588177e-7;
    amrex::Real const a6 = 0.385852041e-9;
    amrex::Real const a7 = 0.146898966e-11;
    amrex::Real const a8 = 0.252751365e-14;

    amrex::Real dtt = t-273.16;
    AMREX_ALWAYS_ASSERT(dtt>-85);

    amrex::Real esati;
    if (dtt > 0.0) {
        esati = 0.0;
    } else {
        esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return esati;
}

// From Flatau et al. (1992):
// https://doi.org/10.1175/1520-0450(1992)031<1507:PFTSVP>2.0.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-85  70] C. Assert we are in this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw (amrex::Real t) {
    amrex::Real const a0 = 6.11239921;
    amrex::Real const a1 = 0.443987641;
    amrex::Real const a2 = 0.142986287e-1;
    amrex::Real const a3 = 0.264847430e-3;
    amrex::Real const a4 = 0.302950461e-5;
    amrex::Real const a5 = 0.206739458e-7;
    amrex::Real const a6 = 0.640689451e-10;
    amrex::Real const a7 = -0.952447341e-13;
    amrex::Real const a8 = -0.976195544e-15;

    amrex::Real dtt_in = t-273.16;
    AMREX_ALWAYS_ASSERT(dtt_in>-85);
    //
    // Rather than assert dt < 70, we impose that the saturation pressure does not change above 70
    //
    amrex::Real dtt = std::min(dtt_in,70.);
    amrex::Real esatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    return esatw;
}

// From Flatau et al. (1992):
// https://doi.org/10.1175/1520-0450(1992)031<1507:PFTSVP>2.0.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-90  0] C. Return 0 if above this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesati (amrex::Real t) {
    amrex::Real const a0 = 0.503223089;
    amrex::Real const a1 = 0.377174432e-1;
    amrex::Real const a2 = 0.126710138e-2;
    amrex::Real const a3 = 0.249065913e-4;
    amrex::Real const a4 = 0.312668753e-6;
    amrex::Real const a5 = 0.255653718e-8;
    amrex::Real const a6 = 0.132073448e-10;
    amrex::Real const a7 = 0.390204672e-13;
    amrex::Real const a8 = 0.497275778e-16;

    amrex::Real dtt = t-273.16;
    AMREX_ALWAYS_ASSERT(dtt>-85);
    amrex::Real dtesati;
    if (dtt > 0.0) {
        dtesati = 0.0;
    } else {
        dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    return dtesati;
}

// From Flatau et al. (1992):
// https://doi.org/10.1175/1520-0450(1992)031<1507:PFTSVP>2.0.CO;2
// Coefficients come from Table 4 and the data is valid over a
// temperature range of [-85  70] C. Assert we are in this temp range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_dtesatw (amrex::Real t) {
    amrex::Real const a0 = 0.443956472;
    amrex::Real const a1 = 0.285976452e-1;
    amrex::Real const a2 = 0.794747212e-3;
    amrex::Real const a3 = 0.121167162e-4;
    amrex::Real const a4 = 0.103167413e-6;
    amrex::Real const a5 = 0.385208005e-9;
    amrex::Real const a6 = -0.604119582e-12;
    amrex::Real const a7 = -0.792933209e-14;
    amrex::Real const a8 = -0.599634321e-17;

    amrex::Real dtt = t-273.16;
    AMREX_ALWAYS_ASSERT(dtt>-85);

    amrex::Real dtesatw;

    //
    // Rather than assert dt < 70, we impose that the saturation pressure does not change above 70,
    //     hence the derivative is 0.
    //
    if (dtt >= 70.) {
        dtesatw = amrex::Real(0.0);
    } else {
        dtesatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }

    return dtesatw;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsati (amrex::Real t, amrex::Real p, amrex::Real &qsati) {
    amrex::Real esati;
    esati = erf_esati(t);
    qsati = Rd_on_Rv*esati/std::max(esati,p-esati);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsatw (amrex::Real t, amrex::Real p, amrex::Real &qsatw) {
    amrex::Real esatw;
    esatw = erf_esatw(t);
    qsatw = Rd_on_Rv*esatw/std::max(esatw,p-esatw);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsati (amrex::Real t, amrex::Real p, amrex::Real &dtqsati) {
    dtqsati = Rd_on_Rv*erf_dtesati(t)/p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsatw (amrex::Real t, amrex::Real p, amrex::Real &dtqsatw) {
    dtqsatw = Rd_on_Rv*erf_dtesatw(t)/p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void z0_est (amrex::Real z, amrex::Real bflx, amrex::Real wnd, amrex::Real ustar, amrex::Real &z0) {
    amrex::Real vonk = 0.4;
    amrex::Real eps = 1.0e-10;
    amrex::Real am = 4.8;
    amrex::Real bm = 19.3;
    amrex::Real c1 = 3.14159/2.0 - 3.0*log(2.0);
    amrex::Real rlmo = -bflx*vonk/(ustar*ustar*ustar+eps);
    amrex::Real zeta = std::min(1.0,z*rlmo);
    amrex::Real x;
    amrex::Real psi1;
    if(zeta >= 0.0) {
        psi1 = -am*zeta;
    }
    else {
        x = std::sqrt(sqrt(1.0-bm*zeta));
        psi1 = 2.0*std::log(1.0+x) + std::log(1.0+x*x) -2.0*std::atan(x) + c1;
    }
    amrex::Real lnz = std::max(0.0, vonk*wnd/(ustar+eps) +psi1);
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
    amrex::Real term_vel = 0.0;
    if(qploc > qp_threshold) {
        amrex::Real omp = std::max(0.0,std::min(1.0,(tabs-tprmin)*a_pr));
        amrex::Real omg = std::max(0.0,std::min(1.0,(tabs-tgrmin)*a_gr));
        amrex::Real qrr = omp*qploc;
        amrex::Real qss = (1.0-omp)*(1.0-omg)*qploc;
        amrex::Real qgg = (1.0-omp)*(omg)*qploc;
        term_vel = omp*vrain*std::pow(rho*qrr,crain)
                 + (1.0-omp)*( (1.0-omg)*vsnow*std::pow(rho*qss,csnow)
                             +      omg *vgrau*std::pow(rho*qgg,cgrau) );
    }
    return term_vel;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pp (amrex::Real y) {
    return std::max(0.0,y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real pn (amrex::Real y) {
    return -std::min(0.0,y);
}
#endif
