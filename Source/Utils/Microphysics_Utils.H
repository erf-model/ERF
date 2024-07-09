/*
 * utility tools for microphysics
 *
 */
#ifndef Microphysics_Utils_H
#define Microphysics_Utils_H

#include <cmath>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <ERF_Constants.H>

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_gammafff (amrex::Real x){
    return std::exp(lgamma(x));
}


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
    amrex::Real esati;
    if(dtt > -80.0) {
        esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    else {
        esati = 0.01*std::exp(9.550426 - 5723.265/t + 3.53068*std::log(t) - 0.00728332*t);
    }
    return esati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real erf_esatw (amrex::Real t) {
#if 1
    amrex::Real const a0 = 6.105851;
    amrex::Real const a1 = 0.4440316;
    amrex::Real const a2 = 0.1430341e-1;
    amrex::Real const a3 = 0.2641412e-3;
    amrex::Real const a4 = 0.2995057e-5;
    amrex::Real const a5 = 0.2031998e-7;
    amrex::Real const a6 = 0.6936113e-10;
    amrex::Real const a7 = 0.2564861e-13;
    amrex::Real const a8 = -0.3704404e-15;

    amrex::Real dtt = t-273.16;

    amrex::Real esatw;
    if(dtt > -80.0) {
        esatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    else {
        esatw = 2.0*0.01*std::exp(9.550426 - 5723.265/t + 3.53068*std::log(t) - 0.00728332*t);
    }
    return esatw;
#else
    amrex::Real a0 = 6.11239921;
    amrex::Real a1 = 0.443987641;
    amrex::Real a2 = 0.142986287e-1;
    amrex::Real a3 = 0.264847430e-3;
    amrex::Real a4 = 0.302950461e-5;
    amrex::Real a5 = 0.206739458e-7;
    amrex::Real a6 = 0.640689451e-10;
    amrex::Real a7 = -0.952447341e-13;
    amrex::Real a8 = -0.976195544e-15;

    amrex::Real a0i = 6.11147274;
    amrex::Real a1i = 0.503160820;
    amrex::Real a2i = 0.188439774e-1;
    amrex::Real a3i = 0.420895665e-3;
    amrex::Real a4i = 0.615021634e-5;
    amrex::Real a5i = 0.602588177e-7;
    amrex::Real a6i = 0.385852041e-9;
    amrex::Real a7i = 0.146898966e-11;
    amrex::Real a8i = 0.252751365e-14;

    amrex::Real Tc = std::max(-80.0, t-273.16);

    amrex::Real esatw;
    if (t>=273.15){
        esatw = (a0  + Tc*(a1 +Tc*(a2 +Tc*(a3 +Tc*(a4 +Tc*(a5 +Tc*(a6 +Tc*(a7 +a8 *Tc))))))))*100.;
    }else{
        esatw = (a0i + Tc*(a1i+Tc*(a2i+Tc*(a3i+Tc*(a4i+Tc*(a5i+Tc*(a6i+Tc*(a7i+a8i*Tc))))))))*100.;
    }
    return esatw;
#endif
}

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
    amrex::Real dtesati;
    if(dtt > -80.0) {
        dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    else {
        dtesati= erf_esati(t+1.0)-erf_esati(t);
    }

    return dtesati;
}

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
    amrex::Real dtesatw;
    if(dtt > -80.0) {
        dtesatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
    }
    else {
        dtesatw = erf_esatw(t+1.0)-erf_esatw(t);
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
#if 1
    qsatw = Rd_on_Rv*esatw/std::max(esatw,p-esatw);
#else
    qsatw = esatw/(R_v*t);
#endif
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
