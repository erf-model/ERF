/*
 * utilitity tools for microphysics are here
 * 
 */
#ifndef Microphysics_Utils_H
#define Microphysics_Utils_H

#include <cmath>
#include <vector>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real erf_gammafff(Real x){
  return exp(lgamma(x));
}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real erf_esati(Real t) {
  Real const a0 = 6.11147274;
  Real const a1 = 0.503160820;
  Real const a2 = 0.188439774e-1;
  Real const a3 = 0.420895665e-3;
  Real const a4 = 0.615021634e-5;
  Real const a5 = 0.602588177e-7;
  Real const a6 = 0.385852041e-9;
  Real const a7 = 0.146898966e-11;
  Real const a8 = 0.252751365e-14;

  Real dtt = t-273.16;
  Real esati;
  if(dtt > -80.0) {
    esati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    esati = 0.01*exp(9.550426 - 5723.265/t + 3.53068*log(t) - 0.00728332*t);
  }
  return esati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real erf_esatw(Real t) {
  Real const a0 = 6.105851;
  Real const a1 = 0.4440316;
  Real const a2 = 0.1430341e-1;
  Real const a3 = 0.2641412e-3;
  Real const a4 = 0.2995057e-5;
  Real const a5 = 0.2031998e-7;
  Real const a6 = 0.6936113e-10;
  Real const a7 = 0.2564861e-13;
  Real const a8 = -0.3704404e-15;

  Real dtt = t-273.16;

  Real esatw;
  if(dtt > -80.0) {
    esatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    esatw = 2.0*0.01*exp(9.550426 - 5723.265/t + 3.53068*log(t) - 0.00728332*t);
  }
  return esatw;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real erf_dtesati(Real t) {
  Real const a0 = 0.503223089;
  Real const a1 = 0.377174432e-1;
  Real const a2 = 0.126710138e-2;
  Real const a3 = 0.249065913e-4;
  Real const a4 = 0.312668753e-6;
  Real const a5 = 0.255653718e-8;
  Real const a6 = 0.132073448e-10;
  Real const a7 = 0.390204672e-13;
  Real const a8 = 0.497275778e-16;

  Real dtt = t-273.16;
  Real dtesati;
  if(dtt > -80.0) {
    dtesati = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    dtesati= erf_esati(t+1.0)-erf_esati(t);
  }

  return dtesati;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real erf_dtesatw(Real t) {
  Real const a0 = 0.443956472;
  Real const a1 = 0.285976452e-1;
  Real const a2 = 0.794747212e-3;
  Real const a3 = 0.121167162e-4;
  Real const a4 = 0.103167413e-6;
  Real const a5 = 0.385208005e-9;
  Real const a6 = -0.604119582e-12;
  Real const a7 = -0.792933209e-14;
  Real const a8 = -0.599634321e-17;

  Real dtt = t-273.16;
  Real dtesatw;
  if(dtt > -80.0) {
    dtesatw = a0 + dtt*(a1+dtt*(a2+dtt*(a3+dtt*(a4+dtt*(a5+dtt*(a6+dtt*(a7+a8*dtt)))))));
  }
  else {
    dtesatw = erf_esatw(t+1.0)-erf_esatw(t);
  }
  return dtesatw;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsati(Real t, Real p, Real &qsati) {
  Real esati;
  esati = erf_esati(t);
  qsati = 0.622*esati/max(esati,p-esati);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_qsatw(Real t, Real p, Real &qsatw) {
  Real esatw;
  esatw = erf_esatw(t);
  qsatw = 0.622*esatw/max(esatw,p-esatw);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsati(Real t, Real p, Real &dtqsati) {
  dtqsati = 0.622*erf_dtesati(t)/p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void erf_dtqsatw(Real t, Real p, Real &dtqsatw) {
  dtqsatw = 0.622*erf_dtesatw(t)/p;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void z0_est(Real z, Real bflx, Real wnd, Real ustar, Real &z0) {
  Real vonk = 0.4;
  Real eps = 1.0e-10;
  Real am = 4.8;
  Real bm = 19.3;
  Real c1 = 3.14159/2.0 - 3.0*log(2.0);
  Real rlmo = -bflx*vonk/(ustar*ustar*ustar+eps);
  Real zeta = min(1.0,z*rlmo);
  Real x;
  Real psi1;
  if(zeta >= 0.0) {
    psi1 = -am*zeta;
  }
  else {
    x = sqrt(sqrt(1.0-bm*zeta));
    psi1 = 2.0*log(1.0+x) + log(1.0+x*x) -2.0*atan(x) + c1;
  }
  Real lnz = max(0.0, vonk*wnd/(ustar+eps) +psi1);
  z0 = z*exp(-lnz);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real term_vel_qp(int i, int j, int k, Real qploc, Real vrain, Real vsnow,
                 Real vgrau, Real crain, Real csnow, Real cgrau, Real rho,
                 Real tabs, Real a_pr, Real a_gr) {
  Real term_vel = 0.0;
  if(qploc > qp_threshold) {
    Real omp = max(0.0,min(1.0,(tabs-tprmin)*a_pr));
    if(omp == 1.0) {
      term_vel = vrain*pow(rho*qploc,crain);
    }
    else if(omp == 0.0) {
      Real omg = max(0.0,min(1.0,(tabs-tgrmin)*a_gr));
      Real qgg=omg*qploc;
      Real qss=qploc-qgg;
      term_vel = (omg*vgrau*pow(rho*qgg,cgrau) + (1.0-omg)*vsnow*pow(rho*qss,csnow));
    }
    else {
      Real omg = max(0.0,min(1.0,(tabs-tgrmin)*a_gr));
      Real qrr=omp*qploc;
      Real qss=qploc-qrr;
      Real qgg=omg*qss;
      qss = qss-qgg;
      term_vel = (omp*vrain*pow(rho*qrr,crain) + (1.0-omp)*(omg*vgrau*pow(rho*qgg,cgrau) + (1.0-omg)*vsnow*pow(rho*qss,csnow)));
    }
  }
  return term_vel;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real pp(Real y) {
  return max(0.0,y);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real pn(Real y) {
  return -min(0.0,y);
}

#endif
