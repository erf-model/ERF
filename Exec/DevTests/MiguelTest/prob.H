#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

struct ProbParm {
  amrex::Real rho_0 = 1.0;
  amrex::Real Theta_0 = 300.0;
  amrex::Real QKE_0   = 0.1;

  amrex::Real Xc_0 = 1000;
  amrex::Real Yc_0 = 1000;
  amrex::Real VMAX = 15;
  amrex::Real RMAX = 100;
  amrex::Real RZERO = 800;
  amrex::Real ZZERO = 2000;

  amrex::Real T_0 = 300.0;
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0;

  // rayleigh damping
  amrex::Real dampcoef = 0.2; // inverse time scale [1/s]
  amrex::Real zdamp = 500.0; // damping depth [m] from model top
      
}; // namespace ProbParm

extern ProbParm parms;

#endif
