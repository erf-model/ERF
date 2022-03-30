#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

struct ProbParm {
  amrex::Real rho_0   = 0.0;
  amrex::Real Theta_0 = 300.0;
  amrex::Real A_0     = 1.0;
  amrex::Real QKE_0     = 0.1;

  amrex::Real U_0 = 10.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0;
  amrex::Real U_0_Pert_Mag = 1.0;
  amrex::Real V_0_Pert_Mag = 1.0;
  amrex::Real W_0_Pert_Mag = 0.0;
}; // namespace ProbParm

extern ProbParm parms;

#endif
