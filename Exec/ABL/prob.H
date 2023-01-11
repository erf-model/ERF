#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

struct ProbParm {
  amrex::Real rho_0   = 0.0;
  amrex::Real Theta_0 = 300.0;
  amrex::Real A_0     = 1.0;
  amrex::Real QKE_0   = 0.1;

  amrex::Real U_0 = 10.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0;

  // random perturbations (legacy code)
  amrex::Real U_0_Pert_Mag = 0.0;
  amrex::Real V_0_Pert_Mag = 0.0;
  amrex::Real W_0_Pert_Mag = 0.0;

  // divergence-free perturbations
  amrex::Real pert_deltaU = 1.0;
  amrex::Real pert_deltaV = 1.0;
  amrex::Real pert_periods_U = 5.0;
  amrex::Real pert_periods_V = 5.0;
  amrex::Real pert_ref_height = 100.0;

  // helper vars
  amrex::Real aval;
  amrex::Real bval;
  amrex::Real ufac;
  amrex::Real vfac;
}; // namespace ProbParm

extern ProbParm parms;

#endif
