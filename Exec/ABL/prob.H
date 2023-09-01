#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

struct ProbParm {
  amrex::Real rho_0 = 0.0;
  amrex::Real T_0   = 0.0;
  amrex::Real A_0   = 1.0;
  amrex::Real QKE_0 = 0.1;

  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0;

  // random initial perturbations (legacy code)
  amrex::Real U_0_Pert_Mag = 0.0;
  amrex::Real V_0_Pert_Mag = 0.0;
  amrex::Real W_0_Pert_Mag = 0.0;
  amrex::Real T_0_Pert_Mag = 0.0; // perturbation to rho*Theta

  // divergence-free initial perturbations
  amrex::Real pert_deltaU = 0.0;
  amrex::Real pert_deltaV = 0.0;
  amrex::Real pert_periods_U = 5.0;
  amrex::Real pert_periods_V = 5.0;
  amrex::Real pert_ref_height = 100.0;

  // rayleigh damping
  amrex::Real dampcoef = 0.2; // inverse time scale [1/s]
  amrex::Real zdamp = 500.0; // damping depth [m] from model top

  // helper vars
  amrex::Real aval;
  amrex::Real bval;
  amrex::Real ufac;
  amrex::Real vfac;
}; // namespace ProbParm

extern ProbParm parms;

#endif
