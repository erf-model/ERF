#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {
  amrex::Real T_0 = 300.0; // surface temperature == mean potential temperature
  amrex::Real U_0 = 10.0;
  amrex::Real Theta_0 = 1.0;
  amrex::Real V_0 = 0.0;
  amrex::Real x_c = 0.0; // center of thermal perturbation
  amrex::Real z_c = 3200.0;
  amrex::Real x_r = 1000.0;
  amrex::Real z_r = 1000.0;
  amrex::Real T_pert = 5.0; // perturbation temperature
  // overridden physical constants
  amrex::Real C_p = 1004.0;
}; // namespace ProbParm

extern ProbParm parms;

#endif

