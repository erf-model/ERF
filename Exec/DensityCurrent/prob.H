#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {
  amrex::Real T_0 = 300.0; // surface temperature == mean potential temperature
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real x_c = 0.0; // center of thermal perturbation
  amrex::Real z_c = 3000.0;
  amrex::Real x_r = 4000.0;
  amrex::Real z_r = 2000.0;
  amrex::Real T_pert = -15.0; // perturbation temperature
  // overridden physical constants
  amrex::Real local_c_p = 1004.0; // to match Straka study
}; // namespace ProbParm

extern ProbParm parms;

#endif

