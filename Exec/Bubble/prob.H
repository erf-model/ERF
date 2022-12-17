#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {
  // background conditions
  // if init_type != "" then these are perturbations
  amrex::Real T_0 = 0.0;
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;

  // center of thermal perturbation
  amrex::Real x_c = 0.0;
  amrex::Real y_c = 0.0;
  amrex::Real z_c = 0.0;

  // radial extent of thermal perturbation
  amrex::Real x_r = 0.0;
  amrex::Real y_r = 0.0;
  amrex::Real z_r = 0.0;

  // perturbation temperature
  amrex::Real T_pert = -15.0;

  // physical constants that may be overwritten
  amrex::Real c_p = 1004.5; 
}; // namespace ProbParm

extern ProbParm parms;

#endif

