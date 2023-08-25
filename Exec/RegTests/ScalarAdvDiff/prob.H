#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

struct ProbParm {
  amrex::Real rho_0   = 0.0;
  amrex::Real Theta_0 = 1.0;
  amrex::Real A_0     = 1.0;
  amrex::Real B_0     = 0.0;
  amrex::Real u_0     = 0.0;
  amrex::Real v_0     = 0.0;
  amrex::Real rad_0   = 0.0;
  amrex::Real z0      = 0.1;  // Surface Roughness
  amrex::Real zRef    = 80.0; // Reference Height
  amrex::Real uRef    = 0.0;  // Reference Wind Speed

  amrex::Real xc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real yc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)
  amrex::Real zc_frac = 0.5; // Location of "center" of scalar (multiplies domain length)

  int prob_type = -1;
}; // namespace ProbParm

extern ProbParm parms;

#endif
