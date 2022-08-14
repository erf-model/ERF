#ifndef _PROB_H_
#define _PROB_H_

#include <AMReX_REAL.H>

AMREX_FORCE_INLINE
amrex::Real
dhdt(int i, int j, int k, amrex::Real time)
{
    return 0.;
}

struct ProbParm {
  amrex::Real rho_0 = 1.0;
  amrex::Real Theta_0 = 300.0;
  amrex::Real V_0 = 1.0;
}; // namespace ProbParm

extern ProbParm parms;

#endif
