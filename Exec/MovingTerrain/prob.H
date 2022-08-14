#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

AMREX_FORCE_INLINE
amrex::Real
dhdt(int i, int j, int k, amrex::Real time)
{
    return 0.;
}

struct ProbParm {
  amrex::Real Ampl = 0.16;
  amrex::Real T_0  = 300.;
}; // namespace ProbParm

extern ProbParm parms;

#endif

