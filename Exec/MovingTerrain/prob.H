#ifndef _PROB_H_
#define _PROB_H_

#include "AMReX_REAL.H"

struct ProbParm {
  amrex::Real Ampl = 0.16;
  amrex::Real T_0  = 300.;
}; // namespace ProbParm

extern ProbParm parms;

#endif

