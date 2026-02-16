#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0;
  amrex::Real   T_0;
};

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* problo, const amrex::Real* probhi);

private:
    ProbParm parms;

public:

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitCustomPert.H"
};

#endif
