#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 1.0;
  amrex::Real M_0   = 1.0;
  amrex::Real V_0   = 1.0;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitCustomPert.H"

protected:

private:
    ProbParm parms;
};

#endif
