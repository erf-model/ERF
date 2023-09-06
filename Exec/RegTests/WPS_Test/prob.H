#ifndef _PROB_H_
#define _PROB_H_

#include <string>

#include <AMReX_REAL.H>

#include "prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 1.0;
  amrex::Real T_0   = 300.0;
  amrex::Real V_0   = 1.0;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/init_constant_density_hse.H"

protected:
    std::string name() override { return "WPS Test"; }

private:
    ProbParm parms;
};

#endif
