#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include <AMReX_REAL.H>

#include "ERF_prob_common.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = 1.2;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_init_constant_density_hse.H"

protected:
    std::string name() override { return "Stokes Second Problem"; }

private:
    ProbParm parms;
};

#endif
