#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  // perturbations to initial conditions
  amrex::Real rho_0 = 0.0;
  amrex::Real T_0 = 0.0;
  amrex::Real U_0 = 0.0; // Needed for Rayleigh damping
  amrex::Real V_0 = 0.0; // Needed for Rayleigh damping
  amrex::Real W_0 = 0.0; // Needed for Rayleigh damping
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitRayleighDamping.H"
#include "Prob/ERF_InitCustomPert.H"

protected:

private:
    ProbParm parms;
};

#endif
