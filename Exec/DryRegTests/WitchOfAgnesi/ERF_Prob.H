#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"
#include "ERF_ProbCommon.H"

struct ProbParmCore : ProbParmDefaults {
  // perturbations to initial conditions
  amrex::Real rho_0 = 0.0;
  amrex::Real T_0 = 0.0;
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real W_0 = 0.0; // needed for rayleigh damping

  // hill parameters
  amrex::Real hmax = 0.0;   // full hill height, can be negative to simulate valleys
  amrex::Real L = 0.0;      // hill length at half-height
  amrex::Real z_offset = 0.0;  // shift in height

  // Default to the x-direction
  int dir = 0;
}; // namespace ProbParm

struct ProbParm : ProbParmCore {
  std::string custom_terrain_type;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem();

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitRayleighDamping.H"
#include "Prob/ERF_InitCustomPert.H"

protected:

private:
    ProbParm parms;
};

#endif
