#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParmCore : ProbParmDefaults {
  amrex::Real T_0 = 300.0; // surface temperature == mean potential temperature
  amrex::Real U_0 = 0.0;
  amrex::Real V_0 = 0.0;
  amrex::Real x_c = 0.0; // center of thermal perturbation
  amrex::Real z_c = 3000.0;
  amrex::Real x_r = 4000.0;
  amrex::Real z_r = 2000.0;
  amrex::Real T_pert = -15.0; // perturbation temperature

  // Default to the x-direction
  int dir = 0;

  // overridden physical constants
  amrex::Real C_p = 1004.0;

};

struct ProbParm : ProbParmCore {
  std::string custom_terrain_type;
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
