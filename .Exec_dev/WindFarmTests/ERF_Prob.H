#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = amrex::Real(0.0);
  amrex::Real T_0   = amrex::Real(0.0);
  amrex::Real A_0   = amrex::Real(1.0);
  amrex::Real KE_0 = amrex::Real(0.1);

  amrex::Real U_0 = amrex::Real(0.0);
  amrex::Real V_0 = amrex::Real(0.0);
  amrex::Real W_0 = amrex::Real(0.0);

  // random initial perturbations (legacy code)
  amrex::Real U_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real V_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real W_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real T_0_Pert_Mag = amrex::Real(0.0); // perturbation to rho*Theta

  // divergence-free initial perturbations
  amrex::Real pert_deltaU = amrex::Real(0.0);
  amrex::Real pert_deltaV = amrex::Real(0.0);
  amrex::Real pert_periods_U = amrex::Real(5.0);
  amrex::Real pert_periods_V = amrex::Real(5.0);
  amrex::Real pert_ref_height = amrex::Real(100.0);

  // helper vars
  amrex::Real aval;
  amrex::Real bval;
  amrex::Real ufac;
  amrex::Real vfac;
}; // namespace ProbParm

class Problem : public ProblemBase
{
public:
    Problem(const amrex::Real* problo, const amrex::Real* probhi);

#include "Prob/ERF_InitDensityHSE.H"
#include "Prob/ERF_InitRayleighDamping.H"
#include "Prob/ERF_InitCustomPert.H"

private:
    ProbParm parms;
};

#endif
