#ifndef ERF_PROB_H_
#define ERF_PROB_H_

#include <string>

#include "AMReX_REAL.H"

#include "ERF_ProbCommon.H"

struct ProbParm : ProbParmDefaults {
  amrex::Real rho_0 = amrex::Real(0.0);
  amrex::Real T_0   = amrex::Real(0.0);
  amrex::Real A_0   = amrex::Real(1.0);
  amrex::Real KE_0  = amrex::Real(0.1);  // initialize to calculated rho_hse times input KE
  amrex::Real rhoKE_0 = amrex::Real(-1); // alternatively, initialize to specified (rho*KE)

  amrex::Real KE_decay_height = amrex::Real(-1);
  amrex::Real KE_decay_order  = amrex::Real( 1);

  amrex::Real U_0 = amrex::Real(0.0);
  amrex::Real V_0 = amrex::Real(0.0);
  amrex::Real W_0 = amrex::Real(0.0);

  // random initial perturbations (legacy code)
  amrex::Real U_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real V_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real W_0_Pert_Mag = amrex::Real(0.0);
  amrex::Real T_0_Pert_Mag = amrex::Real(0.0); // perturbation to rho*Theta
  bool pert_rhotheta = true;

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

#include "Prob/ERF_InitRayleighDamping.H"
#include "Prob/ERF_InitCustomPert.H"

      void update_w_subsidence (
        const amrex::Real& /*time*/,
        amrex::Vector<amrex::Real>& wbar,
        amrex::Gpu::DeviceVector<amrex::Real>& d_wbar,
        const amrex::MultiFab& /*state*/,
        const amrex::Geometry& geom,
        std::unique_ptr<amrex::MultiFab>& z_phys_nd) override;

protected:

private:
    ProbParm parms;
    amrex::Vector<amrex::Real> zlevels;
};

#endif
