#include "prob.H"
#include "prob_common.H"

#include "AMReX_ParmParse.H"
#include "ERF.H"
#include "ABLFieldInit.H"

ProbParm parms;

void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::GeometryData const& geomdata,
                  const int ng_dens_hse)
{
  const int khi = geomdata.Domain().bigEnd()[2];
  for (int k = -ng_dens_hse; k <= khi+ng_dens_hse; k++)
  {
      dens_hse_ptr[k] = parms.rho_0;
  }
}

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& tau,
                  amrex::Vector<amrex::Real>& ubar,
                  amrex::Vector<amrex::Real>& vbar,
                  amrex::Vector<amrex::Real>& thetabar,
                  amrex::GeometryData        const& geomdata)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const auto dx              = geomdata.CellSize();
  const int khi              = geomdata.Domain().bigEnd()[2];

  for (int k = 0; k <= khi; k++)
  {
      const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

      tau[k]      = 1.0;

      ubar[k]     = parms.U_0;
      vbar[k]     = parms.V_0;
      thetabar[k] = parms.Theta_0;
  }
}

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
    ERF::ablinit(bx,state,x_vel,y_vel,z_vel,geomdata);
}

void
amrex_probinit(
  const amrex_real* /**problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
  pp.query("A_0", parms.A_0);

  pp.query("U_0", parms.U_0);
  pp.query("V_0", parms.V_0);
  pp.query("W_0", parms.W_0);
  pp.query("U_0_Pert_Mag", parms.U_0_Pert_Mag);
  pp.query("V_0_Pert_Mag", parms.V_0_Pert_Mag);
  pp.query("W_0_Pert_Mag", parms.W_0_Pert_Mag);
}
