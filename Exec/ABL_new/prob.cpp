#include "prob.H"
#include "ABLFieldInit.H"
#include "ERF.H"

ProbParm parms;

void
erf_init_dens_hse(amrex::Vector<amrex::Real>& dens_hse)
{
  const int klen = dens_hse.size();
  for (int k = 0; k < klen; k++)
  {
      dens_hse[k] = parms.rho_0;
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

AMREX_GPU_DEVICE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const amrex::Real s_int[NVAR],
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata)
{
  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = s_int[n];
  }
}

void
erf_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* init,
  const int* name,
  const int* namelen,
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
  pp.query("A_0", parms.A_0);

  pp.query("U0", parms.U0);
  pp.query("V0", parms.V0);
  pp.query("W0", parms.W0);
  pp.query("U0_Pert_Mag", parms.U0_Pert_Mag);
  pp.query("V0_Pert_Mag", parms.V0_Pert_Mag);
  pp.query("W0_Pert_Mag", parms.W0_Pert_Mag);
}
}
