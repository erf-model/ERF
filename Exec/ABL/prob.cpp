#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_0 = 0.0; // left density (g/cc)
AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_0;
AMREX_GPU_DEVICE_MANAGED amrex::Real Theta_0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real A_0 = 1.0;

AMREX_GPU_DEVICE_MANAGED amrex::Real U0 = 10.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real V0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real W0 = 0.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real U0_Pert_Mag = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real V0_Pert_Mag = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real W0_Pert_Mag = 0.0;
} // namespace ProbParm

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
  pp.query("rho_0", ProbParm::rho_0);
  pp.query("T_0", ProbParm::Theta_0);
  pp.query("A_0", ProbParm::A_0);

  pp.query("U0", ProbParm::U0);
  pp.query("V0", ProbParm::V0);
  pp.query("W0", ProbParm::W0);
  pp.query("U0_Pert_Mag", ProbParm::U0_Pert_Mag);
  pp.query("V0_Pert_Mag", ProbParm::V0_Pert_Mag);
  pp.query("W0_Pert_Mag", ProbParm::W0_Pert_Mag);
}
}
