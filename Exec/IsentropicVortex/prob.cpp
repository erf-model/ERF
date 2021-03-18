#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real M_inf = 0.2; //freestream Mach number [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real alpha = 0.0; //inflow angle, 0 --> x-aligned [rad] 
AMREX_GPU_DEVICE_MANAGED amrex::Real gamma = 1.4; //specific heat ratio [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real beta = 1.0; //non-dimensional max perturbation strength [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real sigma = 1.0; //Gaussian standard deviation, i.e., spreading parameter [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real R = 1.0; //characteristic length scale for grid [m]
} // namespace ProbParm

void
pc_prob_close()
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
  pp.query("M_inf", ProbParm::M_inf);
  pp.query("alpha", ProbParm::alpha);
  pp.query("gamma", ProbParm::gamma);
  pp.query("beta", ProbParm::beta);
  pp.query("sigma", ProbParm::sigma);
  pp.query("R", ProbParm::R);
}
}
