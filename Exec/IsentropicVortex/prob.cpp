#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real p_inf = p_0; //freestream temperature [K]
AMREX_GPU_DEVICE_MANAGED amrex::Real T_inf = 300.0; //freestream temperature [K]
AMREX_GPU_DEVICE_MANAGED amrex::Real M_inf = 0.2; //freestream Mach number [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real alpha = 0.0; //inflow angle, 0 --> x-aligned [rad] 
AMREX_GPU_DEVICE_MANAGED amrex::Real gamma = Gamma; //specific heat ratio [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real beta = 0.01; //non-dimensional max perturbation strength [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real sigma = 2.5; //Gaussian standard deviation, i.e., spreading parameter [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real R = 1.0; //characteristic length scale for grid [m]
AMREX_GPU_DEVICE_MANAGED amrex::Real xc = 0.5; //normalized x-location of vortex center [-]
AMREX_GPU_DEVICE_MANAGED amrex::Real yc = 0.5; //normalized y-location of vortex center [-]
// calculated quantiites
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_inf; //characteristic density [m/s]
AMREX_GPU_DEVICE_MANAGED amrex::Real a_inf; //speed of sound [m/s]
AMREX_GPU_DEVICE_MANAGED amrex::Real inv_gm1; //1/(gamma - 1) [-]
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
  pp.query("p_inf", ProbParm::p_inf);
  pp.query("T_inf", ProbParm::T_inf);
  pp.query("M_inf", ProbParm::M_inf);
  pp.query("alpha", ProbParm::alpha);
  pp.query("gamma", ProbParm::gamma);
  pp.query("beta", ProbParm::beta);
  pp.query("sigma", ProbParm::sigma);
  pp.query("R", ProbParm::R);
  pp.query("xc", ProbParm::xc);
  pp.query("yc", ProbParm::yc);

  ProbParm::xc = problo[0] + ProbParm::xc * (probhi[0] - problo[0]);
  ProbParm::yc = problo[1] + ProbParm::yc * (probhi[1] - problo[1]);
  amrex::Print() << "  vortex initialized at ("
                 << ProbParm::xc << ", "
                 << ProbParm::yc << ")"
                 << std::endl;

  ProbParm::inv_gm1 = 1.0 / (ProbParm::gamma - 1.0);

  amrex::Print() << "  reference pressure = " << ProbParm::p_inf << " Pa" << std::endl;
  amrex::Print() << "  reference temperature = " << ProbParm::T_inf << " K" << std::endl;

  ProbParm::rho_inf = ProbParm::p_inf / (R_d * ProbParm::T_inf);
  amrex::Print() << "  calculated freestream air density = "
                 << ProbParm::rho_inf << " kg/m^3"
                 << std::endl;

  ProbParm::a_inf = std::sqrt(ProbParm::gamma * R_d * ProbParm::T_inf);
  amrex::Print() << "  calculated speed of sound, a = "
                 << ProbParm::a_inf << " m/s"
                 << std::endl;

  amrex::Print() << "  freestream u/a = "
                 << ProbParm::M_inf * std::cos(ProbParm::alpha)
                 << std::endl;
  amrex::Print() << "  freestream v/a = "
                 << ProbParm::M_inf * std::sin(ProbParm::alpha)
                 << std::endl;

}
}
