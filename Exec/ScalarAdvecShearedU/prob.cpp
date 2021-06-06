#include "prob.H"

namespace ProbParm {
AMREX_GPU_DEVICE_MANAGED amrex::Real rho_0 = 0.0; // left density (kg/m^3)
AMREX_GPU_DEVICE_MANAGED amrex::Real T_0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real A_0 = 1.0;
AMREX_GPU_DEVICE_MANAGED amrex::Real z0 = 0.1;    // Surface Roughness
AMREX_GPU_DEVICE_MANAGED amrex::Real zRef = 80.0;  // Reference Height
AMREX_GPU_DEVICE_MANAGED amrex::Real uRef = 8.0;  // Rerefence Wind Speed
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
    pp.query("T_0", ProbParm::T_0);
    pp.query("A_0", ProbParm::A_0);
    pp.query("z0", ProbParm::z0);
    pp.query("zRef", ProbParm::zRef);
    pp.query("uRef", ProbParm::uRef);
  }
}
