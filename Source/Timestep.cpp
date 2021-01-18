#include "EOS.H"
#include "Timestep.H"

// EstDt routines
namespace TimeStep {
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || \
  (__CUDACC_VER_MINOR__ != 2)
AMREX_GPU_DEVICE_MANAGED amrex::Real max_dt =
  std::numeric_limits<amrex::Real>::max();
#else
AMREX_GPU_DEVICE_MANAGED amrex::Real max_dt = 1.e37;
#endif
} // namespace TimeStep

AMREX_GPU_HOST_DEVICE
amrex::Real
pc_estdt_hydro(
  amrex::Box const& bx,
  const amrex::Array4<const amrex::Real>& state,
  const amrex::Real& dx,
  const amrex::Real& dy,
  const amrex::Real& dz) noexcept
{
  amrex::Real dt = TimeStep::max_dt;

  // NOTE THIS CODE IS HACKED TO ONLY USE C, NOT U+C, at the moment

  amrex::Loop(bx, [=, &dt](int i, int j, int k) {
      const amrex::Real rho   = state(i, j, k, Density_comp);
      const amrex::Real theta = state(i, j, k,   Theta_comp);

      amrex::Real pressure = getPgivenRTh(rho,theta);
      amrex::Real c = std::sqrt(Gamma * pressure / rho);
      const amrex::Real dt = dx / c;

#if 0
      const amrex::Real rhoInv = 1.0/rho;
      
      const amrex::Real ux  = u(i, j, k, UMX) * rhoInv;
      const amrex::Real dt1 = dx / (c + amrex::Math::abs(ux));
      dt = amrex::min(dt, dt1);

      const amrex::Real uy = u(i, j, k, UMY) * rhoInv;
      const amrex::Real dt2 = dy / (c + amrex::Math::abs(uy));
      dt = amrex::min(dt, dt2);

      const amrex::Real uz = u(i, j, k, UMZ) * rhoInv;
      const amrex::Real dt3 = dz / (c + amrex::Math::abs(uz));
      dt = amrex::min(dt, dt3););
#endif
  });
  return dt;
}
