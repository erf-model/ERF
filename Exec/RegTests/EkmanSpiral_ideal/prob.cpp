#include "prob.H"
#include "prob_common.H"

#include "ERF_Constants.H"

using namespace amrex;

ProbParm parms;

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  Geometry const&)
{
    Real R0 = parms.rho_0;
    for ( MFIter mfi(rho_hse, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
       Array4<Real> rho_arr = rho_hse.array(mfi);
       const Box& gbx = mfi.growntilebox(1);
       ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
           rho_arr(i,j,k) = R0;
       });
    }
}

void
erf_init_rayleigh(Vector<Real>& /*tau*/,
                  Vector<Real>& /*ubar*/,
                  Vector<Real>& /*vbar*/,
                  Vector<Real>& /*wbar*/,
                  Vector<Real>& /*thetabar*/,
                  Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for WPS tests problem");
}

void
init_custom_prob(
    const Box& /*bx*/,
    const Box& /*xbx*/,
    const Box& /*ybx*/,
    const Box& /*zbx*/,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
#if defined(ERF_USE_MOISTURE)
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
#elif defined(ERF_USE_WARM_NO_PRECIP)
    Array4<Real      > const&,
    Array4<Real      > const&,
#endif
    GeometryData const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
    Array4<Real const> const&,
    const SolverChoice&)
{
  amrex::Print() << "Dummy function..Needed for linking" << std::endl;
}

void
init_custom_terrain (const Geometry& /*geom*/,
                           MultiFab& z_phys_nd,
                     const Real& /*time*/)
{
    // Number of ghost cells
    int ngrow = z_phys_nd.nGrow();

    // Bottom of domain
    int k0 = 0;

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Grown box with no z range
        Box xybx = mfi.growntilebox(ngrow);
        xybx.setRange(2,0);

        Array4<Real> const& z_arr = z_phys_nd.array(mfi);

        ParallelFor(xybx, [=] AMREX_GPU_DEVICE (int i, int j, int) {

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k0) = 0.0;
        });
    }
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  ParmParse pp("prob");
  pp.query("rho_0", parms.rho_0);
  pp.query("T_0", parms.Theta_0);
}
