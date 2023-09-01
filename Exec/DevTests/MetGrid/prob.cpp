#include "prob.H"
#include "prob_common.H"

#include "ERF_Constants.H"

using namespace amrex;

ProbParm parms;

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  amrex::Geometry const&)
{
   amrex::Error("Should never get here for metgrid problems");
}

void
erf_init_rayleigh(Vector<Real>& /*tau*/,
                  Vector<Real>& /*ubar*/,
                  Vector<Real>& /*vbar*/,
                  Vector<Real>& /*wbar*/,
                  Vector<Real>& /*thetabar*/,
                  Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for metgrid problems");
}

void
amrex_probinit(
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
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
   // NOTE: this adds perturbations to what has already been defined; here we add nothing
}

void
init_custom_terrain (const Geometry& /*geom*/,
                           MultiFab& /*z_phys_nd*/,
                     const Real& /*time*/)
{
    amrex::Error("We don't belong in init_custom_terrain!");
}
