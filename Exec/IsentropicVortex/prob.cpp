#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

using namespace amrex;

ProbParm parms;

AMREX_GPU_DEVICE
static
Real
erf_vortex_Gaussian(
  Real x,  Real y,
  Real xc, Real yc,
  Real R,  Real beta,
  Real sigma)
{
  // Evaluate Gaussian function
  const Real r2 = ((x-xc)*(x-xc) + (y-yc)*(y-yc)) / (R*R);
  return beta * std::exp(-r2/(2.*sigma*sigma));
}

void
erf_init_rayleigh(amrex::Vector<Real>& /*tau*/,
                  amrex::Vector<Real>& /*ubar*/,
                  amrex::Vector<Real>& /*vbar*/,
                  amrex::Vector<Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for Isentropic Vortex problem");
}

void
erf_init_dens_hse(MultiFab& rho_hse,
                  std::unique_ptr<MultiFab>&,
                  std::unique_ptr<MultiFab>&,
                  amrex::Geometry const&)
{
    Real rho_inf = parms.rho_inf;
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(rho_hse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(1);
        const Array4<Real> rho_hse_arr = rho_hse[mfi].array();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            rho_hse_arr(i,j,k) = rho_inf;
        });
    }
}

void
init_custom_prob(
  const Box& bx,
  Array4<Real> const& state,
  Array4<Real> const& x_vel,
  Array4<Real> const& y_vel,
  Array4<Real> const& z_vel,
  Array4<Real> const&,
  Array4<Real> const&,
  Array4<Real const> const&,
  Array4<Real const> const&,
  amrex::GeometryData const& geomdata,
  Array4<Real const> const& mf_m,
  Array4<Real const> const& mf_u,
  Array4<Real const> const& mf_v)
{

  Real xc = parms.xc; Real yc = parms.yc;
  Real R  = parms.R ; Real beta = parms.beta;
  Real sigma = parms.sigma;

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
    const Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center

    // Calculate perturbation temperature
    const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);
    const Real deltaT = -(parms.gamma - 1.0)/(2.0*sigma*sigma) * Omg*Omg;

    // Set the density
    const Real rho_norm = std::pow(1.0 + deltaT, parms.inv_gm1);
    state(i, j, k, Rho_comp) = rho_norm * parms.rho_inf;

    // Initial _potential_ temperature
    const Real T = (1.0 + deltaT) * parms.T_inf;
    const Real p = std::pow(rho_norm, Gamma) / Gamma  // isentropic relation
                          * parms.rho_inf*parms.a_inf*parms.a_inf;
    state(i, j, k, RhoTheta_comp) = T * std::pow(p_0 / p, R_d/c_p); // T --> theta
    state(i, j, k, RhoTheta_comp) *= state(i, j, k, Rho_comp);

    // Set scalar = 0 -- unused
    state(i, j, k, RhoScalar_comp) = 0.0;

#ifdef ERF_USE_MOISTURE
    state(i, j, k, RhoQv_comp) = 0.0;
    state(i, j, k, RhoQc_comp) = 0.0;
#endif
  });

  // Construct a box that is on x-faces
  const Box& xbx = amrex::surroundingNodes(bx,0);

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();

      const Real x = prob_lo[0] + i * dx[0]; // face center
      const Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center
      const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

      x_vel(i, j, k) = (parms.M_inf * std::cos(parms.alpha)
                     - (y - parms.yc)/parms.R * Omg) * parms.a_inf;
  });

  // Construct a box that is on y-faces
  const Box& ybx = amrex::surroundingNodes(bx,1);

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      const Real* prob_lo = geomdata.ProbLo();
      const Real* dx = geomdata.CellSize();

      const Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
      const Real y = prob_lo[1] +  j        * dx[1]; // face center
      const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);

      y_vel(i, j, k) = (parms.M_inf * std::sin(parms.alpha)
                     + (x - parms.xc)/parms.R * Omg) * parms.a_inf;
  });

  // Construct a box that is on z-faces
  const Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity

  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });
}

void
init_custom_terrain(const Geometry& geom, MultiFab& z_phys_nd,
                    const Real& /*time*/)
{
    auto dx = geom.CellSizeArray();

    for ( MFIter mfi(z_phys_nd, TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        const Box& gbx = mfi.growntilebox(1);
        Array4<Real> z_arr = z_phys_nd.array(mfi);
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {

            Real z = k * dx[2];

            // Flat terrain with z = 0 at k = 0
            z_arr(i,j,k) = z;
        });
    }
    z_phys_nd.FillBoundary(geom.periodicity());
}

void
amrex_probinit(
  const amrex_real* problo,
  const amrex_real* probhi)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("p_inf", parms.p_inf);
  pp.query("T_inf", parms.T_inf);
  pp.query("M_inf", parms.M_inf);
  pp.query("alpha", parms.alpha);
  pp.query("gamma", parms.gamma);
  pp.query("beta", parms.beta);
  pp.query("sigma", parms.sigma);
  pp.query("R", parms.R);
  pp.query("xc", parms.xc);
  pp.query("yc", parms.yc);

  parms.xc = problo[0] + parms.xc * (probhi[0] - problo[0]);
  parms.yc = problo[1] + parms.yc * (probhi[1] - problo[1]);
  amrex::Print() << "  vortex initialized at ("
                 << parms.xc << ", "
                 << parms.yc << ")"
                 << std::endl;

  parms.inv_gm1 = 1.0 / (parms.gamma - 1.0);

  amrex::Print() << "  reference pressure = " << parms.p_inf << " Pa" << std::endl;
  amrex::Print() << "  reference temperature = " << parms.T_inf << " K" << std::endl;

  parms.rho_inf = parms.p_inf / (R_d * parms.T_inf);
  amrex::Print() << "  calculated freestream air density = "
                 << parms.rho_inf << " kg/m^3"
                 << std::endl;

  parms.a_inf = std::sqrt(parms.gamma * R_d * parms.T_inf);
  amrex::Print() << "  calculated speed of sound, a = "
                 << parms.a_inf << " m/s"
                 << std::endl;

  amrex::Print() << "  freestream u/a = "
                 << parms.M_inf * std::cos(parms.alpha)
                 << std::endl;
  amrex::Print() << "  freestream v/a = "
                 << parms.M_inf * std::sin(parms.alpha)
                 << std::endl;

}

AMREX_GPU_DEVICE
Real
dhdt(int /*i*/, int /*j*/,
     const GpuArray<Real,AMREX_SPACEDIM> /*dx*/,
     const Real /*time_mt*/, const Real /*delta_t*/)
{
    return 0.;
}
