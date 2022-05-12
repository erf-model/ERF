#include "prob.H"
#include "prob_common.H"

#include "IndexDefines.H"
#include "AMReX_ParmParse.H"
#include "AMReX_MultiFab.H"

ProbParm parms;

using namespace amrex;

AMREX_GPU_DEVICE
amrex::Real
erf_vortex_Gaussian(
  amrex::Real x,
  amrex::Real y,
  const ProbParm& parms)
{
  // Evaluate Gaussian function
  const amrex::Real r2 = ((x-parms.xc)*(x-parms.xc) + (y-parms.yc)*(y-parms.yc))
                         / (parms.R*parms.R);
  return parms.beta * std::exp(-r2/(2*parms.sigma*parms.sigma));
}

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
                  amrex::Geometry      const& /*geom*/)
{
   amrex::Error("Should never get here for Isentropic Vortex problem");
}

#ifdef ERF_USE_TERRAIN
void
erf_init_dens_hse(MultiFab& rho_hse,
                  MultiFab const& z_phys_nd,
                  MultiFab const& z_phys_cc,
                  Geometry const& geom)
{
    amrex::Error("This problem not set up to use with terrain");
}

#else
void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::Geometry const& geom,
                  const int ng_dens_hse)
{
  const int khi = geom.Domain().bigEnd()[2];
  for (int k = -ng_dens_hse; k <= khi+ng_dens_hse; k++)
  {
      dens_hse_ptr[k] = parms.rho_inf;
  }
}
#endif

void
init_custom_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
#ifdef ERF_USE_TERRAIN
  amrex::Array4<amrex::Real> const& r_hse,
  amrex::Array4<amrex::Real> const& p_hse,
  amrex::Array4<amrex::Real const> const& z_nd,
  amrex::Array4<amrex::Real const> const& z_cc,
#endif
  amrex::GeometryData const& geomdata)
{
  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
  //  const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center

    // Calculate perturbation temperature
    const amrex::Real Omg = erf_vortex_Gaussian(x,y,parms);
    const amrex::Real deltaT = -(parms.gamma - 1.0)/(2.0*parms.sigma*parms.sigma) * Omg*Omg;
    // Set the density
    const amrex::Real rho_norm = std::pow(1.0 + deltaT, parms.inv_gm1);
    state(i, j, k, Rho_comp) = rho_norm * parms.rho_inf;

    // Initial _potential_ temperature
    const amrex::Real T = (1.0 + deltaT) * parms.T_inf;
    const amrex::Real p = std::pow(rho_norm, Gamma) / Gamma  // isentropic relation
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
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Note that this is called on a box of x-faces
    const amrex::Real* prob_lo = geomdata.ProbLo();
    //  const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + i * dx[0]; // face center
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center
    const amrex::Real Omg = erf_vortex_Gaussian(x,y,parms);

    // Set the x-velocity
    x_vel(i, j, k) = (parms.M_inf * std::cos(parms.alpha)
                   - (y - parms.yc)/parms.R * Omg) * parms.a_inf;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Note that this is called on a box of y-faces
    const amrex::Real* prob_lo = geomdata.ProbLo();
    //  const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
    const amrex::Real y = prob_lo[1] + j * dx[1]; // face center
    const amrex::Real Omg = erf_vortex_Gaussian(x,y,parms);

    // Set the y-velocity
    y_vel(i, j, k) = (parms.M_inf * std::sin(parms.alpha)
                   + (x - parms.xc)/parms.R * Omg) * parms.a_inf;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Note that this is called on a box of z-faces

    // Set the z-velocity
    z_vel(i, j, k) = 0.0;
  });
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
