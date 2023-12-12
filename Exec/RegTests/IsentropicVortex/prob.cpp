#include "prob.H"
#include "EOS.H"
#include "ERF_Constants.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex_real* problo, const amrex_real* probhi)
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
  amrex::Print() << "  reference potential temperature (not used) = " << parms.T_0 << " K" << std::endl;

  parms.rho_0 = parms.p_inf / (R_d * parms.T_inf);
  amrex::Print() << "  calculated freestream air density = "
                 << parms.rho_0 << " kg/m^3"
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

  init_base_parms(parms.rho_0, parms.T_0);
}

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
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real> const& state,
    Array4<Real> const& x_vel,
    Array4<Real> const& y_vel,
    Array4<Real> const& z_vel,
    Array4<Real> const& r_hse,
    Array4<Real> const& p_hse,
    Array4<Real const> const&,
    Array4<Real const> const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    Array4<Real      > const&,
    amrex::GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

  Real xc = parms.xc; Real yc = parms.yc;
  Real R  = parms.R ; Real beta = parms.beta;
  Real sigma = parms.sigma;

  const Real rdOcp = sc.rdOcp;
  //const Real T_0 = parms.T_0;

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real x = prob_lo[0] + (i + 0.5) * dx[0]; // cell center
    const Real y = prob_lo[1] + (j + 0.5) * dx[1]; // cell center

    // Calculate perturbation temperature
    const Real Omg = erf_vortex_Gaussian(x,y,xc,yc,R,beta,sigma);
    const Real deltaT = -(parms.gamma - 1.0)/(2.0*sigma*sigma) * Omg*Omg;

    // Set the perturbation density
    const Real rho_norm = std::pow(1.0 + deltaT, parms.inv_gm1);
    state(i, j, k, Rho_comp) = rho_norm * parms.rho_0 - r_hse(i,j,k);

    // Initial _potential_ temperature
    const Real T = (1.0 + deltaT) * parms.T_inf;
    const Real p = std::pow(rho_norm, Gamma) / Gamma  // isentropic relation
                          * parms.rho_0*parms.a_inf*parms.a_inf;
    const Real rho_theta = parms.rho_0 * rho_norm * (T * std::pow(p_0 / p, rdOcp)); // T --> theta
    state(i, j, k, RhoTheta_comp) = rho_theta - getRhoThetagivenP(p_hse(i,j,k)); // Set the perturbation rho*theta

    // Set scalar = 0 -- unused
    state(i, j, k, RhoScalar_comp) = 0.0;

    if (use_moisture) {
        state(i, j, k, RhoQ1_comp) = 0.0;
        state(i, j, k, RhoQ2_comp) = 0.0;
    }
  });

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

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });
}

void
Problem::init_custom_terrain(
    const Geometry& geom,
    MultiFab& z_phys_nd,
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

#if 0
AMREX_GPU_DEVICE
Real
dhdt(int /*i*/, int /*j*/,
     const GpuArray<Real,AMREX_SPACEDIM> /*dx*/,
     const Real /*time_mt*/, const Real /*delta_t*/)
{
    return 0.;
}
#endif
