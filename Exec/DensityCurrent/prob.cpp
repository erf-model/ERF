#include "prob.H"
#include "EOS.H"

ProbParm parms;

void
init_isentropic_hse(const amrex::Real& r_sfc, const amrex::Real& theta,
                          amrex::Real* r,           amrex::Real* p,
                    const amrex::Real& dz, const amrex::Real&  prob_lo_z,
                    const int& khi)
{
  // r_sfc / p_0 are the density / pressure at the surface
  r[0] = r_sfc;
  p[0] = p_0 - (0.5*dz) * r[0] * CONST_GRAV;

  int MAX_ITER = 10;
  amrex::Real TOL = 1.e-8;

  int k = 0;
  {
      // We do a Newton iteration to satisfy the EOS (with constant theta) and to discretely satisfy HSE
      bool converged_hse = false;

      // Initial guess
      r[k] = r_sfc;
      p[k] = getPgivenRTh(r[k]*theta);

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p_0 -  (0.5*dz) * r_sfc * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          //amrex::Print() << "PHSE PEOS " << p_hse << " " << p_eos << std::endl;

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);

          Real drho = A / (dpdr + 0.5 * dz * CONST_GRAV);

          //amrex::Print() << "DRHO " << drho << std::endl;

          r[k] = std::max(0.9*r_sfc, std::min(r[k] + drho, 1.1*r_sfc));
          p[k] = getPgivenRTh(r[k]*theta);

          //amrex::Print() << "NEW R P " << r[0] << " " << p[0] << std::endl;

          if (std::abs(drho) < TOL * r_sfc)
          {
              converged_hse = true;
              break;
          }
      }

      if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      if (!converged_hse) amrex::Error("Didn't converge the iterations in init");

      //amrex::Print() << " " << std::endl;
      //amrex::Print() << "P AT SFC VS K=0 " << p_0 << " " << p[0] << std::endl;
      //amrex::Print() << "DPDZ VS RHOG   " << k << " " << (p[0]-p_0)/(0.5*dz) << " " << r[0]*CONST_GRAV << std::endl;
      //amrex::Print() << " " << std::endl;
  }

  // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
  for (int k = 1; k <= khi; k++)
  {
      // To get values at k > 0 we do a Newton iteration to satisfy the EOS (with constant theta) and
      // to discretely satisfy HSE -- here we assume spatial_order = 2 -- we can generalize this later if needed
      bool converged_hse = false;

      r[k] = r[k-1];

      Real p_eos = getPgivenRTh(r[k]*theta);
      Real p_hse;

      for (int iter = 0; iter < MAX_ITER && !converged_hse; iter++)
      {
          p_hse = p[k-1] -  dz * 0.5 * (r[k-1]+r[k]) * CONST_GRAV;
          p_eos = getPgivenRTh(r[k]*theta);

          Real A = p_hse - p_eos;

          Real dpdr = getdPdRgivenConstantTheta(r[k],theta);
          // Gamma * p_0 * std::pow( (R_d * theta / p_0), Gamma) * std::pow(r[k], Gamma-1.0) ;

          Real drho = A / (dpdr + 0.5 * dz * CONST_GRAV);

          r[k] = std::max(0.9*r[k-1], std::min(r[k] + drho, 1.1*r[k-1]));
          p[k] = getPgivenRTh(r[k]*theta);

          if (std::abs(drho) < TOL * r[k-1])
          {
              converged_hse = true;
              //amrex::Print() << " converged " << std::endl;
              break;
          }
      }

      if (!converged_hse) amrex::Print() << "DOING ITERATIONS AT K = " << k << std::endl;
      if (!converged_hse) amrex::Error("Didn't converge the iterations in init");

      //amrex::Print() << " " << std::endl;
      //amrex::Print() << "P AT K vs K-1 " << p[k-1] << " " << p[k] << std::endl;
      //amrex::Print() << "DPDZ VS RHOG " << (p[k]-p[k-1])/dz << " " << 0.5*(r[k]+r[k-1])*CONST_GRAV << std::endl;
      //amrex::Print() << " " << std::endl;
  }
}

void
erf_init_dens_hse(amrex::Real* dens_hse_ptr,
                  amrex::GeometryData const& geomdata,
                  const int /*ng_dens_hse*/)
{
  const Real& prob_lo = geomdata.ProbLo()[2];
  const Real& dz      = geomdata.CellSize()[2];
  const int khi       = geomdata.Domain().bigEnd()[2];

  const amrex::Real& rho_sfc   = p_0 / (R_d*parms.T_0);
  const amrex::Real& Thetabar = parms.T_0;

  amrex::Vector<amrex::Real> r;
  amrex::Vector<amrex::Real> p;

  r.resize(khi+1);
  p.resize(khi+1);

  init_isentropic_hse(rho_sfc,Thetabar,r.data(),p.data(),dz,prob_lo,khi);

  for (int k = 0; k <= khi; k++)
  {
      dens_hse_ptr[k] = r[k];
  }

  dens_hse_ptr[   -1] = dens_hse_ptr[  0];
  dens_hse_ptr[khi+1] = dens_hse_ptr[khi];
}

void
erf_init_prob(
  const amrex::Box& bx,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& x_vel,
  amrex::Array4<amrex::Real> const& y_vel,
  amrex::Array4<amrex::Real> const& z_vel,
  amrex::GeometryData const& geomdata)
{
  const auto prob_lo         = geomdata.ProbLo();
  const auto dx              = geomdata.CellSize();
  const int khi              = geomdata.Domain().bigEnd()[2];

  AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

  // This is what we do at k = 0 -- note we assume p = p_0 and T = T_0 at z=0
  const amrex::Real z0 = (0.5) * dx[2] + prob_lo[2];
  amrex::Real Tbar = parms.T_0 - z0 * CONST_GRAV / parms.C_p;
  amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
  amrex::Real rhobar = pbar / (R_d*Tbar);

  amrex::Real theta = parms.T_0;

  amrex::Vector<amrex::Real> r;
  amrex::Vector<amrex::Real> p;

  r.resize(khi+1);
  p.resize(khi+1);

  const Real& rho_sfc   = p_0 / (R_d*parms.T_0);
  const Real& Thetabar  = parms.T_0;
  const Real& dz        = dx[2];
  const Real& prob_lo_z = prob_lo[2];

  init_isentropic_hse(rho_sfc,theta,r.data(),p.data(),dz,prob_lo_z,khi);

  amrex::ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // Geometry
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

//  amrex::Real Tbar = parms.T_0 - z * CONST_GRAV / parms.C_p;
//  amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, R_d/parms.C_p); // from Straka1993
//  amrex::Real pbar = p_0 * std::pow(Tbar/parms.T_0, parms.C_p/R_d); // isentropic relation, consistent with exner pressure def
//  amrex::Real rhobar = pbar / (R_d*Tbar);

    // Temperature that satisfies the EOS given the hydrostatically balanced (r,p)
    const amrex::Real Tbar_hse = p[0] / (R_d * r[0]);

    amrex::Real L = std::sqrt(
        std::pow((x - parms.x_c)/parms.x_r, 2) +
        std::pow((z - parms.z_c)/parms.z_r, 2)
    );
    amrex::Real dT;
    if (L > 1.0) {
        dT = 0.0;
    }
    else {
        dT = parms.T_pert * (std::cos(PI*L) + 1.0)/2.0;
    }

    // Set the density
    state(i, j, k, Rho_comp) = r[k];

    // Note: dT is a perturbation temperature, which should be converted to a delta theta
    state(i, j, k, RhoTheta_comp) = r[k] * (Tbar_hse+dT)*std::pow(p_0/pbar, R_d/parms.C_p);

    // Using this is a test of whether the initial state is in fact hydrostatically stratified
    //state(i, j, k, RhoTheta_comp) = r[k] * parms.T_0;

    // Set scalar = 0 everywhere
    state(i, j, k, RhoScalar_comp) = 0.0;
  });

  // Construct a box that is on x-faces
  const amrex::Box& xbx = amrex::surroundingNodes(bx,0);
  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    x_vel(i, j, k) = parms.U_0;
  });

  // Construct a box that is on y-faces
  const amrex::Box& ybx = amrex::surroundingNodes(bx,1);
  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    y_vel(i, j, k) = 0.0;
  });

  // Construct a box that is on z-faces
  const amrex::Box& zbx = amrex::surroundingNodes(bx,2);
  // Set the z-velocity
  amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    z_vel(i, j, k) = 0.0;
  });
}

void
erf_init_rayleigh(amrex::Vector<amrex::Real>& /*tau*/,
                  amrex::Vector<amrex::Real>& /*ubar*/,
                  amrex::Vector<amrex::Real>& /*vbar*/,
                  amrex::Vector<amrex::Real>& /*thetabar*/,
                  amrex::GeometryData  const& /*geomdata*/)
{
   amrex::Error("Should never get here for DensityCurrent problem");
}

void
erf_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex_real* /*problo*/,
  const amrex_real* /*probhi*/)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);
}
}
