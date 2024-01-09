#include "prob.H"
#include "EOS.H"

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(
    const amrex_real* /*problo*/,
    const amrex_real* /*probhi*/)
{
    return std::make_unique<Problem>();
}

Problem::Problem()
{
  // Parse params
  ParmParse pp("prob");
  pp.query("T_0", parms.T_0);
  pp.query("U_0", parms.U_0);
  pp.query("x_c", parms.x_c);
  pp.query("z_c", parms.z_c);
  pp.query("x_r", parms.x_r);
  pp.query("z_r", parms.z_r);
  pp.query("T_pert", parms.T_pert);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real      > const& state,
    Array4<Real      > const& x_vel,
    Array4<Real      > const& y_vel,
    Array4<Real      > const& z_vel,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const int khi = geomdata.Domain().bigEnd()[2];

    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    AMREX_ALWAYS_ASSERT(bx.length()[2] == khi+1);

    const Real l_x_r = parms.x_r;
    //const Real l_x_r = parms.x_r * mf_u(0,0,0); //used to validate constant msf
    const Real l_z_r = parms.z_r;
    const Real l_x_c = parms.x_c;
    const Real l_z_c = parms.z_c;
    const Real l_Tpt = parms.T_pert;
    const Real rdOcp = sc.rdOcp;

    if (z_cc) {
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo = geomdata.ProbLo();
        const auto dx      = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = z_cc(i,j,k);

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2));
        if (L <= 1.0)
        {
            Real dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
            Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: dT is a perturbation in temperature, theta_perturbed is base state + perturbation
            Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);

            state(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
        }

        // Set scalar = 0 everywhere
        state(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state(i, j, k, RhoQ1_comp) = 0.0;
            state(i, j, k, RhoQ2_comp) = 0.0;
        }
      });
  } else {
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
      {
        // Geometry (note we must include these here to get the data on device)
        const auto prob_lo = geomdata.ProbLo();
        const auto dx      = geomdata.CellSize();

        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        Real L = std::sqrt(
            std::pow((x - l_x_c)/l_x_r, 2) +
            std::pow((z - l_z_c)/l_z_r, 2));
        if (L <= 1.0)
        {
            Real dT = l_Tpt * (std::cos(PI*L) + 1.0)/2.0;
            Real Tbar_hse = p_hse(i,j,k) / (R_d * r_hse(i,j,k));

            // Note: dT is a perturbation in temperature, theta_perturbed is base state + perturbation
            Real theta_perturbed = (Tbar_hse+dT)*std::pow(p_0/p_hse(i,j,k), rdOcp);

            state(i, j, k, Rho_comp) = getRhoThetagivenP(p_hse(i,j,k)) / theta_perturbed - r_hse(i,j,k);
        }

        // Set scalar = 0 everywhere
        state(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state(i, j, k, RhoQ1_comp) = 0.0;
            state(i, j, k, RhoQ2_comp) = 0.0;
        }
      });
  }

  const Real u0 = parms.U_0;

  // Set the x-velocity
  amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      x_vel(i, j, k) = u0;
  });

  // Set the y-velocity
  amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      y_vel(i, j, k) = 0.0;
  });

  // Set the z-velocity
  amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
  {
      z_vel(i, j, k) = 0.0;
  });

  amrex::Gpu::streamSynchronize();
}
