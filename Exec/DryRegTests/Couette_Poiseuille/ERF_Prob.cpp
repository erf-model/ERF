#include "ERF_Prob.H"

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
  pp.query("rho_0", parms.rho_0);
  pp.query("u_0", parms.u_0);
  pp.query("v_0", parms.v_0);
  pp.query("w_0", parms.w_0);
  pp.query("T_0", parms.T_0);

  pp.get("prob_type", parms.prob_type);

  pp.query("u_0_pert_mag", parms.u_0_pert_mag);
  pp.query("v_0_pert_mag", parms.v_0_pert_mag);
  pp.query("w_0_pert_mag", parms.w_0_pert_mag);

  init_base_parms(parms.rho_0, parms.T_0);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& /*state*/,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& /*r_hse*/,
    Array4<Real      > const& /*p_hse*/,
    Array4<Real const> const& /*z_nd*/,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    AMREX_ALWAYS_ASSERT (parms.prob_type == 1 ||
                         parms.prob_type == 10 || parms.prob_type == 11 ||
                         parms.prob_type == 20 || parms.prob_type == 21);

    // Couette flow
    if (parms.prob_type == 1) {

        ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto *const prob_hi  = geomdata.ProbHi();
            const auto *const dx       = geomdata.CellSize();
            const Real z = (k + 0.5) * dx[2];
            x_vel_pert(i, j, k) = parms_d.u_0 * z / prob_hi[2];
        });

        ParallelFor(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const auto *const prob_hi  = geomdata.ProbHi();
            const auto *const dx       = geomdata.CellSize();
            const Real z = (k + 0.5) * dx[2];
            y_vel_pert(i, j, k) = parms_d.v_0 * z / prob_hi[2];
        });

        ParallelFor(zbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            z_vel_pert(i, j, k) = parms_d.w_0;
        });

    // Poiseuille flow
    } else if (parms.prob_type == 10 || parms.prob_type == 11) {

        ParallelFor(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx      = geomdata.CellSize();
            const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

            // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
            if (parms_d.prob_type == 10) {
                x_vel_pert(i, j, k) = 1.0 - z_h * z_h;
            } else {
                x_vel_pert(i, j, k) = 0.0;
            }
        });

        ParallelFor(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* dx      = geomdata.CellSize();
            const Real z_h = prob_lo[2] + (k + 0.5) *  dx[2];

            // Set the x-velocity to be a parabolic profile with max 1 at z = 0 and 0 at z = +/-1
            if (parms_d.prob_type == 11) {
               y_vel_pert(i, j, k) = 1.0 - z_h * z_h;
            } else {
               y_vel_pert(i, j, k) = 0.0;
            }
        });

        ParallelFor(zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            z_vel_pert(i, j, k) = 0.0;
        });

    // plane channel flow initialization
    } else if (parms.prob_type == 20 || parms.prob_type == 21) {

        ParallelForRNG(xbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
        {
            const Real* prob_lo = geomdata.ProbLo();
            const Real* prob_hi = geomdata.ProbHi();
            const Real* dx      = geomdata.CellSize();

            // Normalized values between -1 and 1
            Real y_h = (parms_d.prob_type == 20) ? 2.0 * (j + 0.5) * dx[1]          / (prob_hi[1] - prob_lo[1]) - 1.0
                                                 : 2.0 * (z_cc(i,j,k) - prob_lo[2]) / (prob_hi[2] - prob_lo[2]) - 1.0;

            x_vel_pert(i, j, k) = parms_d.u_0 * (1.0 - y_h * y_h);

            if (parms_d.u_0_pert_mag != 0.0) {
                Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                x_vel_pert(i, j, k) += (rand_double*2.0 - 1.0)*parms_d.u_0_pert_mag;
            }
        });

        ParallelForRNG(ybx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
        {
            if (parms_d.v_0_pert_mag != 0.0) {
                Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                y_vel_pert(i, j, k) = (rand_double*2.0 - 1.0)*parms_d.v_0_pert_mag;
            } else {
                y_vel_pert(i, j, k) = 0.0;
            }
        });

        ParallelForRNG(zbx, [=, parms_d=parms] AMREX_GPU_DEVICE(int i, int j, int k, const amrex::RandomEngine& engine) noexcept
        {
            if (parms_d.w_0_pert_mag != 0.0) {
                Real rand_double = amrex::Random(engine); // Between 0.0 and 1.0
                z_vel_pert(i, j, k) = (rand_double*2.0 - 1.0)*parms_d.w_0_pert_mag;
            } else {
                z_vel_pert(i, j, k) = 0.0;
            }
        });
    } // prob_type
}
