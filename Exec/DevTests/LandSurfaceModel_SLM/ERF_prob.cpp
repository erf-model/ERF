#include "ERF_prob.H"
#include <ERF_SLM.H>

using namespace amrex;

std::unique_ptr<ProblemBase>
amrex_probinit(const amrex_real* problo, const amrex_real* probhi)
{
    return std::make_unique<Problem>(problo, probhi);
}

Problem::Problem(const amrex::Real* /*problo*/, const amrex::Real* /*probhi*/)
{
    // Parse params
    ParmParse pp("slm");
    pp.query("SLM_num_ref_inputs", parms.SLM_num_ref_inputs);
    pp.query("SLM_ref_sounding_file", parms.SLM_ref_sounding_file);
    pp.query("SLM_ref_flux_file", parms.SLM_ref_flux_file);
    pp.query("SLM_ref_sst_file", parms.SLM_ref_sst_file);

    init_base_parms(parms.rho_0, parms.T_0);
}

void Problem::read_SLM_inputs()
{
    // Reads the SLM input sounding file assuming the following fields:
    //  t[day], pres0[mb], tabs[K], q[g/kg], uvel[m/s], vvel[m/s]
    sounding = SLM::read_cols(parms.SLM_ref_sounding_file, 1);

    // Reads the SLM input flux file assuming the following fields:
    //  t[day], swdn[W/m2], lwdn[W/m2], swup[W/m2], lwup[W/m2]
    fluxes = SLM::read_cols(parms.SLM_ref_flux_file, 1);

    // Reads the SLM input SST file assuming the following fields:
    // t[day], sst[K], precip[mm/s]
    sst = SLM::read_cols(parms.SLM_ref_sst_file, 1);
}

void
Problem::init_custom_pert(
    const Box& bx,
    const Box& xbx,
    const Box& ybx,
    const Box& zbx,
    Array4<Real const> const& state,
    Array4<Real      > const& state_pert,
    Array4<Real      > const& x_vel_pert,
    Array4<Real      > const& y_vel_pert,
    Array4<Real      > const& z_vel_pert,
    Array4<Real      > const& r_hse,
    Array4<Real      > const& p_hse,
    Array4<Real const> const& z_nd,
    Array4<Real const> const& z_cc,
    GeometryData const& geomdata,
    Array4<Real const> const& /*mf_m*/,
    Array4<Real const> const& /*mf_u*/,
    Array4<Real const> const& /*mf_v*/,
    const SolverChoice& sc)
{
    const bool use_moisture = (sc.moisture_type != MoistureType::None);

    read_SLM_inputs();

    const Real q_ref = sounding[3][0];
    const Real T_ref = sounding[2][0];
    const Real p_ref = sounding[1][0];
    const Real vx_ref = sounding[4][0];
    const Real vy_ref = sounding[5][0];

    ParallelFor(bx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        // Geometry
        const Real* prob_lo = geomdata.ProbLo();
        const Real* prob_hi = geomdata.ProbHi();
        const Real* dx = geomdata.CellSize();
        const Real x = prob_lo[0] + (i + 0.5) * dx[0];
        const Real y = prob_lo[1] + (j + 0.5) * dx[1];
        const Real z = prob_lo[2] + (k + 0.5) * dx[2];

        const Real qv = q_ref / 1000.0; // qv converted to kg/kg
        const Real pres = p_ref * 100.0; // pressure converted from mbar to Pa
        const Real theta = getThgivenTandP(T_ref, pres, R_d / Cp_d);
        const Real rho = getRhogivenThetaPress(theta, pres, R_d / Cp_d, qv);

        // NOTE: these are pertubations from the initial state
        state_pert(i, j, k, Rho_comp) = rho - 1.0;
        state_pert(i, j, k, RhoTheta_comp) = (rho * theta) - (parms.rho_0 * parms.T_0);
        state_pert(i, j, k, RhoScalar_comp) = 0.0;

        state_pert(i, j, k, RhoKE_comp) = 0.0;

        if (use_moisture) {
            state_pert(i, j, k, RhoQ1_comp) = 0.0;
            state_pert(i, j, k, RhoQ2_comp) = 0.0;
        }
    });

    amrex::ParallelFor(xbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto *const prob_hi  = geomdata.ProbHi();
        const auto *const dx       = geomdata.CellSize();
        const Real z = (k + 0.5) * dx[2];
        x_vel_pert(i, j, k) = vx_ref;
    });

    amrex::ParallelFor(ybx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto *const prob_hi  = geomdata.ProbHi();
        const auto *const dx       = geomdata.CellSize();
        const Real z = (k + 0.5) * dx[2];
        y_vel_pert(i, j, k) = vy_ref;
    });

    amrex::ParallelFor(zbx, [=, parms=parms] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        z_vel_pert(i, j, k) = parms.w_0;
    });

}
