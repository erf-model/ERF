#include <memory>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_GTestSatAdjCommon.H"

// These tests exercise the AMReX-portable SatAdj integration path. Setup and
// error computation use ParallelFor so the same test code runs in CPU and GPU
// builds. Host-side GTest assertions inspect reduced errors after
// synchronization.

using namespace satadj_test;

namespace {

// Exercise the public SatAdj flow used by ERF: initialize microphysics
// storage, copy conserved state into microphysics variables, advance, and copy
// back.
void run_public_flow (SatAdj& satadj,
                      const SolverChoice& sc,
                      const amrex::Geometry& geom,
                      amrex::MultiFab& cons)
{
    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    run_and_sync([&]() {
        satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Advance(amrex::Real(1.0), sc);
        satadj.Copy_Micro_to_State(cons);
    });
}

std::vector<CellState> make_kernel_cases ()
{
    CellState evap_then_recond;
    const bool found = find_evaporation_then_recondensation_state(evap_then_recond);
    AMREX_ALWAYS_ASSERT(found);

    return {
        make_cell_state(amrex::Real(290.0), amrex::Real(900.0),
                        qsat(amrex::Real(290.0), amrex::Real(900.0)) + amrex::Real(8.0e-4),
                        amrex::Real(4.0e-4)),
        make_cell_state(amrex::Real(290.0), amrex::Real(900.0),
                        amrex::Real(4.0e-3), amrex::Real(1.0e-3)),
        evap_then_recond,
        make_cell_state(amrex::Real(288.0), amrex::Real(850.0),
                        qsat(amrex::Real(288.0), amrex::Real(850.0)) + amrex::Real(4.0e-4),
                        -amrex::Real(1.0e-4)),
        make_cell_state(amrex::Real(290.0), amrex::Real(900.0),
                        qsat(amrex::Real(290.0), amrex::Real(900.0)) + amrex::Real(1.0e-6),
                        amrex::Real(0.0)),
        make_cell_state(amrex::Real(295.0), amrex::Real(950.0),
                        amrex::Real(7.0e-3), amrex::Real(0.0)),
        make_cell_state(amrex::Real(285.0), amrex::Real(880.0),
                        qsat(amrex::Real(285.0), amrex::Real(880.0)) + amrex::Real(2.0e-4),
                        amrex::Real(5.0e-5)),
        make_cell_state(amrex::Real(287.0), amrex::Real(870.0),
                        amrex::Real(5.0e-3), amrex::Real(4.0e-4))
    };
}

// Keep the device lambda in a namespace-scope helper. CUDA rejects extended
// device lambdas enclosed by GTest's private TestBody() member function.
void launch_adjust_cell_kernel (const int ncases,
                                const amrex::Real* tabs_in_ptr,
                                const amrex::Real* pres_in_ptr,
                                const amrex::Real* theta_in_ptr,
                                const amrex::Real* qv_in_ptr,
                                const amrex::Real* qc_in_ptr,
                                amrex::Real* tabs_out_ptr,
                                amrex::Real* theta_out_ptr,
                                amrex::Real* qv_out_ptr,
                                amrex::Real* qc_out_ptr)
{
    amrex::ParallelFor(ncases, [=] AMREX_GPU_DEVICE (int idx) noexcept {
        amrex::Real tabs = tabs_in_ptr[idx];
        const amrex::Real pres_mbar = pres_in_ptr[idx];
        amrex::Real theta = theta_in_ptr[idx];
        amrex::Real qv = qv_in_ptr[idx];
        amrex::Real qc = qc_in_ptr[idx];

        SatAdj::AdjustSatAdjCell(kFacCond, kRdOcp, tabs, pres_mbar, theta, qv, qc);

        tabs_out_ptr[idx] = tabs;
        theta_out_ptr[idx] = theta;
        qv_out_ptr[idx] = qv;
        qc_out_ptr[idx] = qc;
    });

    satadj_test::sync();
}

} // namespace

// Motivation: When SHOC owns condensation, SatAdj must be a no-op. This
// protects against double-adjusting vapor and cloud water in SHOC-enabled
// runs.
TEST(SatAdjMultiFab, ShocNoOpKeepsStateUnchangedPortable)
{
    const amrex::Geometry geom = make_geometry(2, 2, 1);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab cons_initial(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab err(ba, dm, NumStateDiffComps, 0);

    cons.setVal(amrex::Real(0.0));
    cons_initial.setVal(amrex::Real(0.0));
    err.setVal(amrex::Real(0.0));
    fill_conserved_state_portable(cons, false);
    amrex::MultiFab::Copy(cons_initial, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(true);
    satadj.Define(sc);
    run_public_flow(satadj, sc, geom, cons);
    compute_state_difference(cons, cons_initial, err);

    EXPECT_LE(err.max(StateDiffRho),
              scaled_tol(std::max(cons.max(Rho_comp), cons_initial.max(Rho_comp)), kStateTolFactor));
    EXPECT_LE(err.max(StateDiffRhoTheta),
              scaled_tol(std::max(cons.max(RhoTheta_comp), cons_initial.max(RhoTheta_comp)), kStateTolFactor));
    EXPECT_LE(err.max(StateDiffQ1),
              scaled_tol(std::max(cons.max(RhoQ1_comp), cons_initial.max(RhoQ1_comp)), kStateTolFactor));
    EXPECT_LE(err.max(StateDiffQ2),
              scaled_tol(std::max(cons.max(RhoQ2_comp), cons_initial.max(RhoQ2_comp)), kStateTolFactor));
}

// Motivation: The scalar tests verify cell physics; this verifies the public
// MultiFab copy-in, advance, and copy-out path preserves the same conserved
// state contract and scalar-reference behavior.
TEST(SatAdjMultiFab, PublicFlowPreservesSatAdjInvariantsPortable)
{
    const amrex::Geometry geom = make_geometry(4, 3, 2);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab cons_initial(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab err(ba, dm, NumInvariantErrComps, 0);

    cons.setVal(amrex::Real(0.0));
    cons_initial.setVal(amrex::Real(0.0));
    err.setVal(amrex::Real(0.0));
    fill_conserved_state_portable(cons, true);
    amrex::MultiFab::Copy(cons_initial, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);
    run_public_flow(satadj, sc, geom, cons);
    compute_satadj_invariant_errors(cons_initial, cons, err);

    const amrex::Real normalized_tol =
#ifdef AMREX_USE_FLOAT
        amrex::Real(4.0);
#else
        amrex::Real(2.0);
#endif

    EXPECT_LE(err.max(InvariantErrRho), normalized_tol);
    EXPECT_LE(err.max(InvariantErrQt), normalized_tol);
    EXPECT_LE(err.max(InvariantErrQcNegative), normalized_tol);
    EXPECT_LE(err.max(InvariantErrRhoThetaRef), normalized_tol);
    EXPECT_LE(err.max(InvariantErrQ1Ref), normalized_tol);
    EXPECT_LE(err.max(InvariantErrQ2Ref), normalized_tol);
}

// Motivation: This is the end-to-end public-flow SatAdj conservation test.
// It exercises Init -> Copy_State_to_Micro -> Advance -> Copy_Micro_to_State
// on physical qv/qc states and checks the quantities SatAdj should conserve:
//   - rho
//   - rho * (qv + qc)
//   - rho * [T + (L/cp) qv] at the pre-adjustment fixed pressure.
//
// SatAdj has no precipitating species, so precipitation-reaching-surface
// conservation is intentionally out of scope here. The reported rho/qt
// residuals are expected to be roundoff-sized because the public path converts
// rho*q to q and back around the cell adjustment.
TEST(SatAdjMultiFab, PublicFlowConservesWaterAndLatentEnergyPortable)
{
    const amrex::Geometry geom = make_geometry(4, 3, 2);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab cons_initial(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab err(ba, dm, NumConservationErrComps, 0);

    cons.setVal(amrex::Real(0.0));
    cons_initial.setVal(amrex::Real(0.0));
    err.setVal(amrex::Real(0.0));

    fill_conserved_state_conservation_portable(cons);
    amrex::MultiFab::Copy(cons_initial, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);

    run_public_flow(satadj, sc, geom, cons);
    compute_satadj_conservation_errors(cons_initial, cons, err);

    const amrex::Real rho_err = err.max(ConservationErrRho);
    const amrex::Real qt_err = err.max(ConservationErrQt);
    const amrex::Real latent_energy_err = err.max(ConservationErrLatentEnergy);
    const amrex::Real initial_qc_coverage = err.max(ConservationErrInitialQcPositive);
    const amrex::Real final_qc_coverage = err.max(ConservationErrFinalQcPositive);
    const amrex::Real final_qc_zero_coverage = err.max(ConservationErrFinalQcZero);

    EXPECT_LE(rho_err, kNormalizedWaterConservationTol)
        << "Observed normalized max rho conservation error: " << rho_err;
    EXPECT_LE(qt_err, kNormalizedWaterConservationTol)
        << "Observed normalized max total-water conservation error: " << qt_err;
    EXPECT_LE(latent_energy_err, kNormalizedLatentEnergyConservationTol)
        << "Observed normalized max latent-energy conservation error: " << latent_energy_err;

    EXPECT_GT(initial_qc_coverage, amrex::Real(0.5));
    EXPECT_GT(final_qc_coverage, amrex::Real(0.5));
    EXPECT_GT(final_qc_zero_coverage, amrex::Real(0.5));
}

// Motivation: Production SatAdj calls AdjustSatAdjCell from ParallelFor. This
// test checks the kernel-compiled path without duplicating CPU/GPU test logic.
TEST(SatAdjKernel, AdjustCellMatchesHostReferencePortable)
{
    // Launch the scalar SatAdj helper through ParallelFor and compare with the
    // host helper result. This verifies the device-compiled path without
    // maintaining a separate GPU-only test implementation.
    const std::vector<CellState> initial_cases = make_kernel_cases();
    std::vector<CellState> host_reference = initial_cases;

    std::vector<amrex::Real> tabs_in(initial_cases.size());
    std::vector<amrex::Real> pres_in(initial_cases.size());
    std::vector<amrex::Real> theta_in(initial_cases.size());
    std::vector<amrex::Real> qv_in(initial_cases.size());
    std::vector<amrex::Real> qc_in(initial_cases.size());

    for (int idx = 0; idx < static_cast<int>(initial_cases.size()); ++idx) {
        tabs_in[idx] = initial_cases[idx].tabs;
        pres_in[idx] = initial_cases[idx].pres_mbar;
        theta_in[idx] = initial_cases[idx].theta;
        qv_in[idx] = initial_cases[idx].qv;
        qc_in[idx] = initial_cases[idx].qc;
        adjust(host_reference[idx]);
    }

    amrex::Gpu::DeviceVector<amrex::Real> d_tabs_in(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_pres_in(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_theta_in(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_qv_in(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_qc_in(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_tabs_out(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_theta_out(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_qv_out(initial_cases.size());
    amrex::Gpu::DeviceVector<amrex::Real> d_qc_out(initial_cases.size());

    amrex::Gpu::copy(amrex::Gpu::hostToDevice, tabs_in.begin(), tabs_in.end(), d_tabs_in.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, pres_in.begin(), pres_in.end(), d_pres_in.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, theta_in.begin(), theta_in.end(), d_theta_in.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, qv_in.begin(), qv_in.end(), d_qv_in.begin());
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, qc_in.begin(), qc_in.end(), d_qc_in.begin());

    const amrex::Real* tabs_in_ptr = d_tabs_in.data();
    const amrex::Real* pres_in_ptr = d_pres_in.data();
    const amrex::Real* theta_in_ptr = d_theta_in.data();
    const amrex::Real* qv_in_ptr = d_qv_in.data();
    const amrex::Real* qc_in_ptr = d_qc_in.data();
    amrex::Real* tabs_out_ptr = d_tabs_out.data();
    amrex::Real* theta_out_ptr = d_theta_out.data();
    amrex::Real* qv_out_ptr = d_qv_out.data();
    amrex::Real* qc_out_ptr = d_qc_out.data();

    launch_adjust_cell_kernel(static_cast<int>(initial_cases.size()),
                              tabs_in_ptr, pres_in_ptr, theta_in_ptr, qv_in_ptr, qc_in_ptr,
                              tabs_out_ptr, theta_out_ptr, qv_out_ptr, qc_out_ptr);

    std::vector<amrex::Real> tabs_out(initial_cases.size());
    std::vector<amrex::Real> theta_out(initial_cases.size());
    std::vector<amrex::Real> qv_out(initial_cases.size());
    std::vector<amrex::Real> qc_out(initial_cases.size());

    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_tabs_out.begin(), d_tabs_out.end(), tabs_out.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_theta_out.begin(), d_theta_out.end(), theta_out.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_qv_out.begin(), d_qv_out.end(), qv_out.begin());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost, d_qc_out.begin(), d_qc_out.end(), qc_out.begin());

    const amrex::Real kernel_tol_factor =
#ifdef AMREX_USE_FLOAT
        amrex::Real(25.0) * kStateTolFactor;
#else
        amrex::Real(10.0) * kStateTolFactor;
#endif

    for (int idx = 0; idx < static_cast<int>(host_reference.size()); ++idx) {
        EXPECT_NEAR(tabs_out[idx], host_reference[idx].tabs,
                    scaled_tol(host_reference[idx].tabs, kernel_tol_factor));
        EXPECT_NEAR(theta_out[idx], host_reference[idx].theta,
                    scaled_tol(host_reference[idx].theta, kernel_tol_factor));
        EXPECT_NEAR(qv_out[idx], host_reference[idx].qv,
                    scaled_tol(host_reference[idx].qv, kernel_tol_factor));
        EXPECT_NEAR(qc_out[idx], host_reference[idx].qc,
                    scaled_tol(host_reference[idx].qc, kernel_tol_factor));
    }
}
