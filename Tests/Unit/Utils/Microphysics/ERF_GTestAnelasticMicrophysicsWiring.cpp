#include <AMReX_Gpu.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>

#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include "ERF_BaseStateHSE.H"
#include "ERF_BaseStateRestart.H"
#include "ERF_GTestMicrophysicsCommon.H"
#include "ERF_IndexDefines.H"
#include "ERF_Kessler.H"
#include "ERF_Morrison.H"
#include "ERF_SAMUtils.H"
#include "ERF_WSM6.H"

using namespace amrex;
using namespace microphysics_test;

namespace {

void initialize_cell_fixture (MultiFab& states, MultiFab& base);

struct CellFixture {
    BoxArray boxes{Box(IntVect(0, 0, 0), IntVect(1, 0, 0))};
    DistributionMapping dm{boxes};
    MultiFab states{boxes, dm, RhoQ11_comp + 1, 1};
    MultiFab base{boxes, dm, BaseState::num_comps, 1};

    CellFixture ()
    {
        states.setVal(Real(0.0));
        base.setVal(Real(-777.0));
        initialize_cell_fixture(states, base);
        Gpu::streamSynchronize();
    }
};

void initialize_cell_fixture (MultiFab& states, MultiFab& base)
{
    for (MFIter mfi(states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto state = states.array(mfi);
        const auto base_array = base.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real rho = (i == 0) ? Real(1.10) : Real(0.90);
            const Real theta = (i == 0) ? Real(303.0) : Real(287.0);
            state(i,j,k,Rho_comp) = rho;
            state(i,j,k,RhoTheta_comp) = rho * theta;
            state(i,j,k,RhoQ1_comp) = rho * ((i == 0) ? Real(0.003) : Real(0.007));
            state(i,j,k,RhoQ2_comp) = rho * ((i == 0) ? Real(0.001) : Real(0.002));
            state(i,j,k,RhoQ3_comp) = rho * ((i == 0) ? Real(0.0004) : Real(0.0008));
            state(i,j,k,RhoQ4_comp) = rho * ((i == 0) ? Real(0.0007) : Real(0.0011));
            state(i,j,k,RhoQ5_comp) = rho * ((i == 0) ? Real(0.0002) : Real(0.0005));
            state(i,j,k,RhoQ6_comp) = rho * ((i == 0) ? Real(0.0003) : Real(0.0006));
            state(i,j,k,RhoQ7_comp) = rho * Real(0.0001);
            state(i,j,k,RhoQ8_comp) = rho * Real(0.0002);
            state(i,j,k,RhoQ9_comp) = rho * Real(0.0003);
            state(i,j,k,RhoQ10_comp) = rho * Real(0.0004);
            state(i,j,k,RhoQ11_comp) = rho * Real(0.0005);

            // Deliberately do not make pi0 the EOS value. The production
            // wiring must consume this stored value rather than rebuild it.
            base_array(i,j,k,BaseState::p0_comp) =
                (i == 0) ? Real(90000.0) : Real(82000.0);
            base_array(i,j,k,BaseState::pi0_comp) =
                (i == 0) ? Real(1.07) : Real(0.93);
        });
    }
}

struct WorkingArrays {
    MultiFab rho;
    MultiFab theta;
    MultiFab qv;
    MultiFab qc;
    MultiFab qi;
    MultiFab qn;
    MultiFab qt;
    MultiFab qpr;
    MultiFab qps;
    MultiFab qpg;
    MultiFab qp;
    MultiFab tabs;
    MultiFab pres;

    WorkingArrays (const BoxArray& boxes, const DistributionMapping& dm)
        : rho(boxes, dm, 1, 1), theta(boxes, dm, 1, 1), qv(boxes, dm, 1, 1),
          qc(boxes, dm, 1, 1), qi(boxes, dm, 1, 1), qn(boxes, dm, 1, 1),
          qt(boxes, dm, 1, 1), qpr(boxes, dm, 1, 1), qps(boxes, dm, 1, 1),
          qpg(boxes, dm, 1, 1), qp(boxes, dm, 1, 1), tabs(boxes, dm, 1, 1),
          pres(boxes, dm, 1, 1)
    {
        rho.setVal(Real(-1.0)); theta.setVal(Real(-1.0)); qv.setVal(Real(-1.0));
        qc.setVal(Real(-1.0)); qi.setVal(Real(-1.0)); qn.setVal(Real(-1.0));
        qt.setVal(Real(-1.0)); qpr.setVal(Real(-1.0)); qps.setVal(Real(-1.0));
        qpg.setVal(Real(-1.0)); qp.setVal(Real(-1.0)); tabs.setVal(Real(-1.0));
        pres.setVal(Real(-1.0));
    }
};

void expect_copy_in_outputs (const WorkingArrays& w, const Real pressure_factor,
                             const int moisture_layout)
{
    MultiFab errors(w.rho.boxArray(), w.rho.DistributionMap(), 1, 0);
    for (MFIter mfi(errors, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto error = errors.array(mfi);
        const auto rho = w.rho.const_array(mfi);
        const auto theta = w.theta.const_array(mfi);
        const auto qv = w.qv.const_array(mfi);
        const auto qc = w.qc.const_array(mfi);
        const auto qi = w.qi.const_array(mfi);
        const auto qn = w.qn.const_array(mfi);
        const auto qt = w.qt.const_array(mfi);
        const auto qpr = w.qpr.const_array(mfi);
        const auto qps = w.qps.const_array(mfi);
        const auto qpg = w.qpg.const_array(mfi);
        const auto qp = w.qp.const_array(mfi);
        const auto tabs = w.tabs.const_array(mfi);
        const auto pres = w.pres.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real expected_rho_value = i == 0 ? Real(1.10) : Real(0.90);
            const Real expected_theta_value = i == 0 ? Real(303.0) : Real(287.0);
            const Real expected_qv_value = i == 0 ? Real(0.003) : Real(0.007);
            const Real expected_qc_value = i == 0 ? Real(0.001) : Real(0.002);
            const Real expected_qi_value = i == 0 ? Real(0.0004) : Real(0.0008);
            const Real expected_qpr_value = i == 0 ? Real(0.0007) : Real(0.0011);
            const Real expected_qps_value = i == 0 ? Real(0.0002) : Real(0.0005);
            const Real expected_qpg_value = i == 0 ? Real(0.0003) : Real(0.0006);
            const Real expected_p0_value = i == 0 ? Real(90000.0) : Real(82000.0);
            const Real expected_pi0_value = i == 0 ? Real(1.07) : Real(0.93);
            const Real expected_qn_value = expected_qc_value + expected_qi_value;
            const Real expected_qt_value = expected_qv_value + expected_qn_value;
            const Real expected_qp_value = expected_qpr_value + expected_qps_value + expected_qpg_value;
            const Real expected_tabs_value = expected_theta_value * expected_pi0_value;
            const Real expected_pres_value = expected_p0_value * pressure_factor;
            Real max_error = normalized_error(rho(i,j,k), expected_rho_value, kValueRelTol);
            max_error = amrex::max(max_error, normalized_error(theta(i,j,k), expected_theta_value, kValueRelTol));
            max_error = amrex::max(max_error, normalized_error(qv(i,j,k), expected_qv_value, kValueRelTol));
            if (moisture_layout == 0) {
                max_error = amrex::max(max_error, normalized_error(qc(i,j,k), expected_qc_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qp(i,j,k), expected_qi_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(
                    qt(i,j,k), expected_qv_value + expected_qc_value, kValueRelTol));
            } else if (moisture_layout == 1) {
                max_error = amrex::max(max_error, normalized_error(qc(i,j,k), expected_qc_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qi(i,j,k), expected_qi_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qn(i,j,k), expected_qn_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qt(i,j,k), expected_qt_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qpr(i,j,k), expected_qpr_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qps(i,j,k), expected_qps_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qpg(i,j,k), expected_qpg_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qp(i,j,k), expected_qp_value, kValueRelTol));
            } else {
                max_error = amrex::max(max_error, normalized_error(qc(i,j,k), expected_qc_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qi(i,j,k), expected_qi_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qpr(i,j,k), expected_qpr_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qps(i,j,k), expected_qps_value, kValueRelTol));
                max_error = amrex::max(max_error, normalized_error(qpg(i,j,k), expected_qpg_value, kValueRelTol));
            }
            max_error = amrex::max(max_error, normalized_error(tabs(i,j,k), expected_tabs_value, kValueRelTol));
            max_error = amrex::max(max_error, normalized_error(pres(i,j,k), expected_pres_value, kValueRelTol));
            error(i,j,k) = max_error;
        });
    }
    Gpu::streamSynchronize();
    EXPECT_LE(errors.max(0, 0, false), Real(2.0))
        << " array-level copy-in mismatch (layout=" << moisture_layout << ")";
    EXPECT_NE(w.tabs.min(0, 0, false), w.tabs.max(0, 0, false));
    EXPECT_NE(w.pres.min(0, 0, false), w.pres.max(0, 0, false));
}

void launch_kessler_copy_in (const CellFixture& fixture, WorkingArrays& working)
{
    for (MFIter mfi(fixture.states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = fixture.states.const_array(mfi);
        const auto base = fixture.base.const_array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta = working.theta.array(mfi);
        const auto qv = working.qv.array(mfi);
        const auto qc = working.qc.array(mfi);
        const auto qp = working.qp.array(mfi);
        const auto qt = working.qt.array(mfi);
        const auto tabs = working.tabs.array(mfi);
        const auto pres = working.pres.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            kessler_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qp, qt, tabs, pres,
                true, i, j, k);
        });
    }
    Gpu::streamSynchronize();
}

void launch_sam_copy_in (const CellFixture& fixture, WorkingArrays& working)
{
    for (MFIter mfi(fixture.states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = fixture.states.const_array(mfi);
        const auto base = fixture.base.const_array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta = working.theta.array(mfi);
        const auto qv = working.qv.array(mfi);
        const auto qc = working.qc.array(mfi);
        const auto qi = working.qi.array(mfi);
        const auto qn = working.qn.array(mfi);
        const auto qt = working.qt.array(mfi);
        const auto qpr = working.qpr.array(mfi);
        const auto qps = working.qps.array(mfi);
        const auto qpg = working.qpg.array(mfi);
        const auto qp = working.qp.array(mfi);
        const auto tabs = working.tabs.array(mfi);
        const auto pres = working.pres.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            sam_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qi, qn, qt, qpr, qps, qpg, qp,
                tabs, pres, true, i, j, k);
        });
    }
    Gpu::streamSynchronize();
}

void launch_morrison_copy_in (const CellFixture& fixture, WorkingArrays& working)
{
    for (MFIter mfi(fixture.states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = fixture.states.const_array(mfi);
        const auto base = fixture.base.const_array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta = working.theta.array(mfi);
        const auto qv = working.qv.array(mfi);
        const auto qc = working.qc.array(mfi);
        const auto qi = working.qi.array(mfi);
        const auto qn = working.qn.array(mfi);
        const auto qt = working.qt.array(mfi);
        const auto qpr = working.qpr.array(mfi);
        const auto qps = working.qps.array(mfi);
        const auto qpg = working.qpg.array(mfi);
        const auto qp = working.qp.array(mfi);
        const auto tabs = working.tabs.array(mfi);
        const auto pres = working.pres.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            morrison_copy_state_to_micro_cell(
                states, base, rho, theta, qv, qc, qi, qn, qt, qpr, qps, qpg, qp,
                tabs, pres, true, i, j, k);
        });
    }
    Gpu::streamSynchronize();
}

void launch_wsm6_copy_in (const CellFixture& fixture, WorkingArrays& working)
{
    for (MFIter mfi(fixture.states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = fixture.states.const_array(mfi);
        const auto base = fixture.base.const_array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta = working.theta.array(mfi);
        const auto tabs = working.tabs.array(mfi);
        const auto pres = working.pres.array(mfi);
        const auto qv = working.qv.array(mfi);
        const auto qc = working.qc.array(mfi);
        const auto qi = working.qi.array(mfi);
        const auto qr = working.qpr.array(mfi);
        const auto qs = working.qps.array(mfi);
        const auto qg = working.qpg.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            wsm6_copy_state_to_micro_cell(
                states, base, rho, theta, tabs, pres, qv, qc, qi, qr, qs, qg,
                true, i, j, k);
        });
    }
    Gpu::streamSynchronize();
}

void prepare_wsm6_anelastic_state (const CellFixture& fixture, WorkingArrays& working)
{
    launch_wsm6_copy_in(fixture, working);
    for (MFIter mfi(fixture.states, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto tabs = working.tabs.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            tabs(i,j,k) += (i == 0) ? Real(2.0) : Real(-1.5);
        });
    }
    Gpu::streamSynchronize();
}

void launch_wsm6_anelastic_copy_out (WorkingArrays& working,
                                     MultiFab& state, MultiFab& errors)
{
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = state.array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta = working.theta.array(mfi);
        const auto tabs = working.tabs.const_array(mfi);
        const auto pres = working.pres.const_array(mfi);
        const auto qv = working.qv.const_array(mfi);
        const auto qc = working.qc.const_array(mfi);
        const auto qi = working.qi.const_array(mfi);
        const auto qr = working.qpr.const_array(mfi);
        const auto qs = working.qps.const_array(mfi);
        const auto qg = working.qpg.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            wsm6_copy_micro_to_state_cell(
                states, rho, theta, tabs, pres, qv, qc, qi, qr, qs, qg,
                true, RdoCp, i, j, k);
        });
        const auto theta_read = working.theta.const_array(mfi);
        const auto error = errors.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real expected_theta_value = getThgivenTandP(
                tabs(i,j,k), pres(i,j,k), RdoCp);
            const Real expected_rho_theta = rho(i,j,k) * expected_theta_value;
            const Real wrong_mode_theta =
                getThgivenRandT(rho(i,j,k), tabs(i,j,k), RdoCp,
                                 i == 0 ? Real(0.003) : Real(0.007));
            const Real expected_qv_value = i == 0 ? Real(0.003) : Real(0.007);
            const Real expected_qc_value = i == 0 ? Real(0.001) : Real(0.002);
            const Real expected_qi_value = i == 0 ? Real(0.0004) : Real(0.0008);
            const Real expected_qr_value = i == 0 ? Real(0.0007) : Real(0.0011);
            const Real expected_qs_value = i == 0 ? Real(0.0002) : Real(0.0005);
            const Real expected_qg_value = i == 0 ? Real(0.0003) : Real(0.0006);
            error(i,j,k,0) = normalized_error(
                theta_read(i,j,k), expected_theta_value, kValueRelTol);
            error(i,j,k,1) = normalized_error(
                states(i,j,k,RhoTheta_comp), expected_rho_theta, kValueRelTol);
            error(i,j,k,2) = normalized_error(
                states(i,j,k,RhoQ1_comp), rho(i,j,k) * expected_qv_value, kValueRelTol);
            error(i,j,k,3) = normalized_error(
                states(i,j,k,RhoQ2_comp), rho(i,j,k) * expected_qc_value, kValueRelTol);
            error(i,j,k,4) = normalized_error(
                states(i,j,k,RhoQ3_comp), rho(i,j,k) * expected_qi_value, kValueRelTol);
            error(i,j,k,5) = normalized_error(
                states(i,j,k,RhoQ4_comp), rho(i,j,k) * expected_qr_value, kValueRelTol);
            error(i,j,k,6) = normalized_error(
                states(i,j,k,RhoQ5_comp), rho(i,j,k) * expected_qs_value, kValueRelTol);
            error(i,j,k,7) = normalized_error(
                states(i,j,k,RhoQ6_comp), rho(i,j,k) * expected_qg_value, kValueRelTol);
            const Real mode_scale = amrex::max(Real(1.0),
                                                amrex::Math::abs(expected_theta_value));
            const bool mode_separated = amrex::Math::abs(
                expected_theta_value - wrong_mode_theta) > Real(1.0e-3) * mode_scale;
            // Keep this logical guard separate from numerical errors: a value
            // of one must fail directly rather than being accepted by the
            // numerical error threshold below.
            error(i,j,k,8) = mode_separated ? Real(0.0) : Real(1.0);
        });
    }
    Gpu::streamSynchronize();
}

void launch_wsm6_compressible_copy_out (WorkingArrays& working,
                                        MultiFab& state, MultiFab& theta,
                                        MultiFab& errors)
{
    for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto states = state.array(mfi);
        const auto rho = working.rho.array(mfi);
        const auto theta_array = theta.array(mfi);
        const auto tabs = working.tabs.const_array(mfi);
        const auto pres = working.pres.const_array(mfi);
        const auto qv = working.qv.const_array(mfi);
        const auto qc = working.qc.const_array(mfi);
        const auto qi = working.qi.const_array(mfi);
        const auto qr = working.qpr.const_array(mfi);
        const auto qs = working.qps.const_array(mfi);
        const auto qg = working.qpg.const_array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            wsm6_copy_micro_to_state_cell(
                states, rho, theta_array, tabs, pres, qv, qc, qi, qr, qs, qg,
                false, RdoCp, i, j, k);
        });
        const auto theta_read = theta.const_array(mfi);
        const auto error = errors.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real expected = getThgivenRandT(
                i == 0 ? Real(1.10) : Real(0.90),
                tabs(i,j,k), RdoCp,
                i == 0 ? Real(0.003) : Real(0.007));
            error(i,j,k) = normalized_error(theta_read(i,j,k), expected, kValueRelTol);
        });
    }
    Gpu::streamSynchronize();
}

void populate_hse_base_state (MultiFab& base)
{
    for (MFIter mfi(base, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto pressure = base.array(mfi, BaseState::p0_comp);
        const auto exner = base.array(mfi, BaseState::pi0_comp);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real density = (k == 0) ? Real(1.10) : Real(0.95);
            const Real pressure_value = Real(100000.0) - density * CONST_GRAV *
                                        (Real(k) + Real(0.5)) * Real(12.0);
            erf_base_state::set_pressure_and_exner(
                pressure, exner, i, j, k, pressure_value, RdoCp);
        });
    }
    base.FillBoundary(Periodicity::NonPeriodic());
    Gpu::streamSynchronize();
}

void check_hse_pressure_exner (const MultiFab& base, MultiFab& errors)
{
    for (MFIter mfi(errors, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto base_array = base.const_array(mfi);
        const auto error = errors.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real pressure = base_array(i,j,k,BaseState::p0_comp);
            const Real expected_exner = getExnergivenP(pressure, RdoCp);
            const Real exner_scale = amrex::max(
                Real(1.0), amrex::Math::abs(expected_exner));
            error(i,j,k) = amrex::Math::abs(
                base_array(i,j,k,BaseState::pi0_comp) - expected_exner) / exner_scale;
        });
    }
    Gpu::streamSynchronize();
}

void populate_legacy_restart_base_state (MultiFab& base, const int legacy_ncomp)
{
    for (MFIter mfi(base, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto arr = base.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k,BaseState::r0_comp) = Real(1.0);
            arr(i,j,k,BaseState::p0_comp) = Real(100000.0) -
                Real(500.0) * i - Real(100.0) * k;
            erf_checkpoint::reconstruct_missing_pi0(
                arr, legacy_ncomp, RdoCp, i, j, k);
        });
    }
    base.FillBoundary(Periodicity::NonPeriodic());
    Gpu::streamSynchronize();
}

void check_legacy_restart_ghosts (const MultiFab& base, MultiFab& diagnostics)
{
    for (MFIter mfi(diagnostics, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto arr = base.const_array(mfi);
        const auto diag = diagnostics.array(mfi);
        const Box valid = mfi.validbox();
        const int ilo = valid.smallEnd(0), ihi = valid.bigEnd(0);
        const int jlo = valid.smallEnd(1), jhi = valid.bigEnd(1);
        const int klo = valid.smallEnd(2), khi = valid.bigEnd(2);
        ParallelFor(mfi.growntilebox(1), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const bool in_domain = i >= 0 && i <= 1 && j == 0 && k >= 0 && k <= 1;
            const bool is_valid = i >= ilo && i <= ihi && j >= jlo && j <= jhi &&
                                  k >= klo && k <= khi;
            diag(i,j,k,0) = in_domain ?
                amrex::Math::abs(arr(i,j,k,BaseState::pi0_comp) -
                                 getExnergivenP(arr(i,j,k,BaseState::p0_comp), RdoCp)) : Real(0.0);
            diag(i,j,k,1) = in_domain && !is_valid ? Real(1.0) : Real(0.0);
        });
    }
    Gpu::streamSynchronize();
}

void populate_current_restart_base_state (MultiFab& current)
{
    for (MFIter mfi(current, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto arr = current.array(mfi);
        ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            arr(i,j,k,BaseState::p0_comp) = Real(88000.0) + Real(100.0) * i;
            arr(i,j,k,BaseState::pi0_comp) = Real(0.88) + Real(0.01) * i;
            erf_checkpoint::reconstruct_missing_pi0(
                arr, BaseState::num_comps, RdoCp, i, j, k);
        });
    }
    Gpu::streamSynchronize();
}

} // namespace

// Motivation: Kessler production copy-in must select both base-state
// components, preserve cell-local values, convert Pa to mbar, and populate its
// actual working arrays rather than only a scalar adapter.
TEST(AnelasticMicrophysicsWiring, KesslerPopulatesProductionArrays)
{
    CellFixture fixture;
    WorkingArrays working(fixture.boxes, fixture.dm);
    launch_kessler_copy_in(fixture, working);
    expect_copy_in_outputs(working, Real(0.01), 0);
}

// Motivation: SAM's production copy-in owns its mbar working pressure while
// taking temperature from theta times the stored base Exner.
TEST(AnelasticMicrophysicsWiring, SAMPopulatesProductionArrays)
{
    CellFixture fixture;
    WorkingArrays working(fixture.boxes, fixture.dm);
    launch_sam_copy_in(fixture, working);
    expect_copy_in_outputs(working, Real(0.01), 1);
}

// Motivation: Morrison's production arrays must retain Pa pressure and its
// existing nonnegative moisture clipping while using cell-local stored Exner.
TEST(AnelasticMicrophysicsWiring, MorrisonPopulatesProductionArrays)
{
    CellFixture fixture;
    WorkingArrays working(fixture.boxes, fixture.dm);
    launch_morrison_copy_in(fixture, working);
    expect_copy_in_outputs(working, Real(1.0), 1);
}

// Motivation: WSM6 copy-in must populate its Pa pressure/temperature arrays
// from the local p0/pi0 pair, not from a neighboring cell or an EOS Exner.
TEST(AnelasticMicrophysicsWiring, WSM6CopyInPopulatesProductionArrays)
{
    CellFixture fixture;
    WorkingArrays working(fixture.boxes, fixture.dm);
    launch_wsm6_copy_in(fixture, working);
    expect_copy_in_outputs(working, Real(1.0), 2);
}

// Motivation: WSM6 copy-out must independently preserve the updated theta,
// RhoTheta, and six moisture components. In anelastic mode theta is rebuilt
// from held pressure, while the compressible branch uses its EOS path; both
// modes must remain distinguishable.
TEST(AnelasticMicrophysicsWiring, WSM6CopyOutWritesConservedComponents)
{
    CellFixture fixture;
    WorkingArrays working(fixture.boxes, fixture.dm);
    MultiFab anelastic_state(fixture.boxes, fixture.dm, RhoQ11_comp + 1, 1);
    MultiFab errors(fixture.boxes, fixture.dm, 9, 0);
    anelastic_state.setVal(Real(0.0));
    prepare_wsm6_anelastic_state(fixture, working);
    launch_wsm6_anelastic_copy_out(working, anelastic_state, errors);
    for (int component = 0; component < 8; ++component) {
        EXPECT_LE(errors.max(component, 0, false), Real(2.0))
            << " WSM6 anelastic copy-out component=" << component;
    }
    EXPECT_EQ(errors.max(8, 0, false), Real(0.0))
        << " WSM6 anelastic and compressible theta paths are not separated";

    MultiFab compressible_state(fixture.boxes, fixture.dm, RhoQ11_comp + 1, 1);
    MultiFab compressible_theta(fixture.boxes, fixture.dm, 1, 1);
    MultiFab compressible_errors(fixture.boxes, fixture.dm, 1, 0);
    MultiFab::Copy(compressible_state, fixture.states, 0, 0, RhoQ11_comp + 1, 1);
    launch_wsm6_compressible_copy_out(
        working, compressible_state, compressible_theta, compressible_errors);
    EXPECT_LE(compressible_errors.max(0, 0, false), Real(2.0));
}

// Motivation: the shared production HSE operation must populate a nontrivial
// p0/pi0 column. The test chooses the pressure profile independently and only
// checks the Exner relation after the production setter has filled both arrays.
TEST(AnelasticBaseState, ProductionHSEPopulationUsesSharedOperation)
{
    const BoxArray boxes(Box(IntVect(0, 0, 0), IntVect(0, 0, 1)));
    const DistributionMapping dm(boxes);
    MultiFab base(boxes, dm, BaseState::num_comps, 1);
    base.setVal(Real(-777.0));
    populate_hse_base_state(base);

    MultiFab hse_errors(boxes, dm, 1, 0);
    hse_errors.setVal(Real(0.0));
    check_hse_pressure_exner(base, hse_errors);

    EXPECT_TRUE(base.is_finite());
    EXPECT_GT(base.min(BaseState::p0_comp, 0, false), Real(0.0));
    EXPECT_GT(base.min(BaseState::pi0_comp, 0, false), Real(0.0));
    EXPECT_NE(base.min(BaseState::p0_comp, 0, false),
              base.max(BaseState::p0_comp, 0, false));
    EXPECT_NE(base.min(BaseState::pi0_comp, 0, false),
              base.max(BaseState::pi0_comp, 0, false));
    EXPECT_LE(hse_errors.max(0, 0, false),
              Real(64.0) * std::numeric_limits<Real>::epsilon())
        << " stored pi0 does not match getExnergivenP(p0, RdoCp)";
}

// Motivation: legacy restart must reconstruct missing pi0, propagate valid
// values through internal MultiFab ghosts, and leave current checkpoint pi0
// untouched. Physical boundary ghosts are intentionally excluded because ERF
// applies their boundary conditions after restart.
TEST(LegacyBaseStateRestart, ReconstructsMissingExnerAndFillsInternalGhosts)
{
    BoxList box_list;
    box_list.push_back(Box(IntVect(0, 0, 0), IntVect(0, 0, 1)));
    box_list.push_back(Box(IntVect(1, 0, 0), IntVect(1, 0, 1)));
    const BoxArray boxes(box_list);
    const DistributionMapping dm(boxes);
    MultiFab base(boxes, dm, BaseState::num_comps, 1);
    base.setVal(Real(777.0));
    const int legacy_ncomp = BaseState::pi0_comp;
    populate_legacy_restart_base_state(base, legacy_ncomp);

    MultiFab diagnostics(boxes, dm, 2, 1);
    diagnostics.setVal(Real(0.0));
    check_legacy_restart_ghosts(base, diagnostics);
    EXPECT_LT(diagnostics.max(0, 1, false), Real(256.0) * std::numeric_limits<Real>::epsilon());
    EXPECT_EQ(diagnostics.max(1, 1, false), Real(1.0));

    MultiFab current(boxes, dm, BaseState::num_comps, 0);
    current.setVal(Real(0.0));
    populate_current_restart_base_state(current);
    EXPECT_EQ(current.min(BaseState::pi0_comp, 0, false), Real(0.88));
    EXPECT_EQ(current.max(BaseState::pi0_comp, 0, false), Real(0.89));
}
