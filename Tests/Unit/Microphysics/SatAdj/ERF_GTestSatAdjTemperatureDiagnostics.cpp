#include <array>
#include <memory>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include <ERF_Derive.H>
#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

namespace {

struct ThermoDiagnosticError {
    amrex::Real normalized_error{amrex::Real(0.0)};
    amrex::Real absolute_error{amrex::Real(0.0)};
    amrex::Real actual{amrex::Real(0.0)};
    amrex::Real expected{amrex::Real(0.0)};
    int i{0};
    int j{0};
    int k{0};
    int box_id{-1};
    const char* quantity{"unknown"};
};

enum DiagnosticQuantity {
    DIAG_RHO = 0,
    DIAG_THETA,
    DIAG_QV,
    DIAG_QC,
    DIAG_PRESSURE,
    DIAG_T_EOS,
    DIAG_T_FIXED_P,
    DIAG_NUM
};

const char* diag_name[DIAG_NUM] = {
    "rho",
    "theta",
    "qv",
    "qc",
    "pressure_mbar",
    "T_eos",
    "T_fixed_p"
};

} // namespace

// Motivation: The observed artifact appears in thermodynamic temperature but
// not in theta. SatAdj writes conserved theta and moisture using the old
// density, so this test separates rho, theta, qv, qc, pressure, EOS-projected
// temperature, and fixed-pressure reconstructed temperature against the scalar
// fixed-pressure reference.
TEST(SatAdjTemperatureDiagnostics, PublicFlowMatchesScalarDiagnosticInputs)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);
    amrex::BoxArray ba = make_boxarray(geom.Domain(), amrex::IntVect(AMREX_D_DECL(7, 5, 5)));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
    amrex::MultiFab before(ba, dm, RhoQ2_comp + 1, 0);
    fill_satadj_active_conserved_state(cons);
    amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), 0);

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);

    run_and_sync([&]() {
        std::unique_ptr<amrex::MultiFab> z_phys_nd;
        std::unique_ptr<amrex::MultiFab> detJ_cc;
        satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Advance(amrex::Real(1.0), sc);
        satadj.Copy_Micro_to_State(cons);
    });

    const amrex::MultiFab initial_host = make_host_mirror(before);
    const amrex::MultiFab final_host = make_host_mirror(cons);

    std::array<ThermoDiagnosticError, DIAG_NUM> worst{};
    for (int idx = 0; idx < DIAG_NUM; ++idx) {
        worst[idx].quantity = diag_name[idx];
    }

    for (amrex::MFIter mfi(final_host); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto initial_arr = initial_host.const_array(mfi);
        const auto final_arr = final_host.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const ConservedState reference = scalar_reference_from_initial_conserved(
                        initial_arr(i,j,k,Rho_comp),
                        initial_arr(i,j,k,RhoTheta_comp),
                        initial_arr(i,j,k,RhoQ1_comp),
                        initial_arr(i,j,k,RhoQ2_comp));

                    const CellState initial_state = make_state_from_conserved(
                        initial_arr(i,j,k,Rho_comp),
                        initial_arr(i,j,k,RhoTheta_comp),
                        initial_arr(i,j,k,RhoQ1_comp),
                        initial_arr(i,j,k,RhoQ2_comp));

                    // Public-flow actuals
                    const amrex::Real rho_public = final_arr(i,j,k,Rho_comp);
                    const amrex::Real theta_public = final_arr(i,j,k,RhoTheta_comp) / rho_public;
                    const amrex::Real qv_public = final_arr(i,j,k,RhoQ1_comp) / rho_public;
                    const amrex::Real qc_public = final_arr(i,j,k,RhoQ2_comp) / rho_public;
                    const amrex::Real pressure_public_mbar = getPgivenRTh(final_arr(i,j,k,RhoTheta_comp), qv_public) * amrex::Real(0.01);
                    const amrex::Real T_eos_public = getTgivenRandRTh(rho_public, final_arr(i,j,k,RhoTheta_comp), qv_public);

                    // Scalar-reference expected values (reconstructed from initial)
                    const amrex::Real rho_ref = reference.rho;
                    const amrex::Real theta_ref = reference.rhotheta / reference.rho;
                    const amrex::Real qv_ref = reference.rhoqv / reference.rho;
                    const amrex::Real qc_ref = reference.rhoqc / reference.rho;
                    const amrex::Real pressure_ref_mbar = getPgivenRTh(reference.rhotheta, qv_ref) * amrex::Real(0.01);
                    const amrex::Real T_eos_ref = getTgivenRandRTh(reference.rho, reference.rhotheta, qv_ref);
                    const amrex::Real T_fixed_p_ref = getTgivenPandTh(amrex::Real(100.0) * initial_state.pres_mbar, theta_ref, kRdOcp);

                    // Compute errors and update worst per diagnostic
                    const amrex::Real abs_err_rho = amrex::Math::abs(rho_public - rho_ref);
                    const amrex::Real norm_err_rho = abs_err_rho / scaled_tol(rho_ref, amrex::Real(20.0) * kStateTolFactor);
                    if (norm_err_rho > worst[DIAG_RHO].normalized_error) {
                        worst[DIAG_RHO] = ThermoDiagnosticError{norm_err_rho, abs_err_rho, rho_public, rho_ref, i, j, k, mfi.index(), "rho"};
                    }

                    const amrex::Real abs_err_theta = amrex::Math::abs(theta_public - theta_ref);
                    const amrex::Real norm_err_theta = abs_err_theta / scaled_tol(theta_ref, amrex::Real(20.0) * kThermoTolFactor);
                    if (norm_err_theta > worst[DIAG_THETA].normalized_error) {
                        worst[DIAG_THETA] = ThermoDiagnosticError{norm_err_theta, abs_err_theta, theta_public, theta_ref, i, j, k, mfi.index(), "theta"};
                    }

                    const amrex::Real abs_err_qv = amrex::Math::abs(qv_public - qv_ref);
                    const amrex::Real norm_err_qv = abs_err_qv / mixing_ratio_tol(qv_ref);
                    if (norm_err_qv > worst[DIAG_QV].normalized_error) {
                        worst[DIAG_QV] = ThermoDiagnosticError{norm_err_qv, abs_err_qv, qv_public, qv_ref, i, j, k, mfi.index(), "qv"};
                    }

                    const amrex::Real abs_err_qc = amrex::Math::abs(qc_public - qc_ref);
                    const amrex::Real norm_err_qc = abs_err_qc / mixing_ratio_tol(qc_ref);
                    if (norm_err_qc > worst[DIAG_QC].normalized_error) {
                        worst[DIAG_QC] = ThermoDiagnosticError{norm_err_qc, abs_err_qc, qc_public, qc_ref, i, j, k, mfi.index(), "qc"};
                    }

                    const amrex::Real abs_err_p = amrex::Math::abs(pressure_public_mbar - pressure_ref_mbar);
                    const amrex::Real norm_err_p = abs_err_p / scaled_tol(pressure_ref_mbar, amrex::Real(20.0) * kThermoTolFactor);
                    if (norm_err_p > worst[DIAG_PRESSURE].normalized_error) {
                        worst[DIAG_PRESSURE] = ThermoDiagnosticError{norm_err_p, abs_err_p, pressure_public_mbar, pressure_ref_mbar, i, j, k, mfi.index(), "pressure_mbar"};
                    }

                    const amrex::Real abs_err_teos = amrex::Math::abs(T_eos_public - T_eos_ref);
                    const amrex::Real norm_err_teos = abs_err_teos / derived_temperature_tol(T_eos_ref);
                    if (norm_err_teos > worst[DIAG_T_EOS].normalized_error) {
                        worst[DIAG_T_EOS] = ThermoDiagnosticError{norm_err_teos, abs_err_teos, T_eos_public, T_eos_ref, i, j, k, mfi.index(), "T_eos"};
                    }

                    const amrex::Real tfix_actual = getTgivenPandTh(amrex::Real(100.0) * initial_state.pres_mbar, theta_public, kRdOcp);
                    const amrex::Real abs_err_tfix = amrex::Math::abs(tfix_actual - T_fixed_p_ref);
                    const amrex::Real norm_err_tfix = abs_err_tfix / derived_temperature_tol(T_fixed_p_ref);
                    if (norm_err_tfix > worst[DIAG_T_FIXED_P].normalized_error) {
                        worst[DIAG_T_FIXED_P] = ThermoDiagnosticError{norm_err_tfix, abs_err_tfix, tfix_actual, T_fixed_p_ref, i, j, k, mfi.index(), "T_fixed_p"};
                    }
                }
            }
        }
    }

    for (int d = 0; d < DIAG_NUM; ++d) {
        EXPECT_LE(worst[d].normalized_error, amrex::Real(1.0))
            << "quantity=" << worst[d].quantity
            << " normalized_error=" << worst[d].normalized_error
            << " absolute_error=" << worst[d].absolute_error
            << " actual=" << worst[d].actual
            << " expected=" << worst[d].expected
            << " i=" << worst[d].i
            << " j=" << worst[d].j
            << " k=" << worst[d].k
            << " box_id=" << worst[d].box_id;
    }
}

// Motivation: The observed signal appears in thermodynamic temperature but not
// theta. The production moist-temperature derived path should use the same EOS
// projection from conserved rho, rhoTheta, and rhoQv that the SatAdj diagnostic
// tests use. This directly guards the plotfile/derived-variable contract.
TEST(SatAdjTemperatureDiagnostics, ProductionMoistDerivedTemperatureMatchesEOSProjection)
{
    const amrex::Geometry geom = make_geometry(4, 4, 2);
    const amrex::BoxArray ba = make_boxarray(geom.Domain(), amrex::IntVect(AMREX_D_DECL(4, 4, 2)));
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
    amrex::MultiFab der(ba, dm, 1, 0);
    amrex::MultiFab zcc(ba, dm, 1, 0);
    fill_conserved_state_portable(cons);
    zcc.setVal(amrex::Real(0.0));

    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        derived::erf_dermoisttemp(bx, der[mfi], 0, 1, cons[mfi], zcc[mfi], geom,
                                  amrex::Real(0.0), nullptr, 0);
    }
    satadj_test::sync();

    const amrex::MultiFab der_host = make_host_mirror(der);
    const amrex::MultiFab cons_host = make_host_mirror(cons);

    for (amrex::MFIter mfi(cons_host); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto der_arr = der_host.const_array(mfi);
        const auto dat = cons_host.const_array(mfi);
        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    const amrex::Real rho = dat(i, j, k, Rho_comp);
                    const amrex::Real rhotheta = dat(i, j, k, RhoTheta_comp);
                    const amrex::Real qv = dat(i, j, k, RhoQ1_comp) / rho;
                    const amrex::Real expected = getTgivenRandRTh(rho, rhotheta, qv);
                    const amrex::Real actual = der_arr(i, j, k);
                    EXPECT_NEAR(actual, expected,
                                derived_temperature_tol(expected))
                        << "i=" << i << " j=" << j << " k=" << k
                        << " rho=" << rho
                        << " rhotheta=" << rhotheta
                        << " qv=" << qv;
                }
            }
        }
    }
}
