#include <memory>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

#ifndef AMREX_USE_GPU
namespace {

CellState make_multifab_cell_state (const int i, const int j, const int k)
{
    const int selector = (i + 2 * j + 3 * k) % 4;

    if (selector == 0) {
        const amrex::Real tabs = amrex::Real(290.0);
        const amrex::Real pres_mbar = amrex::Real(900.0);
        return make_cell_state(tabs, pres_mbar,
                               qsat(tabs, pres_mbar) + amrex::Real(6.0e-4),
                               amrex::Real(3.0e-4));
    }

    if (selector == 1) {
        return make_cell_state(amrex::Real(290.0), amrex::Real(900.0),
                               amrex::Real(4.0e-3), amrex::Real(1.0e-3));
    }

    if (selector == 2) {
        CellState state;
        const bool found = find_evaporation_then_recondensation_state(state);
        AMREX_ALWAYS_ASSERT(found);
        return state;
    }

    const amrex::Real tabs = amrex::Real(288.0);
    const amrex::Real pres_mbar = amrex::Real(850.0);
    return make_cell_state(tabs, pres_mbar,
                           qsat(tabs, pres_mbar) + amrex::Real(4.0e-4),
                           -amrex::Real(1.0e-4));
}

void fill_conserved_state (amrex::MultiFab& cons)
{
    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        auto arr = cons.array(mfi);

        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    const CellState state = make_multifab_cell_state(i, j, k);
                    const ConservedState conserved = make_conserved_state(state);

                    arr(i,j,k,Rho_comp) = conserved.rho;
                    arr(i,j,k,RhoTheta_comp) = conserved.rhotheta;
                    arr(i,j,k,RhoKE_comp) = amrex::Real(0.0);
                    arr(i,j,k,RhoScalar_comp) = amrex::Real(0.0);
                    arr(i,j,k,RhoQ1_comp) = conserved.rhoqv;
                    arr(i,j,k,RhoQ2_comp) = conserved.rhoqc;
                }
            }
        }
    }
}

void run_public_flow (SatAdj& satadj,
                      const SolverChoice& sc,
                      const amrex::Geometry& geom,
                      amrex::MultiFab& cons)
{
    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
    satadj.Copy_State_to_Micro(cons);
    satadj.Advance(amrex::Real(1.0), sc);
    satadj.Copy_Micro_to_State(cons);
}

} // namespace

TEST(SatAdjMultiFab, ShocNoOpKeepsStateUnchanged)
{
    const amrex::Geometry geom = make_geometry(2, 2, 1);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab cons_initial(ba, dm, RhoQ2_comp + 1, 1);

    cons.setVal(amrex::Real(0.0));
    cons_initial.setVal(amrex::Real(0.0));
    fill_conserved_state(cons);
    amrex::MultiFab::Copy(cons_initial, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(true);
    satadj.Define(sc);
    run_public_flow(satadj, sc, geom, cons);

    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        auto before = cons_initial.const_array(mfi);
        auto after = cons.const_array(mfi);

        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    EXPECT_NEAR(after(i,j,k,Rho_comp), before(i,j,k,Rho_comp),
                                scaled_tol(before(i,j,k,Rho_comp), kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoTheta_comp), before(i,j,k,RhoTheta_comp),
                                scaled_tol(before(i,j,k,RhoTheta_comp), kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoQ1_comp), before(i,j,k,RhoQ1_comp),
                                scaled_tol(before(i,j,k,RhoQ1_comp), kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoQ2_comp), before(i,j,k,RhoQ2_comp),
                                scaled_tol(before(i,j,k,RhoQ2_comp), kStateTolFactor));
                }
            }
        }
    }
}

TEST(SatAdjMultiFab, PublicFlowMatchesScalarHelper)
{
    const amrex::Geometry geom = make_geometry(4, 3, 2);
    amrex::BoxArray ba(geom.Domain());
    amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    amrex::MultiFab cons_initial(ba, dm, RhoQ2_comp + 1, 1);

    cons.setVal(amrex::Real(0.0));
    cons_initial.setVal(amrex::Real(0.0));
    fill_conserved_state(cons);
    amrex::MultiFab::Copy(cons_initial, cons, 0, 0, cons.nComp(), cons.nGrowVect());

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);
    satadj.Set_RealWidth(0);
    run_public_flow(satadj, sc, geom, cons);

    for (amrex::MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox();
        auto before = cons_initial.const_array(mfi);
        auto after = cons.const_array(mfi);

        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    CellState state = make_state_from_conserved(before(i,j,k,Rho_comp),
                                                               before(i,j,k,RhoTheta_comp),
                                                               before(i,j,k,RhoQ1_comp),
                                                               before(i,j,k,RhoQ2_comp));
                    adjust(state);

                    const amrex::Real rho = before(i,j,k,Rho_comp);
                    const amrex::Real qt_before = before(i,j,k,RhoQ1_comp) + before(i,j,k,RhoQ2_comp);
                    const amrex::Real qt_after = after(i,j,k,RhoQ1_comp) + after(i,j,k,RhoQ2_comp);

                    EXPECT_NEAR(after(i,j,k,Rho_comp), rho, scaled_tol(rho, kStateTolFactor));
                    EXPECT_NEAR(qt_after, qt_before, scaled_tol(qt_before, amrex::Real(20.0) * kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoTheta_comp), rho * state.theta,
                                scaled_tol(rho * state.theta, amrex::Real(20.0) * kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoQ1_comp), rho * state.qv,
                                scaled_tol(rho * state.qv, amrex::Real(20.0) * kStateTolFactor));
                    EXPECT_NEAR(after(i,j,k,RhoQ2_comp), rho * state.qc,
                                scaled_tol(rho * state.qc, amrex::Real(20.0) * kStateTolFactor));
                }
            }
        }
    }
}
#else
TEST(SatAdjMultiFab, ShocNoOpKeepsStateUnchanged)
{
    GTEST_SKIP() << "Host-side MultiFab verification is not enabled in GPU builds.";
}

TEST(SatAdjMultiFab, PublicFlowMatchesScalarHelper)
{
    GTEST_SKIP() << "Host-side MultiFab verification is not enabled in GPU builds.";
}
#endif