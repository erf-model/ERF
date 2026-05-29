#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

namespace {

std::string component_name (const int component)
{
    switch (component) {
    case Rho_comp: return "Rho_comp";
    case RhoTheta_comp: return "RhoTheta_comp";
    case RhoQ1_comp: return "RhoQ1_comp";
    case RhoQ2_comp: return "RhoQ2_comp";
    case -1: return "EOSProjectedTemperature";
    default: return "unknown";
    }
}

std::string max_size_string (const amrex::IntVect& max_size)
{
    std::ostringstream os;
    os << "(" << max_size[0] << "," << max_size[1] << "," << max_size[2] << ")";
    return os.str();
}

long num_valid_cells (const amrex::Box& domain)
{
    return static_cast<long>(domain.numPts());
}

int wrap_index (int idx,
                const int lo,
                const int hi)
{
    const int extent = hi - lo + 1;
    while (idx < lo) { idx += extent; }
    while (idx > hi) { idx -= extent; }
    return idx;
}

void poison_ghost_cells (amrex::MultiFab& mf,
                         const std::initializer_list<int>& components,
                         const amrex::Real sentinel)
{
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box& fab_box = mfi.fabbox();
        const amrex::Box& valid_box = mfi.validbox();
        auto arr = mf.array(mfi);

        run_and_sync([=] {
            amrex::ParallelFor(fab_box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                if (!valid_box.contains(amrex::IntVect(AMREX_D_DECL(i, j, k)))) {
                    for (const int comp : components) {
                        arr(i,j,k,comp) = sentinel;
                    }
                }
            });
        });
    }
}

amrex::Real expected_full_flow_component (const int component,
                                     const int i,
                                     const int j,
                                     const int k,
                                     const CellState& evap_then_recond)
{
    const CellState initial_state = make_satadj_active_cell_state(i, j, k, evap_then_recond);
    const ConservedState initial_conserved = make_conserved_state(initial_state);
    const ConservedState reference = scalar_reference_from_initial_conserved(
        initial_conserved.rho,
        initial_conserved.rhotheta,
        initial_conserved.rhoqv,
        initial_conserved.rhoqc);
    return conserved_component(reference, component);
}

} // namespace

// Motivation: The serial box tests catch layout dependence on one rank. This
// MPI test exercises the same SatAdj public flow with distributed MultiFabs and
// verifies that rank ownership does not change the result of the local
// fixed-pressure cell adjustment.
TEST(SatAdjParallel, DistributedMultiBoxMatchesScalarReference)
{
    constexpr int nx = 31;
    constexpr int ny = 27;
    constexpr int nz = 6;
    const amrex::Geometry geom = make_geometry(nx, ny, nz);
    const amrex::IntVect max_size(AMREX_D_DECL(7, 5, 3));
    const amrex::BoxArray ba = make_boxarray(geom.Domain(), max_size);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
    amrex::MultiFab before(ba, dm, RhoQ2_comp + 1, 0);
    fill_satadj_active_conserved_state(cons);
    amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), 0);

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);
    run_satadj_public_flow(satadj, sc, geom, cons);

    LocatedError local_worst{};
    ActiveCellCounts local_counts{};
    compare_public_flow_to_scalar_reference(before, cons, local_worst, local_counts);

    long global_checked = local_counts.checked;
    long global_changed = local_counts.changed;
    amrex::ParallelDescriptor::ReduceLongSum(global_checked);
    amrex::ParallelDescriptor::ReduceLongSum(global_changed);
    amrex::Real global_max_error = local_worst.normalized_error;
    amrex::ParallelDescriptor::ReduceRealMax(global_max_error);

    const int rank = amrex::ParallelDescriptor::MyProc();
    const long expected_cells = num_valid_cells(geom.Domain());
    const amrex::Real normalized_tol =
#ifdef AMREX_USE_FLOAT
        amrex::Real(4.0);
#else
        amrex::Real(2.0);
#endif

    EXPECT_EQ(global_checked, expected_cells)
        << "rank=" << rank
        << " max_size=" << max_size_string(max_size)
        << " local_checked=" << local_counts.checked
        << " global_checked=" << global_checked;
    EXPECT_EQ(global_changed, expected_cells)
        << "rank=" << rank
        << " max_size=" << max_size_string(max_size)
        << " local_changed=" << local_counts.changed
        << " global_changed=" << global_changed
        << " local_unexpected_unchanged=" << local_counts.unexpected_unchanged;
    EXPECT_LE(global_max_error, normalized_tol)
        << "rank=" << rank
        << " max_size=" << max_size_string(max_size)
        << " local_component=" << component_name(local_worst.component)
        << " local_i=" << local_worst.i
        << " local_j=" << local_worst.j
        << " local_k=" << local_worst.k
        << " local_box_id=" << local_worst.box_id
        << " local_actual=" << local_worst.actual_value
        << " local_expected=" << local_worst.expected_value
        << " local_absolute_error=" << local_worst.absolute_error
        << " local_normalized_error=" << local_worst.normalized_error
        << " global_max_normalized_error=" << global_max_error
        << " eos_actual_temperature=" << local_worst.eos_projected_actual_temperature
        << " eos_expected_temperature=" << local_worst.eos_projected_expected_temperature;
}

// Motivation: This extends the parallel copy-back ghost test to the full SatAdj
// public flow. After distributed Advance, ghosts are poisoned before
// Copy_Micro_to_State so the test directly verifies that final FillBoundary
// overwrites stale data across rank and box boundaries.
TEST(SatAdjParallel, FullPublicFlowFillBoundaryParallel)
{
    constexpr int nx = 16;
    constexpr int ny = 12;
    constexpr int nz = 4;
    const amrex::Geometry geom = make_geometry(nx, ny, nz);
    const amrex::IntVect max_size(AMREX_D_DECL(5, 4, 2));
    const amrex::BoxArray ba = make_boxarray(geom.Domain(), max_size);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 1);
    fill_satadj_active_conserved_state(cons);

    CellState evap_then_recond;
    ASSERT_TRUE(find_evaporation_then_recondensation_state(evap_then_recond));

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    // Split the public flow so ghosts are poisoned after Advance and before
    // Copy_Micro_to_State. This checks that FillBoundary is responsible for
    // restoring wrapped valid-state values across rank boundaries.
    run_and_sync([&]() {
        satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Advance(amrex::Real(1.0), sc);
    });

    poison_ghost_cells(cons,
                       {Rho_comp, RhoTheta_comp, RhoQ1_comp, RhoQ2_comp},
                       amrex::Real(-999.0));

    run_and_sync([&]() { satadj.Copy_Micro_to_State(cons); });

    const amrex::MultiFab after_host = make_host_mirror(cons);
    const amrex::Box& domain = geom.Domain();
    LocatedError local_worst{};
    long local_checked_ghosts = 0;

    for (amrex::MFIter mfi(after_host); mfi.isValid(); ++mfi) {
        const amrex::Box& fab_box = mfi.fabbox();
        const amrex::Box& valid_box = mfi.validbox();
        const auto after_arr = after_host.const_array(mfi);

        for (int k = fab_box.smallEnd(2); k <= fab_box.bigEnd(2); ++k) {
            for (int j = fab_box.smallEnd(1); j <= fab_box.bigEnd(1); ++j) {
                for (int i = fab_box.smallEnd(0); i <= fab_box.bigEnd(0); ++i) {
                    const amrex::IntVect iv(AMREX_D_DECL(i, j, k));
                    if (valid_box.contains(iv)) {
                        continue;
                    }

                    ++local_checked_ghosts;
                    const int wrapped_i = wrap_index(i, domain.smallEnd(0), domain.bigEnd(0));
                    const int wrapped_j = wrap_index(j, domain.smallEnd(1), domain.bigEnd(1));
                    const int wrapped_k = wrap_index(k, domain.smallEnd(2), domain.bigEnd(2));

                    for (const int comp : {Rho_comp, RhoTheta_comp, RhoQ1_comp, RhoQ2_comp}) {
                        const amrex::Real actual = after_arr(i,j,k,comp);
                        const amrex::Real expected = expected_full_flow_component(comp,
                                                                  wrapped_i,
                                                                  wrapped_j,
                                                                  wrapped_k,
                                                                  evap_then_recond);
                        const amrex::Real absolute_error = amrex::Math::abs(actual - expected);
                        const amrex::Real tol = [&]() -> amrex::Real {
                            if (comp == RhoQ1_comp || comp == RhoQ2_comp) {
                                const amrex::Real rho_exp = expected_full_flow_component(
                                    Rho_comp, wrapped_i, wrapped_j, wrapped_k, evap_then_recond);
                                return conserved_moisture_tol(rho_exp, expected / rho_exp);
                            }
                            return scaled_tol(expected, amrex::Real(10.0) * kStateTolFactor);
                        }();
                        const amrex::Real normalized_error = absolute_error / tol;
                        if (normalized_error > local_worst.normalized_error) {
                            local_worst = LocatedError{normalized_error,
                                               absolute_error,
                                               actual,
                                               expected,
                                               amrex::Real(0.0),
                                               amrex::Real(0.0),
                                               i,
                                               j,
                                               k,
                                               comp,
                                               mfi.index()};
                        }
                    }
                }
            }
        }
    }

    long global_checked_ghosts = local_checked_ghosts;
    amrex::ParallelDescriptor::ReduceLongSum(global_checked_ghosts);
    amrex::Real global_max_error = local_worst.normalized_error;
    amrex::ParallelDescriptor::ReduceRealMax(global_max_error);
    const int rank = amrex::ParallelDescriptor::MyProc();

    EXPECT_GT(global_checked_ghosts, 0)
        << "rank=" << rank
        << " max_size=" << max_size_string(max_size)
        << " local_checked_ghosts=" << local_checked_ghosts;
    EXPECT_LE(global_max_error, amrex::Real(1.0))
        << "rank=" << rank
        << " max_size=" << max_size_string(max_size)
        << " local_component=" << component_name(local_worst.component)
        << " ghost_i=" << local_worst.i
        << " ghost_j=" << local_worst.j
        << " ghost_k=" << local_worst.k
        << " local_box_id=" << local_worst.box_id
        << " local_actual=" << local_worst.actual_value
        << " local_expected=" << local_worst.expected_value
        << " local_absolute_error=" << local_worst.absolute_error
        << " global_max_normalized_error=" << global_max_error;
}
