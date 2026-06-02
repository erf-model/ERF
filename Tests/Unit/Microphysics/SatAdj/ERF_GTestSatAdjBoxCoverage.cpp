#include <array>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include <ERF_Kessler.H>

#include "ERF_GTestSatAdjCommon.H"

using namespace satadj_test;

namespace {

struct LayoutCase {
    std::string name;
    bool single_box;
    amrex::IntVect max_size;
};

std::vector<LayoutCase> make_layout_cases ()
{
    return {
        {"single-box", true, amrex::IntVect(AMREX_D_DECL(23, 19, 5))},
        {"maxsize-7x5x5", false, amrex::IntVect(AMREX_D_DECL(7, 5, 5))},
        {"maxsize-4x11x5", false, amrex::IntVect(AMREX_D_DECL(4, 11, 5))},
        {"maxsize-6x6x2", false, amrex::IntVect(AMREX_D_DECL(6, 6, 2))}
    };
}

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

std::string max_size_string (const LayoutCase& layout)
{
    if (layout.single_box) {
        return "domain";
    }

    std::ostringstream os;
    os << "(" << layout.max_size[0] << "," << layout.max_size[1] << "," << layout.max_size[2] << ")";
    return os.str();
}

amrex::BoxArray make_layout_boxarray (const amrex::Box& domain,
                                      const LayoutCase& layout)
{
    return layout.single_box ? make_boxarray(domain) : make_boxarray(domain, layout.max_size);
}

SolverChoice make_kessler_norain_solver_choice ()
{
    SolverChoice sc{};
    sc.c_p = Cp_d;
    sc.rdOcp = kRdOcp;
    sc.moisture_type = MoistureType::Kessler_NoRain;
    sc.use_shoc = false;
    return sc;
}

void run_kessler_norain_public_flow (Kessler& kessler,
                                     const SolverChoice& sc,
                                     const amrex::Geometry& geom,
                                     amrex::MultiFab& cons,
                                     const amrex::Real dt = amrex::Real(1.0))
{
    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    run_and_sync([&]() {
        kessler.Set_dzmin(geom.CellSize(2));
        kessler.Init(cons, cons.boxArray(), geom, dt, z_phys_nd, detJ_cc);
        kessler.Copy_State_to_Micro(cons);
        kessler.Advance(dt, sc);
        kessler.Copy_Micro_to_State(cons);
    });
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

// Motivation: Current SatAdj advances over each MFIter tilebox with no
// real-width trimming. Every valid cell initialized with an active SatAdj state
// should change. An unchanged active cell indicates skipped execution, likely
// from tilebox, BoxArray, or copy staging rather than valid no-op physics.
TEST(SatAdjBoxCoverage, ActiveRegionCoversEveryValidCell)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);

    for (const LayoutCase& layout : make_layout_cases()) {
        const amrex::BoxArray ba = make_layout_boxarray(geom.Domain(), layout);
        const amrex::DistributionMapping dm(ba);

        for (const int real_width : {0, 1, 2}) {
            amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
            amrex::MultiFab before(ba, dm, RhoQ2_comp + 1, 0);
            fill_satadj_active_conserved_state(cons);
            amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), 0);

            SatAdj satadj;
            SolverChoice sc = make_solver_choice(false);
            satadj.Define(sc);
            satadj.Set_RealWidth(real_width);
            run_satadj_public_flow(satadj, sc, geom, cons);

            const amrex::MultiFab before_host = make_host_mirror(before);
            const amrex::MultiFab after_host = make_host_mirror(cons);

            ActiveCellCounts counts{};
            int first_i = -1;
            int first_j = -1;
            int first_k = -1;
            int first_box = -1;
            std::array<amrex::Real, 3> first_deltas{amrex::Real(0.0), amrex::Real(0.0), amrex::Real(0.0)};

            for (amrex::MFIter mfi(after_host); mfi.isValid(); ++mfi) {
                const amrex::Box& bx = mfi.validbox();
                const auto before_arr = before_host.const_array(mfi);
                const auto after_arr = after_host.const_array(mfi);

                for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
                    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                            ++counts.checked;
                            const amrex::Real delta_rhotheta = amrex::Math::abs(after_arr(i,j,k,RhoTheta_comp) - before_arr(i,j,k,RhoTheta_comp));
                            const amrex::Real delta_rhoq1 = amrex::Math::abs(after_arr(i,j,k,RhoQ1_comp) - before_arr(i,j,k,RhoQ1_comp));
                            const amrex::Real delta_rhoq2 = amrex::Math::abs(after_arr(i,j,k,RhoQ2_comp) - before_arr(i,j,k,RhoQ2_comp));

                            const amrex::Real rho_before = before_arr(i,j,k,Rho_comp);
                            const bool changed =
                                satadj_written_component_changed(RhoTheta_comp,
                                                                  rho_before,
                                                                  before_arr(i,j,k,RhoTheta_comp),
                                                                  after_arr(i,j,k,RhoTheta_comp)) ||
                                satadj_written_component_changed(RhoQ1_comp,
                                                                  rho_before,
                                                                  before_arr(i,j,k,RhoQ1_comp),
                                                                  after_arr(i,j,k,RhoQ1_comp)) ||
                                satadj_written_component_changed(RhoQ2_comp,
                                                                  rho_before,
                                                                  before_arr(i,j,k,RhoQ2_comp),
                                                                  after_arr(i,j,k,RhoQ2_comp));

                            if (changed) {
                                ++counts.changed;
                            } else {
                                ++counts.unexpected_unchanged;
                                if (first_i < 0) {
                                    first_i = i;
                                    first_j = j;
                                    first_k = k;
                                    first_box = mfi.index();
                                    first_deltas = {delta_rhotheta, delta_rhoq1, delta_rhoq2};
                                }
                            }
                        }
                    }
                }
            }

            EXPECT_EQ(counts.changed, geom.Domain().numPts())
                << "layout=" << layout.name
                << " max_size=" << max_size_string(layout)
                << " real_width=" << real_width
                << " rank=0"
                << " checked=" << counts.checked
                << " changed=" << counts.changed
                << " unexpected_unchanged=" << counts.unexpected_unchanged
                << " first_box=" << first_box
                << " first_i=" << first_i
                << " first_j=" << first_j
                << " first_k=" << first_k
                << " delta_rhotheta=" << first_deltas[0]
                << " delta_rhoq1=" << first_deltas[1]
                << " delta_rhoq2=" << first_deltas[2];
        }
    }
}

// Motivation: SatAdj is a local fixed-pressure cell update. Its public
// MultiFab path should match the scalar AdjustSatAdjCell reference regardless
// of how AMReX decomposes the valid domain into boxes. A box-aligned
// temperature signal should appear as a localized conserved-state or
// EOS-projected-temperature mismatch.
TEST(SatAdjBoxCoverage, DecompositionInvariantPublicFlow)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);
    const amrex::Real normalized_tol =
#ifdef AMREX_USE_FLOAT
        amrex::Real(4.0);
#else
        amrex::Real(2.0);
#endif

    for (const LayoutCase& layout : make_layout_cases()) {
        const amrex::BoxArray ba = make_layout_boxarray(geom.Domain(), layout);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
        amrex::MultiFab before(ba, dm, RhoQ2_comp + 1, 0);
        fill_satadj_active_conserved_state(cons);
        amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), 0);

        SatAdj satadj;
        SolverChoice sc = make_solver_choice(false);
        satadj.Define(sc);
        run_satadj_public_flow(satadj, sc, geom, cons);

        LocatedError worst{};
        ActiveCellCounts counts{};
        compare_public_flow_to_scalar_reference(before, cons, worst, counts);

        EXPECT_EQ(counts.checked, geom.Domain().numPts())
            << "layout=" << layout.name
            << " max_size=" << max_size_string(layout);
        EXPECT_EQ(counts.changed, geom.Domain().numPts())
            << "layout=" << layout.name
            << " max_size=" << max_size_string(layout)
            << " unexpected_unchanged=" << counts.unexpected_unchanged;
        EXPECT_LE(worst.normalized_error, normalized_tol)
            << "layout=" << layout.name
            << " max_size=" << max_size_string(layout)
            << " component=" << component_name(worst.component)
            << " i=" << worst.i
            << " j=" << worst.j
            << " k=" << worst.k
            << " box_id=" << worst.box_id
            << " actual=" << worst.actual_value
            << " expected=" << worst.expected_value
            << " absolute_error=" << worst.absolute_error
            << " normalized_error=" << worst.normalized_error
            << " eos_actual_temperature=" << worst.eos_projected_actual_temperature
            << " eos_expected_temperature=" << worst.eos_projected_expected_temperature;
    }
}

// Motivation: This isolates SatAdj staging from the saturation-adjustment
// solve. If a box-aligned signal appears here, the cause is in Init,
// Copy_State_to_Micro, Copy_Micro_to_State, or FillBoundary rather than in
// AdjustSatAdjCell.
TEST(SatAdjBoxCoverage, CopyOnlyRoundTripMultiBox)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);
    const LayoutCase layout{"maxsize-7x5x5", false, amrex::IntVect(AMREX_D_DECL(7, 5, 5))};
    const amrex::BoxArray ba = make_layout_boxarray(geom.Domain(), layout);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 0);
    amrex::MultiFab before(ba, dm, RhoQ2_comp + 1, 0);
    fill_satadj_active_conserved_state(cons);
    amrex::MultiFab::Copy(before, cons, 0, 0, cons.nComp(), 0);

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    run_and_sync([&]() {
        satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Copy_Micro_to_State(cons);
    });

    const amrex::MultiFab before_host = make_host_mirror(before);
    const amrex::MultiFab after_host = make_host_mirror(cons);
    LocatedError worst{};

    for (amrex::MFIter mfi(after_host); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.validbox();
        const auto before_arr = before_host.const_array(mfi);
        const auto after_arr = after_host.const_array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                    for (const int comp : {Rho_comp, RhoTheta_comp, RhoQ1_comp, RhoQ2_comp}) {
                        const amrex::Real actual = after_arr(i,j,k,comp);
                        const amrex::Real expected = before_arr(i,j,k,comp);
                        const amrex::Real absolute_error = amrex::Math::abs(actual - expected);
                        const amrex::Real tol = (comp == RhoQ1_comp || comp == RhoQ2_comp)
                            ? conserved_moisture_tol(before_arr(i,j,k,Rho_comp), expected / before_arr(i,j,k,Rho_comp))
                            : scaled_tol(expected, amrex::Real(10.0) * kStateTolFactor);
                        const amrex::Real normalized_error = absolute_error / tol;
                        if (normalized_error > worst.normalized_error) {
                            worst = LocatedError{normalized_error,
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

    EXPECT_LE(worst.normalized_error, amrex::Real(1.0))
        << "layout=" << layout.name
        << " max_size=" << max_size_string(layout)
        << " component=" << component_name(worst.component)
        << " i=" << worst.i
        << " j=" << worst.j
        << " k=" << worst.k
        << " box_id=" << worst.box_id
        << " actual=" << worst.actual_value
        << " expected=" << worst.expected_value
        << " absolute_error=" << worst.absolute_error
        << " normalized_error=" << worst.normalized_error;
}

// Motivation: The copy-only ghost-fill test isolates staging. This full-flow
// variant poisons ghosts after SatAdj Advance and before Copy_Micro_to_State,
// checking that the final FillBoundary reflects adjusted valid cells rather
// than stale ghost data.
TEST(SatAdjBoxCoverage, PeriodicGhostFillAfterFullPublicFlow)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);
    const LayoutCase layout{"maxsize-6x6x2", false, amrex::IntVect(AMREX_D_DECL(6, 6, 2))};
    const amrex::BoxArray ba = make_layout_boxarray(geom.Domain(), layout);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab cons(ba, dm, RhoQ2_comp + 1, 2);
    fill_satadj_active_conserved_state(cons);

    CellState evap_then_recond;
    ASSERT_TRUE(find_evaporation_then_recondensation_state(evap_then_recond));

    SatAdj satadj;
    SolverChoice sc = make_solver_choice(false);
    satadj.Define(sc);

    std::unique_ptr<amrex::MultiFab> z_phys_nd;
    std::unique_ptr<amrex::MultiFab> detJ_cc;
    // Explicit public-flow split: Init, Copy_State_to_Micro, Advance,
    // poison ghosts, Copy_Micro_to_State. Poisoning after Advance isolates
    // FillBoundary as the sole mechanism restoring valid wrapped values.
    run_and_sync([&]() {
        satadj.Init(cons, cons.boxArray(), geom, amrex::Real(1.0), z_phys_nd, detJ_cc);
        satadj.Copy_State_to_Micro(cons);
        satadj.Advance(amrex::Real(1.0), sc);
    });

    poison_ghost_cells(cons,
                       amrex::GpuArray<int, 4>{Rho_comp, RhoTheta_comp, RhoQ1_comp, RhoQ2_comp},
                       amrex::Real(-999.0));

    run_and_sync([&]() { satadj.Copy_Micro_to_State(cons); });

    const amrex::MultiFab after_host = make_host_mirror(cons);
    const amrex::Box& domain = geom.Domain();
    LocatedError worst{};
    long checked_ghosts = 0;

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

                    ++checked_ghosts;
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
                        if (normalized_error > worst.normalized_error) {
                            worst = LocatedError{normalized_error,
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

    EXPECT_GT(checked_ghosts, 0);
    EXPECT_LE(worst.normalized_error, amrex::Real(1.0))
        << "layout=" << layout.name
        << " max_size=" << max_size_string(layout)
        << " component=" << component_name(worst.component)
        << " ghost_i=" << worst.i
        << " ghost_j=" << worst.j
        << " ghost_k=" << worst.k
        << " box_id=" << worst.box_id
        << " actual=" << worst.actual_value
        << " expected=" << worst.expected_value
        << " absolute_error=" << worst.absolute_error
        << " normalized_error=" << worst.normalized_error;
}

// Motivation: SatAdj and Kessler_NoRain use different adjustment formulas, but
// current source advances both over the same MFIter tilebox region. This test
// compares update masks, not values, to catch future indexing or tilebox-region
// mismatches.
TEST(SatAdjBoxCoverage, ActiveRegionMatchesKesslerNoRain)
{
    const amrex::Geometry geom = make_geometry(23, 19, 5);
    const LayoutCase layout{"maxsize-7x5x5", false, amrex::IntVect(AMREX_D_DECL(7, 5, 5))};
    const amrex::BoxArray ba = make_layout_boxarray(geom.Domain(), layout);
    const amrex::DistributionMapping dm(ba);

    for (const int real_width : {0, 1, 2}) {
        amrex::MultiFab satadj_cons(ba, dm, RhoQ3_comp + 1, 0);
        amrex::MultiFab satadj_before(ba, dm, RhoQ3_comp + 1, 0);
        amrex::MultiFab kessler_cons(ba, dm, RhoQ3_comp + 1, 0);
        amrex::MultiFab kessler_before(ba, dm, RhoQ3_comp + 1, 0);
        fill_satadj_kessler_active_conserved_state(satadj_cons);
        amrex::MultiFab::Copy(kessler_cons, satadj_cons, 0, 0, satadj_cons.nComp(), 0);
        amrex::MultiFab::Copy(satadj_before, satadj_cons, 0, 0, satadj_cons.nComp(), 0);
        amrex::MultiFab::Copy(kessler_before, kessler_cons, 0, 0, kessler_cons.nComp(), 0);

        SatAdj satadj;
        SolverChoice satadj_sc = make_solver_choice(false);
        satadj.Define(satadj_sc);
        satadj.Set_RealWidth(real_width);
        run_satadj_public_flow(satadj, satadj_sc, geom, satadj_cons);

        Kessler kessler;
        SolverChoice kessler_sc = make_kessler_norain_solver_choice();
        kessler.Define(kessler_sc);
        kessler.Set_RealWidth(real_width);
        run_kessler_norain_public_flow(kessler, kessler_sc, geom, kessler_cons);

        const amrex::MultiFab satadj_before_host = make_host_mirror(satadj_before);
        const amrex::MultiFab satadj_after_host = make_host_mirror(satadj_cons);
        const amrex::MultiFab kessler_before_host = make_host_mirror(kessler_before);
        const amrex::MultiFab kessler_after_host = make_host_mirror(kessler_cons);

        long satadj_changed = 0;
        long kessler_changed = 0;
        int mismatch_i = -1;
        int mismatch_j = -1;
        int mismatch_k = -1;
        int mismatch_box = -1;
        bool satadj_mask_value = false;
        bool kessler_mask_value = false;

        for (amrex::MFIter mfi(satadj_after_host); mfi.isValid(); ++mfi) {
            const amrex::Box& bx = mfi.validbox();
            const auto satadj_before_arr = satadj_before_host.const_array(mfi);
            const auto satadj_after_arr = satadj_after_host.const_array(mfi);
            const auto kessler_before_arr = kessler_before_host.const_array(mfi);
            const auto kessler_after_arr = kessler_after_host.const_array(mfi);

            for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        const amrex::Real satadj_rho_before = satadj_before_arr(i,j,k,Rho_comp);
                        const bool satadj_changed_here =
                            satadj_written_component_changed(RhoTheta_comp,
                                                              satadj_rho_before,
                                                              satadj_before_arr(i,j,k,RhoTheta_comp),
                                                              satadj_after_arr(i,j,k,RhoTheta_comp)) ||
                            satadj_written_component_changed(RhoQ1_comp,
                                                              satadj_rho_before,
                                                              satadj_before_arr(i,j,k,RhoQ1_comp),
                                                              satadj_after_arr(i,j,k,RhoQ1_comp)) ||
                            satadj_written_component_changed(RhoQ2_comp,
                                                              satadj_rho_before,
                                                              satadj_before_arr(i,j,k,RhoQ2_comp),
                                                              satadj_after_arr(i,j,k,RhoQ2_comp));
                        const amrex::Real kessler_rho_before = kessler_before_arr(i,j,k,Rho_comp);
                        const bool kessler_changed_here =
                            satadj_written_component_changed(RhoTheta_comp,
                                                              kessler_rho_before,
                                                              kessler_before_arr(i,j,k,RhoTheta_comp),
                                                              kessler_after_arr(i,j,k,RhoTheta_comp)) ||
                            satadj_written_component_changed(RhoQ1_comp,
                                                              kessler_rho_before,
                                                              kessler_before_arr(i,j,k,RhoQ1_comp),
                                                              kessler_after_arr(i,j,k,RhoQ1_comp)) ||
                            satadj_written_component_changed(RhoQ2_comp,
                                                              kessler_rho_before,
                                                              kessler_before_arr(i,j,k,RhoQ2_comp),
                                                              kessler_after_arr(i,j,k,RhoQ2_comp));

                        satadj_changed += satadj_changed_here ? 1 : 0;
                        kessler_changed += kessler_changed_here ? 1 : 0;

                        if (mismatch_i < 0 && satadj_changed_here != kessler_changed_here) {
                            mismatch_i = i;
                            mismatch_j = j;
                            mismatch_k = k;
                            mismatch_box = mfi.index();
                            satadj_mask_value = satadj_changed_here;
                            kessler_mask_value = kessler_changed_here;
                        }
                    }
                }
            }
        }

        EXPECT_EQ(satadj_changed, geom.Domain().numPts())
            << "layout=" << layout.name
            << " real_width=" << real_width
            << " satadj_changed=" << satadj_changed;
        EXPECT_EQ(kessler_changed, geom.Domain().numPts())
            << "layout=" << layout.name
            << " real_width=" << real_width
            << " kessler_changed=" << kessler_changed;
        EXPECT_LT(mismatch_i, 0)
            << "layout=" << layout.name
            << " real_width=" << real_width
            << " box_id=" << mismatch_box
            << " i=" << mismatch_i
            << " j=" << mismatch_j
            << " k=" << mismatch_k
            << " satadj_mask_value=" << satadj_mask_value
            << " kessler_mask_value=" << kessler_mask_value;
    }
}