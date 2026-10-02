#include <gtest/gtest.h>

#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

#include "ERF_IndexDefines.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::Geometry;
using amrex::IntVect;
using amrex::MultiFab;
using amrex::Real;

BoxArray
project_to_xy (const BoxArray& cell_ba)
{
    amrex::BoxList boxes = cell_ba.boxList();
    for (auto& box : boxes)
        box.setRange(2, 0);
    return BoxArray(std::move(boxes));
}

erf_sbm::SBMLayout
make_parallel_layout (const erf_sbm::MomentMode mode)
{
    erf_sbm::SpectralGridSpec grid;
    grid.coordinate_kind = erf_sbm::CoordinateKind::Mass;
    grid.coordinate_units = "kg particle^-1";
    grid.edges = {Real(0.1), Real(0.5), Real(1.0)};
    grid.pivots = {Real(0.3), Real(0.75)};

    erf_sbm::SpectralPopulationSpec population;
    population.population_id = 0;
    population.semantic_id = "liquid";
    population.phase = erf_sbm::PopulationPhase::Liquid;
    population.grid = std::move(grid);
    population.moment_mode = mode;

    erf_sbm::SBMLayoutSpec spec;
    spec.populations.push_back(std::move(population));
    spec.liquid_projection = {0, 1};
    return erf_sbm::SBMLayout(std::move(spec));
}

std::vector<Real>
run_decomposition (const int max_grid_size,
                   const erf_sbm::MomentMode mode)
{
    constexpr int nx = 16;
    constexpr int ny = 4;
    constexpr int nz = 4;
    const Box domain(IntVect(0, 0, 0), IntVect(nx - 1, ny - 1, nz - 1));
    const amrex::RealBox physical({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const int periodic[AMREX_SPACEDIM] = {1, 1, 1};
    const Geometry geom(domain, &physical, amrex::CoordSys::cartesian,
                        periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    const DistributionMapping dm(ba);
    auto layout = make_parallel_layout(mode);

    erf_sbm::SBMStateManager manager(layout, 1);
    manager.define(0, ba, dm);
    erf_sbm::SBMTransport transport(manager.layout(), 1, 1);
    transport.define(0, ba, dm);

    MultiFab detj(ba, dm, 1, 0);
    const BoxArray map_ba = project_to_xy(ba);
    MultiFab mx(map_ba, dm, 1, 0);
    MultiFab my(map_ba, dm, 1, 0);
    detj.setVal(Real(1.0));
    mx.setVal(Real(1.0));
    my.setVal(Real(1.0));
    std::string diagnostic;
    if (!transport.rebuild_static_measure(0, detj, mx, my, diagnostic)) {
        ADD_FAILURE() << diagnostic;
        return {};
    }

    MultiFab conserved_anchor(ba, dm, 3, 0);
    MultiFab conserved_input(ba, dm, 3, 0);
    MultiFab conserved_target(ba, dm, 3, 0);
    conserved_anchor.setVal(Real(0.0));
    conserved_anchor.setVal(Real(1.0), Rho_comp, 1, 0);
    MultiFab::Copy(conserved_input, conserved_anchor, 0, 0, 3, 0);
    MultiFab::Copy(conserved_target, conserved_anchor, 0, 0, 3, 0);

    auto& spectrum = manager.state(0);
    spectrum.setVal(Real(0.0));
    const auto& population = manager.layout().populations().front();
    for (amrex::MFIter mfi(spectrum); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = spectrum.array(mfi);
        const int mass0 = population.mass_offset;
        const int mass1 = mass0 + 1;
        const int number0 = population.number_offset;
        const int number1 = number0 < 0 ? -1 : number0 + 1;
        const bool two_moment = mode == erf_sbm::MomentMode::TwoMoment;
        amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                if (two_moment) {
                    // Deliberately sharp per-bin number and mean-mass changes
                    // exercise endpoint reconstruction and grouped limiting.
                    const int pattern = (17 * i + 31 * j + 13 * k) % 7;
                    const Real n0 = pattern == 0 ? Real(0.01) : Real(0.00001);
                    const Real mean0 = pattern % 2 == 0 ? Real(0.11) : Real(0.45);
                    const Real n1 = pattern == 1 ? Real(0.02) : Real(0.00002);
                    const Real mean1 = pattern % 3 == 0 ? Real(0.51) : Real(0.99);
                    state(i, j, k, mass0) = n0 * mean0;
                    state(i, j, k, mass1) = n1 * mean1;
                    state(i, j, k, number0) = n0;
                    state(i, j, k, number1) = n1;
                } else {
                    const Real phase = static_cast<Real>(i % 8);
                    state(i, j, k, mass0) =
                        Real(0.001) + Real(0.0002) * phase;
                    state(i, j, k, mass1) =
                        Real(0.002) +
                        Real(0.0001) * static_cast<Real>((3 * i + j + k) % 7);
                }
            });
    }

    MultiFab avg_xmom(amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, 1,
                      0);
    MultiFab avg_ymom(amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, 1,
                      0);
    MultiFab avg_zmom(amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, 1,
                      0);
    avg_xmom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.13)
                                                          : Real(0.05));
    avg_ymom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.11)
                                                          : Real(0.0));
    avg_zmom.setVal(mode == erf_sbm::MomentMode::TwoMoment ? Real(0.09)
                                                          : Real(0.0));

    const double dt = mode == erf_sbm::MomentMode::TwoMoment ? 0.25 : 0.01;
    transport.advance_stage_from_host(
        0, erf_auxiliary::HostIntegrator::CompressibleRK3,
                            0, 0.0, 0.0, dt, dt, manager, conserved_anchor,
                            conserved_input, conserved_target, avg_xmom,
                            avg_ymom, avg_zmom, geom, 1, 2);

    const auto& accepted = manager.state(0);
    MultiFab signatures(ba, dm, 2 * layout.ncomp(), 0);
    for (amrex::MFIter mfi(accepted); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto state = accepted.const_array(mfi);
        const auto out = signatures.array(mfi);
        const int ncomp = layout.ncomp();
        amrex::ParallelFor(
            bx, ncomp,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                const Real value = state(i, j, k, n);
                const Real weight = Real(1.0) + Real(0.001) * i +
                                    Real(0.0001) * j + Real(0.00001) * k;
                out(i, j, k, n) = weight * value;
                out(i, j, k, ncomp + n) = value * value;
            });
    }

    std::vector<Real> result;
    for (int component = 0; component < layout.ncomp(); ++component) {
        result.push_back(accepted.sum(component));
        result.push_back(signatures.sum(component));
        result.push_back(signatures.sum(layout.ncomp() + component));
    }
    return result;
}

TEST(SBMTransportParallel, DecompositionInvariant)
{
    const auto one_box = run_decomposition(16, erf_sbm::MomentMode::OneMoment);
    const auto many_boxes = run_decomposition(4, erf_sbm::MomentMode::OneMoment);
    ASSERT_EQ(one_box.size(), many_boxes.size());
    for (std::size_t i = 0; i < one_box.size(); ++i) {
        EXPECT_NEAR(one_box[i], many_boxes[i],
                    Real(512.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(std::abs(one_box[i]),
                                 std::numeric_limits<Real>::min()));
    }
}

TEST(SBMTransportParallel, TwoMomentEndpointGroupLimiterDecompositionInvariant)
{
    const auto one_box = run_decomposition(16, erf_sbm::MomentMode::TwoMoment);
    const auto many_boxes = run_decomposition(4, erf_sbm::MomentMode::TwoMoment);
    ASSERT_EQ(one_box.size(), many_boxes.size());
    for (std::size_t i = 0; i < one_box.size(); ++i) {
        EXPECT_NEAR(one_box[i], many_boxes[i],
                    Real(512.0) * std::numeric_limits<Real>::epsilon() *
                        std::max(std::abs(one_box[i]),
                                 std::numeric_limits<Real>::min()));
    }
}

} // namespace
