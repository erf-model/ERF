#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_ForestDrag.H"
#include "ERF_ForestUtils.H"
#include "../ERF_GTestTempDir.H"

namespace {

void
fill_vertical_coordinates (amrex::MultiFab& field, amrex::Real offset)
{
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const auto z = field.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            z(i, j, k) = (static_cast<amrex::Real>(k) + offset) * amrex::Real(62.5);
        });
    }
}

// Drag and frontal-area fields of the patches in rows, on a 16^3 box of 62.5 m
// cells (cell centres at (k + 1/2) 62.5 m above a flat ground at z = 0).
struct DiscreteForestFields
{
    amrex::MultiFab drag;
    amrex::MultiFab frontal_area;
};

DiscreteForestFields
discrete_forest_fields (const std::string& rows, const std::string& stem)
{
    const std::filesystem::path filename = erf_gtest_temp_path(stem).string() + ".txt";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream output(filename);
        output << rows;
    }
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 15));
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1000.0, 1000.0, 1000.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{1, 1, 0};
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());
    amrex::MultiFab z_cc(ba, dm, 1, 1);
    amrex::MultiFab z_nd(amrex::convert(ba, amrex::IntVect(1, 1, 1)), dm, 1, 1);
    fill_vertical_coordinates(z_cc, amrex::Real(0.5));
    fill_vertical_coordinates(z_nd, amrex::Real(0.0));
    amrex::Gpu::streamSynchronize();

    ForestDrag forest(filename.string());
    forest.define_drag_field(ba, dm, geom, &z_cc, &z_nd, true, 0);
    amrex::Gpu::streamSynchronize();

    // Pinned, so the tests can read the fields on the host in a GPU build too
    const amrex::MFInfo pinned = amrex::MFInfo().SetArena(amrex::The_Pinned_Arena());
    DiscreteForestFields fields{amrex::MultiFab(ba, dm, 1, 1, pinned),
                                amrex::MultiFab(ba, dm, 1, 1, pinned)};
    amrex::MultiFab::Copy(fields.drag, *forest.get_drag_field(), 0, 0, 1, 1);
    amrex::MultiFab::Copy(fields.frontal_area, *forest.get_frontal_area(), 0, 0, 1, 1);
    amrex::Gpu::streamSynchronize();

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::filesystem::remove(filename);
    }
    amrex::ParallelDescriptor::Barrier();
    return fields;
}

// max |a + b - sum| over the valid cells
amrex::Real
max_sum_error (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::MultiFab& sum)
{
    amrex::MultiFab diff(sum.boxArray(), sum.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(diff, sum, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(diff, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(diff, b, 0, 0, 1, 0);
    return diff.norm0(0, 0, false);
}

} // namespace

TEST(ForestUtils, ValidatesLaimaxWithoutAborting)
{
    EXPECT_TRUE(erf_forest_utils::validate_laimax(0.0, "erf.forest_laimax").empty());
    EXPECT_TRUE(erf_forest_utils::validate_laimax(
        std::nextafter(amrex::Real(1.0), amrex::Real(0.0)), "erf.forest_laimax").empty());

    for (const amrex::Real value : {amrex::Real(-0.01), amrex::Real(1.0), amrex::Real(1.01),
                                    amrex::Real(NAN), amrex::Real(INFINITY)}) {
        const auto error = erf_forest_utils::validate_laimax(value, "erf.forest_laimax");
        EXPECT_FALSE(error.empty());
        EXPECT_NE(error.find("erf.forest_laimax"), std::string::npos);
    }
}

TEST(ForestUtils, RejectsRawLonLatCoordinates)
{
    const auto error = erf_forest_utils::validate_cartesian_coordinates(
        false, false, true, true, "forest field 'LAI' in 'canopy.nc'");
    EXPECT_NE(error.find("lon/lat"), std::string::npos);
    EXPECT_NE(error.find("Cartesian"), std::string::npos);
    EXPECT_TRUE(erf_forest_utils::validate_cartesian_coordinates(
        true, true, false, false, "forest field 'LAI'").empty());
}

TEST(ForestDrag, OptionalFrontalAreaStorageWorksForDiscretePatches)
{
    // Per-process scratch name (see ERF_GTestTempDir.H): ctest -j runs this test and
    // erf_unit_tests_shuffle_repeat at once, and a fixed name in the working directory
    // let one rewrite the file the other was reading. Made on every rank.
    const std::string scratch_file = (erf_gtest_temp_path("erf_unit_forest_drag_discrete").string() + ".txt");
    const std::filesystem::path filename = scratch_file;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream output(filename);
        output << "1 1000 1000 1000 2000 0.15 2.0 0.8\n";
    }
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(15, 15, 15));
    const amrex::BoxArray ba(domain);
    const amrex::DistributionMapping dm(ba);
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {2000.0, 2000.0, 1000.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{1, 1, 0};
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());
    amrex::MultiFab z_cc(ba, dm, 1, 1);
    amrex::MultiFab z_nd(amrex::convert(ba, amrex::IntVect(1, 1, 1)), dm, 1, 1);
    fill_vertical_coordinates(z_cc, amrex::Real(0.5));
    fill_vertical_coordinates(z_nd, amrex::Real(0.0));
    amrex::Gpu::streamSynchronize();

    ForestDrag forest(filename.string());
    forest.define_drag_field(ba, dm, geom, &z_cc, &z_nd, false, 0);
    EXPECT_NE(forest.get_drag_field(), nullptr);
    EXPECT_EQ(forest.get_frontal_area(), nullptr);

    forest.define_drag_field(ba, dm, geom, &z_cc, &z_nd, true, 0);
    amrex::Gpu::streamSynchronize();
    ASSERT_NE(forest.get_frontal_area(), nullptr);
    const amrex::Real max_frontal_area = forest.get_frontal_area()->max(0, 1, false);
    EXPECT_TRUE(std::isfinite(static_cast<double>(max_frontal_area)));
    EXPECT_GT(max_frontal_area, amrex::Real(0.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::filesystem::remove(filename);
    }
    amrex::ParallelDescriptor::Barrier();
}

// erf-model/ERF#4632: where discrete patches overlap, each one's Cd a and a add.
// Before the fix the last patch overwrote the earlier ones.
TEST(ForestDrag, OverlappingDiscretePatchesAddTheirDrag)
{
    // A uniform-LAD patch and a Lalic--Mihailovic patch of a different height and
    // drag coefficient, overlapping over part of their footprints
    const std::string row_a = "1 450 500 300 400 0.15 2.0 0.8\n";
    const std::string row_b = "2 600 500 400 400 0.30 3.0 0.6\n";
    const auto a  = discrete_forest_fields(row_a, "erf_unit_forest_overlap_a");
    const auto b  = discrete_forest_fields(row_b, "erf_unit_forest_overlap_b");
    const auto ab = discrete_forest_fields(row_a + row_b, "erf_unit_forest_overlap_ab");

    // The overlap is not empty, so an overwrite cannot pass the sum check
    amrex::MultiFab both(a.drag.boxArray(), a.drag.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(both, a.drag, 0, 0, 1, 0);
    amrex::MultiFab::Multiply(both, b.drag, 0, 0, 1, 0);
    ASSERT_GT(both.max(0, 0, false), amrex::Real(0.0));

    const amrex::Real tol = sizeof(amrex::Real) == 8 ? amrex::Real(1e-14) : amrex::Real(1e-6);
    EXPECT_LE(max_sum_error(a.drag, b.drag, ab.drag), tol * ab.drag.norm0(0, 0, false));
    EXPECT_LE(max_sum_error(a.frontal_area, b.frontal_area, ab.frontal_area),
              tol * ab.frontal_area.norm0(0, 0, false));
}

// erf-model/ERF#4541: a type-2 patch whose top is exactly at a cell centre. The
// profile tends to zero at the tree top; before the fix that cell divided by
// hf - z = 0 and gave NaN drag (sqrt(inf) * exp(-inf)).
TEST(ForestDrag, TypeTwoTreeTopAtCellCentreHasZeroDrag)
{
    // Cell centres are at 31.25, 93.75, 156.25 m: the top is the k = 1 centre
    const auto f = discrete_forest_fields("2 500 500 93.75 400 0.2 2.0 0.4\n",
                                          "erf_unit_forest_tree_top");

    bool all_finite = true;
    amrex::Real drag_below = 0.0;
    amrex::Real drag_top = 0.0;
    amrex::Real area_top = 0.0;
    for (amrex::MFIter mfi(f.drag); mfi.isValid(); ++mfi) {
        const auto drag = f.drag.const_array(mfi);
        const auto area = f.frontal_area.const_array(mfi);
        const amrex::Box& bx = mfi.validbox();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            all_finite = all_finite && std::isfinite(static_cast<double>(drag(i, j, k)))
                                    && std::isfinite(static_cast<double>(area(i, j, k)));
        });
        // Column at the patch centre, x = y = 500 m (between cells 7 and 8)
        drag_below = drag(7, 7, 0);
        drag_top   = drag(7, 7, 1);
        area_top   = area(7, 7, 1);
    }
    // One box: only its owner saw it
    amrex::ParallelDescriptor::ReduceBoolAnd(all_finite);
    amrex::ParallelDescriptor::ReduceRealMax(drag_below);
    amrex::ParallelDescriptor::ReduceRealMax(drag_top);
    amrex::ParallelDescriptor::ReduceRealMax(area_top);
    EXPECT_TRUE(all_finite);
    EXPECT_GT(drag_below, amrex::Real(0.0));
    EXPECT_EQ(drag_top, amrex::Real(0.0));
    EXPECT_EQ(area_top, amrex::Real(0.0));
}
