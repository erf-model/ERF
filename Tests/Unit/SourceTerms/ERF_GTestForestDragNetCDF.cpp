#include <cmath>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>

#include <gtest/gtest.h>

#include "ERF_ForestDrag.H"

#ifdef ERF_USE_NETCDF
namespace {

#ifndef ERF_FOREST_TEST_FIXTURE_DIR
#define ERF_FOREST_TEST_FIXTURE_DIR "../test_files/BellForest"
#endif

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

} // namespace

TEST(ForestDrag, TypeTwoProfileProducesFinitePositiveFrontalArea)
{
    const std::string fixture_dir = ERF_FOREST_TEST_FIXTURE_DIR;
    ForestDrag forest(fixture_dir + "/forest_lai_bell.nc",
                      fixture_dir + "/forest_height_bell.nc",
                      0.15, 2, 0.6);

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

    forest.define_drag_field(ba, dm, geom, &z_cc, &z_nd, true, 0);
    amrex::Gpu::streamSynchronize();
    ASSERT_NE(forest.get_frontal_area(), nullptr);
    const amrex::Real max_frontal_area = forest.get_frontal_area()->max(0, 1, false);
    EXPECT_TRUE(std::isfinite(static_cast<double>(max_frontal_area)));
    EXPECT_GT(max_frontal_area, amrex::Real(0.0));
}

TEST(ForestDrag, ConstantCdModeDoesNotRequireAConstantCdGrid)
{
    const std::string fixture_dir = ERF_FOREST_TEST_FIXTURE_DIR;
    ForestDrag forest(fixture_dir + "/forest_lai_bell.nc",
                      fixture_dir + "/forest_height_bell.nc", 0.15, 1, 0.6);
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
    forest.define_drag_field(ba, dm, geom, &z_cc, &z_nd, false, 0);
    EXPECT_NE(forest.get_drag_field(), nullptr);
    EXPECT_EQ(forest.get_frontal_area(), nullptr);
}

#endif
