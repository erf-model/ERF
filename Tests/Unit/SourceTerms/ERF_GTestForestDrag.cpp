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
    const std::filesystem::path filename = "erf_unit_forest_drag_discrete.txt";
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
