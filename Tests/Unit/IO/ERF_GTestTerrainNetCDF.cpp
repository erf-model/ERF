#ifdef ERF_USE_NETCDF

#include <cstdio>
#include <string>

#include <AMReX_Arena.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_NCInterface.H"
#include "ERF_ProbCommon.H"

namespace {

constexpr const char* terrain_filename = "erf_unit_terrain_endpoints.nc";

void
write_terrain_file ()
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    auto file = ncutils::NCFile::create(terrain_filename, NC_CLOBBER | NC_NETCDF4);
    file.enter_def_mode();
    file.def_dim("x", 3);
    file.def_dim("y", 3);
    file.def_var("x", ncutils::NCDType::Real, {"x"});
    file.def_var("y", ncutils::NCDType::Real, {"y"});
    file.def_var("height", ncutils::NCDType::Real, {"y", "x"});
    file.exit_def_mode();

    const amrex::Vector<amrex::Real> x{0.0, 1.0, 2.0};
    const amrex::Vector<amrex::Real> y{0.0, 1.0, 2.0};
    const amrex::Vector<amrex::Real> height{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0};
    file.var("x").put(x.data());
    file.var("y").put(y.data());
    file.var("height").put(height.data());
    file.close();
}

struct TerrainSamples
{
    amrex::Real x_upper;
    amrex::Real y_upper;
    amrex::Real upper_corner;
    amrex::Real outside_corner;
};

TerrainSamples
read_terrain_samples (int horizontal_cells)
{
    const amrex::Box domain(
        amrex::IntVect(0, 0, 0),
        amrex::IntVect(horizontal_cells - 1, horizontal_cells - 1, 0));
    const amrex::RealBox real_box(
        {0.0, 0.0, 0.0},
        {static_cast<amrex::Real>(horizontal_cells),
         static_cast<amrex::Real>(horizontal_cells), 1.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{0, 0, 0};
    const amrex::Geometry geometry(
        domain, &real_box, amrex::CoordSys::cartesian, periodic.data());
    const amrex::Box node_box = amrex::convert(domain, amrex::IntVect(1, 1, 0));

    amrex::FArrayBox terrain(node_box, 1);
    terrain.setVal(-999.0);
    ProblemBase problem;
    problem.read_terrain_netcdf(terrain_filename, geometry, terrain, 0.0);
    amrex::Gpu::streamSynchronize();

    amrex::FArrayBox host_terrain(node_box, 1, amrex::The_Pinned_Arena());
    amrex::Gpu::copy(amrex::Gpu::deviceToHost,
                     terrain.dataPtr(), terrain.dataPtr() + terrain.size(),
                     host_terrain.dataPtr());
    amrex::Gpu::streamSynchronize();

    const auto values = host_terrain.const_array();
    const int endpoint = 2;
    const int outside = horizontal_cells;
    return TerrainSamples{
        values(endpoint, 0, 0),
        values(0, endpoint, 0),
        values(endpoint, endpoint, 0),
        values(outside, outside, 0)};
}

} // namespace

// Motivation: a terrain file covering the ERF domain used to lose its final
// row and column because the upper coordinate selected an invalid lower
// stencil index. A larger second domain checks that real extrapolation remains
// disabled and produces the documented zero value.
TEST(TerrainNetCDF, PreservesUpperRowColumnAndCornerWithoutExtrapolation)
{
    write_terrain_file();
    amrex::ParallelDescriptor::Barrier();

    const auto endpoint_samples = read_terrain_samples(2);
    EXPECT_EQ(endpoint_samples.x_upper, amrex::Real(3.0));
    EXPECT_EQ(endpoint_samples.y_upper, amrex::Real(7.0));
    EXPECT_EQ(endpoint_samples.upper_corner, amrex::Real(9.0));

    const auto outside_samples = read_terrain_samples(3);
    EXPECT_EQ(outside_samples.upper_corner, amrex::Real(9.0));
    EXPECT_EQ(outside_samples.outside_corner, amrex::Real(0.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(terrain_filename);
    }
}

#endif
