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
constexpr const char* terrain_time_filename = "erf_unit_terrain_time.nc";
constexpr const char* terrain_wps_filename = "erf_unit_terrain_geo_em.nc";

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

void
write_wps_terrain_file ()
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    auto file = ncutils::NCFile::create(terrain_wps_filename, NC_CLOBBER | NC_NETCDF4);
    file.enter_def_mode();
    file.def_dim("west_east", 2);
    file.def_dim("south_north", 2);
    file.def_dim("Time", 2);
    file.def_var("HGT_M", ncutils::NCDType::Real,
                 {"Time", "south_north", "west_east"});
    file.def_var("XLAT_M", ncutils::NCDType::Real,
                 {"Time", "south_north", "west_east"});
    file.def_var("XLONG_M", ncutils::NCDType::Real,
                 {"Time", "south_north", "west_east"});
    file.exit_def_mode();

    // MAP_PROJ = 1 is Lambert conformal, the usual WPS default.  It must be one
    // of the projections the reader supports (1 Lambert, 2 polar, 3 Mercator):
    // projection 6 (latitude/longitude) is deliberately rejected because WPS
    // writes DX/DY in degrees for it, while everything downstream consumes them
    // as metres.  Lambert additionally requires TRUELAT2 alongside TRUELAT1.
    file.put_attr("MAP_PROJ", std::vector<int>{1});
    file.put_attr("CEN_LAT", std::vector<double>{40.0});
    file.put_attr("CEN_LON", std::vector<double>{-105.0});
    file.put_attr("STAND_LON", std::vector<double>{-105.0});
    file.put_attr("TRUELAT1", std::vector<double>{30.0});
    file.put_attr("TRUELAT2", std::vector<double>{60.0});
    file.put_attr("DX", std::vector<double>{1.0});
    file.put_attr("DY", std::vector<double>{1.0});
    file.put_attr("WEST-EAST_GRID_DIMENSION", std::vector<int>{3});
    file.put_attr("SOUTH-NORTH_GRID_DIMENSION", std::vector<int>{3});

    const amrex::Vector<amrex::Real> height{
        1.0, 2.0,
        3.0, 4.0,
        101.0, 102.0,
        103.0, 104.0};
    const amrex::Vector<amrex::Real> lat{
        40.0, 40.0,
        41.0, 41.0,
        40.0, 40.0,
        41.0, 41.0};
    const amrex::Vector<amrex::Real> lon{
        -105.0, -104.0,
        -105.0, -104.0,
        -105.0, -104.0,
        -105.0, -104.0};
    file.var("HGT_M").put(height.data());
    file.var("XLAT_M").put(lat.data());
    file.var("XLONG_M").put(lon.data());
    file.close();
}

void
write_time_leading_terrain_file ()
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    auto file = ncutils::NCFile::create(terrain_time_filename, NC_CLOBBER | NC_NETCDF4);
    file.enter_def_mode();
    file.def_dim("x", 3);
    file.def_dim("y", 3);
    file.def_dim("Time", 2);
    file.def_var("x", ncutils::NCDType::Real, {"x"});
    file.def_var("y", ncutils::NCDType::Real, {"y"});
    file.def_var("height", ncutils::NCDType::Real,
                 {"Time", "y", "x"});
    file.exit_def_mode();

    const amrex::Vector<amrex::Real> x{0.0, 1.0, 2.0};
    const amrex::Vector<amrex::Real> y{0.0, 1.0, 2.0};
    const amrex::Vector<amrex::Real> height{
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 9.0,
        101.0, 102.0, 103.0,
        104.0, 105.0, 106.0,
        107.0, 108.0, 109.0};
    file.var("x").put(x.data());
    file.var("y").put(y.data());
    file.var("height").put(height.data());
    file.close();
}

struct TerrainSamples
{
    amrex::Real lower_corner;
    amrex::Real interior;
    amrex::Real bottom_edge;
    amrex::Real left_edge;
    amrex::Real right_edge;
    amrex::Real top_edge;
    amrex::Real x_upper;
    amrex::Real y_upper;
    amrex::Real upper_corner;
    amrex::Real outside_corner;
};

TerrainSamples
read_terrain_samples (const char* filename, int horizontal_cells)
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
    terrain.template setVal<amrex::RunOn::Device>(amrex::Real(-999.0));
    ProblemBase problem;
    problem.read_terrain_netcdf(filename, geometry, terrain, 0.0);
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
        values(0, 0, 0),
        values(1, 1, 0),
        values(1, 0, 0),
        values(0, 1, 0),
        values(endpoint, 1, 0),
        values(1, endpoint, 0),
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

    const auto endpoint_samples = read_terrain_samples(terrain_filename, 2);
    EXPECT_EQ(endpoint_samples.lower_corner, amrex::Real(1.0));
    EXPECT_EQ(endpoint_samples.x_upper, amrex::Real(3.0));
    EXPECT_EQ(endpoint_samples.y_upper, amrex::Real(7.0));
    EXPECT_EQ(endpoint_samples.upper_corner, amrex::Real(9.0));

    const auto outside_samples = read_terrain_samples(terrain_filename, 3);
    EXPECT_EQ(outside_samples.upper_corner, amrex::Real(9.0));
    EXPECT_EQ(outside_samples.outside_corner, amrex::Real(0.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(terrain_filename);
    }
}

TEST(TerrainNetCDF, ReadsFirstRecordOfExplicitTimeLeadingField)
{
    write_time_leading_terrain_file();
    amrex::ParallelDescriptor::Barrier();

    const auto samples = read_terrain_samples(terrain_time_filename, 2);
    EXPECT_EQ(samples.lower_corner, amrex::Real(1.0));
    EXPECT_EQ(samples.x_upper, amrex::Real(3.0));
    EXPECT_EQ(samples.y_upper, amrex::Real(7.0));
    EXPECT_EQ(samples.upper_corner, amrex::Real(9.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(terrain_time_filename);
    }
}

TEST(TerrainNetCDF, ReadsGenuineWpsGeoEmMassField)
{
    write_wps_terrain_file();
    amrex::ParallelDescriptor::Barrier();

    const auto samples = read_terrain_samples(terrain_wps_filename, 2);
    EXPECT_EQ(samples.lower_corner, amrex::Real(1.0));
    EXPECT_EQ(samples.interior, amrex::Real(2.5));
    EXPECT_EQ(samples.bottom_edge, amrex::Real(1.5));
    EXPECT_EQ(samples.left_edge, amrex::Real(2.0));
    EXPECT_EQ(samples.right_edge, amrex::Real(3.0));
    EXPECT_EQ(samples.top_edge, amrex::Real(3.5));
    EXPECT_EQ(samples.x_upper, amrex::Real(2.0));
    EXPECT_EQ(samples.y_upper, amrex::Real(3.0));
    EXPECT_EQ(samples.upper_corner, amrex::Real(4.0));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(terrain_wps_filename);
    }
}

#endif
