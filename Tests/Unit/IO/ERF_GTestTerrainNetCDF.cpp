#ifdef ERF_USE_NETCDF

#include <cmath>
#include <cstdio>
#include <limits>
#include <string>

#include <AMReX_Arena.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include "ERF_NCInterface.H"
#include "ERF_MetgridUtils.H"
#include "ERF_ProbCommon.H"

namespace {

constexpr const char* terrain_filename = "erf_unit_terrain_endpoints.nc";
constexpr const char* terrain_time_filename = "erf_unit_terrain_time.nc";
constexpr const char* terrain_wps_filename = "erf_unit_terrain_geo_em.nc";
constexpr const char* metgrid_filename_0 = "erf_unit_metgrid_surface_0.nc";
constexpr const char* metgrid_filename_1 = "erf_unit_metgrid_surface_1.nc";
constexpr const char* metgrid_missing_psfc_filename = "erf_unit_metgrid_surface_missing_psfc.nc";
constexpr const char* metgrid_fill_psfc_filename = "erf_unit_metgrid_surface_fill_psfc.nc";
constexpr const char* metgrid_missing_value_psfc_filename =
    "erf_unit_metgrid_surface_missing_value_psfc.nc";

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

void
write_metgrid_surface_file (const char* filename,
                            const amrex::Real psfc,
                            const char* timestamp,
                            const bool include_psfc = true,
                            const bool include_tsk = false,
                            const amrex::Real hgt = 0.0,
                            const bool add_fill_value = false,
                            const bool add_missing_value = false,
                            const bool write_psfc = true)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) {
        return;
    }

    // The regression only needs classic NetCDF data types.  Keeping these
    // fixtures in the portable classic format also avoids requiring HDF5 file
    // locking support from the test environment.
    auto file = ncutils::NCFile::create(filename, NC_CLOBBER);
    file.enter_def_mode();
    file.def_dim("Time", 1);
    file.def_dim("DateStrLen", 19);
    file.def_dim("south_north", 1);
    file.def_dim("west_east", 1);
    file.def_var("Times", NC_CHAR, {"Time", "DateStrLen"});
    file.def_var("SST", ncutils::NCDType::Real,
                 {"Time", "south_north", "west_east"});
    if (include_psfc) {
        file.def_var("PSFC", ncutils::NCDType::Real,
                     {"Time", "south_north", "west_east"});
        if (add_fill_value) {
            file.var("PSFC").put_attr("_FillValue", std::vector<double>{psfc});
        }
        if (add_missing_value) {
            file.var("PSFC").put_attr("missing_value", std::vector<double>{psfc});
        }
    }
    if (include_tsk) {
        file.def_var("SKINTEMP", ncutils::NCDType::Real,
                     {"Time", "south_north", "west_east"});
    }
    file.def_var("HGT_M", ncutils::NCDType::Real,
                 {"Time", "south_north", "west_east"});
    file.put_attr("WEST-EAST_GRID_DIMENSION", std::vector<int>{2});
    file.put_attr("SOUTH-NORTH_GRID_DIMENSION", std::vector<int>{2});
    file.put_attr("DX", std::vector<double>{1.0});
    file.put_attr("DY", std::vector<double>{1.0});
    file.exit_def_mode();

    const amrex::Real sst = 290.0;
    file.var("SST").put(&sst);
    if (include_psfc && write_psfc) {
        file.var("PSFC").put(&psfc);
    }
    if (include_tsk) {
        file.var("SKINTEMP").put(&sst);
    }
    file.var("HGT_M").put(&hgt);

    int times_var = -1;
    ASSERT_EQ(nc_inq_varid(file.ncid, "Times", &times_var), NC_NOERR);
    ASSERT_EQ(nc_put_var_text(file.ncid, times_var, timestamp), NC_NOERR);
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

// Motivation: Metgrid SST/SKINTEMP is normalized with the pressure in the
// same forcing file. Reading PSFC only for itime 0 silently reused stale
// pressure at later forcing times, so two complete files with the same SST
// and different PSFC must produce different potential temperatures.
TEST(MetgridNetCDF, ReadsSurfacePressureForEveryForcingTime)
{
    write_metgrid_surface_file(metgrid_filename_0, 100000.0,
                                "2010-01-01_00:00:00");
    write_metgrid_surface_file(metgrid_filename_1, 90000.0,
                                "2010-01-01_01:00:00");
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{0, 0, 0};
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());

    amrex::FArrayBox xvel, yvel, temp, rhum, pres, ght, hgt, psfc;
    amrex::FArrayBox msfu, msfv, msfm, sst, tsk, lat, lon;
    amrex::IArrayBox lmask;
    std::string date_time;
    double epoch_time = 0.0;
    int flag_psfc = 0, flag_msf = 0, flag_sst = 0, flag_tsk = 0, flag_lmask = 0;
    int nc_nx = 0, nc_ny = 0;
    amrex::Real nc_dx = 0.0, nc_dy = 0.0;

    read_from_metgrid(0, 0, domain, metgrid_filename_0,
                      date_time, epoch_time, flag_psfc, flag_msf,
                      flag_sst, flag_tsk, flag_lmask, nc_nx, nc_ny,
                      nc_dx, nc_dy, xvel, yvel, temp, rhum, pres, ght,
                      hgt, psfc, msfu, msfv, msfm, sst, tsk, lat, lon,
                      lmask, geom);
    ASSERT_EQ(flag_psfc, 1);
    ASSERT_EQ(flag_sst, 1);
    ASSERT_FALSE(psfc.box().isEmpty());
    ASSERT_FALSE(sst.box().isEmpty());

    amrex::FArrayBox psfc_host_0(psfc.box(), 1, amrex::The_Pinned_Arena());
    amrex::FArrayBox sst_host_0(sst.box(), 1, amrex::The_Pinned_Arena());
    psfc_host_0.copy<amrex::RunOn::Device>(psfc);
    sst_host_0.copy<amrex::RunOn::Device>(sst);
    amrex::Gpu::streamSynchronize();
    const amrex::Real psfc_0 = psfc_host_0.const_array()(0, 0, 0);
    const amrex::Real sst_0 = sst_host_0.const_array()(0, 0, 0);

    read_from_metgrid(0, 1, domain, metgrid_filename_1,
                      date_time, epoch_time, flag_psfc, flag_msf,
                      flag_sst, flag_tsk, flag_lmask, nc_nx, nc_ny,
                      nc_dx, nc_dy, xvel, yvel, temp, rhum, pres, ght,
                      hgt, psfc, msfu, msfv, msfm, sst, tsk, lat, lon,
                      lmask, geom);
    ASSERT_EQ(flag_psfc, 1);
    ASSERT_EQ(flag_sst, 1);

    amrex::FArrayBox psfc_host_1(psfc.box(), 1, amrex::The_Pinned_Arena());
    amrex::FArrayBox sst_host_1(sst.box(), 1, amrex::The_Pinned_Arena());
    psfc_host_1.copy<amrex::RunOn::Device>(psfc);
    sst_host_1.copy<amrex::RunOn::Device>(sst);
    amrex::Gpu::streamSynchronize();
    const amrex::Real psfc_1 = psfc_host_1.const_array()(0, 0, 0);
    const amrex::Real sst_1 = sst_host_1.const_array()(0, 0, 0);

    // Independent oracle: theta = T / (p/p0)^(Rd/cp), not the conversion
    // routine used by init_state_from_metgrid.
    const amrex::Real expected_theta_0 = 290.0 *
        std::pow(psfc_0 / p_0, -RdoCp);
    const amrex::Real expected_theta_1 = 290.0 *
        std::pow(psfc_1 / p_0, -RdoCp);
    EXPECT_EQ(sst_0, sst_1);
    EXPECT_NE(psfc_0, psfc_1);
    EXPECT_NE(expected_theta_0, expected_theta_1);
    EXPECT_NEAR(expected_theta_0, 290.0 * std::pow(100000.0 / p_0, -RdoCp), 1.0e-10);
    EXPECT_NEAR(expected_theta_1, 290.0 * std::pow(90000.0 / p_0, -RdoCp), 1.0e-10);

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(metgrid_filename_0);
        std::remove(metgrid_filename_1);
    }
}

// Motivation: SST and SKINTEMP are valid Metgrid inputs even when PSFC is not
// present. The initialization must use the physical terrain height in the
// documented standard-atmosphere fallback, while debug_psfc and file PSFC
// retain their explicit precedence.
TEST(MetgridNetCDF, SurfacePressurePolicyHandlesMissingAndDebugPsfc)
{
    const amrex::Real z_sfc = amrex::Real(100.0);
    const amrex::Real T00 = amrex::Real(290.0);
    const amrex::Real P00 = p_0;
    const amrex::Real TLP = amrex::Real(50.0);
    const amrex::Real toa = T00 / TLP;
    const amrex::Real expected_missing = P00 * std::exp(
        -toa + std::sqrt(toa*toa - amrex::Real(2.0) * CONST_GRAV * z_sfc /
                         (TLP * R_d)));

    amrex::Real pressure = 0.0;
    ASSERT_TRUE(metgrid_surface_pressure(false, 0, 0.0, z_sfc,
                                         P00, T00, TLP, pressure));
    EXPECT_NEAR(pressure, expected_missing, amrex::Real(1.0e-12) * expected_missing);
    ASSERT_TRUE(metgrid_surface_pressure(true, 1, amrex::Real(90000.0), z_sfc,
                                         P00, T00, TLP, pressure));
    EXPECT_EQ(pressure, amrex::Real(100000.0));
    ASSERT_TRUE(metgrid_surface_pressure(false, 1, amrex::Real(90000.0), z_sfc,
                                         P00, T00, TLP, pressure));
    EXPECT_EQ(pressure, amrex::Real(90000.0));

    EXPECT_FALSE(metgrid_surface_pressure(false, 1,
                                          std::numeric_limits<amrex::Real>::quiet_NaN(),
                                          z_sfc, P00, T00, TLP, pressure));
    EXPECT_FALSE(metgrid_surface_pressure(false, 1, amrex::Real(0.0),
                                          z_sfc, P00, T00, TLP, pressure));
    EXPECT_FALSE(metgrid_surface_pressure(false, 1, amrex::Real(-1.0),
                                          z_sfc, P00, T00, TLP, pressure));
    EXPECT_TRUE(metgrid_surface_pressure(true, 1,
                                         std::numeric_limits<amrex::Real>::quiet_NaN(),
                                         z_sfc, P00, T00, TLP, pressure));
    EXPECT_EQ(pressure, amrex::Real(100000.0));

    EXPECT_FALSE(metgrid_surface_pressure(false, 0, 0.0, z_sfc,
                                          P00, T00, amrex::Real(0.0), pressure));
    EXPECT_FALSE(metgrid_surface_pressure(false, 0, 0.0,
                                          std::numeric_limits<amrex::Real>::quiet_NaN(),
                                          P00, T00, TLP, pressure));
    // Motivation: the analytic standard-atmosphere inversion assumes a
    // positive tropospheric lapse rate. A negative value would otherwise
    // silently select the wrong branch of the quadratic pressure relation.
    EXPECT_FALSE(metgrid_surface_pressure(false, 0, 0.0, z_sfc,
                                          P00, T00, amrex::Real(-1.0), pressure));
}

// Motivation: the pressure-selection policy is only useful if the production
// Metgrid normalization applies it to the actual temperature field. Exercise
// the GPU-safe conversion primitive with file, analytic, and debug pressures
// so the theta contract cannot regress to an unconverted absolute temperature.
TEST(MetgridNetCDF, SurfaceTemperatureNormalizationUsesPressurePolicy)
{
    const amrex::Real temperature = amrex::Real(290.0);
    const amrex::Real P00 = p_0;
    const amrex::Real T00 = amrex::Real(290.0);
    const amrex::Real TLP = amrex::Real(50.0);
    const amrex::Real z_sfc = amrex::Real(100.0);
    const amrex::Real tolerance = amrex::Real(64.0) *
        std::numeric_limits<amrex::Real>::epsilon() * temperature;

    // Independent oracle: theta = T * (p0/p_s)^(Rd/cp), evaluated at the
    // pressure values the production policy is expected to select.
    const amrex::Real expected_at_100k = temperature *
        std::pow(P00 / amrex::Real(100000.0), RdoCp);
    const amrex::Real expected_at_90k = temperature *
        std::pow(P00 / amrex::Real(90000.0), RdoCp);
    amrex::Real theta = 0.0;
    ASSERT_TRUE(metgrid_surface_theta(
        temperature, false, 1, amrex::Real(100000.0), z_sfc,
        P00, T00, TLP, RdoCp, theta));
    EXPECT_NEAR(theta, expected_at_100k, tolerance);
    ASSERT_TRUE(metgrid_surface_theta(
        temperature, false, 1, amrex::Real(90000.0), z_sfc,
        P00, T00, TLP, RdoCp, theta));
    EXPECT_NEAR(theta, expected_at_90k, tolerance);
    EXPECT_NE(expected_at_100k, expected_at_90k);

    const amrex::Real toa = T00 / TLP;
    const amrex::Real expected_missing_pressure = P00 * std::exp(
        -toa + std::sqrt(toa*toa - amrex::Real(2.0) * CONST_GRAV * z_sfc /
                         (TLP * R_d)));
    const amrex::Real expected_missing = temperature *
        std::pow(P00 / expected_missing_pressure, RdoCp);
    ASSERT_TRUE(metgrid_surface_theta(
        temperature, false, 0, amrex::Real(0.0), z_sfc,
        P00, T00, TLP, RdoCp, theta));
    EXPECT_NEAR(theta, expected_missing, tolerance);

    // Debug pressure deliberately wins over the file value.
    ASSERT_TRUE(metgrid_surface_theta(
        temperature, true, 1, amrex::Real(90000.0), z_sfc,
        P00, T00, TLP, RdoCp, theta));
    EXPECT_NEAR(theta, expected_at_100k, tolerance);

    EXPECT_FALSE(metgrid_surface_theta(
        std::numeric_limits<amrex::Real>::quiet_NaN(), false, 1,
        amrex::Real(90000.0), z_sfc, P00, T00, TLP, RdoCp, theta));
    EXPECT_FALSE(metgrid_surface_theta(
        temperature, false, 1, amrex::Real(90000.0), z_sfc,
        P00, T00, TLP, amrex::Real(0.0), theta));
}

// Motivation: the old Metgrid reader rejected an SST-only forcing file before
// the analytic pressure fallback could run. Keep this test at the actual
// NetCDF ingestion boundary so missing PSFC remains an accepted input and the
// current SST record is still made available to initialization.
TEST(MetgridNetCDF, AcceptsSurfaceTemperatureWithoutSurfacePressure)
{
    write_metgrid_surface_file(metgrid_missing_psfc_filename, 0.0,
                               "2010-01-01_00:00:00", false, true, 100.0);
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{0, 0, 0};
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());

    amrex::FArrayBox xvel, yvel, temp, rhum, pres, ght, hgt, psfc;
    amrex::FArrayBox msfu, msfv, msfm, sst, tsk, lat, lon;
    amrex::IArrayBox lmask;
    std::string date_time;
    double epoch_time = 0.0;
    int flag_psfc = 0, flag_msf = 0, flag_sst = 0, flag_tsk = 0, flag_lmask = 0;
    int nc_nx = 0, nc_ny = 0;
    amrex::Real nc_dx = 0.0, nc_dy = 0.0;

    read_from_metgrid(0, 0, domain, metgrid_missing_psfc_filename,
                      date_time, epoch_time, flag_psfc, flag_msf,
                      flag_sst, flag_tsk, flag_lmask, nc_nx, nc_ny,
                      nc_dx, nc_dy, xvel, yvel, temp, rhum, pres, ght,
                      hgt, psfc, msfu, msfv, msfm, sst, tsk, lat, lon,
                      lmask, geom);

    EXPECT_EQ(flag_psfc, 0);
    EXPECT_EQ(flag_sst, 1);
    EXPECT_EQ(flag_tsk, 1);
    ASSERT_TRUE(psfc.box().isEmpty());
    ASSERT_FALSE(sst.box().isEmpty());

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(metgrid_missing_psfc_filename);
    }
}

// Motivation: NetCDF may encode a present PSFC field with a declared,
// missing-value, or type-default fill sentinel. Those cells must remain an
// invalid provided pressure so the initialization path fails loudly rather
// than silently selecting the analytic missing-variable fallback.
TEST(MetgridNetCDF, CanonicalizesDeclaredAndDefaultPsfcFillValues)
{
    write_metgrid_surface_file(metgrid_fill_psfc_filename,
                               static_cast<amrex::Real>(NC_FILL_DOUBLE),
                               "2010-01-01_00:00:00", true, false, 0.0,
                               false, false, true);
    write_metgrid_surface_file(metgrid_missing_value_psfc_filename,
                               amrex::Real(123456.0),
                               "2010-01-01_00:00:00", true, false, 0.0,
                               false, true, true);
    amrex::ParallelDescriptor::Barrier();

    const amrex::Box domain(amrex::IntVect(0, 0, 0), amrex::IntVect(0, 0, 0));
    const amrex::RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    const amrex::Array<int, AMREX_SPACEDIM> periodic{0, 0, 0};
    amrex::Geometry geom(domain, &real_box, amrex::CoordSys::cartesian, periodic.data());

    const auto read_psfc = [&] (const char* filename) -> amrex::Real {
        amrex::FArrayBox xvel, yvel, temp, rhum, pres, ght, hgt, psfc;
        amrex::FArrayBox msfu, msfv, msfm, sst, tsk, lat, lon;
        amrex::IArrayBox lmask;
        std::string date_time;
        double epoch_time = 0.0;
        int flag_psfc = 0, flag_msf = 0, flag_sst = 0, flag_tsk = 0, flag_lmask = 0;
        int nc_nx = 0, nc_ny = 0;
        amrex::Real nc_dx = 0.0, nc_dy = 0.0;
        read_from_metgrid(0, 0, domain, filename,
                          date_time, epoch_time, flag_psfc, flag_msf,
                          flag_sst, flag_tsk, flag_lmask, nc_nx, nc_ny,
                          nc_dx, nc_dy, xvel, yvel, temp, rhum, pres, ght,
                          hgt, psfc, msfu, msfv, msfm, sst, tsk, lat, lon,
                          lmask, geom);
        EXPECT_EQ(flag_psfc, 1);
        if (psfc.box().isEmpty()) {
            ADD_FAILURE() << "PSFC unexpectedly has an empty box";
            return std::numeric_limits<amrex::Real>::quiet_NaN();
        }
        amrex::FArrayBox host(psfc.box(), 1, amrex::The_Pinned_Arena());
        host.copy<amrex::RunOn::Device>(psfc);
        amrex::Gpu::streamSynchronize();
        return host.const_array()(0, 0, 0);
    };

    EXPECT_FALSE(std::isfinite(read_psfc(metgrid_fill_psfc_filename)));
    EXPECT_FALSE(std::isfinite(read_psfc(metgrid_missing_value_psfc_filename)));

    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(metgrid_fill_psfc_filename);
        std::remove(metgrid_missing_value_psfc_filename);
    }
}

#endif
