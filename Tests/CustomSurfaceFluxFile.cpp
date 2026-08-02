#include <netcdf.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void check_nc (const int status, const std::string& context)
{
    if (status != NC_NOERR) {
        std::cerr << context << ": " << nc_strerror(status) << '\n';
        std::exit(2);
    }
}

double ustar_value (const int lev, const int it, const int j, const int i)
{
    return 0.20 + 0.01 * lev + 0.05 * it + 0.002 * j + 0.001 * i;
}

double theta_value (const int lev, const int it, const int j, const int i)
{
    return 0.010 * (lev + 1) + 0.004 * it + 0.0002 * j + 0.0001 * i;
}

double qv_value (const int lev, const int it, const int j, const int i)
{
    return 0.00010 * (lev + 1) + 0.00004 * it +
           0.000002 * j + 0.000001 * i;
}

void generate (const std::string& filename, const int nx, const int ny, const int lev)
{
    int ncid;
    check_nc(nc_create(filename.c_str(), NC_CLOBBER | NC_NETCDF4, &ncid),
             "create forcing file");

    const int schema_version = 1;
    check_nc(nc_put_att_int(ncid, NC_GLOBAL, "schema_version", NC_INT, 1,
                            &schema_version), "write schema_version");

    int time_dim, y_dim, x_dim;
    check_nc(nc_def_dim(ncid, "time", 3, &time_dim), "define time");
    check_nc(nc_def_dim(ncid, "y", ny, &y_dim), "define y");
    check_nc(nc_def_dim(ncid, "x", nx, &x_dim), "define x");

    int time_var;
    check_nc(nc_def_var(ncid, "time", NC_DOUBLE, 1, &time_dim, &time_var),
             "define time variable");
    const std::string time_units = "seconds since simulation start";
    check_nc(nc_put_att_text(ncid, time_var, "units", time_units.size(),
                             time_units.c_str()), "write time units");

    const int dims[3] = {time_dim, y_dim, x_dim};
    auto define_field = [&] (const char* name, const std::string& units) {
        int varid;
        check_nc(nc_def_var(ncid, name, NC_DOUBLE, 3, dims, &varid),
                 std::string("define ") + name);
        check_nc(nc_put_att_text(ncid, varid, "units", units.size(), units.c_str()),
                 std::string("write units for ") + name);
        return varid;
    };
    const int ustar_var = define_field("ustar", "m s-1");
    const int theta_var = define_field("theta_flux", "K m s-1");
    const int qv_var = define_field("qv_flux", "kg kg-1 m s-1");
    check_nc(nc_enddef(ncid), "end define mode");

    const double times[3] = {0.0, 2.0, 4.0};
    check_nc(nc_put_var_double(ncid, time_var, times), "write time");
    const std::size_t plane = static_cast<std::size_t>(ny) * nx;
    std::vector<double> ustar(3 * plane);
    std::vector<double> theta(3 * plane);
    std::vector<double> qv(3 * plane);
    for (int it = 0; it < 3; ++it) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const std::size_t n = (static_cast<std::size_t>(it) * ny + j) * nx + i;
                ustar[n] = ustar_value(lev, it, j, i);
                theta[n] = theta_value(lev, it, j, i);
                qv[n] = qv_value(lev, it, j, i);
            }
        }
    }
    check_nc(nc_put_var_double(ncid, ustar_var, ustar.data()), "write ustar");
    check_nc(nc_put_var_double(ncid, theta_var, theta.data()), "write theta_flux");
    check_nc(nc_put_var_double(ncid, qv_var, qv.data()), "write qv_flux");
    check_nc(nc_close(ncid), "close forcing file");
}

std::vector<double> read_var (const int ncid, const char* name, const std::size_t count)
{
    int varid;
    check_nc(nc_inq_varid(ncid, name, &varid), std::string("find ") + name);
    std::vector<double> values(count);
    check_nc(nc_get_var_double(ncid, varid, values.data()), std::string("read ") + name);
    return values;
}

void check_output (const std::string& filename, const int expected_lev,
                   const int expected_nx, const int expected_ny)
{
    int ncid;
    check_nc(nc_open(filename.c_str(), NC_NOWRITE, &ncid), "open output " + filename);
    int level = -1;
    check_nc(nc_get_att_int(ncid, NC_GLOBAL, "CurrentLevel", &level),
             "read CurrentLevel");
    if (level != expected_lev) {
        std::cerr << filename << " has CurrentLevel=" << level
                  << ", expected " << expected_lev << '\n';
        std::exit(3);
    }

    auto dim_length = [&] (const char* name) {
        int dimid;
        std::size_t length;
        check_nc(nc_inq_dimid(ncid, name, &dimid), std::string("find dimension ") + name);
        check_nc(nc_inq_dimlen(ncid, dimid, &length), std::string("read dimension ") + name);
        return length;
    };
    const std::size_t nx = dim_length("nx");
    const std::size_t ny = dim_length("ny");
    const std::size_t nz = dim_length("nz");
    if (nx != static_cast<std::size_t>(expected_nx) ||
        ny != static_cast<std::size_t>(expected_ny) || nz != 1) {
        std::cerr << filename << " has unexpected 2D shape\n";
        std::exit(3);
    }

    const std::size_t count = nx * ny;
    const auto ustar = read_var(ncid, "u_star", count);
    const auto theta = read_var(ncid, "t_star", count);
    const auto qv = read_var(ncid, "q_star", count);
    const auto sensible = read_var(ncid, "sens_flux", count);
    const auto latent = read_var(ncid, "laten_flux", count);

    const double tolerance = 2.0e-12;
    std::size_t active_cells = 0;
    for (int j = 0; j < expected_ny; ++j) {
        for (int i = 0; i < expected_nx; ++i) {
            const std::size_t n = static_cast<std::size_t>(j) * expected_nx + i;
            // Fine-level NetCDF variables span the full level geometry but
            // cells outside the active BoxArray retain the NetCDF fill value.
            if (!std::isfinite(ustar[n]) || std::abs(ustar[n]) > 1.0e20) {
                continue;
            }
            ++active_cells;
            // The output is at t=1 s, halfway between forcing records 0 and 2 s.
            const double expected_u = 0.5 * (ustar_value(expected_lev, 0, j, i) +
                                             ustar_value(expected_lev, 1, j, i));
            const double expected_t = 0.5 * (theta_value(expected_lev, 0, j, i) +
                                             theta_value(expected_lev, 1, j, i));
            const double expected_q = 0.5 * (qv_value(expected_lev, 0, j, i) +
                                             qv_value(expected_lev, 1, j, i));
            if (std::abs(ustar[n] - expected_u) > tolerance ||
                std::abs(theta[n] - expected_t) > tolerance ||
                std::abs(qv[n] - expected_q) > tolerance) {
                std::cerr << filename << " forcing mismatch at (" << i << ',' << j << ")\n";
                std::exit(4);
            }
            if (!std::isfinite(sensible[n]) || !std::isfinite(latent[n]) ||
                std::abs(sensible[n] / theta[n] - latent[n] / qv[n]) > 2.0e-10) {
                std::cerr << filename << " conservative scalar flux mismatch at ("
                          << i << ',' << j << ")\n";
                std::exit(5);
            }
        }
    }
    if ((expected_lev == 0 && active_cells != count) ||
        (expected_lev > 0 && (active_cells == 0 || active_cells == count))) {
        std::cerr << filename << " did not contain the expected partial fine patch\n";
        std::exit(6);
    }
    check_nc(nc_close(ncid), "close output " + filename);
}

} // namespace

int main (int argc, char** argv)
{
    if (argc == 6 && std::string(argv[1]) == "generate") {
        generate(argv[2], std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
        return 0;
    }
    if (argc == 4 && std::string(argv[1]) == "check") {
        check_output(argv[2], 0, 8, 8);
        check_output(argv[3], 1, 16, 16);
        return 0;
    }
    std::cerr << "usage: " << argv[0]
              << " generate FILE NX NY LEVEL | check LEVEL0_NC LEVEL1_NC\n";
    return 1;
}
