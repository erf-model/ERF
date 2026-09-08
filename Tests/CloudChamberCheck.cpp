#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "../Source/Utils/ERF_EOS.H"
#include "../Source/Utils/ERF_MicrophysicsUtils.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::PlotFileData;
using amrex::Real;
using amrex::TilingIfNotGPU;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

int fail (const std::string& message)
{
    std::cerr << "CloudChamberCheck: " << message << "\n";
    return 1;
}

struct BudgetRow {
    std::string scalar;
    std::array<Real, 2 * AMREX_SPACEDIM> faces{};
    Real net_boundary = Real(0.0);
    Real volume_change = Real(0.0);
    Real internal_source = Real(0.0);
    Real residual = Real(0.0);
    Real tolerance = Real(0.0);
    std::string status;
};

bool read_budget (const std::string& path, std::vector<BudgetRow>& rows,
                  std::string& error)
{
    std::ifstream in(path);
    if (!in) {
        error = "cannot open budget file " + path;
        return false;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream fields(line);
        int start_step = 0;
        int end_step = 0;
        double start_time = 0.0;
        double end_time = 0.0;
        std::string units;
        BudgetRow row;
        if (!(fields >> start_step >> start_time >> end_step >> end_time
                    >> row.scalar >> units)) {
            error = "malformed budget row: " + line;
            return false;
        }
        for (auto& value : row.faces) {
            if (!(fields >> value)) {
                error = "malformed budget face data: " + line;
                return false;
            }
        }
        if (!(fields >> row.net_boundary >> row.volume_change >> row.internal_source
                    >> row.residual >> row.tolerance >> row.status)) {
            error = "malformed budget closure data: " + line;
            return false;
        }
        rows.push_back(row);
    }
    return true;
}

Real budget_tolerance (const BudgetRow& row)
{
    return std::max(row.tolerance,
                    Real(128.0) * std::numeric_limits<Real>::epsilon());
}

bool close_to_zero (Real value, const BudgetRow& row)
{
    return std::abs(value) <= budget_tolerance(row);
}

bool is_budget_mode (const std::string& mode)
{
    return mode == "all_dry" || mode == "wet_budget" || mode == "thermal_budget";
}

bool is_checker_mode (const std::string& mode)
{
    return mode == "dry" || mode == "cloudy" || mode == "parity" || is_budget_mode(mode);
}

struct BudgetSummary {
    int rows = 0;
    int thermal_rows = 0;
    int total_rows = 0;
    int vapor_rows = 0;
    int cloud_rows = 0;
    Real max_residual_ratio = Real(0.0);
};

Real integral (PlotFileData& plotfile, const std::string& scalar)
{
    MultiFab rho = plotfile.get(0, "density");
    MultiFab value = plotfile.get(0, scalar == "rhoTheta" ? "theta" : scalar);
    MultiFab::Multiply(value, rho, 0, 0, 1, 0);
    const auto dx = plotfile.cellSize(0);
    const Real volume = dx[0] * dx[1] * dx[2];
    return value.sum(0, true) * volume;
}

bool check_budget_rows (const std::vector<BudgetRow>& rows,
                        const std::string& mode, BudgetSummary& summary,
                        std::string& error)
{
    const bool all_dry = mode == "all_dry";
    const bool wet = mode == "wet_budget";
    const bool thermal = mode == "thermal_budget";
    summary = {};
    summary.rows = static_cast<int>(rows.size());
    int total_rows = 0;
    int vapor_rows = 0;
    int cloud_rows = 0;
    int thermal_rows = 0;
    const auto check_closure = [&](const BudgetRow& row, const char* name) {
        if (!std::isfinite(static_cast<double>(row.residual)) ||
            !std::isfinite(static_cast<double>(row.tolerance)) ||
            row.tolerance < Real(0.0)) {
            error = std::string(name) + " budget closure is non-finite or has negative tolerance";
            return false;
        }
        const Real effective_tolerance = budget_tolerance(row);
        summary.max_residual_ratio = std::max(
            summary.max_residual_ratio, std::abs(row.residual) / effective_tolerance);
        if (row.status != "PASS" || std::abs(row.residual) > effective_tolerance) {
            error = std::string(name) + " budget did not PASS";
            return false;
        }
        return true;
    };
    for (const auto& row : rows) {
        const Real tol = budget_tolerance(row);
        if (row.scalar == "rhoTheta") {
            if (thermal) {
                ++thermal_rows;
                bool finite = true;
                for (const auto value : row.faces) {
                    finite = finite && std::isfinite(static_cast<double>(value));
                }
                finite = finite && std::isfinite(static_cast<double>(row.net_boundary)) &&
                    std::isfinite(static_cast<double>(row.volume_change)) &&
                    std::isfinite(static_cast<double>(row.internal_source)) &&
                    std::isfinite(static_cast<double>(row.residual)) &&
                    std::isfinite(static_cast<double>(row.tolerance));
                if (!finite || row.tolerance < Real(0.0) || row.status != "PASS" ||
                    std::abs(row.residual) > tol) {
                    error = "dry thermal rhoTheta budget did not PASS";
                    return false;
                }
                summary.max_residual_ratio = std::max(
                    summary.max_residual_ratio, std::abs(row.residual) / tol);
            }
        } else if (!thermal && row.scalar == "total_nonprecipitating_water") {
            ++total_rows;
            if (!check_closure(row, "total-water")) { return false; }
            if (all_dry && (!close_to_zero(row.net_boundary, row) ||
                            !close_to_zero(row.volume_change, row) ||
                            !close_to_zero(row.internal_source, row))) {
                error = "all-dry total-water budget is not closed";
                return false;
            }
        } else if (!thermal && row.scalar == "water_vapor") {
            ++vapor_rows;
            if (!check_closure(row, "water-vapor")) { return false; }
            if (all_dry) {
                for (const auto value : row.faces) {
                    if (!close_to_zero(value, row)) {
                        error = "all-dry vapor wall flux is nonzero";
                        return false;
                    }
                }
            } else if (wet) {
                for (int face = 0; face < 4; ++face) {
                    if (!close_to_zero(row.faces[face], row)) {
                        error = "wet-wall side vapor flux is nonzero";
                        return false;
                    }
                }
                if (!std::isfinite(static_cast<double>(row.faces[4])) ||
                    !std::isfinite(static_cast<double>(row.faces[5]))) {
                    error = "wet-wall vapor flux is non-finite";
                    return false;
                }
            }
        } else if (!thermal && row.scalar == "cloud_water") {
            ++cloud_rows;
            if (!check_closure(row, "cloud-water")) { return false; }
            for (const auto value : row.faces) {
                if (!close_to_zero(value, row)) {
                    error = "cloud-water wall flux is nonzero";
                    return false;
                }
            }
        }
    }
    if (thermal) {
        if (thermal_rows < 3) {
            error = "expected at least three dry rhoTheta budget intervals";
            return false;
        }
        summary.thermal_rows = thermal_rows;
        return true;
    }
    if (total_rows < 3 || vapor_rows < 3 || cloud_rows < 3) {
        error = "expected at least three budget intervals for every water scalar";
        return false;
    }
    summary.total_rows = total_rows;
    summary.vapor_rows = vapor_rows;
    summary.cloud_rows = cloud_rows;
    return true;
}

} // namespace

int main (int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : std::string();
    const bool budget_mode = is_budget_mode(mode);
    const int expected_argc = budget_mode ? 5 : 4;
    if (!is_checker_mode(mode) || argc != expected_argc) {
        std::cerr << "usage: checker mode initial_plotfile final_plotfile\n"
                  << "       checker parity budget_off_plotfile budget_on_plotfile\n"
                  << "       checker all_dry|wet_budget|thermal_budget initial_plotfile final_plotfile budget_file\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);
    if (mode == "parity") {
        PlotFileData budget_off(argv[2]);
        PlotFileData budget_on(argv[3]);
        const std::vector<std::string> fields = {
            "density", "theta", "temp", "qv", "qc",
            "x_velocity", "y_velocity", "z_velocity"};
        Real max_difference = Real(0.0);
        for (const auto& name : fields) {
            if (!has_variable(budget_off, name) || !has_variable(budget_on, name)) {
                amrex::Finalize();
                return fail("parity comparison missing variable " + name);
            }
            MultiFab difference = budget_off.get(0, name);
            MultiFab::Subtract(difference, budget_on.get(0, name), 0, 0, 1, 0);
            max_difference = std::max(max_difference, difference.norm0(0, 0, false));
        }
        Real max_integral_difference = Real(0.0);
        for (const auto& scalar : {std::string("rhoTheta"), std::string("qv"), std::string("qc")}) {
            max_integral_difference = std::max(max_integral_difference,
                std::abs(integral(budget_off, scalar) - integral(budget_on, scalar)));
        }
        const Real parity_scale = std::max(Real(1.0), max_integral_difference);
        const Real parity_tolerance = Real(4096.0) *
            std::numeric_limits<Real>::epsilon() * parity_scale;
        std::cout << "budget_on_off_max_difference=" << max_difference
                  << " budget_on_off_max_integral_difference=" << max_integral_difference
                  << " tolerance=" << parity_tolerance << "\n";
        amrex::Finalize();
        return (max_difference <= parity_tolerance &&
                max_integral_difference <= parity_tolerance) ? 0 :
            fail("budget-on/off state mismatch: max difference=" +
                 std::to_string(static_cast<double>(max_difference)));
    }
    PlotFileData initial(argv[2]);
    PlotFileData final(argv[3]);
    const bool cloudy = (mode == "cloudy" || mode == "all_dry" || mode == "wet_budget");
    if (!cloudy && mode != "dry" && mode != "thermal_budget") {
        amrex::Finalize();
        return fail("mode must be dry, cloudy, all_dry, wet_budget, or thermal_budget");
    }

    for (const char* name : {"density", "theta", "temp", "x_velocity",
                             "y_velocity", "z_velocity"}) {
        if (!has_variable(initial, name) || !has_variable(final, name)) {
            amrex::Finalize();
            return fail("missing required variable " + std::string(name));
        }
    }
    if (cloudy) {
        for (const char* name : {"qv", "qc"}) {
            if (!has_variable(initial, name) || !has_variable(final, name)) {
                amrex::Finalize();
                return fail("missing cloudy variable " + std::string(name));
            }
        }
        for (const char* name : {"qsat", "rel_humidity"}) {
            if (!has_variable(initial, name) || !has_variable(final, name)) {
                amrex::Finalize();
                return fail("missing cloudy diagnostic " + std::string(name));
            }
        }
    } else if (has_variable(initial, "qv") || has_variable(initial, "qc")) {
        amrex::Finalize();
        return fail("dry case unexpectedly contains moisture variables");
    }

    MultiFab density = initial.get(0, "density");
    MultiFab theta = initial.get(0, "theta");
    MultiFab temp = initial.get(0, "temp");
    MultiFab pressure = initial.get(0, "pressure");
    MultiFab initial_qv;
    MultiFab initial_qc;
    MultiFab initial_qsat;
    MultiFab initial_rh;
    if (cloudy) {
        initial_qv = initial.get(0, "qv");
        initial_qc = initial.get(0, "qc");
        initial_qsat = initial.get(0, "qsat");
        initial_rh = initial.get(0, "rel_humidity");
    }
    MultiFab x_velocity = initial.get(0, "x_velocity");
    MultiFab y_velocity = initial.get(0, "y_velocity");
    MultiFab z_velocity = initial.get(0, "z_velocity");
    if (!density.is_finite() || !theta.is_finite() || !temp.is_finite() ||
        !pressure.is_finite() ||
        !x_velocity.is_finite() || !y_velocity.is_finite() ||
        !z_velocity.is_finite()) {
        amrex::Finalize();
        return fail("initial required field is non-finite");
    }

    const auto domain = initial.probDomain(0);
    const auto prob_lo = initial.probLo();
    const auto dx = initial.cellSize(0);
    const Real lx = dx[0] * Real(domain.length(0));
    const Real ly = dx[1] * Real(domain.length(1));
    const Real lz = dx[2] * Real(domain.length(2));
    const Real temperature_bottom = Real(300.0);
    const Real temperature_top = Real(284.0);
    const Real amplitude = Real(0.02);
    const Real relative_humidity = Real(0.95);
    MultiFab theta_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab temperature_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab qv_error;
    MultiFab qsat_error;
    MultiFab rh_error;
    if (cloudy) {
        qv_error.define(theta.boxArray(), theta.DistributionMap(), 1, 0);
        qsat_error.define(theta.boxArray(), theta.DistributionMap(), 1, 0);
        rh_error.define(theta.boxArray(), theta.DistributionMap(), 1, 0);
    }
    MultiFab velocity_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    for (amrex::MFIter mfi(theta_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto bx = mfi.tilebox();
        const auto th = theta.const_array(mfi);
        const auto t = temp.const_array(mfi);
        const auto p = pressure.const_array(mfi);
        const auto qv = cloudy ? initial_qv.const_array(mfi) : amrex::Array4<const Real>{};
        const auto qsat = cloudy ? initial_qsat.const_array(mfi) : amrex::Array4<const Real>{};
        const auto rh = cloudy ? initial_rh.const_array(mfi) : amrex::Array4<const Real>{};
        const auto u = x_velocity.const_array(mfi);
        const auto v = y_velocity.const_array(mfi);
        const auto w = z_velocity.const_array(mfi);
        const auto therr = theta_error.array(mfi);
        const auto terr = temperature_error.array(mfi);
        const auto qverr = cloudy ? qv_error.array(mfi) : amrex::Array4<Real>{};
        const auto qsaterr = cloudy ? qsat_error.array(mfi) : amrex::Array4<Real>{};
        const auto rherr = cloudy ? rh_error.array(mfi) : amrex::Array4<Real>{};
        const auto velerr = velocity_error.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real x = prob_lo[0] + (Real(i-domain.smallEnd(0)) + Real(0.5))*dx[0];
            const Real y = prob_lo[1] + (Real(j-domain.smallEnd(1)) + Real(0.5))*dx[1];
            const Real z = prob_lo[2] + (Real(k-domain.smallEnd(2)) + Real(0.5))*dx[2];
            const Real expected_temperature = temperature_bottom +
                (temperature_top-temperature_bottom)*(z-prob_lo[2])/lz +
                amplitude * std::sin(Real(2.0)*Real(3.14159265358979323846)*(x-prob_lo[0])/lx) *
                std::sin(Real(2.0)*Real(3.14159265358979323846)*(y-prob_lo[1])/ly) *
                std::sin(Real(3.14159265358979323846)*(z-prob_lo[2])/lz);
            const Real expected_theta = expected_temperature *
                std::pow(p_0/p(i,j,k), R_d/Cp_d);
            const Real expected_qv = cloudy ?
                (RdoRv * relative_humidity *
                 (Real(100.0)*erf_esatw(expected_temperature)) /
                 (p(i,j,k) - relative_humidity *
                  (Real(100.0)*erf_esatw(expected_temperature)))) : Real(0.0);
            Real expected_qsat = Real(0.0);
            if (cloudy) {
                erf_qsatw(expected_temperature, p(i,j,k)*Real(0.01), expected_qsat);
            }
            const Real expected_rh = cloudy ? relative_humidity : Real(0.0);
            therr(i,j,k) = std::abs(th(i,j,k)-expected_theta);
            terr(i,j,k) = std::abs(t(i,j,k)-expected_temperature);
            if (cloudy) {
                qverr(i,j,k) = std::abs(qv(i,j,k)-expected_qv);
                qsaterr(i,j,k) = std::abs(qsat(i,j,k)-expected_qsat);
                rherr(i,j,k) = std::abs(rh(i,j,k)-expected_rh);
            }
            velerr(i,j,k) = std::max(std::abs(u(i,j,k)),
                                     std::max(std::abs(v(i,j,k)),std::abs(w(i,j,k))));
        });
    }
    amrex::Gpu::streamSynchronize();

    const Real theta_error_max = theta_error.norm0(0, 0, false);
    const Real temperature_error_max = temperature_error.norm0(0, 0, false);
    const Real qv_error_max = cloudy ? qv_error.norm0(0, 0, false) : Real(0.0);
    const Real qsat_error_max = cloudy ? qsat_error.norm0(0, 0, false) : Real(0.0);
    const Real rh_error_max = cloudy ? rh_error.norm0(0, 0, false) : Real(0.0);
    const Real initial_velocity_max = velocity_error.norm0(0, 0, false);
    const Real theta_tolerance = Real(2.0e-10) * std::max(Real(1.0), temperature_bottom);
    if (theta_error_max > theta_tolerance) {
        amrex::Finalize();
        return fail("initial theta profile mismatch: max error=" +
                    std::to_string(static_cast<double>(theta_error_max)));
    }
    const Real temperature_tolerance = Real(2.0e-10) * Real(300.0);
    if (temperature_error_max > temperature_tolerance) {
        amrex::Finalize();
        return fail("initial temperature mismatch: max error=" +
                    std::to_string(static_cast<double>(temperature_error_max)));
    }
    if (cloudy && qv_error_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial qv profile mismatch: max error=" +
                    std::to_string(static_cast<double>(qv_error_max)));
    }
    if (cloudy && qsat_error_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial qsat diagnostic mismatch: max error=" +
                    std::to_string(static_cast<double>(qsat_error_max)));
    }
    if (cloudy && rh_error_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial relative humidity diagnostic mismatch: max error=" +
                    std::to_string(static_cast<double>(rh_error_max)));
    }
    if (initial_velocity_max > Real(2.0e-12)) {
        amrex::Finalize();
        return fail("initial velocity is not zero: max=" +
                    std::to_string(static_cast<double>(initial_velocity_max)));
    }

    for (const char* name : {"density", "theta", "temp", "x_velocity",
                             "y_velocity", "z_velocity"}) {
        if (!final.get(0, name).is_finite()) {
            amrex::Finalize();
            return fail("final field is non-finite: " + std::string(name));
        }
    }

    if (cloudy) {
        MultiFab qv = final.get(0, "qv");
        MultiFab qc = final.get(0, "qc");
        MultiFab qsat = final.get(0, "qsat");
        MultiFab rh = final.get(0, "rel_humidity");
        if (!qv.is_finite() || !qc.is_finite() || !qsat.is_finite() || !rh.is_finite()) {
            amrex::Finalize();
            return fail("cloudy scalar is non-finite");
        }
        if (qv.min(0) < Real(-1.0e-12) || qc.min(0) < Real(-1.0e-12)) {
            amrex::Finalize();
            return fail("cloudy scalar became negative");
        }
        if (initial_qc.max(0) > Real(2.0e-12)) {
            amrex::Finalize();
            return fail("physical initialization did not start with zero cloud water");
        }
    }

    if (is_budget_mode(mode)) {
        std::vector<BudgetRow> rows;
        BudgetSummary summary;
        std::string budget_error;
        if (!read_budget(argv[4], rows, budget_error) ||
            !check_budget_rows(rows, mode, summary, budget_error)) {
            amrex::Finalize();
            return fail(budget_error.empty() ? "budget validation failed" : budget_error);
        }
        std::cout << std::setprecision(17)
                  << "budget_rows=" << summary.rows << " mode=" << mode;
        if (mode == "thermal_budget") {
            std::cout << " thermal_rows=" << summary.thermal_rows;
        } else {
            std::cout << " total_rows=" << summary.total_rows
                      << " vapor_rows=" << summary.vapor_rows
                      << " cloud_rows=" << summary.cloud_rows;
        }
        std::cout << " max_residual_ratio=" << summary.max_residual_ratio << "\n";
    }

    const MultiFab final_u = final.get(0, "x_velocity");
    const MultiFab final_v = final.get(0, "y_velocity");
    const MultiFab final_w = final.get(0, "z_velocity");
    const Real evolved_velocity = std::max(final_u.norm0(0,0,false),
        std::max(final_v.norm0(0,0,false), final_w.norm0(0,0,false)));
    std::cout << "mode=" << mode << " initial_theta_error=" << theta_error_max
              << " initial_temperature_error=" << temperature_error_max
              << " evolved_velocity_max=" << evolved_velocity;
    if (cloudy) {
        std::cout << " final_qv_min=" << final.get(0,"qv").min(0)
                  << " final_qc_max=" << final.get(0,"qc").max(0)
                  << " initial_qsat_error=" << qsat_error_max
                  << " initial_relative_humidity_error=" << rh_error_max;
    }
    std::cout << "\n";

    amrex::Finalize();
    return evolved_velocity > Real(1.0e-12) ? 0 : fail("thermal perturbation produced no buoyant response");
}
