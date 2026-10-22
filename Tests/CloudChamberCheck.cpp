#include <AMReX.H>
#include "ERF_Constants.H"
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include "../Source/DataStructs/ERF_DataStruct.H"
#include "../Source/Prob/ERF_CloudChamber.H"
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
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::Box;
using amrex::MFIter;
using amrex::ParallelFor;
using amrex::PlotFileData;
using amrex::Real;
using amrex::TilingIfNotGPU;

Real scaled_tolerance (Real expected, Real ulps = Real(256.0))
{
    const Real scale = std::max(Real(1.0), std::abs(expected));
    return ulps * std::numeric_limits<Real>::epsilon() * scale;
}

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

bool retained_flux_is_active (Real value)
{
    // This is an activation oracle for the directly retained applied face
    // flux, not a conservation residual. CloudChamberBudget accumulates the
    // production face-flux arrays themselves, so an inactive algebraic-zero
    // path remains exactly zero; any finite nonzero retained value proves that
    // a nonzero wall flux was applied.
    return std::isfinite(static_cast<double>(value)) && value != Real(0.0);
}

bool retained_flux_activation_contract_holds ()
{
    return !retained_flux_is_active(Real(0.0)) &&
        retained_flux_is_active(Real(0.1)) &&
        retained_flux_is_active(Real(1.0e-20)) &&
        !retained_flux_is_active(std::numeric_limits<Real>::quiet_NaN()) &&
        !retained_flux_is_active(std::numeric_limits<Real>::infinity());
}

bool is_budget_mode (const std::string& mode)
{
    return mode == "all_dry" || mode == "wet_budget" || mode == "bulk_wet" ||
        mode == "neutral_wet" || mode == "most_wet";
}

bool is_checker_mode (const std::string& mode)
{
    return mode == "dry" || mode == "cloudy" || mode == "parity" ||
        mode == "bulk" || mode == "bulk_legacy" ||
        mode == "neutral_momentum" || mode == "fixed_momentum" ||
        mode == "most_momentum" ||
        is_budget_mode(mode);
}

bool check_bulk_plotfiles (const std::string& model_name,
                           PlotFileData& initial, PlotFileData& final,
                           std::string& error,
                           const bool physical_initialization)
{
    MoistureType model;
    try {
        model = amrex::getEnumCaseInsensitive<MoistureType>(model_name);
    } catch (const std::exception& exception) {
        error = "invalid bulk checker model: " + std::string(exception.what());
        return false;
    }
    const auto moisture = MoistureComponentIndices::from_moisture_model(model);
    if (!erf_cloud_chamber::supports_moisture_model(model) ||
        model == MoistureType::None || model == MoistureType::SatAdj ||
        model == MoistureType::SBM || model == MoistureType::SuperDroplets) {
        error = "checker requires a supported active non-SatAdj moisture model";
        return false;
    }

    const std::vector<std::string> core = {
        "density", "theta", "temp", "pressure", "qv", "qc", "qt", "qsat",
        "x_velocity", "y_velocity", "z_velocity"};
    for (const auto& name : core) {
        if (!has_variable(initial, name) || !has_variable(final, name)) {
            error = "missing bulk field " + name;
            return false;
        }
    }

    const std::vector<std::string> moisture_candidates = {
        "qv", "qc", "qi", "qrain", "qsnow", "qgraup", "nc", "ni", "nr",
        "ns", "ng", "nn", "rain_accum", "snow_accum", "graup_accum",
        "qt", "qsat", "rel_humidity"};
    std::vector<std::string> active_fields;
    for (const auto& name : moisture_candidates) {
        const auto status = moisture.query_var(name);
        const bool in_initial = has_variable(initial, name);
        const bool in_final = has_variable(final, name);
        if (name == "rel_humidity" && (in_initial || in_final)) {
            error = "rel_humidity is not a common bulk moisture diagnostic";
            return false;
        }
        if (status.governed && status.available) {
            if (!in_initial || !in_final) {
                error = "active moisture field is missing from a bulk plotfile: " + name;
                return false;
            }
            active_fields.push_back(name);
        } else if (status.governed && (in_initial || in_final)) {
            error = "inactive moisture diagnostic was requested for this model: " + name;
            return false;
        }
    }
    if (initial.time() >= final.time()) {
        error = "bulk simulation did not advance beyond its initial plotfile";
        return false;
    }

    MultiFab density = initial.get(0, "density");
    MultiFab theta = initial.get(0, "theta");
    MultiFab temperature = initial.get(0, "temp");
    MultiFab pressure = initial.get(0, "pressure");
    MultiFab qv = initial.get(0, "qv");
    MultiFab qc = initial.get(0, "qc");
    if (!density.is_finite() || !theta.is_finite() || !temperature.is_finite() ||
        !pressure.is_finite() || !qv.is_finite() || !qc.is_finite()) {
        error = "initial bulk thermodynamic field is non-finite";
        return false;
    }

    const auto domain = initial.probDomain(0);
    const auto prob_lo = initial.probLo();
    const auto dx = initial.cellSize(0);
    const Real lx = dx[0] * Real(domain.length(0));
    const Real ly = dx[1] * Real(domain.length(1));
    const Real lz = dx[2] * Real(domain.length(2));
    const Real bottom = Real(300.0);
    const Real top = Real(284.0);
    const Real amplitude = Real(0.02);
    const Real relative_humidity = Real(0.95);
    const Real legacy_qv_bottom = Real(0.004);
    const Real legacy_qv_top = Real(0.006);
    MultiFab error_field(theta.boxArray(), theta.DistributionMap(), 3, 0);
    for (MFIter mfi(error_field, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto th = theta.const_array(mfi);
        const auto temp = temperature.const_array(mfi);
        const auto p = pressure.const_array(mfi);
        const auto vapor = qv.const_array(mfi);
        const auto cloud = qc.const_array(mfi);
        const auto out = error_field.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const Real x = prob_lo[0] + (Real(i - domain.smallEnd(0)) + Real(0.5))*dx[0];
            const Real y = prob_lo[1] + (Real(j - domain.smallEnd(1)) + Real(0.5))*dx[1];
            const Real z = prob_lo[2] + (Real(k - domain.smallEnd(2)) + Real(0.5))*dx[2];
            const Real profile_value = bottom + (top - bottom) *
                (z - prob_lo[2]) / lz + amplitude *
                std::sin(Real(2.0)*Real(3.14159265358979323846)*(x-prob_lo[0])/lx) *
                std::sin(Real(2.0)*Real(3.14159265358979323846)*(y-prob_lo[1])/ly) *
                std::sin(Real(3.14159265358979323846)*(z-prob_lo[2])/lz);
            if (physical_initialization) {
                const Real expected_temperature = profile_value;
                const Real expected_theta = expected_temperature *
                    std::pow(p_0 / p(i,j,k), R_d / Cp_d);
                const Real vapor_pressure = relative_humidity *
                    (Real(100.0) * erf_esatw(expected_temperature));
                const Real expected_qv = RdoRv * vapor_pressure /
                    (p(i,j,k) - vapor_pressure);
                out(i,j,k,0) = std::abs(temp(i,j,k) - expected_temperature);
                out(i,j,k,1) = std::abs(th(i,j,k) - expected_theta);
                out(i,j,k,2) = std::abs(vapor(i,j,k) - expected_qv);
            } else {
                const Real expected_qv = legacy_qv_bottom +
                    (legacy_qv_top - legacy_qv_bottom) *
                    (z - prob_lo[2]) / lz;
                out(i,j,k,0) = std::abs(th(i,j,k) - profile_value);
                out(i,j,k,1) = std::abs(vapor(i,j,k) - expected_qv);
                out(i,j,k,2) = cloud(i,j,k) == Real(0.0) ? Real(0.0) :
                    std::numeric_limits<Real>::infinity();
            }
            if (cloud(i,j,k) != Real(0.0)) {
                out(i,j,k,0) = std::numeric_limits<Real>::infinity();
            }
        });
    }
    amrex::Gpu::streamSynchronize();
    const Real first_tolerance = scaled_tolerance(Real(300.0));
    const Real second_tolerance = physical_initialization ?
        scaled_tolerance(Real(300.0)) : scaled_tolerance(Real(1.0));
    const Real third_tolerance = physical_initialization ?
        scaled_tolerance(Real(1.0)) : Real(0.0);
    if (error_field.norm0(0, 0, false) > first_tolerance ||
        error_field.norm0(1, 0, false) > second_tolerance ||
        error_field.norm0(2, 0, false) > third_tolerance) {
        error = physical_initialization ?
            "bulk cold-start temperature, theta, qv, or zero-qc contract failed" :
            "legacy bulk cold-start theta, qv, or zero-qc contract failed";
        return false;
    }

    const std::vector<std::string> cold_start_zero = {
        "qc", "qi", "qrain", "qsnow", "qgraup", "nc", "ni", "nr", "ns",
        "ng", "nn", "rain_accum", "snow_accum", "graup_accum"};
    for (const auto& name : cold_start_zero) {
        if (moisture.query_var(name).available &&
            initial.get(0, name).norm0(0, 0, false) > scaled_tolerance(Real(1.0))) {
            error = "bulk cold start did not leave " + name + " at zero";
            return false;
        }
    }
    for (const auto& name : active_fields) {
        MultiFab initial_field = initial.get(0, name);
        MultiFab final_field = final.get(0, name);
        if (!initial_field.is_finite() || !final_field.is_finite()) {
            error = "bulk moisture field is non-finite: " + name;
            return false;
        }
        if (name != "rel_humidity" &&
            final_field.min(0) < -scaled_tolerance(Real(1.0))) {
            error = "bulk moisture field became negative: " + name;
            return false;
        }
    }
    if (model == MoistureType::WDM6 &&
        final.get(0, "nn").max(0) <= Real(0.0)) {
        error = "WDM6 native microphysics did not populate its CCN number reservoir";
        return false;
    }
    for (const auto& name : core) {
        if (!final.get(0, name).is_finite()) {
            error = "final bulk field is non-finite: " + name;
            return false;
        }
    }
    return true;
}

struct BudgetSummary {
    int rows = 0;
    int thermal_rows = 0;
    int total_rows = 0;
    int vapor_rows = 0;
    int cloud_rows = 0;
    bool heat_face_active = false;
    bool vapor_face_active = false;
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
    const bool wet = mode == "wet_budget" || mode == "bulk_wet" ||
        mode == "neutral_wet" || mode == "most_wet";
    const bool cloudy = is_budget_mode(mode);
    const bool require_transfer_activation = mode == "bulk_wet" ||
        mode == "neutral_wet" || mode == "most_wet";
    summary = {};
    summary.rows = static_cast<int>(rows.size());
    const auto row_is_finite = [](const BudgetRow& row) {
        if (!std::isfinite(static_cast<double>(row.net_boundary)) ||
            !std::isfinite(static_cast<double>(row.volume_change)) ||
            !std::isfinite(static_cast<double>(row.internal_source)) ||
            !std::isfinite(static_cast<double>(row.residual)) ||
            !std::isfinite(static_cast<double>(row.tolerance)) ||
            row.tolerance < Real(0.0)) {
            return false;
        }
        for (const auto value : row.faces) {
            if (!std::isfinite(static_cast<double>(value))) { return false; }
        }
        return true;
    };
    const auto check_closure = [&](const BudgetRow& row, const char* name) {
        if (!row_is_finite(row)) {
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
    const auto check_unsupported_source = [&](const BudgetRow& row, const char* name) {
        if (!row_is_finite(row)) {
            error = std::string(name) + " budget diagnostics are non-finite or have negative tolerance";
            return false;
        }
        if (row.status != "UNSUPPORTED_SOURCE") {
            error = std::string(name) + " budget must be UNSUPPORTED_SOURCE";
            return false;
        }
        return true;
    };
    for (const auto& row : rows) {
        if (row.scalar == "rhoTheta") {
            ++summary.thermal_rows;
            if (cloudy) {
                if (!check_unsupported_source(row, "rhoTheta")) { return false; }
            } else if (!check_closure(row, "rhoTheta")) {
                return false;
            }
            if (require_transfer_activation) {
                for (int face = 2 * (AMREX_SPACEDIM - 1);
                     face < 2 * AMREX_SPACEDIM; ++face) {
                    summary.heat_face_active = summary.heat_face_active ||
                        retained_flux_is_active(row.faces[face]);
                }
            }
        } else if (row.scalar == "total_nonprecipitating_water") {
            ++summary.total_rows;
            if (!check_closure(row, "total-water")) { return false; }
            if (all_dry && (!close_to_zero(row.net_boundary, row) ||
                            !close_to_zero(row.volume_change, row) ||
                            !close_to_zero(row.internal_source, row))) {
                error = "all-dry total-water budget is not closed";
                return false;
            }
        } else if (row.scalar == "water_vapor") {
            ++summary.vapor_rows;
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
                if (require_transfer_activation) {
                    for (int face = 2 * (AMREX_SPACEDIM - 1);
                         face < 2 * AMREX_SPACEDIM; ++face) {
                        summary.vapor_face_active = summary.vapor_face_active ||
                            retained_flux_is_active(row.faces[face]);
                    }
                }
            }
        } else if (row.scalar == "cloud_water") {
            ++summary.cloud_rows;
            if (!check_closure(row, "cloud-water")) { return false; }
            for (int face = 0; face < 2 * AMREX_SPACEDIM; ++face) {
                if (row.faces[face] != Real(0.0)) {
                    error = "cloud-water wall flux is not algebraically zero on face " +
                        std::to_string(face);
                    return false;
                }
            }
        }
    }
    if (summary.thermal_rows < 2 || summary.total_rows < 2 ||
        summary.vapor_rows < 2 || summary.cloud_rows < 2) {
        error = "expected at least two budget intervals for rhoTheta and every water scalar";
        return false;
    }
    if (require_transfer_activation) {
        if (!summary.heat_face_active) {
            error = "active wet-wall heat oracle found no retained rhoTheta wall flux above its activation threshold";
            return false;
        }
        if (!summary.vapor_face_active) {
            error = "active wet-wall vapor oracle found no retained water_vapor wall flux above its activation threshold";
            return false;
        }
    }
    return true;
}

} // namespace

int main (int argc, char** argv)
{
    const std::string mode = argc > 1 ? argv[1] : std::string();
    const bool budget_mode = is_budget_mode(mode);
    const bool bulk_mode = mode == "bulk" || mode == "bulk_legacy";
    const bool activation_mode = mode == "neutral_momentum" ||
        mode == "fixed_momentum" || mode == "most_momentum";
    const int expected_argc = (budget_mode || activation_mode || bulk_mode) ? 5 : 4;
    if (!is_checker_mode(mode) || argc != expected_argc) {
        std::cerr << "usage: checker mode initial_plotfile final_plotfile\n"
                  << "       checker parity budget_off_plotfile budget_on_plotfile\n"
                  << "       checker all_dry|wet_budget|bulk_wet|neutral_wet|most_wet initial_plotfile final_plotfile budget_file\n"
                  << "       checker neutral_momentum|fixed_momentum|most_momentum initial_plotfile final_a final_b\n";
        return 2;
    }
    if (!retained_flux_activation_contract_holds()) {
        return fail("retained-flux activation oracle contract is broken");
    }

    amrex::Initialize(argc, argv, false);
    if (bulk_mode) {
        PlotFileData initial(argv[3]);
        PlotFileData final(argv[4]);
        std::string error;
        const bool valid = check_bulk_plotfiles(
            argv[2], initial, final, error, mode == "bulk");
        if (!valid) {
            amrex::Finalize();
            return fail(error);
        }
        std::cout << "mode=" << mode << " model=" << argv[2]
                  << " active_fields=" << initial.varNames().size()
                  << " initial_time=" << initial.time()
                  << " final_time=" << final.time() << "\n";
        amrex::Finalize();
        return 0;
    }
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
            MultiFab on_field = budget_on.get(0, name);
            if (!difference.is_finite() || !on_field.is_finite()) {
                amrex::Finalize();
                return fail("parity comparison found non-finite variable " + name);
            }
            MultiFab::Subtract(difference, on_field, 0, 0, 1, 0);
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
    std::unique_ptr<PlotFileData> alternate;
    if (activation_mode) {
        alternate = std::make_unique<PlotFileData>(argv[4]);
    }
    const bool cloudy = (mode == "cloudy" || mode == "all_dry" ||
                         mode == "wet_budget" || mode == "bulk_wet" ||
                         mode == "neutral_wet" || mode == "most_wet" ||
                         mode == "most_momentum");
    if (!cloudy && mode != "dry" && mode != "neutral_momentum" &&
        mode != "fixed_momentum") {
        amrex::Finalize();
        return fail("mode must be dry, cloudy, all_dry, wet_budget, bulk_wet, neutral_wet, most_wet, neutral_momentum, fixed_momentum, or most_momentum");
    }

    for (const char* name : {"density", "theta", "temp", "x_velocity",
                             "y_velocity", "z_velocity"}) {
        if (!has_variable(initial, name) || !has_variable(final, name) ||
            !initial.get(0, name).is_finite() || !final.get(0, name).is_finite() ||
            (activation_mode && (!has_variable(*alternate, name) ||
                                 !alternate->get(0, name).is_finite()))) {
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
    const Real theta_tolerance = scaled_tolerance(temperature_bottom);
    if (!budget_mode && theta_error_max > theta_tolerance) {
        amrex::Finalize();
        return fail("initial theta profile mismatch: max error=" +
                    std::to_string(static_cast<double>(theta_error_max)));
    }
    const Real temperature_tolerance = scaled_tolerance(temperature_bottom);
    if (!budget_mode && temperature_error_max > temperature_tolerance) {
        amrex::Finalize();
        return fail("initial temperature mismatch: max error=" +
                    std::to_string(static_cast<double>(temperature_error_max)));
    }
    if (!budget_mode && cloudy && qv_error_max > scaled_tolerance(Real(1.0))) {
        amrex::Finalize();
        return fail("initial qv profile mismatch: max error=" +
                    std::to_string(static_cast<double>(qv_error_max)));
    }
    if (!budget_mode && cloudy && qsat_error_max > scaled_tolerance(Real(1.0))) {
        amrex::Finalize();
        return fail("initial qsat diagnostic mismatch: max error=" +
                    std::to_string(static_cast<double>(qsat_error_max)));
    }
    if (!budget_mode && cloudy && rh_error_max > scaled_tolerance(relative_humidity)) {
        amrex::Finalize();
        return fail("initial relative humidity diagnostic mismatch: max error=" +
                    std::to_string(static_cast<double>(rh_error_max)));
    }
    if (!budget_mode && !activation_mode && initial_velocity_max > scaled_tolerance(Real(1.0))) {
        amrex::Finalize();
        return fail("initial velocity is not zero: max=" +
                    std::to_string(static_cast<double>(initial_velocity_max)));
    }

    for (const char* name : {"density", "theta", "temp", "x_velocity",
                             "y_velocity", "z_velocity"}) {
        if (!final.get(0, name).is_finite() ||
            (activation_mode && !alternate->get(0, name).is_finite())) {
            amrex::Finalize();
            return fail("final field is non-finite: " + std::string(name));
        }
    }

    if (cloudy) {
        MultiFab qv = final.get(0, "qv");
        MultiFab qc = final.get(0, "qc");
        MultiFab qsat = final.get(0, "qsat");
        MultiFab rh = final.get(0, "rel_humidity");
        if (!initial_qv.is_finite() || !initial_qc.is_finite() ||
            !initial_qsat.is_finite() || !initial_rh.is_finite() ||
            !qv.is_finite() || !qc.is_finite() || !qsat.is_finite() || !rh.is_finite()) {
            amrex::Finalize();
            return fail("cloudy scalar is non-finite");
        }
        const Real nonnegative_tolerance = scaled_tolerance(Real(1.0));
        if (qv.min(0) < -nonnegative_tolerance ||
            qc.min(0) < -nonnegative_tolerance) {
            amrex::Finalize();
            return fail("cloudy scalar became negative");
        }
        if (!budget_mode && initial_qc.max(0) > scaled_tolerance(Real(1.0))) {
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
        std::cout << " thermal_rows=" << summary.thermal_rows
                  << " total_rows=" << summary.total_rows
                  << " vapor_rows=" << summary.vapor_rows
                  << " cloud_rows=" << summary.cloud_rows
                  << " heat_face_active=" << summary.heat_face_active
                  << " vapor_face_active=" << summary.vapor_face_active;
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
    Real activation_difference = Real(0.0);
    if (activation_mode) {
        for (const auto& name : {std::string("x_velocity"),
                                 std::string("y_velocity"),
                                 std::string("z_velocity")}) {
            MultiFab difference = final.get(0, name);
            MultiFab::Subtract(difference, alternate->get(0, name), 0, 0, 1, 0);
            activation_difference = std::max(activation_difference,
                                              difference.norm0(0, 0, false));
        }
        std::cout << " alternate_velocity_difference=" << activation_difference
                  << " z0_m_response=nonzero";
    }
    if (cloudy) {
        std::cout << " final_qv_min=" << final.get(0,"qv").min(0)
                  << " final_qc_max=" << final.get(0,"qc").max(0)
                  << " initial_qsat_error=" << qsat_error_max
                  << " initial_relative_humidity_error=" << rh_error_max;
    }
    std::cout << "\n";

    amrex::Finalize();
    if (activation_mode &&
        activation_difference <= scaled_tolerance(Real(1.0), Real(16.0))) {
        return fail("momentum wall parameter change produced no production response");
    }
    // This is a nonzero-response guard, not a magnitude assertion. A few
    // ulps are sufficient because the short regression intentionally starts
    // from a nearly motionless chamber.
    if (budget_mode) { return 0; }
    return evolved_velocity > scaled_tolerance(Real(1.0), Real(4.0)) ? 0 :
        fail("thermal perturbation produced no buoyant response");
}
