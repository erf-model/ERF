#include <ERF_InputPreflight.H>
#include <ERF_DescriptionWriter.H>
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace erf_preflight {
namespace {
struct Diagnostic {
    std::string code, parameter, message;
};

struct Report {
    std::string input;
    erf_description::Format format = erf_description::Format::JSON;
    std::vector<Diagnostic> diagnostics;
    std::vector<std::string> completed;
    int exit_code = 0;

    void error (const std::string& code, const std::string& parameter,
                const std::string& message, int status = 2)
    {
        diagnostics.push_back({code, parameter, message});
        exit_code = std::max(exit_code, status);
    }

    void write () const
    {
        erf_description::Writer writer(std::cout, format);
        writer.value("schema_version", 1);
        writer.value("application", "ERF");
        writer.value("scope", "basic-input-v1");
        writer.value("status", exit_code == 0 ? "passed" : exit_code == 2 ? "invalid" : "error");
        writer.value("exit_code", exit_code);
        writer.value("input", input.c_str());
        writer.value("simulation_validated", false);
        writer.value("device_initialized", false);
        writer.strings("checks_completed", completed);
        writer.strings("checks_deferred", {"unrecognized-parameters", "problem-initializer",
            "physics-data-and-compatibility", "amr-layout", "boundary-conditions",
            "checkpoint-content-and-compatibility", "device-and-mpi-runtime", "numerical-stability"});
        writer.objects("diagnostics", diagnostics, [](auto& item, const Diagnostic& diagnostic) {
            item.value("code", diagnostic.code.c_str());
            item.value("severity", "error");
            item.value("parameter", diagnostic.parameter.empty() ? nullptr : diagnostic.parameter.c_str());
            item.value("message", diagnostic.message.c_str());
            item.value("location", nullptr); // ParmParse has no source locations.
        });
        writer.finish();
    }
};

// Rank zero parses alone on MPI_COMM_SELF with a throwing error handler, so a bad
// file becomes a diagnostic instead of stranding the other ranks in a collective.
class ParserSession {
public:
    ParserSession () : old_out(amrex::system::osout), old_err(amrex::system::oserr),
        old_handler(amrex::system::error_handler)
    {
        int argc = 0;
        char** argv = nullptr;
#ifdef AMREX_USE_MPI
        amrex::ParallelDescriptor::StartParallel(&argc, &argv, MPI_COMM_SELF);
#else
        amrex::ParallelDescriptor::StartParallel(&argc, &argv, amrex::MPI_COMM_WORLD);
#endif
        amrex::system::osout = &messages;
        amrex::system::oserr = &messages;
        amrex::SetErrorHandler([](const char* message) {
            throw std::runtime_error(message ? message : "AMReX input processing failed");
        });
    }

    ~ParserSession ()
    {
        amrex::ParmParse::SetVerbose(0);
        amrex::ParmParse::Finalize();
        amrex::ParallelDescriptor::EndParallel();
        amrex::SetErrorHandler(old_handler);
        amrex::system::osout = old_out;
        amrex::system::oserr = old_err;
    }

    // What AMReX printed while parsing, as one line.
    std::string output () const
    {
        std::string text;
        for (const char c : messages.str()) {
            if (std::isspace(static_cast<unsigned char>(c))) {
                if (!text.empty() && text.back() != ' ') text += ' ';
            } else {
                text += c;
            }
        }
        if (!text.empty() && text.back() == ' ') text.pop_back();
        return text;
    }

private:
    std::ostringstream messages;
    std::ostream *old_out, *old_err;
    amrex::ErrorHandler old_handler;
};

template <class T>
std::optional<std::vector<T>> values (Report& report, const std::string& name, int count = 1)
{
    amrex::ParmParse pp;
    if (!pp.contains(name)) return std::nullopt;
    if (pp.countval(name) != count) {
        report.error("ERF_INPUT_VALUE_COUNT", name, "Expected " + std::to_string(count) + " value(s).");
        return std::nullopt;
    }
    try {
        std::vector<T> value;
        if constexpr (std::is_same_v<T, bool>) {
            bool scalar = false;
            pp.get(name, scalar);
            value.push_back(scalar);
        } else {
            pp.getarr(name, value);
        }
        return value;
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::exception&) {
        const char* expected = std::is_same_v<T, bool> ? "boolean" :
            std::is_integral_v<T> ? "integer" : std::is_floating_point_v<T> ? "real" : "string";
        report.error("ERF_INPUT_VALUE_TYPE", name, std::string("Expected ") + expected + " value(s) accepted by AMReX.");
        return std::nullopt;
    }
}

std::optional<double> real (Report& report, const std::string& name)
{
    const auto value = values<double>(report, name);
    if (!value) return std::nullopt;
    if (!std::isfinite((*value)[0])) {
        report.error("ERF_INPUT_NONFINITE", name, "Expected a finite number.");
        return std::nullopt;
    }
    return (*value)[0];
}

std::string exclusive_alias (Report& report, const std::string& current, const std::string& legacy)
{
    amrex::ParmParse pp;
    if (pp.contains(current) && pp.contains(legacy)) {
        report.error("ERF_INPUT_ALIAS_CONFLICT", current, "Specify either " + current + " or " + legacy + ", not both.");
    }
    return pp.contains(current) ? current : legacy;
}

void mesh (Report& report)
{
    amrex::ParmParse pp;
    const auto cell_name = exclusive_alias(report, "erf.n_cell", "amr.n_cell");
    if (!pp.contains(cell_name)) report.error("ERF_INPUT_REQUIRED", cell_name, "Cell counts are required.");
    if (const auto cells = values<int>(report, cell_name, AMREX_SPACEDIM)) {
        if (std::any_of(cells->begin(), cells->end(), [](int n) { return n <= 0; }))
            report.error("ERF_INPUT_CELL_COUNT", cell_name, "Every cell count must be positive.");
    }
    const auto lo_name = exclusive_alias(report, "erf.prob_lo", "geometry.prob_lo");
    const bool erf_hi = pp.contains("erf.prob_hi"), erf_extent = pp.contains("erf.prob_extent");
    const bool legacy_hi = pp.contains("geometry.prob_hi"), legacy_extent = pp.contains("geometry.prob_extent");
    if ((erf_hi || erf_extent) && !pp.contains("erf.prob_lo"))
        report.error("ERF_INPUT_REQUIRED", "erf.prob_lo", "ERF-prefixed bounds require erf.prob_lo.");
    if ((erf_hi || erf_extent) && (legacy_hi || legacy_extent))
        report.error("ERF_INPUT_ALIAS_CONFLICT", "erf.prob_hi", "Do not mix ERF and legacy upper-bound or extent parameters.");
    if (legacy_hi && legacy_extent)
        report.error("ERF_INPUT_DOMAIN_CONFLICT", "geometry.prob_hi", "Specify only one of geometry.prob_hi and geometry.prob_extent.");
    const std::string hi_name = erf_hi || erf_extent ? "erf.prob_hi" : "geometry.prob_hi";
    const std::string extent_name = erf_hi || erf_extent ? "erf.prob_extent" : "geometry.prob_extent";
    auto lo = values<amrex::Real>(report, lo_name, AMREX_SPACEDIM);
    if (!pp.contains(lo_name)) lo = std::vector<amrex::Real>(AMREX_SPACEDIM, 0);
    auto hi = values<amrex::Real>(report, hi_name, AMREX_SPACEDIM);
    const auto extent = values<amrex::Real>(report, extent_name, AMREX_SPACEDIM);
    if (!pp.contains(hi_name) && !pp.contains(extent_name))
        report.error("ERF_INPUT_REQUIRED", hi_name, "Specify the domain upper bound or extent.");
    if (extent) {
        if (std::any_of(extent->begin(), extent->end(), [](amrex::Real e) { return !std::isfinite(e) || e <= 0; }))
            report.error("ERF_INPUT_DOMAIN_EXTENT", extent_name, "Every extent must be finite and positive.");
    }
    if (erf_hi && erf_extent && lo && extent && hi) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            const auto derived = (*lo)[i] + (*extent)[i];
            const auto scale = std::max({std::abs((*lo)[i]), std::abs((*extent)[i]), std::abs((*hi)[i])});
            if (std::abs(derived - (*hi)[i]) > 100 * std::numeric_limits<amrex::Real>::epsilon() * scale) {
                report.error("ERF_INPUT_DOMAIN_CONFLICT", hi_name, "Upper bound disagrees with lower bound plus extent.");
                break;
            }
        }
    } else if (lo && extent && !pp.contains(hi_name)) {
        hi = *lo;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) (*hi)[i] += (*extent)[i];
    }
    if (lo && hi) {
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            if (!std::isfinite((*lo)[i]) || !std::isfinite((*hi)[i]) || (*hi)[i] <= (*lo)[i]) {
                report.error("ERF_INPUT_DOMAIN_BOUNDS", hi_name, "Bounds must be finite, with upper greater than lower in every dimension.");
                break;
            }
        }
    }
    const auto periodic_name = exclusive_alias(report, "erf.is_periodic", "geometry.is_periodic");
    if (const auto periodic = values<int>(report, periodic_name, AMREX_SPACEDIM)) {
        if (std::any_of(periodic->begin(), periodic->end(), [](int n) { return n != 0 && n != 1; }))
            report.error("ERF_INPUT_PERIODICITY", periodic_name, "Periodicity values must be 0 or 1.");
    }
    report.completed.push_back("mesh-bounds-and-cell-counts");
}

void controls (Report& report)
{
    values<int>(report, "max_step"); // Negative means no step limit.
    real(report, "stop_time");
    const auto dt = real(report, "erf.fixed_dt");
    const auto cfl = real(report, "erf.cfl");
    if (cfl && *cfl <= 0 && (!dt || *dt <= 0))
        report.error("ERF_INPUT_CFL", "erf.cfl", "A positive CFL is required when fixed_dt is not positive.");
    for (const std::string field : {"erf.plot_int_1", "erf.plot_int_2", "erf.plot2d_int_1", "erf.plot2d_int_2"}) {
        const auto interval = values<int>(report, field);
        auto period_name = field;
        period_name.replace(period_name.find("_int"), 4, "_per");
        const auto period = real(report, period_name);
        if (interval && period && (*interval)[0] > 0 && *period > 0)
            report.error("ERF_INPUT_OUTPUT_CONFLICT", field, "Choose a positive step interval or time period, not both.");
    }
    for (const auto* field : {"erf.check_int", "amr.check_int"}) values<int>(report, field);
    for (const auto* field : {"erf.check_per", "amr.check_per"}) real(report, field);
    report.completed.push_back("stop-timestep-and-output-controls");
}

void features_and_paths (Report& report)
{
    amrex::ParmParse pp;
    if (const auto fft = values<bool>(report, "erf.use_fft")) {
#ifndef ERF_USE_FFT
        if ((*fft)[0]) report.error("ERF_INPUT_FEATURE_UNAVAILABLE", "erf.use_fft", "This ERF executable was built without FFT support.");
#endif
    }
    report.completed.push_back("compiled-fft-requirement");
    // amr.restart overrides erf.restart, as in ReadParameters, even when empty.
    const auto erf_restart = values<std::string>(report, "erf.restart");
    const auto amr_restart = values<std::string>(report, "amr.restart");
    const auto& restart = pp.contains("amr.restart") ? amr_restart : erf_restart;
    const std::string parameter = pp.contains("amr.restart") ? "amr.restart" : "erf.restart";
    if (restart && !(*restart)[0].empty()) {
        std::error_code error;
        const bool directory = std::filesystem::is_directory((*restart)[0], error);
        if (error && error != std::errc::no_such_file_or_directory)
            report.error("ERF_INPUT_INSPECTION_FAILED", parameter, "Cannot inspect the checkpoint path: " + error.message(), 3);
        else if (!directory)
            report.error("ERF_INPUT_RESTART_PATH", parameter, "Checkpoint path must name an existing directory; content checks are deferred.");
    }
    report.completed.push_back("restart-path-presence");
}

int inspect (int argc, char* argv[])
{
    Report report;
    std::vector<char*> overrides;
    for (int i = 1; i < argc; ++i) {
        const std::string option(argv[i]);
        if (option == "--validate-input") continue;
        if (option == "--format=json") report.format = erf_description::Format::JSON;
        else if (option == "--format=yaml") report.format = erf_description::Format::YAML;
        else if (option.rfind("--", 0) == 0)
            report.error("ERF_PREFLIGHT_USAGE", "", "Unknown preflight option: " + option);
        else if (report.input.empty()) report.input = option;
        else overrides.push_back(argv[i]);
    }
    if (report.input.empty()) report.error("ERF_PREFLIGHT_USAGE", "", "Supply an input file after --validate-input.");
    if (report.exit_code) { report.write(); return report.exit_code; }
    try {
        auto path = report.input;
        if (const char* prefix = std::getenv("AMREX_INPUTS_FILE_PREFIX")) {
            if (!*prefix) {
                report.error("ERF_INPUT_INSPECTION_FAILED", "", "AMREX_INPUTS_FILE_PREFIX must be unset or nonempty.", 3);
                report.write(); return report.exit_code;
            }
            path = std::string(prefix) + (std::string(prefix).back() == '/' ? "" : "/") + path;
        }
        std::error_code error;
        const bool regular = std::filesystem::is_regular_file(path, error);
        if (error && error != std::errc::no_such_file_or_directory)
            report.error("ERF_INPUT_INSPECTION_FAILED", "", "Cannot inspect the input file: " + error.message(), 3);
        else if (!regular)
            report.error("ERF_INPUT_FILE", "", "Input path must name an existing regular file.");
        else if (!std::ifstream(path).good())
            report.error("ERF_INPUT_INSPECTION_FAILED", "", "Cannot read the input file.", 3);
        if (!report.exit_code) {
            ParserSession parser;
            try {
                amrex::ParmParse::Initialize(static_cast<int>(overrides.size()), overrides.data(), report.input.c_str());
                report.completed.push_back("amrex-input-parsing");
                mesh(report);
                controls(report);
                features_and_paths(report);
            } catch (const std::bad_alloc&) {
                throw;
            } catch (const std::exception& exception) {
                std::string message = "AMReX input processing failed";
                const std::string what = exception.what(), printed = parser.output();
                if (!what.empty() && what != message) message += ": " + what;
                if (!printed.empty()) message += " (" + printed + ")";
                report.error("ERF_INPUT_PARSE_FAILED", "", message);
            }
        }
    } catch (const std::exception&) {
        report.error("ERF_INPUT_INSPECTION_FAILED", "", "Preflight could not complete input inspection.", 3);
    }
    report.write();
    return report.exit_code;
}
} // namespace

int run (int argc, char* argv[])
{
    int rank = 0, result = 0;
#ifdef AMREX_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    if (rank == 0) result = inspect(argc, argv);
#ifdef AMREX_USE_MPI
    MPI_Bcast(&result, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    return result;
}
} // namespace erf_preflight
