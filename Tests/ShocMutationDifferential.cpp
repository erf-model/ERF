#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

struct Arguments {
    std::string field;
    std::string baseline_initial;
    std::string mutant_initial;
    std::string baseline_final;
    std::string mutant_final;
    Real min_final_difference = 0.0;
    Real min_baseline_evolution = 0.0;
    Real max_mutant_to_baseline_evolution_ratio = 0.0;
};

// Paired mutation tests compare only fields that are valid in the initial
// plotfile.  Native SHOC diagnostics are deliberately excluded because they
// are not diagnosed until after the first saved state.
const std::vector<std::string>& host_state_fields ()
{
    static const std::vector<std::string> fields {
        "density", "temp", "theta", "pressure", "qv", "qc", "rhoKE",
        "x_velocity", "y_velocity"
    };
    return fields;
}

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    return std::find(plotfile.varNames().begin(), plotfile.varNames().end(), name) !=
           plotfile.varNames().end();
}

Real maximum_difference (const MultiFab& lhs, const MultiFab& rhs)
{
    MultiFab difference(lhs.boxArray(), lhs.DistributionMap(), lhs.nComp(), 0);
    MultiFab::Copy(difference, lhs, 0, 0, 1, 0);
    MultiFab::Subtract(difference, rhs, 0, 0, 1, 0);
    return difference.norm0();
}

bool parse_real (const std::string& text, Real& value)
{
    try {
        std::size_t consumed = 0;
        const double parsed = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(parsed)) {
            return false;
        }
        value = static_cast<Real>(parsed);
    } catch (...) {
        return false;
    }
    return true;
}

bool read_arguments (int argc, char** argv, Arguments& args, std::string& error)
{
    bool have_min_final_difference = false;
    bool have_min_baseline_evolution = false;
    bool have_max_ratio = false;
    for (int i = 1; i < argc; ++i) {
        const std::string option(argv[i]);
        if (i + 1 >= argc) {
            error = "missing value after " + option;
            return false;
        }
        const std::string value(argv[++i]);
        if (option == "--field") {
            args.field = value;
        } else if (option == "--min-final-difference") {
            have_min_final_difference = parse_real(value, args.min_final_difference);
            if (!have_min_final_difference) {
                error = "invalid --min-final-difference value '" + value + "'";
                return false;
            }
        } else if (option == "--min-baseline-evolution") {
            have_min_baseline_evolution = parse_real(value, args.min_baseline_evolution);
            if (!have_min_baseline_evolution) {
                error = "invalid --min-baseline-evolution value '" + value + "'";
                return false;
            }
        } else if (option == "--max-mutant-to-baseline-evolution-ratio") {
            have_max_ratio = parse_real(value, args.max_mutant_to_baseline_evolution_ratio);
            if (!have_max_ratio) {
                error = "invalid --max-mutant-to-baseline-evolution-ratio value '" + value + "'";
                return false;
            }
        } else if (option == "--baseline-initial") {
            args.baseline_initial = value;
        } else if (option == "--mutant-initial") {
            args.mutant_initial = value;
        } else if (option == "--baseline-final") {
            args.baseline_final = value;
        } else if (option == "--mutant-final") {
            args.mutant_final = value;
        } else {
            error = "unknown option '" + option + "'";
            return false;
        }
    }
    if (args.field.empty() || args.baseline_initial.empty() || args.mutant_initial.empty() ||
        args.baseline_final.empty() || args.mutant_final.empty() ||
        !have_min_final_difference || !have_min_baseline_evolution || !have_max_ratio) {
        error = "missing required mutation argument";
        return false;
    }
    if (args.min_final_difference <= 0.0) {
        error = "--min-final-difference must be positive";
        return false;
    }
    if (args.min_baseline_evolution <= 0.0) {
        error = "--min-baseline-evolution must be positive";
        return false;
    }
    if (args.max_mutant_to_baseline_evolution_ratio <= 0.0 ||
        args.max_mutant_to_baseline_evolution_ratio >= 1.0) {
        error = "--max-mutant-to-baseline-evolution-ratio must be between zero and one";
        return false;
    }
    return true;
}

bool check_finite_fields (PlotFileData& plotfile,
                          const std::vector<std::string>& fields,
                          const std::string& label)
{
    for (const auto& field : fields) {
        if (!has_variable(plotfile, field)) {
            std::cerr << label << " is missing required field '" << field << "'\n";
            return false;
        }
        if (!plotfile.get(0, field).is_finite(0, 1, 0, false)) {
            std::cerr << label << " field '" << field << "' is non-finite\n";
            return false;
        }
    }
    return true;
}

bool check_physical_final_state (PlotFileData& plotfile, const std::string& label)
{
    constexpr Real allowance = 1.0e-12;
    for (const auto& field : {"density", "temp", "theta", "pressure"}) {
        const Real minimum = plotfile.get(0, field).min(0, 0, false);
        if (minimum <= 0.0) {
            std::cerr << label << " field '" << field
                      << "' has non-positive minimum " << minimum << "\n";
            return false;
        }
    }
    for (const auto& field : {"qv", "qc", "rhoKE"}) {
        const Real minimum = plotfile.get(0, field).min(0, 0, false);
        if (minimum < -allowance) {
            std::cerr << label << " field '" << field
                      << "' is below the allowed roundoff floor: minimum="
                      << minimum << ", allowance=" << allowance << "\n";
            return false;
        }
    }
    return true;
}

bool check_compatible_layout (const MultiFab& lhs,
                              const MultiFab& rhs,
                              const std::string& context)
{
    if (lhs.nComp() < 1 || rhs.nComp() < 1) {
        std::cerr << context << " has a field with no components\n";
        return false;
    }
    if (lhs.nComp() != rhs.nComp()) {
        std::cerr << context << " has component-count mismatch: "
                  << lhs.nComp() << " versus " << rhs.nComp() << "\n";
        return false;
    }
    if (lhs.boxArray() != rhs.boxArray()) {
        std::cerr << context << " has BoxArray mismatch\n";
        return false;
    }
    if (lhs.DistributionMap() != rhs.DistributionMap()) {
        std::cerr << context << " has DistributionMapping mismatch\n";
        return false;
    }
    return true;
}

} // namespace

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);

    Arguments args;
    std::string argument_error;
    if (!read_arguments(argc, argv, args, argument_error)) {
        std::cerr << "SHOC mutation checker argument error: " << argument_error << "\n"
                  << "usage: --field FIELD --min-final-difference VALUE "
                     "--min-baseline-evolution VALUE "
                     "--max-mutant-to-baseline-evolution-ratio VALUE "
                     "--baseline-initial PF --mutant-initial PF "
                     "--baseline-final PF --mutant-final PF\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    PlotFileData baseline_initial(args.baseline_initial);
    PlotFileData mutant_initial(args.mutant_initial);
    PlotFileData baseline_final(args.baseline_final);
    PlotFileData mutant_final(args.mutant_final);

    if (!check_finite_fields(baseline_initial, host_state_fields(), "baseline initial") ||
        !check_finite_fields(mutant_initial, host_state_fields(), "mutant initial") ||
        !check_finite_fields(baseline_final, host_state_fields(), "baseline final") ||
        !check_finite_fields(mutant_final, host_state_fields(), "mutant final") ||
        !check_physical_final_state(baseline_final, "baseline final") ||
        !check_physical_final_state(mutant_final, "mutant final")) {
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    constexpr Real initial_tolerance = 2.0e-10;
    Real maximum_initial_difference = 0.0;
    for (const auto& field : host_state_fields()) {
        if (!check_compatible_layout(
                baseline_initial.get(0, field), mutant_initial.get(0, field),
                "baseline/mutant initial field '" + field + "'") ||
            !check_compatible_layout(
                baseline_initial.get(0, field), baseline_final.get(0, field),
                "baseline initial/final field '" + field + "'") ||
            !check_compatible_layout(
                mutant_initial.get(0, field), mutant_final.get(0, field),
                "mutant initial/final field '" + field + "'")) {
            amrex::Finalize();
            return EXIT_FAILURE;
        }
        maximum_initial_difference = std::max(
            maximum_initial_difference,
            maximum_difference(baseline_initial.get(0, field), mutant_initial.get(0, field)));
    }
    if (!check_compatible_layout(
            baseline_final.get(0, args.field), mutant_final.get(0, args.field),
            "baseline/mutant final target field '" + args.field + "'")) {
        amrex::Finalize();
        return EXIT_FAILURE;
    }
    if (maximum_initial_difference > initial_tolerance) {
        std::cerr << "paired mutation initial states differ by "
                  << maximum_initial_difference << " > " << initial_tolerance << "\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    const Real final_difference = maximum_difference(
        baseline_final.get(0, args.field), mutant_final.get(0, args.field));
    const Real baseline_evolution = maximum_difference(
        baseline_initial.get(0, args.field), baseline_final.get(0, args.field));
    const Real mutant_evolution = maximum_difference(
        mutant_initial.get(0, args.field), mutant_final.get(0, args.field));

    amrex::Print() << "SHOC paired mutation: field=" << args.field
                   << " initial_difference=" << maximum_initial_difference
                   << " final_difference=" << final_difference
                   << " baseline_evolution=" << baseline_evolution
                   << " mutant_evolution=" << mutant_evolution
                   << " min_final_difference=" << args.min_final_difference
                   << " min_baseline_evolution=" << args.min_baseline_evolution
                   << " max_mutant_to_baseline_evolution_ratio="
                   << args.max_mutant_to_baseline_evolution_ratio << "\n";

    if (final_difference < args.min_final_difference) {
        std::cerr << "debug mutation final difference for '" << args.field
                  << "' is below the configured minimum: observed=" << final_difference
                  << ", minimum=" << args.min_final_difference << "\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }
    if (baseline_evolution < args.min_baseline_evolution) {
        std::cerr << "baseline evolution for '" << args.field
                  << "' is below the configured minimum: observed=" << baseline_evolution
                  << ", minimum=" << args.min_baseline_evolution << "\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }
    if (mutant_evolution >
        args.max_mutant_to_baseline_evolution_ratio * baseline_evolution) {
        std::cerr << "debug mutation mutant evolution for '" << args.field
                  << "' exceeds the configured baseline ratio: mutant=" << mutant_evolution
                  << ", baseline=" << baseline_evolution
                  << ", ratio_limit="
                  << args.max_mutant_to_baseline_evolution_ratio << "\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    amrex::Finalize();
    return EXIT_SUCCESS;
}
