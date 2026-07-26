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
    Real threshold = 0.0;
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

bool read_arguments (int argc, char** argv, Arguments& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string option(argv[i]);
        if (i + 1 >= argc) return false;
        const std::string value(argv[++i]);
        if (option == "--field") {
            args.field = value;
        } else if (option == "--threshold") {
            args.threshold = std::stod(value);
        } else if (option == "--baseline-initial") {
            args.baseline_initial = value;
        } else if (option == "--mutant-initial") {
            args.mutant_initial = value;
        } else if (option == "--baseline-final") {
            args.baseline_final = value;
        } else if (option == "--mutant-final") {
            args.mutant_final = value;
        } else {
            return false;
        }
    }
    return !args.field.empty() && args.threshold > 0.0 &&
           !args.baseline_initial.empty() && !args.mutant_initial.empty() &&
           !args.baseline_final.empty() && !args.mutant_final.empty();
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

} // namespace

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);

    Arguments args;
    if (!read_arguments(argc, argv, args)) {
        std::cerr << "usage: --field FIELD --threshold VALUE "
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
        !check_finite_fields(baseline_final, {args.field}, "baseline final") ||
        !check_finite_fields(mutant_final, {args.field}, "mutant final")) {
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    constexpr Real initial_tolerance = 2.0e-10;
    constexpr Real reduction_factor = 0.5;
    Real maximum_initial_difference = 0.0;
    for (const auto& field : host_state_fields()) {
        maximum_initial_difference = std::max(
            maximum_initial_difference,
            maximum_difference(baseline_initial.get(0, field), mutant_initial.get(0, field)));
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
                   << " threshold=" << args.threshold << "\n";

    if (final_difference <= args.threshold) {
        std::cerr << "debug mutation did not materially change " << args.field << "\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }
    if (mutant_evolution >= reduction_factor * baseline_evolution) {
        std::cerr << "debug mutation did not reduce " << args.field
                  << " evolution by the required factor\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    amrex::Finalize();
    return EXIT_SUCCESS;
}
