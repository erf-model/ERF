#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

// Stable exit codes let CTest distinguish a broken simulation from a
// deliberately broken state-update contract.
enum CheckCode : int {
    kSuccess = 0,
    kUsage = 2,
    kMissingField = 3,
    kInvalidField = 4,
    kTransportContract = 5,
    kCloudContract = 6
};

// The checker is intentionally driven by three plotfile snapshots so that
// every oracle can inspect both the initial state and the resulting evolution.
struct Arguments {
    std::string mode;
    std::string initial;
    std::string midpoint;
    std::string final;
};

// Return whether a named field is present without asking PlotFileData to throw
// for an absent diagnostic.  This keeps missing-output failures actionable.
bool
has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

// Parse the small command-line interface used by RunShocRegression.cmake.
// Keeping parsing here makes the checker independently runnable when debugging
// a failed plotfile contract.
bool
read_arguments (int argc, char** argv, Arguments& args)
{
    for (int i = 1; i < argc; ++i) {
        const std::string option(argv[i]);
        if (i + 1 >= argc) {
            std::cerr << "missing value for " << option << "\n";
            return false;
        }
        if (option == "--mode") {
            args.mode = argv[++i];
        } else if (option == "--initial") {
            args.initial = argv[++i];
        } else if (option == "--midpoint") {
            args.midpoint = argv[++i];
        } else if (option == "--final") {
            args.final = argv[++i];
        } else {
            std::cerr << "unknown option: " << option << "\n";
            return false;
        }
    }

    return !args.mode.empty() && !args.initial.empty() &&
           !args.midpoint.empty() && !args.final.empty();
}

// Verify the fields shared by every SHOC snapshot exist and contain finite
// values.  Physical bounds are checked separately so each failure can report
// the most relevant contract category.
bool
check_finite (PlotFileData& plotfile,
              const std::vector<std::string>& required,
              CheckCode& code,
              std::string& error)
{
    for (const auto& name : required) {
        if (!has_variable(plotfile, name)) {
            code = kMissingField;
            error = "required field '" + name + "' is missing";
            return false;
        }

        const MultiFab field = plotfile.get(0, name);
        if (!field.is_finite(0, 1, 0, false)) {
            code = kInvalidField;
            error = "field '" + name + "' contains a non-finite value";
            return false;
        }
    }

    return true;
}

// Compute the cellwise L-infinity difference between two same-layout fields.
// This is used for both evolution checks and mutation regressions.
Real
maximum_difference (const MultiFab& first, const MultiFab& second)
{
    MultiFab difference(first.boxArray(), first.DistributionMap(), 1, 0);
    MultiFab::Copy(difference, first, 0, 0, 1, 0);
    MultiFab::Subtract(difference, second, 0, 0, 1, 0);
    return difference.norm0();
}

// Enforce non-negativity for mass, energy, and diffusivity diagnostics while
// allowing the tiny negative roundoff produced by floating-point arithmetic.
bool
check_nonnegative (PlotFileData& plotfile,
                   const std::string& name,
                   const Real allowance,
                   CheckCode& code,
                   std::string& error)
{
    const Real minimum = plotfile.get(0, name).min(0, 0, false);
    if (minimum < -allowance) {
        code = kInvalidField;
        error = "field '" + name + "' has minimum " + std::to_string(minimum) +
                ", below allowed negative roundoff " + std::to_string(allowance);
        return false;
    }
    return true;
}

// Check hydrometeor fields only when the selected plotfile actually contains
// them.  The Kessler and WSM6 fixtures request their supported species, so
// this path is real coverage rather than a vacuous list of optional names.
bool
check_optional_hydrometeors (PlotFileData& plotfile,
                             CheckCode& code,
                             std::string& error)
{
    std::string checked;
    for (const auto& name : {"qi", "qrain", "qsnow", "qgraup"}) {
        if (has_variable(plotfile, name)) {
            const std::vector<std::string> field{name};
            if (!check_finite(plotfile, field, code, error) ||
                !check_nonnegative(plotfile, name, 1.0e-12, code, error)) {
                return false;
            }
            if (!checked.empty()) { checked += ", "; }
            checked += name;
        }
    }

    if (!checked.empty()) {
        amrex::Print() << "SHOC property check: verified hydrometeor fields "
                       << checked << "\n";
    }
    return true;
}

// Scheme-specific fixtures must name their hydrometeor inventory explicitly;
// otherwise a future fixture/gold edit could make an optional loop vacuous.
bool
check_required_hydrometeors (PlotFileData& plotfile,
                             const std::vector<std::string>& required,
                             const std::string& label,
                             CheckCode& code,
                             std::string& error)
{
    for (const auto& name : required) {
        if (!has_variable(plotfile, name)) {
            code = kMissingField;
            error = label + " is missing required hydrometeor field '" + name + "'";
            return false;
        }
        const MultiFab field = plotfile.get(0, name);
        if (!field.is_finite(0, 1, 0, false)) {
            code = kInvalidField;
            error = label + " hydrometeor field '" + name + "' contains a non-finite value";
            return false;
        }
        const Real minimum = plotfile.get(0, name).min(0, 0, false);
        if (minimum < -1.0e-12) {
            code = kInvalidField;
            error = label + " hydrometeor field '" + name + "' has minimum " +
                    std::to_string(minimum) + ", below allowed negative roundoff";
            return false;
        }
    }
    return true;
}

std::vector<std::string>
required_hydrometeors_for_mode (const std::string& mode)
{
    if (mode == "unstable_cloud_kessler") {
        return {"qrain"};
    }
    if (mode == "unstable_cloud_wsm6") {
        return {"qi", "qrain", "qsnow", "qgraup"};
    }
    return {};
}

// Check the host state that is meaningful at initialization.  Native SHOC
// diagnostics are intentionally absent here: the first plotfile is written
// before the SHOC driver has diagnosed its PDF and transport fields.
bool
check_initial_state (PlotFileData& plotfile,
                     const std::string& label,
                     const std::vector<std::string>& required_hydrometeors,
                     CheckCode& code,
                     std::string& error)
{
    const std::vector<std::string> required {
        "density", "temp", "theta", "pressure", "qv", "qc", "rhoKE",
        "x_velocity", "y_velocity"
    };
    if (!check_finite(plotfile, required, code, error)) {
        return false;
    }

    for (const auto& name : {"density", "temp", "theta", "pressure"}) {
        if (plotfile.get(0, name).min(0, 0, false) <= 0.0) {
            code = kInvalidField;
            error = "field '" + std::string(name) + "' is not strictly positive";
            return false;
        }
    }

    for (const auto& name : {"density", "qv", "qc", "rhoKE"}) {
        if (!check_nonnegative(plotfile, name, 1.0e-12, code, error)) {
            return false;
        }
    }

    if (!required_hydrometeors.empty()) {
        return check_required_hydrometeors(
            plotfile, required_hydrometeors, label, code, error);
    }
    return check_optional_hydrometeors(plotfile, code, error);
}

// A native diagnostic is unavailable when the writer has no SHOC source for
// it, and ERF represents that state with exactly -999.  Keep this sentinel
// distinct from a diagnosed physical zero at saved times after SHOC runs.
bool
check_not_unavailable (PlotFileData& plotfile,
                       const std::string& name,
                       CheckCode& code,
                       std::string& error)
{
    const Real minimum = plotfile.get(0, name).min(0, 0, false);
    if (std::abs(minimum + 999.0) <= 1.0e-12) {
        code = kMissingField;
        error = "SHOC diagnostic '" + name + "' retains ERF unavailable sentinel -999";
        return false;
    }
    return true;
}

// Require the diagnostics that are guaranteed after the midpoint and final
// SHOC calls.  This is the post-lifecycle counterpart to check_initial_state.
bool
check_available_shoc_diagnostics (PlotFileData& plotfile,
                                  CheckCode& code,
                                  std::string& error)
{
    const std::vector<std::string> diagnostics {
        "Kmv", "Khv", "Lturb", "shoc_cldfrac", "shoc_ql", "shoc_cond",
        "brunt", "shear_prod", "buoy_prod", "diss_tke"
    };
    if (!check_finite(plotfile, diagnostics, code, error)) {
        return false;
    }
    for (const auto& name : diagnostics) {
        if (!check_not_unavailable(plotfile, name, code, error)) {
            return false;
        }
    }

    for (const auto& name : {"Kmv", "Khv", "Lturb", "shoc_ql", "shoc_cond"}) {
        if (!check_nonnegative(plotfile, name, 1.0e-12, code, error)) {
            return false;
        }
    }

    const auto& cldfrac = plotfile.get(0, "shoc_cldfrac");
    const Real cldfrac_min = cldfrac.min(0, 0, false);
    const Real cldfrac_max = cldfrac.max(0, 0, false);
    if (cldfrac_min < -1.0e-12 || cldfrac_max > 1.0 + 1.0e-12) {
        code = kInvalidField;
        error = "field 'shoc_cldfrac' is outside [0,1]: min=" +
                std::to_string(cldfrac_min) + ", max=" + std::to_string(cldfrac_max);
        return false;
    }

    return true;
}

// Check a saved state after SHOC has run: host-state validity, available
// native diagnostics, and any selected scheme-specific hydrometeors.
bool
check_post_shoc_state (PlotFileData& plotfile,
                       const std::string& label,
                       const std::vector<std::string>& required_hydrometeors,
                       CheckCode& code,
                       std::string& error)
{
    return check_initial_state(
               plotfile, label, required_hydrometeors, code, error) &&
           check_available_shoc_diagnostics(plotfile, code, error);
}

// Confirm that the run exercised the native SHOC transport path rather than
// merely producing a numerically valid but unchanged plotfile.
bool
check_active_transport (PlotFileData& initial,
                        PlotFileData& final,
                        CheckCode& code,
                        std::string& error)
{
    const Real momentum_change = maximum_difference(initial.get(0, "x_velocity"),
                                                    final.get(0, "x_velocity"));
    const Real tke_change = maximum_difference(initial.get(0, "rhoKE"),
                                               final.get(0, "rhoKE"));
    const Real thermal_change = maximum_difference(initial.get(0, "theta"),
                                                   final.get(0, "theta"));
    if (momentum_change <= 1.0e-12 || thermal_change <= 1.0e-12) {
        code = kTransportContract;
        error = "native SHOC transport did not evolve momentum and thermal state: "
                "max |delta x_velocity|=" + std::to_string(momentum_change) +
                ", max |delta theta|=" + std::to_string(thermal_change);
        return false;
    }
    if (tke_change <= 1.0e-12) {
        code = kTransportContract;
        error = "TKE did not change (max |delta rhoKE| = " +
                std::to_string(tke_change) + ")";
        return false;
    }
    if (final.get(0, "Kmv").max(0, 0, false) <= 0.0 ||
        final.get(0, "Khv").max(0, 0, false) <= 0.0 ||
        final.get(0, "Lturb").max(0, 0, false) <= 0.0) {
        code = kTransportContract;
        error = "native SHOC diffusivity or mixing-length diagnostics are trivial";
        return false;
    }
    return true;
}

// Clear fixtures are expected to remain free of condensed water at every
// checkpoint.  This oracle intentionally uses qc only because SHOC PDF
// diagnostics are not valid at initialization.
bool
check_clear_snapshot (PlotFileData& plotfile, const std::string& label,
                      CheckCode& code, std::string& error)
{
    const Real qc = plotfile.get(0, "qc").max(0, 0, false);
    if (qc > 1.0e-10) {
        code = kCloudContract;
        error = label + " snapshot is not clear: max qc=" + std::to_string(qc);
        return false;
    }
    return true;
}

// Apply the mode-specific oracle after common state and transport checks have
// passed.  Cloud modes require formation and evolution, clear modes require
// the requested stratification sign and the corresponding production sign.
bool
check_regime (const Arguments& args,
              PlotFileData& initial,
              PlotFileData& midpoint,
              PlotFileData& final,
              CheckCode& code,
              std::string& error)
{
    if (!check_active_transport(initial, final, code, error)) {
        return false;
    }

    if (args.mode == "stable_clear" || args.mode == "unstable_clear") {
        if (!check_clear_snapshot(initial, "initial", code, error) ||
            !check_clear_snapshot(midpoint, "midpoint", code, error) ||
            !check_clear_snapshot(final, "final", code, error)) {
            return false;
        }
        const Real brunt_min = final.get(0, "brunt").min(0, 0, false);
        const Real brunt_max = final.get(0, "brunt").max(0, 0, false);
        const Real buoy_min = final.get(0, "buoy_prod").min(0, 0, false);
        const Real buoy_max = final.get(0, "buoy_prod").max(0, 0, false);
        if (args.mode == "stable_clear" && (brunt_max <= 0.0 || buoy_min >= 0.0)) {
            code = kCloudContract;
            error = "stable case lacks positive stratification and negative buoyancy production: "
                    "brunt_max=" + std::to_string(brunt_max) +
                    ", buoy_min=" + std::to_string(buoy_min);
            return false;
        }
        if (args.mode == "unstable_clear" && (brunt_min >= 0.0 || buoy_max <= 0.0)) {
            code = kCloudContract;
            error = "unstable case lacks negative stratification and positive buoyancy production: "
                    "brunt_min=" + std::to_string(brunt_min) +
                    ", buoy_max=" + std::to_string(buoy_max);
            return false;
        }
        amrex::Print() << "SHOC property check: mode=" << args.mode
                       << " clear snapshots verified at initial/midpoint/final\n";
        return true;
    }

    if (args.mode == "stable_cloud" || args.mode == "unstable_cloud" ||
        args.mode == "unstable_cloud_nocond" ||
        args.mode == "unstable_cloud_kessler" ||
        args.mode == "unstable_cloud_wsm6") {
        const Real qc_initial = initial.get(0, "qc").max(0, 0, false);
        const Real qc_midpoint = midpoint.get(0, "qc").max(0, 0, false);
        const Real qc_final = final.get(0, "qc").max(0, 0, false);
        const Real cldfrac_final = final.get(0, "shoc_cldfrac").max(0, 0, false);
        const Real shoc_ql_final = final.get(0, "shoc_ql").max(0, 0, false);
        const Real qc_change = maximum_difference(midpoint.get(0, "qc"), final.get(0, "qc"));
        if (qc_initial > 1.0e-12) {
            code = kCloudContract;
            error = "cloud-formation fixture is not initialized clear: max qc(t=0) = " +
                    std::to_string(qc_initial);
            return false;
        }
        if (qc_final <= 1.0e-8 || shoc_ql_final <= 1.0e-8 || cldfrac_final <= 1.0e-6) {
            code = kCloudContract;
            error = "cloud-formation oracle failed: max qc=" + std::to_string(qc_final) +
                    ", max shoc_ql=" + std::to_string(shoc_ql_final) +
                    ", max shoc_cldfrac=" + std::to_string(cldfrac_final);
            return false;
        }
        if (qc_change <= 1.0e-10 || qc_midpoint <= 1.0e-8) {
            code = kCloudContract;
            error = "cloud water did not form and evolve through midpoint/final outputs: "
                    "max qc(midpoint)=" + std::to_string(qc_midpoint) +
                    ", max delta qc=" + std::to_string(qc_change);
            return false;
        }
        amrex::Print() << "SHOC property check: mode=" << args.mode
                       << " qc_initial=" << qc_initial
                       << " qc_midpoint=" << qc_midpoint
                       << " qc_final=" << qc_final << "\n";
        return true;
    }

    code = kUsage;
    error = "unsupported SHOC property-check mode '" + args.mode + "'";
    return false;
}

} // namespace

int
main (int argc, char** argv)
{
    // This executable is a contract checker, not a statistical comparison:
    // it validates physical invariants and state evolution before fcompare
    // performs the narrow numerical gold-file comparison.
    amrex::Initialize(argc, argv, false);

    Arguments args;
    CheckCode code = kUsage;
    std::string error;
    bool ok = read_arguments(argc, argv, args);
    if (ok) {
        const auto required_hydrometeors = required_hydrometeors_for_mode(args.mode);
        PlotFileData initial(args.initial);
        PlotFileData midpoint(args.midpoint);
        PlotFileData final(args.final);
        ok = check_initial_state(
                 initial, "initial", required_hydrometeors, code, error) &&
             check_post_shoc_state(
                 midpoint, "midpoint", required_hydrometeors, code, error) &&
             check_post_shoc_state(
                 final, "final", required_hydrometeors, code, error) &&
             check_regime(args, initial, midpoint, final, code, error);
        if (ok) code = kSuccess;
    }

    if (!ok) {
        std::cerr << "SHOC property check failed (code " << static_cast<int>(code) << "): "
                  << (error.empty() ? "usage: --mode MODE --initial PF --midpoint PF --final PF" : error)
                  << "\n";
    }

    amrex::Finalize();
    return static_cast<int>(code);
}
