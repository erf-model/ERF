#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

struct Arguments {
    std::string mode;
    std::string initial;
    std::string midpoint;
    std::string final;
};

bool
has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

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

bool
check_finite_and_bounds (PlotFileData& plotfile,
                         const std::vector<std::string>& required,
                         std::string& error)
{
    for (const auto& name : required) {
        if (!has_variable(plotfile, name)) {
            error = "required field '" + name + "' is missing from " +
                    std::to_string(plotfile.time());
            return false;
        }

        const MultiFab field = plotfile.get(0, name);
        if (!field.is_finite(0, 1, 0, false)) {
            error = "field '" + name + "' contains a non-finite value";
            return false;
        }
    }

    return true;
}

Real
maximum_difference (const MultiFab& first, const MultiFab& second)
{
    Real result = 0.0;
    for (amrex::MFIter mfi(first); mfi.isValid(); ++mfi) {
        const auto& a = first.const_array(mfi);
        const auto& b = second.const_array(mfi);
        const auto& box = mfi.validbox();
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    result = std::max(result, std::abs(a(i,j,k) - b(i,j,k)));
                }
            }
        }
    }
    return result;
}

bool
check_nonnegative (PlotFileData& plotfile,
                   const std::string& name,
                   const Real allowance,
                   std::string& error)
{
    const MultiFab field = plotfile.get(0, name);
    const Real minimum = field.min(0, 0, false);
    if (minimum < -allowance) {
        error = "field '" + name + "' has minimum " + std::to_string(minimum) +
                ", below allowed negative roundoff " + std::to_string(allowance);
        return false;
    }
    return true;
}

bool
check_common_state (PlotFileData& plotfile, bool require_qc, std::string& error)
{
    const std::vector<std::string> required {
        "density", "temp", "theta", "pressure", "qv", "rhoKE",
        "x_velocity", "y_velocity", "Kmv", "Khv", "Lturb", "shoc_cldfrac",
        "shoc_ql", "shoc_cond", "brunt", "shear_prod", "buoy_prod", "diss_tke"
    };
    if (!check_finite_and_bounds(plotfile, required, error)) {
        return false;
    }
    if (require_qc) {
        if (std::find(plotfile.varNames().begin(), plotfile.varNames().end(), "qc") ==
            plotfile.varNames().end()) {
            error = "required field 'qc' is missing";
            return false;
        }
    }

    const std::vector<std::string> nonnegative {
        "density", "qv", "rhoKE", "Kmv", "Khv", "Lturb"
    };
    for (const auto& name : nonnegative) {
        if (!check_nonnegative(plotfile, name, 1.0e-12, error)) {
            return false;
        }
    }

    return true;
}

bool
check_active_transport (PlotFileData& initial,
                        PlotFileData& final,
                        std::string& error)
{
    const Real momentum_change = maximum_difference(initial.get(0, "x_velocity"),
                                                    final.get(0, "x_velocity"));
    const Real tke_change = maximum_difference(initial.get(0, "rhoKE"),
                                               final.get(0, "rhoKE"));
    const Real thermal_change = maximum_difference(initial.get(0, "theta"),
                                                   final.get(0, "theta"));
    if (momentum_change <= 1.0e-12) {
        error = "horizontal momentum did not change (max |delta x_velocity| = " +
                std::to_string(momentum_change) + ")";
        return false;
    }
    if (tke_change <= 1.0e-12) {
        error = "TKE did not change (max |delta rhoKE| = " +
                std::to_string(tke_change) + ")";
        return false;
    }
    if (thermal_change <= 1.0e-12) {
        error = "potential temperature did not change (max |delta theta| = " +
                std::to_string(thermal_change) + ")";
        return false;
    }
    if (final.get(0, "Kmv").max(0, 0, false) <= 0.0 ||
        final.get(0, "Khv").max(0, 0, false) <= 0.0 ||
        final.get(0, "Lturb").max(0, 0, false) <= 0.0) {
        error = "native SHOC diffusivity or mixing-length diagnostics are trivial";
        return false;
    }
    return true;
}

bool
check_regime (const Arguments& args,
              PlotFileData& initial,
              PlotFileData& midpoint,
              PlotFileData& final,
              std::string& error)
{
    if (!check_active_transport(initial, final, error)) {
        return false;
    }

    if (args.mode == "unstable_cloud_nocond") {
        amrex::Print() << "SHOC property check: mode=" << args.mode
                       << " qv_initial=" << initial.get(0, "qv").max(0, 0, false)
                       << " qv_final=" << final.get(0, "qv").max(0, 0, false) << "\n";
        return true;
    }

    if (args.mode == "negative_tke") {
        const Real tke_change = maximum_difference(initial.get(0, "rhoKE"),
                                                   final.get(0, "rhoKE"));
        if (tke_change <= 1.0e-3) {
            error = "negative control detected insufficient TKE state evolution: max |delta rhoKE| = " +
                    std::to_string(tke_change);
            return false;
        }
    }
    if (args.mode == "negative_theta") {
        const Real thermal_change = maximum_difference(initial.get(0, "theta"),
                                                       final.get(0, "theta"));
        if (thermal_change <= 1.0e-3) {
            error = "negative control detected insufficient thermal state evolution: max |delta theta| = " +
                    std::to_string(thermal_change);
            return false;
        }
    }

    const Real qc_initial = initial.get(0, "qc").max(0, 0, false);
    const Real qc_midpoint = midpoint.get(0, "qc").max(0, 0, false);
    const Real qc_final = final.get(0, "qc").max(0, 0, false);
    const Real cldfrac_final = final.get(0, "shoc_cldfrac").max(0, 0, false);
    const Real shoc_ql_final = final.get(0, "shoc_ql").max(0, 0, false);
    const Real condensation_midpoint = midpoint.get(0, "shoc_cond").max(0, 0, false);
    const Real condensation_final = final.get(0, "shoc_cond").max(0, 0, false);
    const Real qc_change = maximum_difference(midpoint.get(0, "qc"), final.get(0, "qc"));

    if (args.mode == "stable_clear" || args.mode == "unstable_clear") {
        if (args.mode == "stable_clear" && qc_final > 1.0e-10) {
            error = "clear-case cloud water reached " + std::to_string(qc_final) +
                    ", above 1.0e-10 tolerance";
            return false;
        }
        const Real brunt_min = final.get(0, "brunt").min(0, 0, false);
        const Real brunt_max = final.get(0, "brunt").max(0, 0, false);
        const Real buoy_min = final.get(0, "buoy_prod").min(0, 0, false);
        const Real buoy_max = final.get(0, "buoy_prod").max(0, 0, false);
        if (args.mode == "stable_clear" && (brunt_max <= 0.0 || buoy_min >= 0.0)) {
            error = "stable case lacks positive stratification and negative buoyancy production";
            return false;
        }
        if (args.mode == "unstable_clear" && (brunt_min >= 0.0 || buoy_max <= 0.0)) {
            error = "unstable case lacks negative stratification and positive buoyancy production";
            return false;
        }
        amrex::Print() << "SHOC property check: mode=" << args.mode
                       << " qc_initial=" << qc_initial
                       << " qc_midpoint=" << qc_midpoint
                       << " qc_final=" << qc_final << "\n";
        return true;
    }

    if (args.mode == "stable_cloud" || args.mode == "unstable_cloud") {
        if (qc_initial > 1.0e-12) {
            error = "cloud-formation fixture is not initialized clear: max qc(t=0) = " +
                    std::to_string(qc_initial);
            return false;
        }
        if (qc_final <= 1.0e-8 || shoc_ql_final <= 1.0e-8 || cldfrac_final <= 1.0e-6) {
            error = "cloud-formation oracle failed: max qc=" + std::to_string(qc_final) +
                    ", max shoc_ql=" + std::to_string(shoc_ql_final) +
                    ", max shoc_cldfrac=" + std::to_string(cldfrac_final) +
                    ", max shoc_cond(mid/final)=" + std::to_string(condensation_midpoint) +
                    "/" + std::to_string(condensation_final);
            return false;
        }
        if (qc_change <= 1.0e-10) {
            error = "cloud water did not evolve between midpoint and final output: max delta qc=" +
                    std::to_string(qc_change);
            return false;
        }
        amrex::Print() << "SHOC property check: mode=" << args.mode
                       << " qc_initial=" << qc_initial
                       << " qc_midpoint=" << qc_midpoint
                       << " qc_final=" << qc_final
                       << " shoc_ql_final=" << shoc_ql_final
                       << " shoc_cond_midpoint/final=" << condensation_midpoint
                       << "/" << condensation_final << "\n";
        return true;
    }

    error = "unsupported SHOC property-check mode '" + args.mode + "'";
    return false;
}

} // namespace

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);

    Arguments args;
    bool ok = read_arguments(argc, argv, args);
    std::string error;
    if (ok) {
        PlotFileData initial(args.initial);
        PlotFileData midpoint(args.midpoint);
        PlotFileData final(args.final);
        const bool require_qc = args.mode != "unstable_cloud_nocond";
        ok = check_common_state(initial, require_qc, error) &&
             check_common_state(midpoint, require_qc, error) &&
             check_common_state(final, require_qc, error) &&
             check_regime(args, initial, midpoint, final, error);
    }

    if (!ok) {
        std::cerr << "SHOC property check failed: "
                  << (error.empty() ? "usage: --mode MODE --initial PF --midpoint PF --final PF" : error)
                  << "\n";
    }

    amrex::Finalize();
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
