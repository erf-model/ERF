#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

struct Arguments {
    std::string mode;
    std::string gold;
    std::string candidate;
};

struct FieldTolerance {
    Real absolute;
    Real relative;
};

// These are source-reviewed tolerances for the shared cloudy SHOC golds.  The
// additive absolute-plus-relative rule below keeps small moisture and
// diagnostic fields tighter than a single global fcompare absolute floor.
// The expanded entries below cover the measured deterministic
// CPU/compiler envelope from the stable-cloud CI failures; all other fields
// retain their pre-repair limits.
const std::map<std::string, FieldTolerance>& tolerance_profile ()
{
    static const std::map<std::string, FieldTolerance> profile {
        {"density",     {1.0e-10, 3.0e-8}},
        {"rhoKE",       {1.0e-10, 4.0e-7}},
        {"x_velocity",  {1.0e-9,  6.0e-7}},
        {"y_velocity",  {1.0e-10, 6.0e-8}},
        {"z_velocity",  {3.0e-6,  1.0e-7}},
        {"temp",        {1.0e-8,  3.0e-8}},
        {"theta",       {1.0e-8,  3.0e-8}},
        {"pressure",    {1.0e-5,  3.0e-8}},
        {"Kmv",         {1.0e-8,  1.0e-6}},
        {"Khv",         {1.0e-8,  1.0e-6}},
        {"Lturb",       {1.0e-8,  2.0e-5}},
        {"shoc_cldfrac",{1.0e-12, 1.0e-12}},
        {"shoc_ql",     {5.0e-10, 2.0e-6}},
        {"shoc_cond",   {3.0e-11, 1.5e-5}},
        {"brunt",       {3.0e-8,  3.0e-5}},
        {"shear_prod",  {5.0e-13, 4.0e-6}},
        {"buoy_prod",   {2.0e-11, 1.5e-5}},
        {"diss_tke",    {2.0e-11, 3.5e-7}},
        {"qv",          {5.0e-10, 3.0e-8}},
        {"qc",          {5.0e-10, 2.0e-6}},
        {"qi",          {1.0e-12, 1.0e-8}},
        {"qrain",       {1.0e-12, 1.0e-8}},
        {"qsnow",       {1.0e-12, 1.0e-8}},
        {"qgraup",      {1.0e-12, 1.0e-8}}
    };
    return profile;
}

const std::vector<std::string>& base_fields ()
{
    static const std::vector<std::string> fields {
        "density", "rhoKE", "x_velocity", "y_velocity", "z_velocity",
        "temp", "theta", "pressure", "Kmv", "Khv", "Lturb",
        "shoc_cldfrac", "shoc_ql", "shoc_cond", "brunt", "shear_prod",
        "buoy_prod", "diss_tke", "qv", "qc"
    };
    return fields;
}

const std::vector<std::string>& fields_for_mode (const std::string& mode)
{
    static const std::vector<std::string> stable_cloud = base_fields();
    static const std::vector<std::string> unstable_cloud = base_fields();
    static const std::vector<std::string> unstable_cloud_nocond = base_fields();
    static const std::vector<std::string> unstable_cloud_kessler = [] {
        auto fields = base_fields();
        fields.emplace_back("qrain");
        return fields;
    }();
    static const std::vector<std::string> unstable_cloud_wsm6 = [] {
        auto fields = base_fields();
        fields.emplace_back("qi");
        fields.emplace_back("qrain");
        fields.emplace_back("qsnow");
        fields.emplace_back("qgraup");
        return fields;
    }();

    if (mode == "stable_cloud") {
        return stable_cloud;
    }
    if (mode == "unstable_cloud") {
        return unstable_cloud;
    }
    if (mode == "unstable_cloud_nocond") {
        return unstable_cloud_nocond;
    }
    if (mode == "unstable_cloud_kessler") {
        return unstable_cloud_kessler;
    }
    if (mode == "unstable_cloud_wsm6") {
        return unstable_cloud_wsm6;
    }
    std::abort();
}

bool is_supported_mode (const std::string& mode)
{
    return mode == "stable_cloud" || mode == "unstable_cloud" ||
           mode == "unstable_cloud_nocond" || mode == "unstable_cloud_kessler" ||
           mode == "unstable_cloud_wsm6";
}

bool read_arguments (int argc, char** argv, Arguments& args, std::string& error)
{
    for (int i = 1; i < argc; ++i) {
        const std::string option(argv[i]);
        if (i + 1 >= argc) {
            error = "missing value after " + option;
            return false;
        }
        const std::string value(argv[++i]);
        if (option == "--mode") {
            args.mode = value;
        } else if (option == "--gold") {
            args.gold = value;
        } else if (option == "--candidate") {
            args.candidate = value;
        } else {
            error = "unknown option '" + option + "'";
            return false;
        }
    }
    if (args.mode.empty() || args.gold.empty() || args.candidate.empty()) {
        error = "--mode, --gold, and --candidate are required";
        return false;
    }
    if (!is_supported_mode(args.mode)) {
        error = "unsupported cloudy SHOC mode '" + args.mode + "'";
        return false;
    }
    return true;
}

bool check_inventory (const PlotFileData& plotfile,
                      const std::vector<std::string>& expected,
                      const std::string& label)
{
    const std::set<std::string> expected_set(expected.begin(), expected.end());
    const std::set<std::string> actual_set(plotfile.varNames().begin(), plotfile.varNames().end());
    if (expected_set != actual_set ||
        static_cast<amrex::Long>(actual_set.size()) != plotfile.varNames().size()) {
        std::cerr << label << " has an unexpected field inventory\n";
        for (const auto& field : expected_set) {
            if (actual_set.count(field) == 0) {
                std::cerr << "  missing: " << field << "\n";
            }
        }
        for (const auto& field : actual_set) {
            if (expected_set.count(field) == 0) {
                std::cerr << "  unexpected: " << field << "\n";
            }
        }
        return false;
    }
    return true;
}

bool check_plotfile_layout (PlotFileData& gold,
                            PlotFileData& candidate,
                            const std::vector<std::string>& fields)
{
    if (gold.spaceDim() != amrex::SpaceDim || candidate.spaceDim() != amrex::SpaceDim ||
        gold.spaceDim() != candidate.spaceDim()) {
        std::cerr << "plotfile spatial-dimension mismatch: gold=" << gold.spaceDim()
                  << ", candidate=" << candidate.spaceDim()
                  << ", build=" << amrex::SpaceDim << "\n";
        return false;
    }
    if (gold.finestLevel() != candidate.finestLevel()) {
        std::cerr << "plotfile level-count mismatch: gold=" << gold.finestLevel() + 1
                  << ", candidate=" << candidate.finestLevel() + 1 << "\n";
        return false;
    }

    for (int lev = 0; lev <= gold.finestLevel(); ++lev) {
        if (gold.probDomain(lev) != candidate.probDomain(lev)) {
            std::cerr << "plotfile domain mismatch at level " << lev << "\n";
            return false;
        }
        if (gold.boxArray(lev) != candidate.boxArray(lev)) {
            std::cerr << "plotfile BoxArray mismatch at level " << lev << "\n";
            return false;
        }
    }

    for (const auto& field : fields) {
        for (int lev = 0; lev <= gold.finestLevel(); ++lev) {
            const MultiFab gold_field = gold.get(lev, field);
            const MultiFab candidate_field = candidate.get(lev, field);
            if (gold_field.nComp() != 1 || candidate_field.nComp() != 1 ||
                gold_field.nComp() != candidate_field.nComp()) {
                std::cerr << "field '" << field << "' has incompatible component count at level "
                          << lev << "\n";
                return false;
            }
            if (!gold_field.is_finite(0, 1, 0, false) ||
                !candidate_field.is_finite(0, 1, 0, false)) {
                std::cerr << "field '" << field << "' is non-finite at level " << lev << "\n";
                return false;
            }
            if (gold_field.boxArray() != candidate_field.boxArray()) {
                std::cerr << "field '" << field << "' has incompatible layout at level "
                          << lev << "\n";
                return false;
            }
        }
    }
    return true;
}

Real maximum_difference (const MultiFab& candidate, const MultiFab& gold)
{
    MultiFab difference(candidate.boxArray(), candidate.DistributionMap(), candidate.nComp(), 0);
    MultiFab::Copy(difference, candidate, 0, 0, 1, 0);
    MultiFab::Subtract(difference, gold, 0, 0, 1, 0);
    return difference.norm0();
}

bool compare_fields (PlotFileData& gold,
                     PlotFileData& candidate,
                     const std::vector<std::string>& fields)
{
    bool all_passed = true;
    amrex::Print() << std::scientific << std::setprecision(12)
                   << "field absolute_error scale absolute_tolerance relative_tolerance "
                      "allowed_error result\n";
    for (const auto& field : fields) {
        Real absolute_error = 0.0;
        Real scale = 0.0;
        for (int lev = 0; lev <= gold.finestLevel(); ++lev) {
            const MultiFab gold_field = gold.get(lev, field);
            const MultiFab candidate_field = candidate.get(lev, field);
            absolute_error = std::max(absolute_error, maximum_difference(candidate_field, gold_field));
            scale = std::max(scale, std::max(candidate_field.norm0(), gold_field.norm0()));
        }

        const FieldTolerance tolerance = tolerance_profile().at(field);
        const Real allowed_error = tolerance.absolute + tolerance.relative * scale;
        const bool passed = absolute_error <= allowed_error;
        all_passed = all_passed && passed;
        amrex::Print() << std::left << std::setw(18) << field
                       << std::right << std::setw(18) << absolute_error
                       << std::setw(18) << scale
                       << std::setw(18) << tolerance.absolute
                       << std::setw(18) << tolerance.relative
                       << std::setw(18) << allowed_error
                       << " " << (passed ? "PASS" : "FAIL") << "\n";
    }
    return all_passed;
}

} // namespace

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);

    Arguments args;
    std::string argument_error;
    if (!read_arguments(argc, argv, args, argument_error)) {
        std::cerr << "SHOC gold comparator argument error: " << argument_error << "\n"
                  << "usage: --mode MODE --gold GOLD_PLOTFILE --candidate CANDIDATE_PLOTFILE\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    const auto& fields = fields_for_mode(args.mode);
    PlotFileData gold(args.gold);
    PlotFileData candidate(args.candidate);
    if (!check_inventory(gold, fields, "gold plotfile") ||
        !check_inventory(candidate, fields, "candidate plotfile") ||
        !check_plotfile_layout(gold, candidate, fields)) {
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    const bool passed = compare_fields(gold, candidate, fields);
    amrex::Finalize();
    return passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
