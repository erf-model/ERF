#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    return std::find(plotfile.varNames().begin(), plotfile.varNames().end(), name) !=
           plotfile.varNames().end();
}

Real maximum_difference (PlotFileData& lhs,
                         PlotFileData& rhs,
                         const std::string& name)
{
    MultiFab difference(lhs.get(0, name));
    MultiFab::Subtract(difference, rhs.get(0, name), 0, 0, 1, 0);
    return difference.norm0();
}

}

int
main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);

    if (argc != 3) {
        std::cerr << "usage: erf_shoc_microphysics_differential SATADJ_PLOTFILE NOCOND_PLOTFILE\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    PlotFileData satadj(argv[1]);
    PlotFileData nocond(argv[2]);

    const std::vector<std::string> fields {"qv", "shoc_ql"};
    Real max_qv_difference = 0.0;
    Real max_ql_difference = 0.0;
    for (const auto& field : fields) {
        if (!has_variable(satadj, field) || !has_variable(nocond, field)) {
            std::cerr << "required common field '" << field << "' is missing\n";
            amrex::Finalize();
            return EXIT_FAILURE;
        }
        const Real difference = maximum_difference(satadj, nocond, field);
        if (field == "qv") {
            max_qv_difference = difference;
        } else {
            max_ql_difference = difference;
        }
    }

    amrex::Print() << "SHOC microphysics differential: max |delta qv|="
                   << max_qv_difference << " max |delta shoc_ql|="
                   << max_ql_difference << "\n";

    const bool distinct = max_qv_difference > 1.0e-12 || max_ql_difference > 1.0e-12;
    if (!distinct) {
        std::cerr << "SatAdj and MoistNoCondensation final states are indistinguishable\n";
    }

    amrex::Finalize();
    return distinct ? EXIT_SUCCESS : EXIT_FAILURE;
}
