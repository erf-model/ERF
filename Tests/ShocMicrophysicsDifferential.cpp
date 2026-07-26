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

// Keep field-presence failures separate from numerical differences so a
// changed plotfile inventory cannot be mistaken for equivalent physics.
bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    return std::find(plotfile.varNames().begin(), plotfile.varNames().end(), name) !=
           plotfile.varNames().end();
}

// Return the largest absolute cellwise difference for one field.  The
// comparison deliberately uses common state and SHOC diagnostics only; the
// two microphysics models may legitimately expose different precipitating
// species outside this shared contract.
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
    // SatAdj with condensation disabled and MoistNoCondensation should take
    // the same native SHOC state-update path.  This test therefore requires
    // every common final-state field to agree to roundoff, while allowing
    // model-specific precipitation fields to differ.
    amrex::Initialize(argc, argv, false);

    if (argc != 3) {
        std::cerr << "usage: erf_shoc_microphysics_differential SATADJ_PLOTFILE NOCOND_PLOTFILE\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    PlotFileData satadj(argv[1]);
    PlotFileData nocond(argv[2]);

    // These fields define the common-state contract.  Keep this list explicit
    // so adding a new diagnostic requires an intentional differential-test
    // decision rather than silently changing the comparison surface.
    const std::vector<std::string> fields {
        "density", "temp", "theta", "pressure", "qv", "qc", "rhoKE",
        "x_velocity", "y_velocity", "Kmv", "Khv", "Lturb", "shoc_cldfrac",
        "shoc_ql", "shoc_cond", "brunt", "shear_prod", "buoy_prod", "diss_tke"
    };
    Real maximum = 0.0;
#ifdef AMREX_USE_FLOAT
    constexpr Real tolerance = 1.0e-6;
#else
    constexpr Real tolerance = 1.0e-12;
#endif
    for (const auto& field : fields) {
        if (!has_variable(satadj, field) || !has_variable(nocond, field)) {
            std::cerr << "required common field '" << field << "' is missing\n";
            amrex::Finalize();
            return EXIT_FAILURE;
        }
        const Real difference = maximum_difference(satadj, nocond, field);
        maximum = std::max(maximum, difference);
        amrex::Print() << "SHOC microphysics differential: " << field
                       << " max |delta|=" << difference << "\n";
    }

    amrex::Print() << "SHOC microphysics differential: max common-field difference="
                   << maximum << " tolerance=" << tolerance << "\n";
    if (maximum > tolerance) {
        std::cerr << "SatAdj and MoistNoCondensation common state fields diverged\n";
        amrex::Finalize();
        return EXIT_FAILURE;
    }

    amrex::Finalize();
    return EXIT_SUCCESS;
}
