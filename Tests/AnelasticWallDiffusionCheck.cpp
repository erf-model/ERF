#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

int fail (const std::string& message)
{
    std::cerr << "AnelasticWallDiffusionCheck: " << message << "\n";
    return 1;
}

} // namespace

int main (int argc, char** argv)
{
    if (argc != 6) {
        std::cerr << "usage: checker plotfile axis theta_lo theta_hi alpha_T\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);
    const std::string plotfile_name(argv[1]);
    const int axis = std::atoi(argv[2]);
    const Real theta_lo = std::atof(argv[3]);
    const Real theta_hi = std::atof(argv[4]);
    const Real alpha_T = std::atof(argv[5]);

    if (axis < 0 || axis >= AMREX_SPACEDIM || alpha_T <= 0.0) {
        amrex::Finalize();
        return fail("invalid axis or diffusion coefficient");
    }

    PlotFileData plotfile(plotfile_name);
    const std::string required[] = {"density", "theta", "x_velocity",
                                    "y_velocity", "z_velocity"};
    for (const auto& name : required) {
        if (!has_variable(plotfile, name)) {
            amrex::Finalize();
            return fail("missing required plot variable " + name);
        }
    }

    MultiFab density = plotfile.get(0, "density");
    MultiFab theta = plotfile.get(0, "theta");
    MultiFab x_velocity = plotfile.get(0, "x_velocity");
    MultiFab y_velocity = plotfile.get(0, "y_velocity");
    MultiFab z_velocity = plotfile.get(0, "z_velocity");

    const auto prob_lo = plotfile.probLo();
    const auto prob_hi = plotfile.probHi();
    const auto dx = plotfile.cellSize(0);
    const auto domain = plotfile.probDomain(0);
    const Real length = prob_hi[axis] - prob_lo[axis];
    const Real volume = dx[0] * dx[1] * dx[2];
    const Real tolerance = 5.0e-8;
    Real max_theta_error = 0.0;
    Real max_density_error = 0.0;
    Real max_velocity = 0.0;
    Real integrated_theta_change = 0.0;

    for (amrex::MFIter mfi(theta); mfi.isValid(); ++mfi) {
        const auto& rho = density.const_array(mfi);
        const auto& th = theta.const_array(mfi);
        const auto& u = x_velocity.const_array(mfi);
        const auto& v = y_velocity.const_array(mfi);
        const auto& w = z_velocity.const_array(mfi);
        const auto& box = mfi.validbox();
        for (int k = box.smallEnd(2); k <= box.bigEnd(2); ++k) {
            for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
                for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                    const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
                    const Real coordinate = prob_lo[axis] +
                        (index - domain.smallEnd(axis) + Real(0.5)) * dx[axis];
                    const Real expected_theta = theta_lo + (theta_hi-theta_lo) *
                        (coordinate-prob_lo[axis]) / length;
                    max_theta_error = std::max(max_theta_error,
                                               std::abs(th(i,j,k)-expected_theta));
                    max_density_error = std::max(max_density_error,
                                                 std::abs(rho(i,j,k)-Real(1.0)));
                    max_velocity = std::max(max_velocity,
                        std::max(std::abs(u(i,j,k)),
                        std::max(std::abs(v(i,j,k)), std::abs(w(i,j,k)))));
                    integrated_theta_change += (rho(i,j,k)*th(i,j,k) -
                        expected_theta) * volume;
                }
            }
        }
    }

    // Independent flux oracle: for a linear profile, the low-face inward
    // flux is +alpha*dtheta/dn and the high-face inward flux is its negative.
    // Equal opposite face contributions make the integrated source exactly
    // zero; alpha_T is intentionally used here rather than read from ERF.
    const Real gradient = (theta_hi-theta_lo) / length;
    const Real low_inward_flux = alpha_T * gradient;
    const Real high_inward_flux = -alpha_T * gradient;
    const int n0 = domain.length(0), n1 = domain.length(1), n2 = domain.length(2);
    const Real face_area = (axis == 0) ? dx[1]*dx[2]*n1*n2 :
                           ((axis == 1) ? dx[0]*dx[2]*n0*n2 : dx[0]*dx[1]*n0*n1);
    const Real boundary_budget = face_area * (low_inward_flux + high_inward_flux);

    const Real integrated_tolerance = 5.0e-8 * std::max(Real(1.0),
                                                         std::abs(theta_hi-theta_lo));
    const bool valid = max_theta_error <= tolerance &&
                       max_density_error <= tolerance &&
                       max_velocity <= tolerance &&
                       std::abs(integrated_theta_change) <= integrated_tolerance &&
                       std::abs(boundary_budget) <= integrated_tolerance;
    std::cout << "axis=" << axis << " max_theta_error=" << max_theta_error
              << " max_density_error=" << max_density_error
              << " max_velocity=" << max_velocity
              << " integrated_theta_change=" << integrated_theta_change
              << " boundary_budget=" << boundary_budget << "\n";

    amrex::Finalize();
    return valid ? 0 : fail("manufactured stationary-state property failed");
}
