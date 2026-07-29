#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

using amrex::GpuArray;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::PlotFileData;
using amrex::Real;
using amrex::TilingIfNotGPU;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

Real scaled_tolerance (Real scale, Real factor)
{
    return factor * std::numeric_limits<Real>::epsilon() *
           std::max(Real(1.0), std::abs(scale));
}

int fail (const std::string& message)
{
    std::cerr << "AnelasticWallDiffusionCheck: " << message << "\n";
    return 1;
}

} // namespace

int main (int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "usage: checker plotfile axis theta_lo theta_hi\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);
    const std::string plotfile_name(argv[1]);
    const int axis = std::atoi(argv[2]);
    const Real theta_lo = std::atof(argv[3]);
    const Real theta_hi = std::atof(argv[4]);

    if (axis < 0 || axis >= AMREX_SPACEDIM) {
        amrex::Finalize();
        return fail("invalid axis");
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
    if (!density.is_finite() || !theta.is_finite() ||
        !x_velocity.is_finite() || !y_velocity.is_finite() ||
        !z_velocity.is_finite()) {
        amrex::Finalize();
        return fail("plotfile contains a non-finite required field");
    }

    const auto prob_lo_vector = plotfile.probLo();
    const auto dx_vector = plotfile.cellSize(0);
    const auto domain = plotfile.probDomain(0);
    const GpuArray<Real, AMREX_SPACEDIM> prob_lo = {
        prob_lo_vector[0], prob_lo_vector[1], prob_lo_vector[2]};
    const GpuArray<Real, AMREX_SPACEDIM> dx = {
        dx_vector[0], dx_vector[1], dx_vector[2]};
    const Real length = dx[axis] * Real(domain.length(axis));
    const Real volume = dx[0] * dx[1] * dx[2];

    MultiFab theta_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab density_error(density.boxArray(), density.DistributionMap(), 1, 0);
    MultiFab velocity_error(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab abs_state_change(theta.boxArray(), theta.DistributionMap(), 1, 0);
    MultiFab state_change(theta.boxArray(), theta.DistributionMap(), 1, 0);

    for (MFIter mfi(theta_error, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto bx = mfi.tilebox();
        const auto rho = density.const_array(mfi);
        const auto th = theta.const_array(mfi);
        const auto u = x_velocity.const_array(mfi);
        const auto v = y_velocity.const_array(mfi);
        const auto w = z_velocity.const_array(mfi);
        const auto theta_err = theta_error.array(mfi);
        const auto density_err = density_error.array(mfi);
        const auto velocity_err = velocity_error.array(mfi);
        const auto abs_change = abs_state_change.array(mfi);
        const auto change = state_change.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const int index = (axis == 0) ? i : ((axis == 1) ? j : k);
            const Real coordinate = prob_lo[axis] +
                (Real(index - domain.smallEnd(axis)) + Real(0.5)) * dx[axis];
            const Real expected_theta = theta_lo + (theta_hi-theta_lo) *
                (coordinate-prob_lo[axis]) / length;
            theta_err(i,j,k) = std::abs(th(i,j,k) - expected_theta);
            density_err(i,j,k) = std::abs(rho(i,j,k) - Real(1.0));
            velocity_err(i,j,k) = std::max(std::abs(u(i,j,k)),
                std::max(std::abs(v(i,j,k)), std::abs(w(i,j,k))));
            change(i,j,k) = rho(i,j,k) * th(i,j,k) - expected_theta;
            abs_change(i,j,k) = std::abs(change(i,j,k));
        });
    }
    amrex::Gpu::streamSynchronize();

    const Real max_theta_error = theta_error.norm0(0, 0, false);
    const Real max_density_error = density_error.norm0(0, 0, false);
    const Real max_velocity_error = velocity_error.norm0(0, 0, false);
    const auto cell_count = domain.numPts();
    const Real mean_abs_state_change =
        abs_state_change.sum(0) / Real(cell_count);
    const Real integrated_state_change = state_change.sum(0) * volume;
    const Real state_scale = std::max(std::abs(theta_lo), std::abs(theta_hi));
    const Real theta_tolerance = scaled_tolerance(state_scale, Real(1.0));
    const Real density_tolerance = scaled_tolerance(Real(1.0), Real(128.0));
    const Real velocity_tolerance = scaled_tolerance(Real(1.0), Real(256.0));
    const Real mean_abs_tolerance = scaled_tolerance(state_scale, Real(1.0));
    const Real integrated_tolerance =
        Real(8.0) * mean_abs_tolerance * volume * Real(cell_count);
    const bool valid = max_theta_error <= theta_tolerance &&
                       max_density_error <= density_tolerance &&
                       max_velocity_error <= velocity_tolerance &&
                       mean_abs_state_change <= mean_abs_tolerance &&
                       std::abs(integrated_state_change) <= integrated_tolerance;

    std::cout << "axis=" << axis << " max_theta_error=" << max_theta_error
              << " max_density_error=" << max_density_error
              << " max_velocity_error=" << max_velocity_error
              << " mean_abs_state_change=" << mean_abs_state_change
              << " integrated_state_change=" << integrated_state_change
              << " theta_tolerance=" << theta_tolerance
              << " mean_abs_tolerance=" << mean_abs_tolerance
              << " integrated_tolerance=" << integrated_tolerance << "\n";

    // This full-run check is deliberately a stationary manufactured-state
    // property test. The direct unit test owns the production face-flux/RHS
    // oracle and its signed integrated boundary balance.
    amrex::Finalize();
    return valid ? 0 : fail("manufactured stationary-state property failed");
}
