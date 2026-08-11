//
// Property oracle for horizontal advection of a z-independent passive scalar on
// a vertically stretched mesh.
//
// The companion input deck advects a Gaussian blob -- initialized with no
// z-dependence -- with a uniform horizontal velocity and no gravity, diffusion
// or terrain slope.  Every k-level therefore sees an identical 1-D advection
// problem and the exact solution stays independent of z forever.
//
// That property is exactly what pins the metric/map-factor weighting of the
// time-averaged mass flux (avg_xmom/avg_ymom) that acoustic substepping hands
// to scalar advection: the consumer re-applies mfsq/detJ to the flux
// divergence, so the flux itself must carry h_zeta = dz_k/dz and 1/mf.  Drop
// h_zeta and level k advects at u/h_zeta(k) instead of u, which shears the blob
// apart in the vertical -- an O(1) error on a stretched mesh, versus roundoff
// when the weighting is present.
//
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using amrex::Box;
using amrex::BoxArray;
using amrex::DistributionMapping;
using amrex::MFIter;
using amrex::MultiFab;
using amrex::ParallelFor;
using amrex::PlotFileData;
using amrex::Real;

bool has_variable (const PlotFileData& plotfile, const std::string& name)
{
    const auto& names = plotfile.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

int fail (const std::string& message)
{
    std::cerr << "ColumnAdvectionCheck: " << message << "\n";
    return 1;
}

// Gather a plotfile variable onto one domain-covering box so the column
// reduction below sees every k of a column at once regardless of how the run
// that wrote the plotfile was decomposed.
MultiFab gather_to_domain (const PlotFileData& plotfile,
                           const Box& domain,
                           const BoxArray& ba,
                           const DistributionMapping& dm,
                           const std::string& name)
{
    PlotFileData& pf = const_cast<PlotFileData&>(plotfile);
    MultiFab src = pf.get(0, name);
    MultiFab dst(ba, dm, 1, 0);
    dst.setVal(Real(0.0));
    dst.ParallelCopy(src, 0, 0, 1);
    amrex::ignore_unused(domain);
    return dst;
}

// Largest max-min spread down any single (i,j) column.
Real column_spread (const MultiFab& mf, const Box& domain)
{
    const Box slab = amrex::makeSlab(domain, 2, domain.smallEnd(2));
    MultiFab spread(BoxArray(slab), mf.DistributionMap(), 1, 0);
    spread.setVal(Real(0.0));

    const int klo = domain.smallEnd(2);
    const int khi = domain.bigEnd(2);

    for (MFIter mfi(spread); mfi.isValid(); ++mfi) {
        const Box bx = mfi.validbox();
        const auto s = mf.const_array(mfi);
        const auto out = spread.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            amrex::ignore_unused(k);
            Real vmin = s(i,j,klo);
            Real vmax = s(i,j,klo);
            for (int kk = klo+1; kk <= khi; ++kk) {
                vmin = amrex::min(vmin, s(i,j,kk));
                vmax = amrex::max(vmax, s(i,j,kk));
            }
            out(i,j,klo) = vmax - vmin;
        });
    }
    amrex::Gpu::streamSynchronize();

    return spread.norm0(0, 0, false);
}

} // namespace

int main (int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: checker plotfile rel_column_tolerance min_amplitude\n";
        return 2;
    }

    amrex::Initialize(argc, argv, false);

    const std::string plotfile_name(argv[1]);
    const Real rel_tolerance = static_cast<Real>(std::atof(argv[2]));
    const Real min_amplitude = static_cast<Real>(std::atof(argv[3]));

    PlotFileData plotfile(plotfile_name);

    const std::string required[] = {"scalar", "density", "x_velocity",
                                    "y_velocity", "z_velocity"};
    for (const auto& name : required) {
        if (!has_variable(plotfile, name)) {
            amrex::Finalize();
            return fail("missing required plot variable " + name);
        }
    }

    const Box domain = plotfile.probDomain(0);
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);

    MultiFab scalar     = gather_to_domain(plotfile, domain, ba, dm, "scalar");
    MultiFab density    = gather_to_domain(plotfile, domain, ba, dm, "density");
    MultiFab x_velocity = gather_to_domain(plotfile, domain, ba, dm, "x_velocity");
    MultiFab y_velocity = gather_to_domain(plotfile, domain, ba, dm, "y_velocity");
    MultiFab z_velocity = gather_to_domain(plotfile, domain, ba, dm, "z_velocity");

    if (!scalar.is_finite() || !density.is_finite() || !x_velocity.is_finite() ||
        !y_velocity.is_finite() || !z_velocity.is_finite()) {
        amrex::Finalize();
        return fail("plotfile contains a non-finite required field");
    }

    // The blob has to still be a blob -- a field that decayed (or was never
    // initialized) is trivially z-independent and must not pass.
    const Real amplitude = scalar.max(0) - scalar.min(0);

    // Uniform initial flow with no forcing stays uniform, so these spreads
    // bound the whole solution, not just the passive scalar.  They are the
    // guard that a shear-free background actually held for the whole run.
    const Real scalar_spread  = column_spread(scalar,  domain);
    const Real density_spread = density.max(0) - density.min(0);
    const Real u_spread       = x_velocity.max(0) - x_velocity.min(0);
    const Real v_spread       = y_velocity.max(0) - y_velocity.min(0);
    const Real w_max          = std::max(std::abs(z_velocity.max(0)),
                                         std::abs(z_velocity.min(0)));

    const Real scalar_tolerance = rel_tolerance * std::max(amplitude, Real(1.0));
    const Real flow_tolerance   = rel_tolerance * Real(10.0);

    const bool valid = amplitude          >= min_amplitude     &&
                       scalar_spread      <= scalar_tolerance  &&
                       density_spread     <= flow_tolerance    &&
                       u_spread           <= flow_tolerance    &&
                       v_spread           <= flow_tolerance    &&
                       w_max              <= flow_tolerance;

    std::cout << "amplitude=" << amplitude
              << " scalar_column_spread=" << scalar_spread
              << " scalar_tolerance=" << scalar_tolerance
              << " density_spread=" << density_spread
              << " u_spread=" << u_spread
              << " v_spread=" << v_spread
              << " max_abs_w=" << w_max
              << " flow_tolerance=" << flow_tolerance << "\n";

    amrex::Finalize();

    if (!valid) {
        return fail("z-independence of the advected scalar was not preserved");
    }
    return 0;
}
