#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Compares a run over an immersed plane that lies on a cell face with the same run on flat
// ground (the surface layer at the bottom of the domain):
//
//   checker <plotfile_flat> <plotfile_plane> <k_offset> <tol> [k_count]
//
// For x_velocity, y_velocity, theta, and (when present) KE and qv, the horizontally averaged profile of
// the plane run, from its cell k_offset up, must match the flat run's profile from its cell 0 up
// (the lowest k_count cells, all by default), within tol relative to the flat run's largest
// |value| of that variable over those cells (theta against its departure from its own bottom value). With the fraction-stress wall law a face-aligned surface
// is the surface-layer law of flat ground, so the two agree; the legacy law leaves the fluid cell
// on the solid unforced and must fail.
// Exit 0: passed; 1: a profile differs by more than tol; 3: the plotfiles cannot be compared
// (sizes, missing variable, non-finite values).

namespace {

using namespace amrex;

int fail (const std::string& message, int code = 3)
{
    std::cerr << "ImmersedFlatEquivalenceCheck: " << message << "\n";
    return code;
}

bool has_variable (const PlotFileData& pf, const std::string& name)
{
    const auto& names = pf.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

// horizontally averaged profile of a level-0 variable
std::vector<Real> profile (PlotFileData& pf, const std::string& name)
{
    const MultiFab m = pf.get(0, name);
    const Box dom = pf.probDomain(0);
    const int nk = dom.length(2), k0 = dom.smallEnd(2);
    std::vector<Real> sum(nk, Real(0.)), cnt(nk, Real(0.));
    for (MFIter mfi(m); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const a = m.const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) { sum[k-k0] += a(i,j,k); cnt[k-k0] += Real(1.); });
    }
    ParallelDescriptor::ReduceRealSum(sum.data(), nk);
    ParallelDescriptor::ReduceRealSum(cnt.data(), nk);
    for (int k = 0; k < nk; ++k) { sum[k] /= cnt[k]; }
    return sum;
}

int compare (PlotFileData& pf_flat, PlotFileData& pf_plane, int k_off, Real tol, int k_count)
{
    const int nk_flat  = pf_flat.probDomain(0).length(2);
    const int nk_plane = pf_plane.probDomain(0).length(2);
    if (nk_plane != nk_flat + k_off) { return fail("the plane run must have k_offset more cells in z than the flat run"); }
    const int nk_cmp = (k_count > 0) ? std::min(k_count, nk_flat) : nk_flat;
    int rc = 0;
    for (const char* name : {"x_velocity", "y_velocity", "theta", "KE", "qv"}) {
        const bool in_flat = has_variable(pf_flat, name), in_plane = has_variable(pf_plane, name);
        if (!in_flat || !in_plane) {
            // KE and qv are compared when both runs write them (the k-eqn RANS, a moisture model)
            const bool optional = (std::string(name) == "KE" || std::string(name) == "qv");
            if (optional && !in_flat && !in_plane) { continue; }
            return fail(std::string("plotfile has no ") + name);
        }
        const std::vector<Real> a = profile(pf_flat, name);
        const std::vector<Real> b = profile(pf_plane, name);
        const bool is_theta = (std::string(name) == "theta");
        Real scale = Real(0.), worst = Real(0.);
        int k_worst = -1;
        for (int k = 0; k < nk_cmp; ++k) {
            const Real ak = is_theta ? a[k] - a[0] : a[k];
            scale = std::max(scale, std::abs(ak));
        }
        if (is_theta) { scale = std::max(scale, Real(1.)); }   // a uniform theta compares in kelvin
        for (int k = 0; k < nk_cmp; ++k) {
            const Real ak = is_theta ? a[k] - a[0] : a[k];
            const Real bk = is_theta ? b[k + k_off] - a[0] : b[k + k_off];
            if (!std::isfinite(ak) || !std::isfinite(bk)) { return fail(std::string("non-finite ") + name); }
            const Real diff = std::abs(ak - bk);
            if (diff > worst) { worst = diff; k_worst = k; }
        }
        const Real rel = worst / scale;
        Print() << name << ": largest |flat - plane| over the lowest " << nk_cmp << " cells = " << worst
                << " at flat cell " << k_worst << " (relative to " << scale << ": " << rel << ", tol " << tol << ")\n";
        if (rel > tol) { rc = 1; }
    }
    if (rc) { return fail("the plane run differs from flat ground by more than tol", 1); }
    return 0;
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int rc = 0;
    {
        if (argc < 5) {
            rc = fail("usage: checker <plotfile_flat> <plotfile_plane> <k_offset> <tol> [k_count]");
        } else {
            PlotFileData pf_flat(argv[1]);
            PlotFileData pf_plane(argv[2]);
            const int k_count = (argc > 5) ? std::atoi(argv[5]) : 0;
            rc = compare(pf_flat, pf_plane, std::atoi(argv[3]), static_cast<Real>(std::atof(argv[4])), k_count);
        }
    }
    amrex::Finalize();
    return rc;
}
