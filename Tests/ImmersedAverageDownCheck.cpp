#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

// Checks the average down of a two-level immersed-forcing plotfile:
//
//   checker <plotfile> <tol>
//
// Every coarse cell covered by level 1 that is not solid on its own level must hold the
// fluid-mass-weighted average of the fine cells over it,
//     theta_c = sum (1 - beta_f) rho_f theta_f / sum (1 - beta_f) rho_f,
// with beta the solid fraction (terrain_IB_mask). The check also requires that at least one of
// those coarse cells lies over a solid fine cell, since a plain average and the weighted one agree
// everywhere else. The weight has no cell volume: the check is for constant-dz decks without map
// factors, where the volume is uniform.
// Exit 0: passed; 1: the coarse cells are not the weighted average; 3: the plotfile cannot be
// checked (missing level or variable, non-finite theta, no coarse cell over a solid fine cell).

namespace {

using namespace amrex;

// exit 1: the average is not the fluid-mass-weighted one; exit 3: the plotfile cannot be checked
int fail (const std::string& message, int code = 3)
{
    std::cerr << "ImmersedAverageDownCheck: " << message << "\n";
    return code;
}

bool has_variable (const PlotFileData& pf, const std::string& name)
{
    const auto& names = pf.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

int check_fluid_average (PlotFileData& pf, Real tol)
{
    if (pf.finestLevel() < 1) { return fail("needs two levels"); }
    for (const char* name : {"theta", "density", "terrain_IB_mask"}) {
        if (!has_variable(pf, name)) { return fail(std::string("plotfile has no ") + name); }
    }
    const MultiFab th_c = pf.get(0, "theta");
    const MultiFab b_c  = pf.get(0, "terrain_IB_mask");
    const MultiFab th_f = pf.get(1, "theta");
    const MultiFab r_f  = pf.get(1, "density");
    const MultiFab b_f  = pf.get(1, "terrain_IB_mask");
    if (!th_c.is_finite() || !th_f.is_finite()) { return fail("non-finite theta"); }

    // the coarse values on the coarsened fine layout
    const IntVect rr(pf.refRatio(0));
    const BoxArray cba = amrex::coarsen(th_f.boxArray(), rr);
    MultiFab tc(cba, th_f.DistributionMap(), 1, 0), bc(cba, th_f.DistributionMap(), 1, 0);
    tc.ParallelCopy(th_c, 0, 0, 1);
    bc.ParallelCopy(b_c, 0, 0, 1);

    ReduceOps<ReduceOpMax, ReduceOpSum> op;
    ReduceData<Real, int> data(op);
    using Tuple = typename decltype(data)::Type;
    for (MFIter mfi(tc); mfi.isValid(); ++mfi) {
        auto const t  = tc.const_array(mfi);
        auto const bb = bc.const_array(mfi);
        auto const tf = th_f.const_array(mfi);
        auto const rf = r_f.const_array(mfi);
        auto const bf = b_f.const_array(mfi);
        const int rx = rr[0], ry = rr[1], rz = rr[2];
        op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Tuple
        {
            if (bb(i,j,k) >= Real(1.)) { return {Real(0.), 0}; }
            Real num = Real(0.), den = Real(0.);
            int nsolid = 0;
            for (int kk = k*rz; kk < (k+1)*rz; ++kk) {
            for (int jj = j*ry; jj < (j+1)*ry; ++jj) {
            for (int ii = i*rx; ii < (i+1)*rx; ++ii) {
                const Real w = amrex::max(Real(0.), Real(1.) - bf(ii,jj,kk)) * rf(ii,jj,kk);
                num += w * tf(ii,jj,kk);
                den += w;
                if (bf(ii,jj,kk) >= Real(1.)) { ++nsolid; }
            }}}
            if (den <= Real(0.)) { return {Real(0.), 0}; }
            return {std::abs(t(i,j,k) - num / den), (nsolid > 0) ? 1 : 0};
        });
    }
    auto v = data.value(op);
    Real err = amrex::get<0>(v);
    int nmixed = amrex::get<1>(v);
    ParallelDescriptor::ReduceRealMax(err);
    ParallelDescriptor::ReduceIntSum(nmixed);
    Print() << "max |coarse theta - fluid-mass-weighted fine average| = " << err
            << " over the covered coarse cells (" << nmixed
            << " of them over a solid fine cell; tol " << tol << ")\n";
    if (nmixed == 0) { return fail("no covered coarse cell lies over a solid fine cell"); }
    if (err > tol) { return fail("the coarse cells are not the fluid-mass-weighted fine average", 1); }
    return 0;
}

} // namespace

int main (int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: checker plotfile tol\n";
        return 2;
    }
    amrex::Initialize(argc, argv, false);
    int result = 2;
    {
        PlotFileData pf(argv[1]);
        const Real tol = static_cast<Real>(std::atof(argv[2]));
        result = check_fluid_average(pf, tol);
    }
    amrex::Finalize();
    return result;
}
