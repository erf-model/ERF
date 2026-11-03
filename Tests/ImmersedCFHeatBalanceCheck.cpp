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

// Checks the heat balance of a two-level immersed-forcing plotfile:
//
//   checker <plotfile> <theta0> <tol>
//
// The fluid-mass-weighted mean theta of the domain, taken from the finest level that holds each
// cell (level 1 where it covers level 0, level 0 elsewhere), must be within tol of theta0: with a
// uniform theta0 and no heat entering, any drift is a spurious source. The weight is
// (1 - beta) rho with beta the solid fraction (terrain_IB_mask); no cell volume enters, so the
// check is for constant-dz decks without map factors, and the fine cells are weighted by their
// share of a coarse volume. Taking the fine level where it exists keeps the check independent of
// how the covered coarse cells are averaged down (erf.if_average_down).
// Exit 0: passed; 1: the mean drifted by more than tol; 3: the plotfile cannot be checked
// (missing level or variable, non-finite theta, no fluid mass).

namespace {

using namespace amrex;

int fail (const std::string& message, int code = 3)
{
    std::cerr << "ImmersedCFHeatBalanceCheck: " << message << "\n";
    return code;
}

bool has_variable (const PlotFileData& pf, const std::string& name)
{
    const auto& names = pf.varNames();
    return std::find(names.begin(), names.end(), name) != names.end();
}

// fluid-mass-weighted sums of (theta - theta0) and of the weight over the cells of a level where
// the mask (null: every cell) is zero; weight w_cell per cell
void weighted_sums (const MultiFab& th, const MultiFab& rho, const MultiFab& b, const iMultiFab* covered,
                    Real theta0, Real w_cell, Real& num, Real& den)
{
    ReduceOps<ReduceOpSum, ReduceOpSum> op;
    ReduceData<Real, Real> data(op);
    using Tuple = typename decltype(data)::Type;
    for (MFIter mfi(th); mfi.isValid(); ++mfi) {
        auto const t  = th.const_array(mfi);
        auto const r  = rho.const_array(mfi);
        auto const bb = b.const_array(mfi);
        auto const cv = covered ? covered->const_array(mfi) : Array4<int const>{};
        const bool use_cv = (covered != nullptr);
        op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> Tuple
        {
            if (use_cv && cv(i,j,k) != 0) { return {Real(0.), Real(0.)}; }
            const Real fluid = amrex::max(Real(0.), Real(1.) - amrex::min(bb(i,j,k), Real(1.)));
            const Real w = w_cell * fluid * r(i,j,k);
            return {w * (t(i,j,k) - theta0), w};
        });
    }
    auto sums = data.value(op);
    num += amrex::get<0>(sums);
    den += amrex::get<1>(sums);
}

int check_mean (PlotFileData& pf, Real theta0, Real tol)
{
    if (pf.finestLevel() < 1) { return fail("needs two levels"); }
    for (const char* name : {"theta", "density", "terrain_IB_mask"}) {
        if (!has_variable(pf, name)) { return fail(std::string("plotfile has no ") + name); }
    }
    const MultiFab th_c = pf.get(0, "theta");
    const MultiFab r_c  = pf.get(0, "density");
    const MultiFab b_c  = pf.get(0, "terrain_IB_mask");
    const MultiFab th_f = pf.get(1, "theta");
    const MultiFab r_f  = pf.get(1, "density");
    const MultiFab b_f  = pf.get(1, "terrain_IB_mask");
    if (!th_c.is_finite() || !th_f.is_finite()) { return fail("non-finite theta"); }

    const IntVect rr(pf.refRatio(0));
    const Geometry cgeom(pf.probDomain(0), RealBox(pf.probLo(), pf.probHi()), pf.coordSys(), {AMREX_D_DECL(0,0,0)});
    const iMultiFab covered = amrex::makeFineMask(th_c.boxArray(), th_c.DistributionMap(), IntVect(0),
                                                  th_f.boxArray(), rr, cgeom.periodicity(), 0, 1);
    const Real w_fine = Real(1.) / static_cast<Real>(rr[0] * rr[1] * rr[2]);

    Real num = Real(0.), den = Real(0.);
    weighted_sums(th_c, r_c, b_c, &covered, theta0, Real(1.), num, den);
    weighted_sums(th_f, r_f, b_f, nullptr,  theta0, w_fine,   num, den);
    ParallelDescriptor::ReduceRealSum(num);
    ParallelDescriptor::ReduceRealSum(den);
    if (!(den > Real(0.))) { return fail("no fluid mass"); }

    const Real drift = num / den;
    Print() << "fluid-mass-weighted mean theta - theta0 = " << drift << " K (tol " << tol << " K)\n";
    if (std::abs(drift) > tol) { return fail("the mean theta drifted by more than tol", 1); }
    return 0;
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int rc = 0;
    {
        if (argc < 4) {
            rc = fail("usage: checker <plotfile> <theta0> <tol>");
        } else {
            PlotFileData pf(argv[1]);
            rc = check_mean(pf, static_cast<Real>(std::atof(argv[2])), static_cast<Real>(std::atof(argv[3])));
        }
    }
    amrex::Finalize();
    return rc;
}
