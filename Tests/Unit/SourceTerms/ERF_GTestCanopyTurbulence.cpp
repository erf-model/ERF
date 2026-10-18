// Unit tests for the canopy terms of the k-eqn RANS model
// (ERF_CanopyTurbulence.H): the inputs, the length scale cap, the TKE source
// formula and the kernel that stages it in the RhoKE source.

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "ERF_CanopyTurbulence.H"

using amrex::Real;

namespace {

Real tol (Real factor = Real(100))
{
    return factor * std::numeric_limits<Real>::epsilon();
}

} // namespace

// Defaults match the ForestDrag.canopy_* inputs of kynema-sgf (f/winsent):
// off, beta_p = 1, beta_d = 4, alpha = 0.04
TEST(CanopyTurbulence, Defaults)
{
    const CanopyTurbulence canopy;
    EXPECT_FALSE(canopy.enabled);
    EXPECT_EQ(canopy.beta_p, Real(1.0));
    EXPECT_EQ(canopy.beta_d, Real(4.0));
    EXPECT_EQ(canopy.length_alpha, Real(0.04));
}

TEST(CanopyTurbulence, ReadsInputs)
{
    amrex::ParmParse pp("canopy_gtest_reads");
    pp.add("forest_canopy_tke", true);
    pp.add("forest_canopy_beta_p", Real(0.8));
    pp.add("forest_canopy_beta_d", Real(5.1));
    pp.add("forest_canopy_length_alpha", Real(0.054));
    CanopyTurbulence canopy;
    canopy.init_params(pp);
    EXPECT_TRUE(canopy.enabled);
    EXPECT_EQ(canopy.beta_p, Real(0.8));
    EXPECT_EQ(canopy.beta_d, Real(5.1));
    EXPECT_EQ(canopy.length_alpha, Real(0.054));
}

// L <= alpha / f where f > 0: alpha / f = 2 m caps L = 10 m
TEST(CanopyTurbulence, LengthCappedInsideCanopy)
{
    const Real f = Real(0.02);
    const Real alpha = Real(0.04);
    EXPECT_NEAR(canopy_turbulence::limited_length(Real(10), f, alpha), alpha / f, tol());
}

// A sparse canopy, alpha / f = 20 m > 10 m, leaves L as it is
TEST(CanopyTurbulence, LengthUnchangedInSparseCanopy)
{
    EXPECT_EQ(canopy_turbulence::limited_length(Real(10), Real(0.002), Real(0.04)), Real(10));
}

// No forest (f = 0, f < 0), a tiny f and alpha = 0 all leave L as it is
TEST(CanopyTurbulence, LengthUnchangedWithoutCap)
{
    EXPECT_EQ(canopy_turbulence::limited_length(Real(10), Real(0),  Real(0.04)), Real(10));
    EXPECT_EQ(canopy_turbulence::limited_length(Real(10), Real(-1), Real(0.04)), Real(10));
    EXPECT_EQ(canopy_turbulence::limited_length(Real(10), Real(1.0e-30), Real(0.04)), Real(10));
    EXPECT_EQ(canopy_turbulence::limited_length(Real(10), Real(0.02), Real(0)), Real(10));
}

// The cap is continuous where it starts to act, alpha / f = L
TEST(CanopyTurbulence, LengthCapContinuous)
{
    const Real L = Real(10);
    const Real alpha = Real(0.04);
    const Real f_edge = alpha / L;
    EXPECT_NEAR(canopy_turbulence::limited_length(L, f_edge * (1 + tol()), alpha), L, L * tol(Real(1000)));
    EXPECT_EQ(canopy_turbulence::limited_length(L, f_edge * (1 - tol()), alpha), L);
}

// Outside the forest the source is exactly zero
TEST(CanopyTurbulence, SourceZeroWithoutForest)
{
    EXPECT_EQ(canopy_turbulence::tke_source(Real(1.2), Real(0.5), Real(0), Real(5), Real(1), Real(4), Real(0.1)),
              Real(0));
}

// With beta_d = 0 only the wake production rho beta_p f |U|^3 is left
TEST(CanopyTurbulence, SourceProductionOnly)
{
    const Real rho = Real(1.2), f = Real(0.05), wind = Real(5), beta_p = Real(0.7);
    EXPECT_NEAR(canopy_turbulence::tke_source(rho, Real(0.5), f, wind, beta_p, Real(0), Real(0.1)),
                rho * beta_p * f * wind * wind * wind, tol(Real(1000)));
}

// For beta_d f |U| dt << 1 the sink tends to the explicit -rho beta_d f |U| k
TEST(CanopyTurbulence, SinkSmallTimeStep)
{
    const Real rho = Real(1.2), tke = Real(0.5), f = Real(0.01), wind = Real(2), beta_d = Real(4);
    const Real dt = Real(1.0e-4);
    const Real explicit_sink = -rho * beta_d * f * wind * tke;
    const Real sink = canopy_turbulence::tke_source(rho, tke, f, wind, Real(0), beta_d, dt);
    EXPECT_NEAR(sink, explicit_sink, Real(1.0e-5) * std::abs(explicit_sink));
}

// A dense canopy and a long step: the sink tends to -rho k / dt and never
// removes more than the local k in one step
TEST(CanopyTurbulence, SinkBoundedByLocalTKE)
{
    const Real rho = Real(1.2), tke = Real(0.5), dt = Real(10);
    const Real sink = canopy_turbulence::tke_source(rho, tke, Real(1), Real(20), Real(0), Real(4), dt);
    EXPECT_GE(sink, -rho * tke / dt);
    EXPECT_NEAR(sink, -rho * tke / dt, tol(Real(1000)));
}

// The staged kernel: cell 0 is in the forest, cell 1 is not. The wind is the
// average of the two faces of each cell, and only RhoKE changes.
TEST(CanopyTurbulence, KernelStagesSourceInCanopyCells)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const int ncomp = RhoKE_comp + 1;

    MultiFab source(ba, dm, ncomp, 0);
    MultiFab state(ba, dm, ncomp, 0);
    MultiFab forest(ba, dm, 1, 0);
    MultiFab xvel(convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    MultiFab yvel(convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    MultiFab zvel(convert(ba, IntVect(0, 0, 1)), dm, 1, 0);

    const Real rho = Real(1.2), tke = Real(0.5), f = Real(0.05), dt = Real(0.1);
    source.setVal(Real(0));
    state.setVal(Real(0));
    state.setVal(rho, Rho_comp, 1);
    state.setVal(rho * tke, RhoKE_comp, 1);
    yvel.setVal(Real(1));
    zvel.setVal(Real(0));
    for (MFIter mfi(xvel); mfi.isValid(); ++mfi) {
        const auto& u = xvel.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i, j, k) = Real(2) + Real(2) * Real(i);   // faces 2, 4, 6
        });
    }
    for (MFIter mfi(forest); mfi.isValid(); ++mfi) {
        const auto& fd = forest.array(mfi);
        ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fd(i, j, k) = (i == 0) ? f : Real(0);
        });
    }

    CanopyTurbulence canopy;
    canopy.enabled = true;
    AddCanopyTKESources(source, state, xvel, yvel, zvel, forest, canopy, dt);

    const Real wind0 = std::sqrt(Real(3) * Real(3) + Real(1));   // u = (2 + 4) / 2, v = 1
    const Real expected = canopy_turbulence::tke_source(rho, tke, f, wind0, canopy.beta_p,
                                                        canopy.beta_d, dt);
    // Reductions rather than host reads of the arrays, so the test also runs on a GPU
    // build: the canopy cell holds the (positive) source, the other cell zero, and no
    // other component is touched
    ASSERT_GT(expected, Real(0));
    EXPECT_NEAR(source.max(RhoKE_comp), expected, std::abs(expected) * tol());
    EXPECT_EQ(source.min(RhoKE_comp), Real(0));
    for (int n = 0; n < RhoKE_comp; ++n) {
        EXPECT_EQ(source.norm0(n), Real(0)) << "component " << n;
    }
}
