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

// Kernels live outside the TEST bodies: nvcc rejects extended device lambdas in
// the private TestBody() that TEST generates
void fill_face_velocity (amrex::MultiFab& xvel)
{
    for (amrex::MFIter mfi(xvel); mfi.isValid(); ++mfi) {
        const auto& u = xvel.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            u(i, j, k) = Real(2) + Real(2) * Real(i);   // faces 2, 4, 6
        });
    }
}

void fill_forest_first_cell (amrex::MultiFab& forest, Real f)
{
    for (amrex::MFIter mfi(forest); mfi.isValid(); ++mfi) {
        const auto& fd = forest.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            fd(i, j, k) = (i == 0) ? f : Real(0);
        });
    }
}

// Stage the canopy source box by box, as make_sources does; blank may be null
void add_canopy_sources (amrex::MultiFab& source, const amrex::MultiFab& state,
                         const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
                         const amrex::MultiFab& zvel, const amrex::MultiFab& forest,
                         const amrex::MultiFab* blank, const CanopyTurbulence& canopy, Real dt)
{
    for (amrex::MFIter mfi(source); mfi.isValid(); ++mfi) {
        AddCanopyTKESources(mfi.validbox(), source.array(mfi), state.const_array(mfi),
                            state.const_array(mfi), xvel.const_array(mfi), yvel.const_array(mfi),
                            zvel.const_array(mfi), forest.const_array(mfi),
                            blank ? blank->const_array(mfi) : amrex::Array4<const Real>{},
                            canopy, dt);
    }
}

} // namespace

// Defaults match the ForestDrag.canopy_* inputs of kynema-sgf (f/winsent):
// off, beta_p = 1, beta_d = 4, alpha = 0.04
TEST(CanopyTurbulence, Defaults)
{
    const CanopyTurbulence canopy{};
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
    EXPECT_EQ(canopy_turbulence::tke_source(Real(1.2), Real(0.6), Real(0), Real(5), Real(1), Real(4), Real(0.1)),
              Real(0));
}

// With beta_d = 0 only the wake production rho beta_p f |U|^3 is left (phi1 = 1)
TEST(CanopyTurbulence, SourceProductionOnly)
{
    const Real rho = Real(1.2), f = Real(0.05), wind = Real(5), beta_p = Real(0.7);
    EXPECT_NEAR(canopy_turbulence::tke_source(rho, Real(0.6), f, wind, beta_p, Real(0), Real(0.1)),
                rho * beta_p * f * wind * wind * wind, tol(Real(1000)));
}

// For beta_d f |U| dt << 1 the source tends to the explicit
// rho beta_p f |U|^3 - beta_d f |U| (rho k)
TEST(CanopyTurbulence, SourceSmallTimeStep)
{
    const Real rho = Real(1.2), rhok = Real(0.6), f = Real(0.01), wind = Real(2);
    const Real dt = Real(1.0e-4);
    const Real explicit_src = rho * f * wind * wind * wind - Real(4) * f * wind * rhok;
    const Real src = canopy_turbulence::tke_source(rho, rhok, f, wind, Real(1), Real(4), dt);
    EXPECT_NEAR(src, explicit_src, Real(1.0e-5) * std::abs(explicit_src));
}

// The stage update (rho k)^n + dt S is the exact solution with production and
// sink frozen, (rho k)^n e^-x + rho (P / c) (1 - e^-x), for small and large x
TEST(CanopyTurbulence, StageUpdateIsExact)
{
    const Real rho = Real(1.2), rhok = Real(0.6), f = Real(0.2), wind = Real(5);
    const Real beta_p = Real(1), beta_d = Real(4);
    const Real P = beta_p * f * wind * wind * wind;
    const Real c = beta_d * f * wind;
    for (const Real dt : {Real(1.0e-3), Real(0.1), Real(1), Real(50)}) {
        const Real x = c * dt;
        const Real exact = rhok * std::exp(-x) + rho * (P / c) * (-std::expm1(-x));
        const Real updated = rhok + dt * canopy_turbulence::tke_source(rho, rhok, f, wind,
                                                                       beta_p, beta_d, dt);
        EXPECT_NEAR(updated, exact, std::abs(exact) * tol(Real(1000))) << "dt = " << dt;
    }
}

// The exact balance rho P / c is a steady state at any time step: the source is
// zero there, so an RK stage never moves it (the per-stage form it replaces
// settled at 1 + c dt / 2 times the balance)
TEST(CanopyTurbulence, BalanceIsSteadyAtAnyTimeStep)
{
    const Real rho = Real(1.2), f = Real(0.2), wind = Real(5);
    const Real rhok_eq = rho * (f * wind * wind * wind) / (Real(4) * f * wind);
    for (const Real dt : {Real(0.01), Real(0.5), Real(2), Real(20)}) {
        const Real src = canopy_turbulence::tke_source(rho, rhok_eq, f, wind, Real(1), Real(4), dt);
        EXPECT_NEAR(src, Real(0), rho * f * wind * wind * wind * tol(Real(1000))) << "dt = " << dt;
    }
}

// Without production, a dense canopy and a long step: the stage update tends to
// zero and never removes more than (rho k)^n, to round-off
TEST(CanopyTurbulence, SinkBoundedByStartOfStepTKE)
{
    const Real rho = Real(1.2), rhok = Real(0.6);
    for (const Real dt : {Real(0.1), Real(10), Real(1000)}) {
        const Real src = canopy_turbulence::tke_source(rho, rhok, Real(1), Real(20), Real(0), Real(4), dt);
        // (rho k)^n - ((rho k)^n / dt) dt is zero only to round-off (about 5e-9 in
        // single precision for dt = 1000)
        EXPECT_GE(rhok + dt * src, -rhok * tol(Real(10))) << "dt = " << dt;
    }
    const Real src = canopy_turbulence::tke_source(rho, rhok, Real(1), Real(20), Real(0), Real(4), Real(10));
    EXPECT_NEAR(rhok + Real(10) * src, Real(0), rhok * tol(Real(1000)));
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
    fill_face_velocity(xvel);
    fill_forest_first_cell(forest, f);

    CanopyTurbulence canopy;
    canopy.enabled = true;
    add_canopy_sources(source, state, xvel, yvel, zvel, forest, nullptr, canopy, dt);

    const Real wind0 = std::sqrt(Real(3) * Real(3) + Real(1));   // u = (2 + 4) / 2, v = 1
    const Real expected = canopy_turbulence::tke_source(rho, rho * tke, f, wind0, canopy.beta_p,
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

// The same box with the forest cell fully immersed (terrain_blank = 1): no source
// anywhere. A partly immersed cell (0 < terrain_blank < 1) keeps its source.
TEST(CanopyTurbulence, KernelSkipsImmersedCells)
{
    using namespace amrex;

    const Box domain(IntVect(0, 0, 0), IntVect(1, 0, 0));
    const BoxArray ba(domain);
    const DistributionMapping dm(ba);
    const int ncomp = RhoKE_comp + 1;

    MultiFab source(ba, dm, ncomp, 0);
    MultiFab state(ba, dm, ncomp, 0);
    MultiFab forest(ba, dm, 1, 0);
    MultiFab blank(ba, dm, 1, 0);
    MultiFab xvel(convert(ba, IntVect(1, 0, 0)), dm, 1, 0);
    MultiFab yvel(convert(ba, IntVect(0, 1, 0)), dm, 1, 0);
    MultiFab zvel(convert(ba, IntVect(0, 0, 1)), dm, 1, 0);

    const Real rho = Real(1.2), tke = Real(0.5), f = Real(0.05), dt = Real(0.1);
    state.setVal(Real(0));
    state.setVal(rho, Rho_comp, 1);
    state.setVal(rho * tke, RhoKE_comp, 1);
    yvel.setVal(Real(1));
    zvel.setVal(Real(0));
    fill_face_velocity(xvel);
    fill_forest_first_cell(forest, f);

    CanopyTurbulence canopy;
    canopy.enabled = true;

    source.setVal(Real(0));
    blank.setVal(Real(1));
    add_canopy_sources(source, state, xvel, yvel, zvel, forest, &blank, canopy, dt);
    for (int n = 0; n < ncomp; ++n) {
        EXPECT_EQ(source.norm0(n), Real(0)) << "component " << n;
    }

    source.setVal(Real(0));
    blank.setVal(Real(0.5));
    add_canopy_sources(source, state, xvel, yvel, zvel, forest, &blank, canopy, dt);
    const Real wind0 = std::sqrt(Real(3) * Real(3) + Real(1));
    const Real expected = canopy_turbulence::tke_source(rho, rho * tke, f, wind0, canopy.beta_p,
                                                        canopy.beta_d, dt);
    EXPECT_NEAR(source.max(RhoKE_comp), expected, std::abs(expected) * tol());
}
