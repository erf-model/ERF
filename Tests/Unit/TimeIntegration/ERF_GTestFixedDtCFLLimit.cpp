// Unit tests for the CFL limit of a fixed time step (ERF_FixedDtCFLLimit.H): erf.fixed_dt must be
// lowered exactly when its Courant number is above min(erf.cfl, 0.8), and then to that Courant
// number (erf-model/ERF#5238).

#include <ERF_FixedDtCFLLimit.H>

#include <gtest/gtest.h>

TEST(FixedDtCFLLimit, WithinLimit)
{
    // Courant number 0.5 and 0.8 (at the limit) with erf.cfl = 0.8
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.25, 2.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.4,  2.0, 0.8), 1.0);
    // no fixed dt, no constraining rate, no CFL limit
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(-1.0, 2.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.5,  0.0, 0.8), 1.0);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.5,  2.0, 0.0), 1.0);
}

TEST(FixedDtCFLLimit, LoweredToCappedCFL)
{
    // erf.cfl = 1 is capped at 0.8: Courant number 1.2 -> 0.8
    EXPECT_DOUBLE_EQ(0.6 * FixedDtCFLScale(0.6, 2.0, 1.0) * 2.0, 0.8);
    // the cap also sets when to lower: with erf.cfl = 1.5, Courant numbers 0.9 and 1.49 are
    // lowered to 0.8 (no step above 0.8 is kept), while 0.8 is kept
    EXPECT_DOUBLE_EQ(0.45  * FixedDtCFLScale(0.45,  2.0, 1.5) * 2.0, 0.8);
    EXPECT_DOUBLE_EQ(0.745 * FixedDtCFLScale(0.745, 2.0, 1.5) * 2.0, 0.8);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(0.4, 2.0, 1.5), 1.0);
    // the default erf.cfl = 0.8 is used as given: Courant number 1 -> 0.8
    EXPECT_DOUBLE_EQ(0.5 * FixedDtCFLScale(0.5, 2.0, 0.8) * 2.0, 0.8);
    // a smaller erf.cfl is used as given: Courant number 1 -> 0.5
    EXPECT_DOUBLE_EQ(0.5 * FixedDtCFLScale(0.5, 2.0, 0.5) * 2.0, 0.5);
    // a step set to the limit by the scaling is not scaled again
    const double dt = 0.5 * FixedDtCFLScale(0.5, 2.0, 0.5);
    EXPECT_DOUBLE_EQ(FixedDtCFLScale(dt, 2.0, 0.5), 1.0);
}

TEST(FixedDtCFLLimit, SameCellsWhateverTheLayout)
{
    // A level re-chopped over the same cells is not a new set of grids for the check (a restart
    // with a smaller amr.max_grid_size); one that covers other cells is
    const amrex::Box domain(amrex::IntVect(0), amrex::IntVect(15));
    amrex::BoxArray two(domain);
    two.maxSize(amrex::IntVect(AMREX_D_DECL(8,16,16)));
    amrex::BoxArray eight(domain);
    eight.maxSize(8);
    EXPECT_FALSE(two == eight);
    EXPECT_TRUE(FixedDtSameCells(two, eight));
    EXPECT_TRUE(FixedDtSameCells(eight, two));

    // half the cells, then the other half with as many cells
    amrex::BoxArray lo(amrex::Box(amrex::IntVect(0), amrex::IntVect(AMREX_D_DECL(7,15,15))));
    amrex::BoxArray hi(amrex::Box(amrex::IntVect(AMREX_D_DECL(8,0,0)), amrex::IntVect(15)));
    EXPECT_FALSE(FixedDtSameCells(two, lo));
    EXPECT_FALSE(FixedDtSameCells(lo, hi));

    // a level not checked yet
    EXPECT_FALSE(FixedDtSameCells(amrex::BoxArray(), two));
}
