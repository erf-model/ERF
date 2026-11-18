// Contract of WriteBndryPlanes::rows_to_keep, which decides how much of an existing time.dat a run
// continues from. A boundary-plane series is a time.dat of "step time" rows, one per plane, read
// back by erf.input_bndry_planes, which stops on a step or time that does not increase. A restart
// keeps the rows up to and including its own step (the run before it wrote those planes) and drops
// the later ones, which the restart writes again; rows that cannot belong to the run being
// continued (not increasing, later than the restart, or at the restart step at another time) end
// the kept rows. A fresh start keeps none. The CTests
// ABL_BndryPlanes_* check the same behaviour end to end against a run without a restart.

#include <gtest/gtest.h>

#include <ERF_WriteBndryPlanes.H>

namespace {

using Steps = amrex::Vector<int>;
using Times = amrex::Vector<double>;

// planes every 2 steps of 0.02 s, as in the CTest deck
Times times_of (const Steps& steps)
{
    Times t;
    for (int s : steps) { t.push_back(0.02 * s); }
    return t;
}

int keep (const Steps& steps, int start_step, bool restarting)
{
    return WriteBndryPlanes::rows_to_keep(steps, times_of(steps), start_step, 0.02 * start_step, restarting);
}

TEST(WriteBndryPlanes, FreshStartKeepsNoRows)
{
    // a run in a directory an earlier run wrote into begins a new series
    EXPECT_EQ(keep(Steps{0, 2, 4, 6}, 0, false), 0);
    EXPECT_EQ(keep(Steps{}, 0, false), 0);
}

TEST(WriteBndryPlanes, RestartOnAnOutputStepKeepsItsRow)
{
    // checkpoint at step 4, planes every 2 steps: the step-4 plane is there already
    EXPECT_EQ(keep(Steps{0, 2, 4}, 4, true), 3);
}

TEST(WriteBndryPlanes, RestartBetweenOutputStepsKeepsTheEarlierRows)
{
    // checkpoint at step 3: rows 0 and 2 stay, the next plane comes at step 4
    EXPECT_EQ(keep(Steps{0, 2}, 3, true), 2);
}

TEST(WriteBndryPlanes, RestartDropsTheRowsItReplays)
{
    // the run went on to step 8 after its step-3 checkpoint; the restart writes 4, 6 and 8 again
    EXPECT_EQ(keep(Steps{0, 2, 4, 6, 8}, 3, true), 2);
    EXPECT_EQ(keep(Steps{0, 2, 4, 6, 8}, 4, true), 3);
}

TEST(WriteBndryPlanes, RestartIntoAnEmptySeriesKeepsNoRows)
{
    // output switched on at the restart
    EXPECT_EQ(keep(Steps{}, 4, true), 0);
    // output on, but its start time came after the checkpoint
    EXPECT_EQ(keep(Steps{4, 6, 8}, 3, true), 0);
}

TEST(WriteBndryPlanes, RestartStopsAtAStepThatDoesNotIncrease)
{
    // the "0 t" row an earlier restart appended: nothing after it is kept
    const Steps steps{0, 2, 4, 0, 6, 8};
    const Times times{0.0, 0.04, 0.08, 0.08, 0.12, 0.16};
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(steps, times, 8, 0.16, true), 3);
    EXPECT_EQ(keep(Steps{0, 2, 2, 4}, 4, true), 2);
}

TEST(WriteBndryPlanes, RestartStopsAtATimeThatDoesNotIncreaseOrIsLate)
{
    // steps increase, a time does not
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, Times{0.0, 0.04, 0.04}, 4, 0.08, true), 2);
    // a row from another run with a larger time step, later than the restart
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, Times{0.0, 0.1, 0.2}, 4, 0.08, true), 1);
    // the row at the restart step is at another time
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, Times{0.0, 0.02, 0.04}, 4, 0.08, true), 2);
    // absolute times from a start_datetime keep their rows
    const double t0 = 1.5778368e9;
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, Times{t0, t0 + 0.04, t0 + 0.08}, 4, t0 + 0.08, true), 3);
    // ... but a row at the restart step 0.1 ms off is another run's, also on a calendar clock
    // (a tolerance relative to the clock, 1e-12 of 1.6e9 s, would have taken it as this run's)
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, Times{t0, t0 + 0.04, t0 + 0.0801}, 4, t0 + 0.08, true), 2);
}

} // namespace
