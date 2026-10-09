// Contract of WriteBndryPlanes::rows_to_keep, which decides how much of an existing time.dat a run
// continues from. A boundary-plane series is a time.dat of "step time" rows, one per plane, read
// back by erf.input_bndry_planes, which stops on a step that does not increase. A restart keeps
// the rows up to and including its own step (the run before it wrote those planes) and drops the
// later ones, which the restart writes again; a fresh start keeps none. The CTests
// ABL_BndryPlanes_* check the same behaviour end to end against a run without a restart.

#include <gtest/gtest.h>

#include <ERF_WriteBndryPlanes.H>

namespace {

using Steps = amrex::Vector<int>;

TEST(WriteBndryPlanes, FreshStartKeepsNoRows)
{
    // a run in a directory an earlier run wrote into begins a new series
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4, 6}, 0, false), 0);
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{}, 0, false), 0);
}

TEST(WriteBndryPlanes, RestartOnAnOutputStepKeepsItsRow)
{
    // checkpoint at step 4, planes every 2 steps: the step-4 plane is there already
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4}, 4, true), 3);
}

TEST(WriteBndryPlanes, RestartBetweenOutputStepsKeepsTheEarlierRows)
{
    // checkpoint at step 3: rows 0 and 2 stay, the next plane comes at step 4
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2}, 3, true), 2);
}

TEST(WriteBndryPlanes, RestartDropsTheRowsItReplays)
{
    // the run went on to step 8 after its step-3 checkpoint; the restart writes 4, 6 and 8 again
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4, 6, 8}, 3, true), 2);
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4, 6, 8}, 4, true), 3);
}

TEST(WriteBndryPlanes, RestartIntoAnEmptySeriesKeepsNoRows)
{
    // output switched on at the restart: the start-up plane begins the series
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{}, 4, true), 0);
    // a series that starts after the restart step belongs to a later run
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{6, 8}, 4, true), 0);
}

TEST(WriteBndryPlanes, RestartStopsAtAStepThatDoesNotIncrease)
{
    // the "0 t" row an earlier restart appended: nothing after it is kept
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 4, 0, 6, 8}, 8, true), 3);
    EXPECT_EQ(WriteBndryPlanes::rows_to_keep(Steps{0, 2, 2, 4}, 4, true), 2);
}

} // namespace
