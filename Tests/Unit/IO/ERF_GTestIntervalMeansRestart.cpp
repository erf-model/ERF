#include <sstream>
#include <string>

#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>

#include <gtest/gtest.h>

#include "ERF_IntervalMeansCheckpoint.H"

namespace {

std::string
make_metadata (const amrex::Vector<amrex::BoxArray>& boxes,
               const amrex::Vector<double>& counts,
               const amrex::Vector<int>& reset_flags)
{
    std::ostringstream output;
    output << "ERF interval means checkpoint v1\n";
    output << boxes.size() << " 10\n";
    for (int lev = 0; lev < static_cast<int>(boxes.size()); ++lev) {
        output << lev << " " << counts[lev] << " " << reset_flags[lev] << "\n";
        boxes[lev].writeOn(output);
        output << "\n";
        output << boxes[lev].size();
        for (int ibox = 0; ibox < boxes[lev].size(); ++ibox) {
            output << " " << (ibox + lev);
        }
        output << "\n";
    }
    return output.str();
}

amrex::Vector<amrex::BoxArray>
make_boxes ()
{
    amrex::Vector<amrex::BoxArray> boxes;
    boxes.emplace_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(3, 3, 1)));
    boxes.emplace_back(amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(7, 7, 3)));
    boxes[0].maxSize(2);
    boxes[1].maxSize(4);
    return boxes;
}

} // namespace

TEST(IntervalMeansRestart, ParsesVersionCountsBoxesAndProcessorMaps)
{
    const auto boxes = make_boxes();
    const std::string text = make_metadata(boxes, {3.5, 8.25}, {0, 1});
    erf_interval_means::Metadata metadata;
    std::string error;
    std::istringstream input(text);
    ASSERT_TRUE(erf_interval_means::parse_metadata(input, metadata, error)) << error;
    EXPECT_EQ(metadata.levels, 2);
    EXPECT_EQ(metadata.components, 10);
    EXPECT_DOUBLE_EQ(metadata.level[0].accumulation_count, 3.5);
    EXPECT_DOUBLE_EQ(metadata.level[1].accumulation_count, 8.25);
    EXPECT_EQ(metadata.level[0].processor_map[0], 0);
    EXPECT_EQ(metadata.level[1].processor_map[0], 1);
    EXPECT_EQ(erf_interval_means::global_reset_done(metadata), 1);
    EXPECT_TRUE(erf_interval_means::validate_metadata(metadata, 2, 10, boxes).empty());
}

TEST(IntervalMeansRestart, BoxArrayMismatchIsRejectedButProcessorMapIsNotAConstraint)
{
    const auto boxes = make_boxes();
    erf_interval_means::Metadata metadata;
    std::string error;
    std::istringstream input(make_metadata(boxes, {1.0, 2.0}, {0, 0}));
    ASSERT_TRUE(erf_interval_means::parse_metadata(input, metadata, error)) << error;

    auto changed_boxes = boxes;
    changed_boxes[0] = amrex::BoxArray(
        amrex::Box(amrex::IntVect(0, 0, 0), amrex::IntVect(4, 3, 1)));
    const auto mismatch = erf_interval_means::validate_metadata(
        metadata, 2, 10, changed_boxes);
    EXPECT_NE(mismatch.find("BoxArray"), std::string::npos);
    EXPECT_TRUE(erf_interval_means::validate_metadata(metadata, 2, 10, boxes).empty());
}

TEST(IntervalMeansRestart, RejectsMalformedAndTruncatedHeaders)
{
    for (const std::string& text : {
             std::string("ERF interval means checkpoint v0\n1 10\n"),
             std::string("ERF interval means checkpoint v1\n1 10\n0 1.0 0\n"),
             std::string("ERF interval means checkpoint v1\n1 10\n0 -1.0 0\n")}) {
        erf_interval_means::Metadata metadata;
        std::string error;
        std::istringstream input(text);
        EXPECT_FALSE(erf_interval_means::parse_metadata(input, metadata, error));
        EXPECT_FALSE(error.empty());
    }
}

TEST(IntervalMeansRestart, LegacyCheckpointFallbackPreservesTimeResetSemantics)
{
    EXPECT_EQ(erf_interval_means::legacy_reset_done(true, 10.0, 10.0), 1);
    EXPECT_EQ(erf_interval_means::legacy_reset_done(true, 9.0, 10.0), 0);
    EXPECT_EQ(erf_interval_means::legacy_reset_done(false, 20.0, 10.0), 0);
}

TEST(IntervalMeansRestart, InitializationPreservesRestoredWindow)
{
    EXPECT_TRUE(erf_interval_means::initialization_accumulates_state(false));
    EXPECT_FALSE(erf_interval_means::initialization_accumulates_state(true));
    EXPECT_TRUE(erf_interval_means::initialization_plot_consumes_interval(false, true));
    EXPECT_FALSE(erf_interval_means::initialization_plot_consumes_interval(true, true));
    EXPECT_FALSE(erf_interval_means::initialization_plot_consumes_interval(false, false));
}
