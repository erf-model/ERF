// Running statistics of a body's diagnostics: mean, root mean square, minimum and maximum of
// the samples since the averaging start, exact for a known sequence, written as a CSV, carried
// across a checkpoint round trip; a non-finite sample or a malformed checkpoint aborts, naming
// the statistics or the file; a quantity renamed since the checkpoint (an older span checkpoint's
// drag_y) continues when it is the same quantity and starts afresh when it is not.

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <AMReX_REAL.H>

#include <gtest/gtest.h>

#include "../ERF_GTestTempDir.H"
#include "ERF_GTestThrowOnAbort.H"
#include "ERF_RunningStats.H"

namespace {

// a scratch root drawn once per test process (ERF_GTestTempDir.H): the fixed names below it are this
// process's alone, so ctest -j and the shuffled rerun never share them; removed when the process exits
const std::filesystem::path& gtest_scratch_root ()
{
    struct Root {
        std::filesystem::path p;
        ~Root () { std::error_code ec; std::filesystem::remove_all(p, ec); }
    };
    static const Root root{[] {
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_runningstats");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}
using amrex::Real;
using erf_actuator::RunningStats;
Real tol () { return (sizeof(Real) == 8) ? Real(1.0e-12) : Real(1.0e-5); }
}

TEST(RunningStats, MeanRmsMinMaxOfAKnownSequence)
{
    RunningStats s("T1", "out/T1", {"a", "b"});
    EXPECT_EQ(s.num_samples(), 0);
    EXPECT_EQ(s.mean(0), 0.0);
    // a: 1, 2, 3, 4  -> mean 2.5, rms sqrt(7.5), min 1, max 4; b: constant 7
    for (int k = 1; k <= 4; ++k) { s.accumulate(10.0 * k, {static_cast<Real>(k), Real(7.0)}); }
    EXPECT_EQ(s.num_samples(), 4);
    EXPECT_DOUBLE_EQ(s.t_first(), 10.0);
    EXPECT_DOUBLE_EQ(s.t_last(), 40.0);
    EXPECT_NEAR(s.mean(0), 2.5, tol());
    EXPECT_NEAR(s.rms(0), std::sqrt(7.5), tol());
    EXPECT_EQ(s.min(0), 1.0);
    EXPECT_EQ(s.max(0), 4.0);
    EXPECT_NEAR(s.mean(1), 7.0, tol());
    EXPECT_NEAR(s.rms(1), 7.0, tol());
    EXPECT_EQ(s.min(1), 7.0);
    EXPECT_EQ(s.max(1), 7.0);
}

TEST(RunningStats, StateRoundTripsThroughACheckpointAndTheFileHasTheRows)
{
    const auto dir = gtest_scratch_root() / "erf_gtest_stats";
    std::filesystem::create_directories(dir);
    RunningStats a("T1", (dir / "T1").string(), {"thrust", "power"});
    a.accumulate(1.0, {Real(2.0e6), Real(1.0e7)});
    a.accumulate(2.0, {Real(3.0e6), Real(1.5e7)});
    a.write_state(dir.string());
    ASSERT_TRUE(std::filesystem::exists(dir / "T1_stats.dat"));

    RunningStats b("T1", (dir / "T1").string(), {"thrust", "power"});
    EXPECT_TRUE(b.read_state(dir.string()));
    EXPECT_EQ(b.num_samples(), 2);
    EXPECT_EQ(b.mean(0), a.mean(0));
    EXPECT_EQ(b.rms(1), a.rms(1));
    EXPECT_EQ(b.min(0), a.min(0));
    EXPECT_EQ(b.max(1), a.max(1));
    EXPECT_DOUBLE_EQ(b.t_first(), 1.0);
    // a third sample continues the same statistics
    a.accumulate(3.0, {Real(4.0e6), Real(2.0e7)});
    b.accumulate(3.0, {Real(4.0e6), Real(2.0e7)});
    EXPECT_EQ(b.mean(0), a.mean(0));
    EXPECT_EQ(b.rms(0), a.rms(0));
    EXPECT_NEAR(a.mean(0), 3.0e6, tol() * 3.0e6);
    // no state for another body
    RunningStats c("T2", (dir / "T2").string(), {"thrust", "power"});
    EXPECT_FALSE(c.read_state(dir.string()));

    a.write();
    std::ifstream csv(dir / "T1_stats.csv");
    std::string line;
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line, "samples,t_first,t_last,quantity,mean,rms,min,max");
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line.rfind("3,1,3,thrust,3000000,", 0), 0u) << line;
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line.rfind("3,1,3,power,", 0), 0u) << line;
}

TEST(RunningStats, ANonFiniteSampleOrAMalformedCheckpointIsRefused)
{
    RunningStats s("L1_span1", "out/L1_span1", {"swing", "tension"});
    s.accumulate(1.0, {Real(3.0), Real(1.0e4)});
    std::string msg = erf_gtest::abort_message([&] { s.accumulate(2.0, {Real(3.0), std::numeric_limits<Real>::quiet_NaN()}); });
    EXPECT_NE(msg.find("RunningStats L1_span1"), std::string::npos) << msg;
    EXPECT_NE(msg.find("tension"), std::string::npos) << msg;
    EXPECT_EQ(s.num_samples(), 1) << "a refused sample is not counted";

    const auto dir = gtest_scratch_root() / "erf_gtest_stats_bad";
    std::filesystem::create_directories(dir);
    s.write_state(dir.string());
    std::string good;
    {
        std::ifstream in(dir / "L1_span1_stats.dat");
        std::stringstream ss;
        ss << in.rdbuf();
        good = ss.str();
    }
    auto refused = [&] (const std::string& from, const std::string& to) {
        std::string text = good;
        const auto at = text.find(from);
        EXPECT_NE(at, std::string::npos) << from;
        if (at != std::string::npos) { text.replace(at, from.size(), to); }
        std::ofstream(dir / "L1_span1_stats.dat", std::ios::trunc) << text;
        RunningStats b("L1_span1", "out/L1_span1", {"swing", "tension"});
        const std::string m = erf_gtest::abort_message([&] { b.read_state(dir.string()); });
        EXPECT_EQ(b.num_samples(), 0) << "a refused checkpoint changes nothing";
        EXPECT_EQ(b.mean(0), 0.0);
        return m;
    };
    EXPECT_NE(refused("count = 1", "count = -1").find("malformed statistics checkpoint"), std::string::npos);
    EXPECT_NE(refused("t_last = 1", "t_last = 0.5").find("t_first <= t_last"), std::string::npos);
    EXPECT_NE(refused("tension", "power").find("does not list tension"), std::string::npos);
    EXPECT_NE(refused("size = 2", "size = 3").find("holds 3 quantities"), std::string::npos);
}

// A checkpoint may list a span's ERF-frame y drag (drag_y) where its drag across the span (drag_normal)
// stands, the same quantity only for a span along +x: a checkpoint holding drag_y continues as
// drag_normal where the renaming is accepted as the same quantity, starts afresh where it is accepted as
// another, and is refused where it is not accepted, as is any other renaming.
TEST(RunningStats, ARenamedQuantityContinuesOnlyWhenItIsTheSame)
{
    const auto dir = erf_gtest_temp_path("erf_gtest_stats_drag");
    std::filesystem::create_directories(dir);
    RunningStats old_names("S1", (dir / "S1").string(), {"swing_deg", "drag_y"});
    old_names.accumulate(1.0, {Real(10.0), Real(250.0)});
    old_names.accumulate(2.0, {Real(12.0), Real(270.0)});
    old_names.write_state(dir.string());
    RunningStats same("S1", (dir / "S1").string(), {"swing_deg", "drag_normal"});
    same.renamed_from("drag_y", "drag_normal", true);
    EXPECT_TRUE(same.read_state(dir.string()));
    EXPECT_EQ(same.num_samples(), 2);
    EXPECT_EQ(same.mean(1), old_names.mean(1));
    RunningStats another("S1", (dir / "S1").string(), {"swing_deg", "drag_normal"});
    another.renamed_from("drag_y", "drag_normal", false);
    EXPECT_TRUE(another.read_state(dir.string()));
    EXPECT_EQ(another.num_samples(), 0) << "drag_y is another quantity here: the statistics start afresh";
    RunningStats unaccepted("S1", (dir / "S1").string(), {"swing_deg", "drag_normal"});
    std::string msg = erf_gtest::abort_message([&] { unaccepted.read_state(dir.string()); });
    EXPECT_NE(msg.find("does not list drag_normal where expected"), std::string::npos) << msg;
    RunningStats other("S1", (dir / "S1").string(), {"swing_deg", "drag_x"});
    other.renamed_from("drag_y", "drag_normal", true);
    msg = erf_gtest::abort_message([&] { other.read_state(dir.string()); });
    EXPECT_NE(msg.find("does not list drag_x where expected"), std::string::npos) << msg;
    std::filesystem::remove_all(dir);
}
