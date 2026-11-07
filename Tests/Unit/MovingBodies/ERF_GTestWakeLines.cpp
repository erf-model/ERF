// Wake sampling lines behind a rotor: the points sit on the lateral and vertical lines at the
// requested x/D behind the hub along the rotor axis, the running average is the mean of the
// samples (an analytic Gaussian wake plus a zero-mean perturbation comes back exactly), the
// running-average state survives a checkpoint round trip, and a malformed checkpoint, an empty
// list of distances or a non-finite sample aborts, naming the file or the rotor.

#include <array>
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
#include "ERF_WakeLines.H"

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
        const std::filesystem::path p = erf_gtest_temp_path("erf_gtest_wakelines");
        std::filesystem::create_directories(p);
        return p;
    }()};
    return root.p;
}

using amrex::Real;
using erf_actuator::WakeLines;

Real tol () { return (sizeof(Real) == 8) ? Real(1.0e-10) : Real(1.0e-4); }

// a Gaussian wake deficit behind the rotor: u = U (1 - a exp(-r^2 / (2 sigma^2))) along the axis
std::vector<Real> gaussian_wake (const WakeLines& w, Real U, Real a, Real sigma)
{
    const auto& pos = w.positions();
    const auto& hub = w.hub();
    const auto& n = w.axis();
    std::vector<Real> vel(pos.size(), 0.0);
    for (std::size_t p = 0; p < pos.size() / 3; ++p) {
        std::array<Real,3> rv;
        for (int d = 0; d < 3; ++d) { rv[d] = pos[3*p+d] - hub[d]; }
        const Real along = rv[0]*n[0] + rv[1]*n[1] + rv[2]*n[2];
        Real r2 = 0.0;
        for (int d = 0; d < 3; ++d) { const Real t = rv[d] - along * n[d]; r2 += t * t; }
        const Real u = U * (1.0 - a * std::exp(-r2 / (2.0 * sigma * sigma)));
        for (int d = 0; d < 3; ++d) { vel[3*p+d] = u * n[d]; }
    }
    return vel;
}

} // namespace

TEST(WakeLines, PointsLieOnTheLateralAndVerticalLinesBehindTheHub)
{
    const std::array<Real,3> hub{{900.0, 800.0, 150.0}};
    // a yawed and tilted axis: the lines follow its horizontal projection (a wake follows the
    // wind at hub height, not the shaft tilt); the lateral line is horizontal and normal to
    // it, the vertical line is z
    const std::array<Real,3> shaft{{static_cast<Real>(std::cos(0.3) * std::cos(0.1)),
                                    static_cast<Real>(std::sin(0.3) * std::cos(0.1)),
                                    static_cast<Real>(-std::sin(0.1))}};
    const std::array<Real,3> axis{{static_cast<Real>(std::cos(0.3)), static_cast<Real>(std::sin(0.3)), Real(0.0)}};
    const Real D = 240.0;
    const std::vector<Real> xD{2.0, 4.0};
    const int npts = 5;
    WakeLines w("T1", "out/T1", hub, shaft, D, xD, 1.5, npts, Real(-1.0e30));
    ASSERT_EQ(w.num_points(), 2 * 2 * npts);
    for (int d = 0; d < 3; ++d) { EXPECT_NEAR(w.axis()[d], axis[d], tol()); }
    const auto& pos = w.positions();
    const std::array<Real,3> lat{{static_cast<Real>(-std::sin(0.3)), static_cast<Real>(std::cos(0.3)), Real(0.0)}};
    int p = 0;
    for (const Real x : xD) {
        for (int line = 0; line < 2; ++line) {
            for (int i = 0; i < npts; ++i, ++p) {
                const Real s = (-1.5 + 3.0 * i / (npts - 1)) * D;
                for (int d = 0; d < 3; ++d) {
                    const Real e = (line == 0) ? lat[d] : ((d == 2) ? Real(1.0) : Real(0.0));
                    EXPECT_NEAR(pos[3*p+d], hub[d] + x * D * axis[d] + s * e, tol() * D) << "point " << p << " dir " << d;
                }
            }
        }
    }
    // the middle point of each line is on the axis
    for (int line = 0; line < 4; ++line) {
        const int mid = line * npts + npts / 2;
        const Real x = xD[line / 2];
        for (int d = 0; d < 3; ++d) { EXPECT_NEAR(pos[3*mid+d], hub[d] + x * D * axis[d], tol() * D); }
    }
}

TEST(WakeLines, VerticalLineIsClippedAtTheGround)
{
    // hub 150 m, D 240 m, +-1.5 D: the lateral line keeps its full span, the vertical line
    // starts at the ground (s = -150/240) and ends at +1.5 D
    const std::array<Real,3> hub{{500.0, 500.0, 150.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    const int npts = 7;
    WakeLines w("T1", "out/T1", hub, axis, 240.0, {2.0}, 1.5, npts, Real(0.0));
    const auto& pos = w.positions();
    EXPECT_NEAR(pos[3*0 + 1], 500.0 - 1.5 * 240.0, tol() * 240.0);             // lateral start
    EXPECT_NEAR(pos[3*(npts-1) + 1], 500.0 + 1.5 * 240.0, tol() * 240.0);      // lateral end
    EXPECT_NEAR(pos[3*npts + 2], 0.0, tol() * 240.0);                          // vertical start on the ground
    EXPECT_NEAR(pos[3*(2*npts-1) + 2], 150.0 + 1.5 * 240.0, tol() * 240.0);    // vertical end
    for (int i = 0; i < 2 * npts; ++i) { EXPECT_GE(pos[3*i + 2], -tol()) << "point " << i << " below ground"; }
    // the file records the clipped offset of the first vertical point
    const auto dir = gtest_scratch_root() / "erf_gtest_wake_clip";
    std::filesystem::create_directories(dir);
    WakeLines wf("T1", (dir / "T1").string(), hub, axis, 240.0, {2.0}, 1.5, npts, Real(0.0));
    wf.write_average();
    std::ifstream csv(dir / "T1_wake_avg.csv");
    std::string line;
    std::getline(csv, line);
    for (int i = 0; i < npts; ++i) { std::getline(csv, line); }
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line.rfind("0,2,vertical,-0.625,", 0), 0u) << line;
}

TEST(WakeLines, EachVerticalLineIsClippedAtItsOwnGround)
{
    // terrain rising downstream: 100 m under the 2 D line, 200 m under the 4 D line; each vertical
    // line starts at its own ground, the lateral lines are untouched, and a hub-height ground
    // (the scalar form) would have put the 4 D line's lowest points inside the hill
    const std::array<Real,3> hub{{500.0, 500.0, 250.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    const int npts = 7;
    WakeLines w("T1", "out/T1", hub, axis, 240.0, {2.0, 4.0}, 1.5, npts, std::vector<Real>{100.0, 200.0});
    const auto& pos = w.positions();
    ASSERT_EQ(pos.size(), std::size_t(3 * 4 * npts));
    EXPECT_NEAR(pos[3*npts + 2], 100.0, tol() * 240.0);                           // 2 D vertical start
    EXPECT_NEAR(pos[3*(3*npts) + 2], 200.0, tol() * 240.0);                       // 4 D vertical start
    EXPECT_NEAR(pos[3*(2*npts) + 1], 500.0 - 1.5 * 240.0, tol() * 240.0);         // 4 D lateral start, full span
    EXPECT_NEAR(pos[3*(4*npts-1) + 2], 250.0 + 1.5 * 240.0, tol() * 240.0);       // 4 D vertical end
    for (int i = 0; i < 4 * npts; ++i) {
        const Real ground = (i < 2 * npts) ? Real(100.0) : Real(200.0);
        if ((i / npts) % 2 == 1) { EXPECT_GE(pos[3*i + 2], ground - tol()) << "point " << i << " below its ground"; }
    }
    WakeLines s("T1", "out/T1", hub, axis, 240.0, {2.0, 4.0}, 1.5, npts, Real(100.0));
    EXPECT_NEAR(s.positions()[3*(3*npts) + 2], 100.0, tol() * 240.0);             // the scalar form: the hub's ground
    erf_gtest::ThrowOnAbort guard;                                                           // one ground per distance
    EXPECT_THROW(WakeLines("T1", "out/T1", hub, axis, 240.0, {2.0, 4.0}, 1.5, npts, std::vector<Real>{100.0}),
                 std::runtime_error);
}

TEST(WakeLines, VerticalLineIsClippedAtTheTerrainUnderTheHub)
{
    // the same rotor on a 100 m hill: the vertical line starts at the terrain surface, not at z = 0
    const std::array<Real,3> hub{{500.0, 500.0, 250.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    const int npts = 7;
    WakeLines w("T1", "out/T1", hub, axis, 240.0, {2.0}, 1.5, npts, Real(100.0));
    const auto& pos = w.positions();
    EXPECT_NEAR(pos[3*npts + 2], 100.0, tol() * 240.0);                        // vertical start on the hill
    EXPECT_NEAR(pos[3*(2*npts-1) + 2], 250.0 + 1.5 * 240.0, tol() * 240.0);    // vertical end unchanged
    for (int i = npts; i < 2 * npts; ++i) { EXPECT_GE(pos[3*i + 2], 100.0 - tol()) << "point " << i << " below the hill"; }
}

TEST(WakeLines, RunningAverageRecoversTheAnalyticWake)
{
    const std::array<Real,3> hub{{0.0, 0.0, 150.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    WakeLines w("T1", "out/T1", hub, axis, 240.0, {2.0, 7.0}, 1.0, 21, Real(-1.0e30));
    const std::vector<Real> exact = gaussian_wake(w, 10.0, 0.4, 60.0);
    EXPECT_EQ(w.num_samples(), 0);
    for (const Real v : w.average()) { EXPECT_EQ(v, 0.0); }
    // three samples with perturbations -1, 0, +1 times a per-point amplitude: the mean is exact
    for (int k = -1; k <= 1; ++k) {
        std::vector<Real> vel = exact;
        for (std::size_t i = 0; i < vel.size(); ++i) { vel[i] += k * Real(0.1) * (i % 7); }
        w.accumulate(vel);
    }
    EXPECT_EQ(w.num_samples(), 3);
    const auto avg = w.average();
    for (std::size_t i = 0; i < avg.size(); ++i) { EXPECT_NEAR(avg[i], exact[i], 30.0 * tol()) << "entry " << i; }
    // the deficit at the axis point of the 2 D line is a = 0.4, and 0.4 at 7 D too (this wake
    // does not recover): the axis point is the middle of each line
    const int mid2 = 21 / 2, mid7 = 2 * 21 + 21 / 2;
    EXPECT_NEAR(avg[3*mid2], 6.0, 30.0 * tol());
    EXPECT_NEAR(avg[3*mid7], 6.0, 30.0 * tol());
}

TEST(WakeLines, StateRoundTripsThroughACheckpoint)
{
    const auto dir = gtest_scratch_root() / "erf_gtest_wake";
    std::filesystem::create_directories(dir);
    const std::array<Real,3> hub{{0.0, 0.0, 150.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    WakeLines a("T1", (dir / "T1").string(), hub, axis, 240.0, {2.0}, 1.0, 7, Real(-1.0e30));
    const std::vector<Real> v1 = gaussian_wake(a, 10.0, 0.3, 50.0);
    std::vector<Real> v2 = v1;
    for (auto& v : v2) { v *= Real(0.5); }
    a.accumulate(v1);
    a.accumulate(v2);
    a.write_state(dir.string());
    ASSERT_TRUE(std::filesystem::exists(dir / "T1_wake_avg.dat"));

    // b is built from a slightly different hub and diameter (as a restarted turbine's float
    // state gives): the checkpoint brings back a's geometry exactly
    const std::array<Real,3> hub_b{{0.0, Real(1.0e-4), Real(150.0 + 2.0e-5)}};
    WakeLines b("T1", (dir / "T1").string(), hub_b, axis, Real(240.0 - 1.0e-4), {2.0}, 1.0, 7, Real(-1.0e30));
    EXPECT_NE(b.positions()[0], a.positions()[0]);
    EXPECT_TRUE(b.read_state(dir.string()));
    EXPECT_EQ(b.num_samples(), 2);
    EXPECT_EQ(b.diameter(), a.diameter());
    for (std::size_t i = 0; i < a.positions().size(); ++i) { EXPECT_EQ(b.positions()[i], a.positions()[i]) << "position " << i; }
    const auto aa = a.average(), ba = b.average();
    for (std::size_t i = 0; i < aa.size(); ++i) { EXPECT_EQ(ba[i], aa[i]) << "entry " << i; }
    // a third sample continues the same average
    b.accumulate(v1);
    a.accumulate(v1);
    const auto aa3 = a.average(), ba3 = b.average();
    for (std::size_t i = 0; i < aa3.size(); ++i) { EXPECT_EQ(ba3[i], aa3[i]) << "entry " << i; }
    // no state for another body
    WakeLines c("T2", (dir / "T2").string(), hub, axis, 240.0, {2.0}, 1.0, 7, Real(-1.0e30));
    EXPECT_FALSE(c.read_state(dir.string()));

    // the average file holds the sample count and the mean
    a.write_average();
    std::ifstream csv(dir / "T1_wake_avg.csv");
    std::string line;
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line, "samples,xD,line,s,x,y,z,u,v,w");
    ASSERT_TRUE(std::getline(csv, line));
    EXPECT_EQ(line.rfind("3,2,lateral,-1,", 0), 0u) << line;
}

TEST(WakeLines, AMalformedCheckpointOrInputIsRefusedNamingIt)
{
    const auto dir = gtest_scratch_root() / "erf_gtest_wake_bad";
    std::filesystem::create_directories(dir);
    const std::array<Real,3> hub{{0.0, 0.0, 150.0}};
    const std::array<Real,3> axis{{1.0, 0.0, 0.0}};
    WakeLines a("T1", (dir / "T1").string(), hub, axis, 240.0, {2.0}, 1.0, 7, Real(-1.0e30));
    a.accumulate(gaussian_wake(a, 10.0, 0.3, 50.0));
    a.write_state(dir.string());
    std::string good;
    {
        std::ifstream in(dir / "T1_wake_avg.dat");
        std::stringstream ss;
        ss << in.rdbuf();
        good = ss.str();
    }
    // rewrite the checkpoint with one entry changed and read it back into a fresh object
    auto refused = [&] (const std::string& from, const std::string& to) {
        std::string text = good;
        const auto at = text.find(from);
        EXPECT_NE(at, std::string::npos) << from;
        if (at != std::string::npos) { text.replace(at, from.size(), to); }
        std::ofstream(dir / "T1_wake_avg.dat", std::ios::trunc) << text;
        WakeLines b("T1", (dir / "T1").string(), hub, axis, 240.0, {2.0}, 1.0, 7, Real(-1.0e30));
        const std::string msg = erf_gtest::abort_message([&] { b.read_state(dir.string()); });
        EXPECT_EQ(b.num_samples(), 0) << "a refused checkpoint changes nothing";
        return msg;
    };
    EXPECT_NE(refused("count = 1", "count = -1").find("malformed wake-average checkpoint"), std::string::npos);
    EXPECT_NE(refused("diameter = 240", "diameter = -240").find("diameter positive"), std::string::npos);
    EXPECT_NE(refused("axis = 1 0 0", "axis = 0.6 0 0.8").find("horizontal unit vector"), std::string::npos);
    EXPECT_NE(refused("size = 42", "size = 41").find("41 values"), std::string::npos);
    // a truncated file
    {
        std::ofstream(dir / "T1_wake_avg.dat", std::ios::trunc) << good.substr(0, good.size() / 2);
        WakeLines b("T1", (dir / "T1").string(), hub, axis, 240.0, {2.0}, 1.0, 7, Real(-1.0e30));
        EXPECT_NE(erf_gtest::abort_message([&] { b.read_state(dir.string()); }).find("truncated"), std::string::npos);
    }
    // the inputs: no distances, and a non-finite sample
    EXPECT_NE(erf_gtest::abort_message([&] { WakeLines c("T3", "out/T3", hub, axis, 240.0, {}, 1.0, 7, Real(0.0)); }).find("WakeLines T3"),
              std::string::npos);
    std::vector<Real> vel(a.positions().size(), Real(1.0));
    vel[5] = std::numeric_limits<Real>::quiet_NaN();
    const std::string msg = erf_gtest::abort_message([&] { a.accumulate(vel); });
    EXPECT_NE(msg.find("point 1 (0-based)"), std::string::npos) << msg;
    EXPECT_EQ(a.num_samples(), 1) << "a refused sample is not counted";
}
