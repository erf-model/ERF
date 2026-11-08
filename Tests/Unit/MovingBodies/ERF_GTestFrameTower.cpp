// A lattice tower that bends as its frame model (FrameTower), and the rigid links that tie the
// tower's drag nodes and attachments to the frame. The frame is case T of
// Tests/test_files/FrameSubDynTower: a 30 m lattice tower with a 12 m cross-arm in tower-local axes.
//
// RigidLink.KeepsTheForceAndItsMomentAndDoesNoWork: a force shared among a link's nodes keeps its
//   sum and its moment, the link's motion is the load sharing's adjoint (F . u_p = sum F_i . u_i), and
//   a point on a node is that node alone.
// FrameTower.RingsAtSubDynsFirstFrequency: the frame's first natural frequency is SubDyn's.
// FrameTower.SettlesAtTheFramesStaticResponse: a damped tower turned 30 degrees, under steady drag
//   and a line's pull, comes to rest where the frame's static solution of the linked loads puts its
//   cross-arm and drag nodes.
// FrameTower.ItsFootingsCarryTheLoadsAndTheWeight: at rest, the footings' force and moment are the
//   applied loads' (ERF axes, about the base centre), the legs carry the weight and the vertical
//   loads, and the legs on the downwind side are the more compressed.
// FrameTower.BeforeItsFirstStepTheFootingsTakeTheStaticLoads: at rest, the footings take the tower's
//   present drag and line pull exactly, as Tower::foundation() reckons a tower at rest.
// FrameTower.OnSpringFootingsItsLegsCarryTheSpringForces: on vertical soil springs (with a spring
//   mass), a damped tower settled under steady loads has the static solution's footing loads, the
//   springs' force counted in each leg's reaction.
// FrameTower.HeldByAStiffSpanItsSwayDecays: coupled to a spring three times its own stiffness at
//   the cross-arm, iterated each step to the spring's mean pull over the step (the cross-arm moving
//   at a constant velocity) as the coupling with MoorDyn is, the tower's oscillation decays.
// FrameTower.AnInsulatorStringOnItsCrossArmDoesNotRing: a 60 kg string of two segments (1e7 N
//   axial stiffness, MoorDyn's damping ratio 0.5) with 30 kg of conductor at its foot hangs from the
//   cross-arm, integrated finely while the cross-arm moves at a constant velocity over each coupling
//   step. At coupling steps of 0.9 s divided by 76, 80 and 84 (a conductor run's), the string's
//   ring decays when the frame steps under the string's mean pull held over the step
//   (step_between); stepped by Newmark's method to the pull at the step's end instead, it grows at
//   one of these steps at least.
// FrameTower.TheStateRestoresTheSameMotion: a state saved and restored continues bit for bit.
// FrameTower.RefusesAFrameThatDoesNotFitItsTowerOrItsType: a frame whose cross-arm is far from the
//   lines' attachment, or 1.5 m above or 7.5 m below the tower's cross-arm, or with two interface joints,
//   aborts naming the frame file; a type giving a frame and a frequency, a weight or a foundation
//   stiffness is refused.
// FrameTower.ItsMembersAreCheckedUnderTheFramesForces: on a generated frame with design data, the
//   members' checks before the first step and once settled under steady loads are those of the
//   frame's static solution under the linked loads and the weight; without design data, none.
// FrameTower.HeatingKeepsItsPlaceThenSettlesOnTheSofterFrame: heated to 600 C while settled, the
//   cross-arm stays put, the first frequency drops by sqrt(k_E), and the tower settles where the hot
//   frame's static solution puts it, its checks those of the hot steel; bad temperatures are refused.
// FrameTower.ItsStateCarriesTheTemperatures: a cold tower restored from a heated one's state heats up
//   and continues bit for bit; a state without temperatures keeps the present ones; a restart
//   (restart_mismatch) refuses temperatures other than the inputs'.
// FrameTower.AStateFromBeforeTheFirstStepKeepsTheStaticFootings: a step-0 state restores the static
//   footings under the present loads.
// FrameTower.ItsCrossArmDisplacementIsItsCentresWhenItTwists: twisted by opposite pulls at the arm's
//   ends, the cross-arm's centre stays put while its ends swing.
// TowerType.RefusesFrameKeysThatDoNotFit: frame_panels with frame_file, frame keys without their frame,
//   angles without two values in order, an unknown bracing, a bad yield strength or temperature.

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_FrameTower.H"
#include "ERF_LatticeFrame.H"
#include "ERF_GTestThrowOnAbort.H"

using amrex::Real;
using namespace erf_towers;
using P3 = std::array<Real,3>;

namespace {

constexpr double g = 9.80665;
// loads reach the tower as Real; the frame works in double
constexpr double tol = (std::is_same<Real, float>::value) ? 2.0e-5 : 1.0e-9;

std::string frame_file () { return std::string(ERF_FRAME_TEST_FILES) + "/caseT/towerT.dat"; }

std::shared_ptr<const Frame> frame_t ()
{
    FrameInputs in;
    const std::string rerr = read_subdyn(frame_file(), in);
    EXPECT_TRUE(rerr.empty()) << rerr;
    std::string err;
    std::shared_ptr<const Frame> f = Frame::create(in, err);
    EXPECT_TRUE(f) << err;
    return f;
}

TowerType framed (double zeta = 0.02)
{
    TowerType t;
    t.name = "lattice";
    t.base_width = 6.0; t.top_width = 1.5; t.solidity = 0.2; t.arm_length = 12.0; t.arm_depth = 1.2;
    t.damping_ratio = static_cast<Real>(zeta);
    t.frame_file = frame_file();
    return t;
}

/** A tower of case T's geometry at (100, 200, 50), its cross-arm turned 30 degrees from y, one line on the cross-arm's centre. */
Tower turned (const TowerType& t, Real arm_height = 30.0)
{
    const double c = std::cos(0.5235987755982988), s = std::sin(0.5235987755982988);
    Tower tw("t1", t, P3{{100.0, 200.0, 50.0}}, arm_height, P3{{static_cast<Real>(-s), static_cast<Real>(c), 0.0}});
    tw.add_attachment(P3{{100.0, 200.0, static_cast<Real>(50.0 + arm_height)}});
    return tw;
}

std::vector<double> subdyn_frequencies ()
{
    std::ifstream f(std::string(ERF_FRAME_TEST_FILES) + "/caseT/towerT_frequencies.txt");
    std::vector<double> v;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream is(line);
        double x = 0.0;
        while (is >> x) { v.push_back(x); }
    }
    return v;
}

/** Steady drag on every drag node (N, ERF axes) and a line's pull at the attachment. */
std::vector<Real> drag (const Tower& tw)
{
    std::vector<Real> f(3 * tw.nodes().size());
    for (std::size_t i = 0; i < tw.nodes().size(); ++i) { f[3 * i] = 150.0; f[3 * i + 1] = -60.0; f[3 * i + 2] = 5.0; }
    return f;
}
const std::array<Real,3> pull{{Real(1.2e4), Real(-3.0e3), Real(-2.0e4)}};

} // namespace

TEST(RigidLink, KeepsTheForceAndItsMomentAndDoesNoWork)
{
    const auto f = frame_t();
    ASSERT_TRUE(f);
    const std::array<double,3> p{{0.3, -0.2, 11.0}}, F{{2.0e3, -7.0e2, 4.0e2}};
    const RigidLink link(*f, p);
    ASSERT_EQ(link.nodes().size(), 4u);
    std::vector<double> loads(f->num_dofs(), 0.0);
    link.add_load(F, loads);
    std::array<double,3> sum{{0, 0, 0}}, mom{{0, 0, 0}};
    for (std::size_t n = 0; n < f->num_nodes(); ++n) {
        const auto& x = f->node_position(n);
        const double* q = &loads[6 * n];
        for (int d = 0; d < 3; ++d) { sum[static_cast<std::size_t>(d)] += q[d]; }
        mom[0] += x[1] * q[2] - x[2] * q[1];
        mom[1] += x[2] * q[0] - x[0] * q[2];
        mom[2] += x[0] * q[1] - x[1] * q[0];
    }
    const std::array<double,3> pm{{p[1] * F[2] - p[2] * F[1], p[2] * F[0] - p[0] * F[2], p[0] * F[1] - p[1] * F[0]}};
    for (std::size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(sum[d], F[d], 1e-9 * 2.0e3);
        EXPECT_NEAR(mom[d], pm[d], 1e-9 * 2.0e4);
    }
    // the motion is the adjoint of the load sharing
    std::vector<double> u(f->num_dofs());
    for (std::size_t i = 0; i < u.size(); ++i) { u[i] = std::sin(0.37 * static_cast<double>(i) + 0.1); }
    const auto up = link.motion(u);
    double work_point = 0.0, work_nodes = 0.0;
    for (std::size_t d = 0; d < 3; ++d) { work_point += F[d] * up[d]; }
    for (std::size_t i = 0; i < u.size(); ++i) { work_nodes += loads[i] * u[i]; }
    EXPECT_NEAR(work_point, work_nodes, 1e-9 * std::abs(work_nodes) + 1e-9);
    // a point on a node is that node
    const RigidLink on(*f, f->node_position(f->node_of_joint(21)));
    ASSERT_EQ(on.nodes().size(), 1u);
    EXPECT_EQ(on.nodes()[0], f->node_of_joint(21));
    const auto uo = on.motion(u);
    for (std::size_t d = 0; d < 3; ++d) { EXPECT_DOUBLE_EQ(uo[d], u[6 * f->node_of_joint(21) + d]); }
}

TEST(FrameTower, RingsAtSubDynsFirstFrequency)
{
    const Tower tw = turned(framed());
    const FrameTower ft(tw, frame_t(), g);
    const std::vector<double> theirs = subdyn_frequencies();
    ASSERT_FALSE(theirs.empty());
    EXPECT_NEAR(static_cast<double>(ft.frequency()), theirs[0], 1e-6 * theirs[0]);
    EXPECT_NEAR(ft.mass(), 6800.110, 0.01);   // SubDyn's total mass of case T (kg)
}

TEST(FrameTower, SettlesAtTheFramesStaticResponse)
{
    const Tower tw = turned(framed(0.3));
    const auto frame = frame_t();
    FrameTower ft(tw, frame, g);
    const std::vector<Real> fd = drag(tw);
    const double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 1500; ++n) { ft.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    // the frame's static solution of the same linked loads, in tower-local axes
    const double c = std::cos(0.5235987755982988), s = std::sin(0.5235987755982988);
    const std::array<double,3> along{{c, s, 0.0}}, across{{-s, c, 0.0}};
    auto local = [&] (double x, double y, double z) { return std::array<double,3>{{x * along[0] + y * along[1], x * across[0] + y * across[1], z}}; };
    std::vector<double> loads(frame->num_dofs(), 0.0);
    for (std::size_t i = 0; i < tw.nodes().size(); ++i) {
        const auto& p = tw.nodes()[i].pos;
        RigidLink(*frame, local(p[0] - 100.0, p[1] - 200.0, p[2] - 50.0)).add_load(local(fd[3 * i], fd[3 * i + 1], fd[3 * i + 2]), loads);
    }
    const RigidLink arm(*frame, {{0.0, 0.0, 30.0}});
    arm.add_load(local(pull[0], pull[1], pull[2]), loads);
    const FrameSolution st = frame->solve(loads, 0.0);
    const auto ul = arm.motion(st.displacement);
    const std::array<double,3> expected{{ul[0] * along[0] + ul[1] * across[0], ul[0] * along[1] + ul[1] * across[1], ul[2]}};
    const auto got = ft.attachment_displacement(0);
    const double scale = std::sqrt(expected[0] * expected[0] + expected[1] * expected[1] + expected[2] * expected[2]);
    ASSERT_GT(scale, 1.0e-4) << "the cross-arm must have moved: the check is not vacuous";
    for (std::size_t d = 0; d < 3; ++d) { EXPECT_NEAR(static_cast<double>(got[d]), expected[d], 1e-6 * scale + tol * scale) << d; }
    // and at rest: the velocity is a small part of the sway's amplitude times its angular frequency
    const double vscale = scale * 2.0 * 3.14159265358979323846 * static_cast<double>(ft.frequency());
    for (std::size_t d = 0; d < 3; ++d) { EXPECT_NEAR(static_cast<double>(ft.attachment_velocity(0)[d]), 0.0, 1e-4 * vscale); }
}

TEST(FrameTower, ItsFootingsCarryTheLoadsAndTheWeight)
{
    const Tower tw = turned(framed(0.3));
    FrameTower ft(tw, frame_t(), g);
    const std::vector<Real> fd = drag(tw);
    const double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 1500; ++n) { ft.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    FoundationLoad L;
    ASSERT_TRUE(ft.foundation(tw, L));
    // the applied loads and their moment about the base centre
    std::array<double,3> F{{0, 0, 0}}, M{{0, 0, 0}};
    auto add = [&] (const std::array<double,3>& r, const std::array<double,3>& f) {
        for (std::size_t d = 0; d < 3; ++d) { F[d] += f[d]; }
        M[0] += r[1] * f[2] - r[2] * f[1];
        M[1] += r[2] * f[0] - r[0] * f[2];
        M[2] += r[0] * f[1] - r[1] * f[0];
    };
    for (std::size_t i = 0; i < tw.nodes().size(); ++i) {
        const auto& p = tw.nodes()[i].pos;
        add({{p[0] - 100.0, p[1] - 200.0, p[2] - 50.0}}, {{fd[3 * i], fd[3 * i + 1], fd[3 * i + 2]}});
    }
    add({{0.0, 0.0, 30.0}}, {{pull[0], pull[1], pull[2]}});
    for (std::size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(static_cast<double>(L.force[d]), F[d], 1e-5 * 2.0e4) << d;
        EXPECT_NEAR(static_cast<double>(L.moment[d]), M[d], 1e-5 * 6.0e5) << d;
    }
    const double weight = ft.mass() * g;
    EXPECT_NEAR(static_cast<double>(L.vertical), weight - F[2], 1e-5 * weight);
    // the load leans along +along (cos 30 * 150 - sin 30 * 60 > 0 per node, and the pull): legs 0 and 1 are downwind
    EXPECT_GT(std::min(L.legs[0], L.legs[1]), std::max(L.legs[2], L.legs[3]));
}

TEST(FrameTower, BeforeItsFirstStepTheFootingsTakeTheStaticLoads)
{
    Tower tw = turned(framed());
    const FrameTower ft(tw, frame_t(), g);
    const std::vector<Real> fd = drag(tw);
    tw.set_loads(fd);
    tw.set_line_loads({pull}, {tw.attachments()[0]});
    FoundationLoad L;
    ASSERT_TRUE(ft.foundation(tw, L));
    std::array<double,3> F{{static_cast<double>(pull[0]), static_cast<double>(pull[1]), static_cast<double>(pull[2])}};
    std::array<double,3> M{{-30.0 * static_cast<double>(pull[1]), 30.0 * static_cast<double>(pull[0]), 0.0}};
    for (std::size_t i = 0; i < tw.nodes().size(); ++i) {
        const auto& p = tw.nodes()[i].pos;
        const std::array<double,3> r{{p[0] - 100.0, p[1] - 200.0, p[2] - 50.0}}, f{{fd[3 * i], fd[3 * i + 1], fd[3 * i + 2]}};
        for (std::size_t d = 0; d < 3; ++d) { F[d] += f[d]; }
        M[0] += r[1] * f[2] - r[2] * f[1];
        M[1] += r[2] * f[0] - r[0] * f[2];
        M[2] += r[0] * f[1] - r[1] * f[0];
    }
    for (std::size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(static_cast<double>(L.force[d]), F[d], tol * 2.0e4 + 1e-9 * 2.0e4) << d;
        EXPECT_NEAR(static_cast<double>(L.moment[d]), M[d], tol * 6.0e5 + 1e-9 * 6.0e5) << d;
    }
    EXPECT_NEAR(static_cast<double>(L.vertical), ft.mass() * g - F[2], 1e-6 * ft.mass() * g);
}

TEST(FrameTower, OnSpringFootingsItsLegsCarryTheSpringForces)
{
    // case T with every leg free to move vertically on a soil spring (Kzz, and a spring mass Mzz)
    FrameInputs in;
    const std::string rerr = read_subdyn(frame_file(), in);
    ASSERT_TRUE(rerr.empty()) << rerr;
    for (auto& sp : in.supports) {
        sp.fixed[2] = false;
        sp.stiffness[5] = 2.0e8;   // Kzz (N/m), the 6th of SubDyn's upper-triangle order
        sp.mass[5] = 500.0;        // Mzz (kg)
    }
    std::string err;
    const std::shared_ptr<const Frame> f = Frame::create(in, err);
    ASSERT_TRUE(f) << err;
    Tower tw = turned(framed(0.3));
    FrameTower stepped(tw, f, g);
    const FrameTower at_rest(tw, f, g);
    const std::vector<Real> fd = drag(tw);
    const double h = 1.0 / (20.0 * static_cast<double>(stepped.frequency()));
    for (int n = 0; n < 1500; ++n) { stepped.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    tw.set_loads(fd);
    tw.set_line_loads({pull}, {tw.attachments()[0]});
    FoundationLoad dyn, stat;
    ASSERT_TRUE(stepped.foundation(tw, dyn));
    ASSERT_TRUE(at_rest.foundation(tw, stat));
    // the legs' loads differ by about 30 kN from leg to leg here; the springs alone would lose that
    EXPECT_GT(std::max({stat.legs[0], stat.legs[1], stat.legs[2], stat.legs[3]}) -
              std::min({stat.legs[0], stat.legs[1], stat.legs[2], stat.legs[3]}), Real(1.0e4));
    for (std::size_t l = 0; l < 4; ++l) { EXPECT_NEAR(dyn.legs[l], stat.legs[l], 1e-5 * 6.0e4) << "leg " << l; }
    for (std::size_t d = 0; d < 3; ++d) {
        EXPECT_NEAR(dyn.force[d], stat.force[d], 1e-5 * 2.0e4) << d;
        EXPECT_NEAR(dyn.moment[d], stat.moment[d], 1e-5 * 6.0e5) << d;
    }
}

TEST(FrameTower, HeldByAStiffSpanItsSwayDecays)
{
    const Tower tw = turned(framed());
    const std::vector<Real> none(3 * tw.nodes().size(), Real(0.0));
    // the cross-arm's stiffness along x: a steady pull, settled on a damped copy
    double K = 0.0;
    {
        const Tower td = turned(framed(0.5));
        FrameTower settle(td, frame_t(), g);
        const double h = 1.0 / (20.0 * static_cast<double>(settle.frequency()));
        const std::array<Real,3> P{{Real(1.0e4), Real(0.0), Real(0.0)}};
        for (int n = 0; n < 1500; ++n) { settle.step(static_cast<Real>(h), none, std::vector<std::array<Real,3>>{P}); }
        K = 1.0e4 / static_cast<double>(settle.attachment_displacement()[0]);
    }
    ASSERT_GT(K, 0.0);
    // a span three times as stiff pulls the cross-arm toward x0; its mean pull over each coupling step,
    // the cross-arm moving at a constant velocity, is found by fixed-point iteration from the pull at
    // the start, the tower's state restored each time
    const double k = 3.0 * K, x0 = 0.01;
    FrameTower ft(tw, frame_t(), g);
    const double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    std::vector<double> x;
    for (int n = 0; n < 600; ++n) {
        const std::vector<double> s0 = ft.state();
        const double xs = static_cast<double>(ft.attachment_displacement()[0]);
        const std::array<Real,3> F0{{static_cast<Real>(-k * (xs - x0)), Real(0.0), Real(0.0)}};
        std::array<Real,3> F = F0;
        for (int it = 0; it < 50; ++it) {
            ASSERT_TRUE(ft.set_state(s0));
            ft.step_between(static_cast<Real>(h), none, {F}, {F});
            const double xe = static_cast<double>(ft.attachment_displacement()[0]);
            const std::array<Real,3> Fn{{static_cast<Real>(-k * (0.5 * (xs + xe) - x0)), Real(0.0), Real(0.0)}};
            const bool done = std::abs(static_cast<double>(Fn[0] - F[0])) <= 1e-6 * k * x0;
            F = Fn;
            if (done) { break; }
        }
        x.push_back(static_cast<double>(ft.attachment_displacement()[0]));
    }
    // the swing about the equilibrium k x0 / (k + K) over the first and the last 100 steps
    const double xe = k * x0 / (k + K);
    double first = 0.0, last = 0.0;
    for (std::size_t n = 0; n < 100; ++n) {
        first = std::max(first, std::abs(x[n] - xe));
        last = std::max(last, std::abs(x[x.size() - 1 - n] - xe));
    }
    EXPECT_LT(last, 0.5 * first) << "first " << first << " m, last " << last << " m";
}

TEST(FrameTower, AnInsulatorStringOnItsCrossArmDoesNotRing)
{
    const Tower tw = turned(framed());
    const std::vector<Real> none(3 * tw.nodes().size(), Real(0.0));
    // the string's two segments (MoorDyn's lumping: half of each segment's mass at its ends), the conductor at its
    // foot, displacements along z from equilibrium; c = zeta sqrt(EA m_l) per segment (MoorDyn's -zeta)
    const double EA = 1.0e7, L = 2.5, mass = 60.0;
    const double k = EA / (0.5 * L), c = 0.5 * std::sqrt(EA * mass / L);
    const double m1 = 0.5 * mass, m2 = 0.25 * mass + 30.0;
    struct String { double y1{0.0}, v1{0.0}, y2{0.0}, v2{0.0}; };
    // the string's pull on the cross-arm (up) with the cross-arm at z u moving at ud
    auto pull = [&] (const String& s, double u, double ud) { return k * (s.y1 - u) + c * (s.v1 - ud); };
    auto rates = [&] (const String& s, double u, double ud) {
        const double f1 = -pull(s, u, ud) + k * (s.y2 - s.y1) + c * (s.v2 - s.v1);
        const double f2 = -k * (s.y2 - s.y1) - c * (s.v2 - s.v1);
        return std::array<double,4>{{s.v1, f1 / m1, s.v2, f2 / m2}};
    };
    // the string over a coupling step, the cross-arm moving from u0 at a constant velocity ud: the pull at the
    // end and its mean over the step (trapezoidal rule over the fine steps), by fourth-order Runge-Kutta
    auto advance = [&] (String& s, double u0, double ud, double h, double& end, double& mean) {
        const int nf = 200;
        const double d = h / nf;
        mean = 0.5 * pull(s, u0, ud);
        for (int i = 0; i < nf; ++i) {
            const double t = i * d;
            auto at = [&] (const String& b, const std::array<double,4>& r, double w) {
                return String{b.y1 + w * r[0], b.v1 + w * r[1], b.y2 + w * r[2], b.v2 + w * r[3]};
            };
            const auto r1 = rates(s, u0 + ud * t, ud);
            const auto r2 = rates(at(s, r1, 0.5 * d), u0 + ud * (t + 0.5 * d), ud);
            const auto r3 = rates(at(s, r2, 0.5 * d), u0 + ud * (t + 0.5 * d), ud);
            const auto r4 = rates(at(s, r3, d), u0 + ud * (t + d), ud);
            s.y1 += d / 6.0 * (r1[0] + 2.0 * r2[0] + 2.0 * r3[0] + r4[0]);
            s.v1 += d / 6.0 * (r1[1] + 2.0 * r2[1] + 2.0 * r3[1] + r4[1]);
            s.y2 += d / 6.0 * (r1[2] + 2.0 * r2[2] + 2.0 * r3[2] + r4[2]);
            s.v2 += d / 6.0 * (r1[3] + 2.0 * r2[3] + 2.0 * r3[3] + r4[3]);
            mean += (i + 1 < nf ? 1.0 : 0.5) * pull(s, u0 + ud * (t + d), ud);
        }
        mean /= nf;
        end = pull(s, u0 + ud * h, ud);
    };
    // the string's ring (the largest swing of the cross-arm's z) over the first and the last 50 coupling steps
    // of 400 after the conductor is pulled 1 mm down, held_mean the coupling with MoorDyn, else Newmark's end pull
    auto ring = [&] (double h, bool held_mean, double& first, double& last) {
        FrameTower ft(tw, frame_t(), g);
        String s;
        s.y2 = -1.0e-3;
        double fe = pull(s, 0.0, 0.0);
        std::vector<double> z;
        for (int n = 0; n < 400; ++n) {
            const std::vector<double> s0 = ft.state();
            const String str0 = s;
            const double u0 = static_cast<double>(ft.attachment_displacement()[2]);
            double F = fe, end = 0.0, mean = 0.0;
            for (int it = 0; it < 60; ++it) {
                if (!ft.set_state(s0)) { ADD_FAILURE() << "state"; return; }
                s = str0;
                const std::array<Real,3> P{{Real(0.0), Real(0.0), static_cast<Real>(F)}};
                if (held_mean) { ft.step_between(static_cast<Real>(h), none, {P}, {P}); }
                else { ft.step(static_cast<Real>(h), none, std::vector<std::array<Real,3>>{P}); }
                const double u1 = static_cast<double>(ft.attachment_displacement()[2]);
                advance(s, u0, (u1 - u0) / h, h, end, mean);
                const double Fn = held_mean ? mean : end;
                const bool done = std::abs(Fn - F) <= 1e-10 * k * 1.0e-3;
                F = Fn;
                if (done) { break; }
            }
            fe = end;
            z.push_back(static_cast<double>(ft.attachment_displacement()[2]));
        }
        double mz = 0.0;
        for (const double v : z) { mz += v / static_cast<double>(z.size()); }
        first = last = 0.0;
        for (std::size_t n = 0; n < 50; ++n) {
            first = std::max(first, std::abs(z[n] - mz));
            last = std::max(last, std::abs(z[z.size() - 1 - n] - mz));
        }
    };
    double grew = 0.0;
    for (const int steps : {76, 80, 84}) {
        const double h = 0.9 / steps;
        double first = 0.0, last = 0.0;
        ring(h, true, first, last);
        EXPECT_LT(last, 0.5 * first) << steps << " coupling steps per 0.9 s, mean pull held: first " << first << " m, last " << last << " m";
        ring(h, false, first, last);
        grew = std::max(grew, last / first);
    }
    // the coupling with Newmark's end pull grows at one of these steps at least: the check is not vacuous
    EXPECT_GT(grew, 2.0);
}

TEST(FrameTower, TheStateRestoresTheSameMotion)
{
    const Tower tw = turned(framed());
    const auto frame = frame_t();
    FrameTower a(tw, frame, g), b(tw, frame, g);
    const std::vector<Real> fd = drag(tw);
    for (int n = 0; n < 10; ++n) { a.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull}); }
    ASSERT_TRUE(b.set_state(a.state()));
    for (int n = 0; n < 10; ++n) {
        a.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull});
        b.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull});
    }
    EXPECT_EQ(a.state(), b.state());
    std::vector<double> bad = a.state();
    bad.pop_back();
    EXPECT_FALSE(b.set_state(bad));
}

TEST(FrameTower, RefusesAFrameThatDoesNotFitItsTowerOrItsType)
{
    // the line hangs 60 m up, the frame's cross-arm is at 30 m
    const Tower tall = turned(framed(), Real(60.0));
    const std::string msg = erf_gtest::abort_message([&] { FrameTower ft(tall, frame_t(), g); });
    EXPECT_NE(msg.find("frame_file"), std::string::npos) << msg;
    EXPECT_NE(msg.find("from the nearest node of"), std::string::npos) << msg;
    // 1.5 m above the frame's cross-arm (and its top): close enough to tie, but built for another height
    const Tower raised = turned(framed(), Real(31.5));
    const std::string higher = erf_gtest::abort_message([&] { FrameTower ft(raised, frame_t(), g); });
    EXPECT_NE(higher.find("the frame's interface joint must be its cross-arm's centre"), std::string::npos) << higher;
    // 22.5 m, a panel joint's height: close to frame nodes, but not to the frame's cross-arm
    const Tower panel = turned(framed(), Real(22.5));
    const std::string lower = erf_gtest::abort_message([&] { FrameTower ft(panel, frame_t(), g); });
    EXPECT_NE(lower.find("a frame is built for one cross-arm height"), std::string::npos) << lower;
    {
        // a second interface joint: the lines' pulls are tied at one, the cross-arm's centre
        FrameInputs in;
        ASSERT_TRUE(read_subdyn(frame_file(), in).empty());
        in.interface_joints.push_back(in.joints.front().id);
        std::string err;
        const std::shared_ptr<const Frame> two = Frame::create(in, err);
        ASSERT_TRUE(two) << err;
        const std::string twice = erf_gtest::abort_message([&] { FrameTower ft(turned(framed(), Real(30.0)), two, g); });
        EXPECT_NE(twice.find("has 2 interface joints"), std::string::npos) << twice;
    }
    TowerType t = framed();
    EXPECT_TRUE(t.validate().empty()) << t.validate();
    t.frequency = 2.0;
    EXPECT_NE(t.validate().find("frequency is not given with"), std::string::npos) << t.validate();
    t = framed();
    t.weight = 6.0e4;
    EXPECT_NE(t.validate().find("weight is not given with"), std::string::npos) << t.validate();
    t = framed();
    t.foundation_rotational_stiffness = 1.0e9;
    EXPECT_NE(t.validate().find("foundation_rotational_stiffness is not given with"), std::string::npos) << t.validate();
}

namespace {

/** The generated case G (a 30 m tower with a 3 m peak) and the type that stands it. */
TowerType generated_type (double zeta = 0.02)
{
    TowerType t;
    t.name = "generated";
    t.base_width = 6.0; t.top_width = 1.5; t.solidity = 0.2; t.arm_length = 12.0; t.arm_depth = 1.5; t.peak = 3.0;
    t.damping_ratio = static_cast<Real>(zeta);
    t.frame_panels = 6;
    t.leg_angle = {Real(0.15), Real(0.012)};
    t.brace_angle = {Real(0.09), Real(0.007)};
    return t;
}

struct Generated
{
    std::shared_ptr<const Frame> frame;
    std::shared_ptr<const std::vector<MemberDesign>> designs;
    std::vector<std::size_t> links;     //!< the nodes loads may be tied to: all but the diagonals' crossings
};

Generated generated (double theta = 20.0)
{
    LatticeSpec s;
    s.base_width = 6.0; s.top_width = 1.5; s.arm_height = 30.0; s.arm_length = 12.0; s.arm_depth = 1.5; s.peak = 3.0;
    s.panels = 6; s.leg_b = 0.15; s.leg_t = 0.012; s.brace_b = 0.09; s.brace_t = 0.007;
    FrameInputs in;
    std::vector<MemberDesign> d;
    std::vector<int> load_joints;
    const std::string gerr = lattice_frame(s, in, d, &load_joints);
    EXPECT_TRUE(gerr.empty()) << gerr;
    if (theta != 20.0) { in.temperature.assign(in.members.size(), theta); }
    std::string err;
    Generated gen;
    gen.frame = Frame::create(in, err);
    EXPECT_TRUE(gen.frame) << err;
    gen.designs = std::make_shared<const std::vector<MemberDesign>>(d);
    for (const int id : load_joints) { gen.links.push_back(gen.frame->node_of_joint(id)); }
    return gen;
}

/**
 * The frame loads of the tower's drag fd and the line's pull at the cross-arm's centre, as FrameTower links them
 * to the nodes links (tower-local axes).
 */
std::vector<double> linked_loads (const Frame& frame, const Tower& tw, const std::vector<Real>& fd, const std::vector<std::size_t>& links)
{
    const double c = std::cos(0.5235987755982988), s = std::sin(0.5235987755982988);
    const std::array<double,3> along{{c, s, 0.0}}, across{{-s, c, 0.0}};
    auto local = [&] (double x, double y, double z) {
        return std::array<double,3>{{x * along[0] + y * along[1], x * across[0] + y * across[1], z}};
    };
    std::vector<double> loads(frame.num_dofs(), 0.0);
    for (std::size_t i = 0; i < tw.nodes().size(); ++i) {
        const auto& p = tw.nodes()[i].pos;
        RigidLink(frame, local(p[0] - 100.0, p[1] - 200.0, p[2] - 50.0), 4, &links)
            .add_load(local(fd[3 * i], fd[3 * i + 1], fd[3 * i + 2]), loads);
    }
    RigidLink(frame, {{0.0, 0.0, 30.0}}, 4, &links).add_load(local(pull[0], pull[1], pull[2]), loads);
    return loads;
}

} // namespace

TEST(FrameTower, ItsMembersAreCheckedUnderTheFramesForces)
{
    Tower tw = turned(generated_type(0.3));
    const Generated gen = generated();
    FrameTower ft(tw, gen.frame, g, gen.designs, "case G", gen.links);
    ASSERT_TRUE(ft.has_member_checks());
    const std::vector<Real> fd = drag(tw);
    // before the first step: the static forces under the tower's present loads and the weight
    tw.set_loads(fd);
    tw.set_line_loads({pull}, {tw.attachments()[0]});
    const std::vector<double> loads = linked_loads(*gen.frame, tw, fd, gen.links);
    const FrameSolution st = gen.frame->solve(loads, g);
    const auto expected = check_members(*gen.frame, *gen.designs, st.element_force, {});
    auto at_rest = ft.member_checks(tw);
    ASSERT_EQ(at_rest.size(), expected.size());
    double umax = 0.0;
    for (std::size_t m = 0; m < expected.size(); ++m) {
        EXPECT_NEAR(at_rest[m].utilisation, expected[m].utilisation, tol * (expected[m].utilisation + 1e-3)) << "member " << m + 1;
        umax = std::max(umax, expected[m].utilisation);
    }
    ASSERT_GT(umax, 0.01) << "the loads must work the members: the check is not vacuous";
    // settled under the same loads, the dynamic forces about the weight's equilibrium add up to the same
    const double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 1500; ++n) { ft.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    const auto settled = ft.member_checks(tw);
    for (std::size_t m = 0; m < expected.size(); ++m) {
        EXPECT_NEAR(settled[m].utilisation, expected[m].utilisation, 1e-5 * umax + tol * umax) << "member " << m + 1;
    }
    // a frame without design data checks nothing
    const FrameTower plain(tw, gen.frame, g);
    EXPECT_FALSE(plain.has_member_checks());
    EXPECT_TRUE(plain.member_checks(tw).empty());
}

TEST(FrameTower, HeatingKeepsItsPlaceThenSettlesOnTheSofterFrame)
{
    const Tower tw = turned(generated_type(0.3));
    const Generated cold = generated(), hot = generated(600.0);
    FrameTower ft(tw, cold.frame, g, cold.designs, "case G", cold.links);
    const std::vector<Real> fd = drag(tw);
    double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 1500; ++n) { ft.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    const auto before = ft.attachment_displacement(0);
    const double f_cold = static_cast<double>(ft.frequency());
    const std::string err = ft.set_temperature(std::vector<double>(cold.frame->inputs().members.size(), 600.0));
    ASSERT_TRUE(err.empty()) << err;
    // the same place at the switch
    const auto after = ft.attachment_displacement(0);
    for (std::size_t d = 0; d < 3; ++d) { EXPECT_NEAR(static_cast<double>(after[d]), static_cast<double>(before[d]), tol * 1.0e-2) << d; }
    // every member softened by k_E = 0.31: the frequency by its square root
    EXPECT_NEAR(static_cast<double>(ft.frequency()), f_cold * std::sqrt(0.31), 1e-6 * f_cold);
    // the Rayleigh damping is the cold frame's, lighter on some of the hot frame's modes: settle twice as long
    h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 3000; ++n) { ft.step(static_cast<Real>(h), fd, std::vector<std::array<Real,3>>{pull}); }
    // settled: the hot frame's static response to the loads and its weight, measured from the cold frame's sag
    const std::vector<double> loads = linked_loads(*cold.frame, tw, fd, cold.links);
    const FrameSolution s_hot = hot.frame->solve(loads, g), sag_cold = cold.frame->solve({}, g);
    const RigidLink arm(*hot.frame, {{0.0, 0.0, 30.0}}, 4, &cold.links);
    const auto u_hot = arm.motion(s_hot.displacement), u_sag = arm.motion(sag_cold.displacement);
    const double c = std::cos(0.5235987755982988), s = std::sin(0.5235987755982988);
    const std::array<double,3> ul{{u_hot[0] - u_sag[0], u_hot[1] - u_sag[1], u_hot[2] - u_sag[2]}};
    const std::array<double,3> expected{{ul[0] * c - ul[1] * s, ul[0] * s + ul[1] * c, ul[2]}};
    const auto got = ft.attachment_displacement(0);
    const double scale = std::sqrt(expected[0] * expected[0] + expected[1] * expected[1] + expected[2] * expected[2]);
    for (std::size_t d = 0; d < 3; ++d) { EXPECT_NEAR(static_cast<double>(got[d]), expected[d], 1e-5 * scale + tol * scale) << d; }
    // and the members' checks see the hot steel: the same forces, 0.47 of the tension strength
    const auto checks = ft.member_checks(tw);
    const auto expected_checks = check_members(*hot.frame, *hot.designs, s_hot.element_force, hot.frame->inputs().temperature);
    double umax = 0.0;
    for (const auto& ck : expected_checks) { umax = std::max(umax, ck.utilisation); }
    for (std::size_t m = 0; m < checks.size(); ++m) {
        EXPECT_NEAR(checks[m].utilisation, expected_checks[m].utilisation, 1e-5 * umax + tol * umax) << m + 1;
    }
    // a temperature for each member is needed, below 1200 C; nothing changes otherwise
    EXPECT_NE(ft.set_temperature({600.0}).find("temperatures for"), std::string::npos);
    EXPECT_NE(ft.set_temperature(std::vector<double>(checks.size(), 1250.0)).find("below 1200"), std::string::npos);
    EXPECT_EQ(ft.temperature().front(), 600.0);
}

TEST(FrameTower, ItsStateCarriesTheTemperatures)
{
    const Tower tw = turned(generated_type());
    const Generated gen = generated();
    FrameTower a(tw, gen.frame, g, gen.designs, "case G", gen.links), b(tw, gen.frame, g, gen.designs, "case G", gen.links);
    const std::vector<Real> fd = drag(tw);
    for (int n = 0; n < 10; ++n) { a.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull}); }
    std::vector<double> theta(gen.frame->inputs().members.size(), 20.0);
    for (std::size_t m = 0; m < theta.size(); ++m) { theta[m] = 20.0 + 5.0 * static_cast<double>(m % 100); }
    ASSERT_TRUE(a.set_temperature(theta).empty());
    for (int n = 0; n < 10; ++n) { a.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull}); }
    // a cold tower restored from the hot one's state heats up and continues bit for bit
    ASSERT_TRUE(b.set_state(a.state()));
    EXPECT_EQ(b.temperature(), theta);
    for (int n = 0; n < 10; ++n) {
        a.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull});
        b.step(Real(0.01), fd, std::vector<std::array<Real,3>>{pull});
    }
    EXPECT_EQ(a.state(), b.state());
    EXPECT_EQ(a.attachment_displacement(0), b.attachment_displacement(0));
    // a state without the temperatures keeps the present ones
    std::vector<double> short_state = a.state();
    short_state.resize(short_state.size() - theta.size());
    ASSERT_TRUE(b.set_state(short_state));
    EXPECT_EQ(b.temperature(), theta);
    // a temperature that is not finite is refused
    std::vector<double> bad = a.state();
    bad.back() = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(b.set_state(bad));
    // a restart, though, cannot change the steel's temperature: a cold run refuses the heated checkpoint
    FrameTower cold(tw, gen.frame, g, gen.designs, "case G", gen.links);
    EXPECT_NE(cold.restart_mismatch(a.state()).find("a restart cannot change the steel's temperature"), std::string::npos);
    EXPECT_TRUE(cold.restart_mismatch(cold.state()).empty());
    EXPECT_TRUE(cold.restart_mismatch(short_state).empty()) << "a state without temperatures is not compared";
}

TEST(FrameTower, AStateFromBeforeTheFirstStepKeepsTheStaticFootings)
{
    // a checkpoint at step 0 holds the frame at rest: restored, its footings are the static ones under
    // the present loads, as the run that wrote it had them, not the zero dynamic reactions
    Tower tw = turned(framed());
    const FrameTower first(tw, frame_t(), g);
    FrameTower restored(tw, frame_t(), g);
    ASSERT_TRUE(restored.set_state(first.state()));
    tw.set_loads(drag(tw));
    tw.set_line_loads({pull}, {tw.attachments()[0]});
    FoundationLoad a, b;
    ASSERT_TRUE(first.foundation(tw, a));
    ASSERT_TRUE(restored.foundation(tw, b));
    EXPECT_GT(a.shear, Real(1.0e3));
    for (std::size_t l = 0; l < 4; ++l) { EXPECT_EQ(a.legs[l], b.legs[l]) << "leg " << l; }
    EXPECT_EQ(a.shear, b.shear);
}

TEST(FrameTower, ItsCrossArmDisplacementIsItsCentresWhenItTwists)
{
    // two lines pulling opposite ways along the line at the arm's ends twist the frame: the arm's centre
    // stays nearly put while the arm's ends swing; the cross-arm's displacement in the logs is the centre's
    const double c = std::cos(0.5235987755982988), s = std::sin(0.5235987755982988);
    const P3 across{{static_cast<Real>(-s), static_cast<Real>(c), 0.0}}, along{{static_cast<Real>(c), static_cast<Real>(s), 0.0}};
    Tower tw("t1", framed(0.3), P3{{100.0, 200.0, 50.0}}, Real(30.0), across);
    for (const double side : {-6.0, 6.0}) {
        tw.add_attachment(P3{{static_cast<Real>(100.0 + side * across[0]), static_cast<Real>(200.0 + side * across[1]), Real(80.0)}});
    }
    FrameTower ft(tw, frame_t(), g);
    const std::vector<Real> none(3 * tw.nodes().size(), Real(0.0));
    const Real F = 5.0e3;
    const std::vector<std::array<Real,3>> twist{{{F * along[0], F * along[1], 0.0}}, {{-F * along[0], -F * along[1], 0.0}}};
    const double h = 1.0 / (20.0 * static_cast<double>(ft.frequency()));
    for (int n = 0; n < 1500; ++n) { ft.step(static_cast<Real>(h), none, twist); }
    std::array<Real,3> centre{};
    ASSERT_TRUE(ft.arm_centre_displacement(centre));
    const auto end = ft.attachment_displacement(0);
    const double moved = std::hypot(static_cast<double>(end[0]), static_cast<double>(end[1]));
    EXPECT_GT(moved, 1.0e-4) << "the arm's end must swing for the check to mean anything";
    EXPECT_LT(std::hypot(static_cast<double>(centre[0]), static_cast<double>(centre[1])), 0.05 * moved);
}

TEST(TowerType, RefusesFrameKeysThatDoNotFit)
{
    auto refused = [] (const TowerType& t, const std::string& part) {
        const std::string e = t.validate();
        EXPECT_NE(e.find(part), std::string::npos) << "expected '" << part << "' in: " << e;
    };
    EXPECT_TRUE(generated_type().validate().empty()) << generated_type().validate();
    { TowerType t = generated_type(); t.frame_file = "x.dat"; refused(t, "both give the tower's frame"); }
    { TowerType t = generated_type(); t.weight = 1.0e4; refused(t, "weight is not given with erf.conductors.generated.frame_panels"); }
    { TowerType t = generated_type(); t.leg_angle = {Real(0.1)}; refused(t, "leg_angle needs two values"); }
    { TowerType t = generated_type(); t.brace_angle = {Real(0.01), Real(0.02)}; refused(t, "brace_angle needs two values"); }
    { TowerType t = generated_type(); t.frame_panels = 0; refused(t, "leg_angle needs erf.conductors.generated.frame_panels"); }
    { TowerType t = generated_type(); t.frame_panels = 201; refused(t, "frame_panels must be in [0, 200]"); }
    { TowerType t = generated_type(); t.bracing = "k"; refused(t, "bracing must be crossed or single"); }
    { TowerType t = generated_type(); t.yield_strength = Real(-1.0); refused(t, "yield_strength must be finite"); }
    { TowerType t = generated_type(); t.steel_temperature = Real(1200.0); refused(t, "below 1200"); }
    { TowerType t = generated_type(); t.steel_temperature = std::numeric_limits<Real>::quiet_NaN(); refused(t, "below 1200"); }
    TowerType bare;
    bare.name = "bare";
    bare.base_width = 6.0; bare.top_width = 1.5; bare.solidity = 0.2; bare.arm_length = 12.0;
    ASSERT_TRUE(bare.validate().empty());
    { TowerType t = bare; t.member_file = "m.dat"; refused(t, "member_file needs erf.conductors.bare.frame_file"); }
    { TowerType t = bare; t.bracing = "single"; refused(t, "bracing needs erf.conductors.bare.frame_panels"); }
    { TowerType t = bare; t.yield_strength = Real(3.0e8); refused(t, "yield_strength needs erf.conductors.bare.frame_panels"); }
    { TowerType t = bare; t.steel_temperature = Real(500.0); refused(t, "steel_temperature needs a frame model"); }
    EXPECT_TRUE(generated_type().moves());
    EXPECT_TRUE(generated_type().has_frame());
    EXPECT_FALSE(bare.moves());
    EXPECT_EQ(generated_type().yield(), Real(3.45e8));
}
