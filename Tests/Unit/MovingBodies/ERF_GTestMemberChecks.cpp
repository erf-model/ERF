// The strength checks of a lattice tower's members (ASCE 10-15), the steel's reduction with
// temperature (EN 1993-1-2), the member design file, and the generated lattice tower.
//
// MemberChecks.AnglePropertiesAreTheTwoLegs: an equal-leg angle's area, centroid, second moments,
//   principal moments and torsion constant equal the two rectangles' by hand, and its least radius of
//   gyration is within 3 % of the tabulated L100x100x10 (whose root fillet adds a little area).
// MemberChecks.EffectiveSlendernessFollowsTheSixCurves: legs take L/r; bracing and redundant members
//   take curves 1 to 3 up to L/r = 120 and curves 4 to 6 beyond, and every curve meets at 120.
// MemberChecks.LocalBucklingAndColumnCurvesAreContinuous: the local-buckling stress is Fy up to
//   80/sqrt(Fy ksi) at E = 29000 ksi and joins its next two pieces continuously at both limits; the
//   column curve gives Fy at KL/r = 0 and Fy/2 from both sides of Cc.
// MemberChecks.SteelReductionIsEurocodeTable31: the table's rows, a linear value between them, 1
//   below 20 C and 0 at 1200 C.
// MemberChecks.TensionAndCompressionStrengths: an angle bolted by one leg takes 0.9 Fy on its net
//   area, by both legs Fy; its compression strength is Fa A at its KL/r; the flags for slenderness and
//   w/t; a heated member's yield and modulus scale by k_y and k_E.
// MemberChecks.TheDesignFileRoundTripsAndItsErrorsNameTheLine: written and read back, the rows are
//   the same; a bad role, a short row and a bad keyword name the line; match_designs orders the rows,
//   and refuses a missing member, a stranger, a repeated row and an angle whose area is not its section's.
// MemberChecks.AColumnCarriesItsLoadAxially: a vertical member under a load along it is in
//   compression (or tension) of exactly that load, element by element.
// LatticeFrame.HasTheTowersGeometry: the generated tower's corners, taper, crossings, cross-arm tips,
//   centre joint and member count, for crossed and single bracing; loads are tied to every joint but the
//   crossings, so a link from the shaft's centre takes leg joints.
// LatticeFrame.ItsMembersBalanceTheLoadAtEveryLevel: under a lateral load at the cross-arm, the end
//   forces of the elements cut by a horizontal plane balance the load above it (forces and moments),
//   and the legs carry most of the overturning moment (over three quarters), as a truss.
// LatticeFrame.AFileFramesLoadsAvoidTheCrossingsAsAGeneratedOnes: square_joints() of a frame without roles
//   finds the joints the generator lists for loads, every one but the diagonals' crossings.
// LatticeFrame.AFrameTooLargeForTheDenseSolverIsRefused: 150 panels, over Frame::max_free_dofs, are refused.
// LatticeFrame.TheWrittenSubDynFileReadsBackTheSameFrame: write_subdyn then read_subdyn gives the same
//   joints, members and sections, so the same stiffness.
// LatticeFrame.TheAnglesLieOnTheirPrincipalAxes: every shaft leg's major principal axis (its axis of symmetry)
//   points to its corner and every shaft face brace's lies at 45 degrees to its (tapering) face; on the
//   cross-arms, whose bottom and sides slope to the tips, every chord's points between its two faces and
//   every face member's lies at 45 degrees to its face; without principal_axes every MSpin is 0.
// LatticeFrame.TheStiffnessAtTheCrossArmIsSubDyns: the generated case G has the file's joints (to 1e-12 m; the
//   file was written on another machine) and, as SubDyn reads that file, SubDyn's KBBt at the cross-arm's
//   centre and its lowest natural frequencies.
// HeatedFrame.DeflectsByOneOverKE: a frame at a uniform 600 C moves 1/k_E as far under the same loads.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_Frame.H"
#include "ERF_FrameDynamics.H"
#include "../ERF_GTestTempDir.H"
#include "ERF_FrameTower.H"
#include "ERF_LatticeFrame.H"
#include "ERF_MemberChecks.H"

using namespace erf_towers;

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double ksi = 6.894757e6;   // Pa

/** Case G of Tests/test_files/FrameSubDynTower: a generated 30 m tower like the GeneratedTowers deck's. */
LatticeSpec case_g (bool crossed = true)
{
    LatticeSpec s;
    s.base_width = 6.0; s.top_width = 1.5; s.arm_height = 30.0; s.arm_length = 12.0; s.arm_depth = 1.5; s.peak = 3.0;
    s.panels = 6; s.crossed = crossed;
    s.leg_b = 0.15; s.leg_t = 0.012; s.brace_b = 0.09; s.brace_t = 0.007;
    return s;
}

std::unique_ptr<Frame> make (const FrameInputs& in)
{
    std::string err;
    auto f = Frame::create(in, err);
    EXPECT_TRUE(f) << err;
    return f;
}

std::vector<double> numbers_in (const std::string& path)
{
    std::ifstream f(path);
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

/** The frame's stiffness condensed to a joint: the inverse of its 6 x 6 flexibility. */
std::array<double,36> condensed (const Frame& f, int id)
{
    const std::size_t node = f.node_of_joint(id);
    std::vector<double> flex(36, 0.0);
    for (std::size_t k = 0; k < 6; ++k) {
        std::vector<double> load(f.num_dofs(), 0.0);
        load[6 * node + k] = 1.0;
        const FrameSolution s = f.solve(load, 0.0);
        for (std::size_t i = 0; i < 6; ++i) { flex[6 * i + k] = s.displacement[6 * node + i]; }
    }
    DenseCholesky c;
    EXPECT_LT(c.factor(flex, 6), 0);
    std::array<double,36> k{};
    for (std::size_t j = 0; j < 6; ++j) {
        std::vector<double> e(6, 0.0);
        e[j] = 1.0;
        c.solve(e);
        for (std::size_t i = 0; i < 6; ++i) { k[6 * i + j] = e[i]; }
    }
    return k;
}

} // namespace

TEST(MemberChecks, AnglePropertiesAreTheTwoLegs)
{
    const double b = 0.1, t = 0.01;
    const AngleProperties a = angle_properties(b, t);
    EXPECT_NEAR(a.A, t * (2.0 * b - t), 1e-15);
    const double c = (b * b + b * t - t * t) / (2.0 * (2.0 * b - t));
    EXPECT_NEAR(a.centroid, c, 1e-14);
    const double ig = (t * std::pow(b - c, 3) + b * std::pow(c, 3) - (b - t) * std::pow(c - t, 3)) / 3.0;
    EXPECT_NEAR(a.Ig, ig, 1e-12 * ig);
    // the principal moments straddle Ig equally (equal legs: the principal axes are the diagonals)
    EXPECT_NEAR(a.Iu + a.Iv, 2.0 * a.Ig, 1e-12 * ig);
    // the minor axis through the centroid along the diagonal: each leg as a thin strip at 45 degrees
    double iv = 0.0;
    {
        // rectangles [0,b]x[0,t] and [0,t]x[t,b]; v = (x + y)/sqrt(2) - (c + c)/sqrt(2)
        const int n = 2000;
        const double dx = b / n;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                const double x = (i + 0.5) * dx, y = (j + 0.5) * dx;
                if (!(y < t || x < t)) { continue; }
                const double v = (x + y - 2.0 * c) / std::sqrt(2.0);
                iv += v * v * dx * dx;
            }
        }
    }
    EXPECT_NEAR(a.Iv, iv, 1e-4 * iv);
    EXPECT_NEAR(a.J, (2.0 * b - t) * t * t * t / 3.0, 1e-18);
    EXPECT_NEAR(a.rv, 0.0195, 0.03 * 0.0195);     // L100x100x10 tables: r_v = 19.5 mm
}

TEST(MemberChecks, EffectiveSlendernessFollowsTheSixCurves)
{
    using R = MemberRole;
    using E = EndLoading;
    using X = EndRestraint;
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Leg, E::BothEccentric, X::BothEnds, 90.0), 90.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Bracing, E::Concentric, X::None, 90.0), 90.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Bracing, E::OneEccentric, X::None, 90.0), 30.0 + 0.75 * 90.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Redundant, E::BothEccentric, X::None, 90.0), 60.0 + 0.5 * 90.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Bracing, E::BothEccentric, X::None, 160.0), 160.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Bracing, E::BothEccentric, X::OneEnd, 160.0), 28.6 + 0.762 * 160.0);
    EXPECT_DOUBLE_EQ(effective_slenderness(R::Bracing, E::BothEccentric, X::BothEnds, 160.0), 46.2 + 0.615 * 160.0);
    for (const E e : {E::Concentric, E::OneEccentric, E::BothEccentric}) {
        for (const X x : {X::None, X::OneEnd, X::BothEnds}) {
            EXPECT_NEAR(effective_slenderness(R::Bracing, e, x, 120.0), 120.0, 0.05);
            EXPECT_NEAR(effective_slenderness(R::Bracing, e, x, 120.0 + 1e-9), 120.0, 0.05);
        }
    }
    EXPECT_EQ(slenderness_limit(R::Leg), 150.0);
    EXPECT_EQ(slenderness_limit(R::Bracing), 200.0);
    EXPECT_EQ(slenderness_limit(R::Redundant), 250.0);
}

TEST(MemberChecks, LocalBucklingAndColumnCurvesAreContinuous)
{
    const double E = 29000.0 * ksi, fy = 50.0 * ksi;
    const double lim1 = 80.0 / std::sqrt(50.0), lim2 = 144.0 / std::sqrt(50.0);
    EXPECT_DOUBLE_EQ(local_buckling_stress(lim1 - 1e-9, fy, E), fy);
    EXPECT_NEAR(local_buckling_stress(lim1 + 1e-9, fy, E), fy, 1e-6 * fy);
    const double below = (1.677 - 0.677 * lim2 / lim1) * fy;
    const double above = 0.0332 * pi * pi * E / (lim2 * lim2);
    EXPECT_NEAR(below, above, 2e-3 * fy);   // the standard's rounded constants meet to 0.1 %
    EXPECT_NEAR(local_buckling_stress(lim2 - 1e-9, fy, E), below, 1e-6 * fy);
    EXPECT_NEAR(local_buckling_stress(lim2 + 1e-9, fy, E), above, 1e-6 * fy);
    EXPECT_LT(local_buckling_stress(20.0, fy, E), local_buckling_stress(12.0, fy, E));
    const double cc = pi * std::sqrt(2.0 * E / fy);
    EXPECT_DOUBLE_EQ(compression_stress(0.0, fy, E), fy);
    EXPECT_NEAR(compression_stress(cc - 1e-9, fy, E), 0.5 * fy, 1e-6 * fy);
    EXPECT_NEAR(compression_stress(cc + 1e-9, fy, E), 0.5 * fy, 1e-6 * fy);
    EXPECT_NEAR(compression_stress(200.0, fy, E), pi * pi * E / 4.0e4, 1e-9 * fy);
}

TEST(MemberChecks, SteelReductionIsEurocodeTable31)
{
    double ky = 0.0, kE = 0.0;
    const double rows[][3] = {{20.0, 1.0, 1.0}, {100.0, 1.0, 1.0}, {200.0, 1.0, 0.9}, {400.0, 1.0, 0.7}, {500.0, 0.78, 0.6},
                              {600.0, 0.47, 0.31}, {700.0, 0.23, 0.13}, {800.0, 0.11, 0.09}, {1100.0, 0.02, 0.0225}};
    for (const auto& r : rows) {
        steel_reduction(r[0], ky, kE);
        EXPECT_NEAR(ky, r[1], 1e-15) << r[0];
        EXPECT_NEAR(kE, r[2], 1e-15) << r[0];
    }
    steel_reduction(550.0, ky, kE);
    EXPECT_NEAR(ky, 0.625, 1e-15);
    EXPECT_NEAR(kE, 0.455, 1e-15);
    steel_reduction(-40.0, ky, kE);
    EXPECT_EQ(ky, 1.0);
    EXPECT_EQ(kE, 1.0);
    steel_reduction(1200.0, ky, kE);
    EXPECT_EQ(ky, 0.0);
    EXPECT_EQ(kE, 0.0);
}

TEST(MemberChecks, TensionAndCompressionStrengths)
{
    MemberDesign d;
    d.member = 7;
    d.role = MemberRole::Bracing;
    d.yield = 3.45e8;
    d.b = 0.09; d.t = 0.007;
    d.net_area = 0.85;
    const AngleProperties a = angle_properties(d.b, d.t);
    const double L = 2.0, E = 2.0e11;
    MemberCheck c = check_member(d, L, 0.0, 0.0, E, 20.0, 1.0e5, 0.0);
    EXPECT_NEAR(c.tension_capacity, 0.9 * d.yield * 0.85 * a.A, 1e-9 * c.tension_capacity);
    EXPECT_NEAR(c.utilisation, 1.0e5 / c.tension_capacity, 1e-12);
    d.one_leg = false;
    c = check_member(d, L, 0.0, 0.0, E, 20.0, 1.0e5, 0.0);
    EXPECT_NEAR(c.tension_capacity, d.yield * 0.85 * a.A, 1e-9 * c.tension_capacity);
    // compression: curve 3 (both ends eccentric) at L/r <= 120
    c = check_member(d, L, 0.0, 0.0, E, 20.0, 0.0, 5.0e4);
    const double lr = L / a.rv;
    ASSERT_LT(lr, 120.0);
    EXPECT_NEAR(c.slenderness, lr, 1e-12 * lr);
    EXPECT_NEAR(c.effective, 60.0 + 0.5 * lr, 1e-12 * lr);
    EXPECT_NEAR(c.w_t, (0.09 - 0.007) / 0.007, 1e-12);
    const double fy = local_buckling_stress(c.w_t, d.yield, E);
    EXPECT_NEAR(c.compression_capacity, compression_stress(c.effective, fy, E) * a.A, 1e-9 * c.compression_capacity);
    EXPECT_NEAR(c.utilisation, 5.0e4 / c.compression_capacity, 1e-12);
    // the same member by hand (an independent script, not this code's functions): r_v 17.8006 mm,
    // KL/r 116.178, w/t 11.857 over (w/t)_1 11.311 so F_cr 333.72 MPa, C_c 108.765 below KL/r, so the
    // elastic F_a = pi^2 E/(KL/r)^2 = 146.246 MPa on 1211 mm^2: 177.104 kN
    EXPECT_NEAR(a.rv, 0.0178006, 1e-6);
    EXPECT_NEAR(c.effective, 116.1778, 1e-3);
    EXPECT_NEAR(fy, 3.337199e8, 1e3);
    EXPECT_NEAR(c.compression_capacity, 1.771036e5, 1.0);
    EXPECT_FALSE(c.slender);
    EXPECT_FALSE(c.thin);
    // a long member is slender, a thin one has w/t over 25
    EXPECT_TRUE(check_member(d, 4.5, 0.0, 0.0, E, 20.0, 0.0, 1.0).slender);
    MemberDesign thin = d;
    thin.t = 0.003;
    EXPECT_TRUE(check_member(thin, L, 0.0, 0.0, E, 20.0, 0.0, 1.0).thin);
    // heated to 600 C: yield by 0.47, E by 0.31
    const MemberCheck hot = check_member(d, L, 0.0, 0.0, E, 600.0, 0.0, 5.0e4);
    EXPECT_NEAR(hot.yield, 0.47 * d.yield, 1e-6);
    EXPECT_NEAR(hot.E, 0.31 * E, 1e-3);
    EXPECT_GT(hot.utilisation, c.utilisation);
    // not an angle: the section's area and radius
    MemberDesign pipe = d;
    pipe.b = 0.0; pipe.t = 0.0;
    const MemberCheck p = check_member(pipe, L, 2.0e-3, 0.05, E, 20.0, 1.0e5, 0.0);
    EXPECT_NEAR(p.r, 0.05, 1e-15);
    EXPECT_EQ(p.w_t, 0.0);
    EXPECT_NEAR(p.tension_capacity, d.yield * 0.85 * 2.0e-3, 1e-6);
}

TEST(MemberChecks, TheDesignFileRoundTripsAndItsErrorsNameTheLine)
{
    LatticeSpec s = case_g();
    FrameInputs in;
    std::vector<MemberDesign> designs;
    ASSERT_TRUE(lattice_frame(s, in, designs).empty());
    designs[3].net_area = 0.8;
    designs[5].ends = EndLoading::OneEccentric;
    designs[6].restraint = EndRestraint::BothEnds;
    designs[7].role = MemberRole::Redundant;
    const std::filesystem::path dir = erf_gtest_temp_path("member_checks_designs");
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "designs.dat").string();
    ASSERT_TRUE(write_member_designs(path, designs, "test"));
    std::vector<MemberDesign> back;
    const std::string err = read_member_designs(path, back);
    ASSERT_TRUE(err.empty()) << err;
    ASSERT_EQ(back.size(), designs.size());
    for (std::size_t i = 0; i < back.size(); ++i) {
        EXPECT_EQ(back[i].member, designs[i].member);
        EXPECT_EQ(back[i].role, designs[i].role);
        EXPECT_EQ(back[i].yield, designs[i].yield);
        EXPECT_EQ(back[i].b, designs[i].b);
        EXPECT_EQ(back[i].t, designs[i].t);
        EXPECT_EQ(back[i].net_area, designs[i].net_area);
        EXPECT_EQ(back[i].one_leg, designs[i].one_leg);
        EXPECT_EQ(back[i].ends, designs[i].ends);
        EXPECT_EQ(back[i].restraint, designs[i].restraint);
    }
    // the rows in another order come back in the members' order
    std::vector<MemberDesign> shuffled(back.rbegin(), back.rend());
    ASSERT_TRUE(match_designs(in, shuffled, path).empty());
    for (std::size_t i = 0; i < shuffled.size(); ++i) { EXPECT_EQ(shuffled[i].member, in.members[i].id); }
    auto refused = [&] (std::vector<MemberDesign> d, const std::string& part) {
        const std::string e = match_designs(in, d, path);
        EXPECT_NE(e.find(part), std::string::npos) << "expected '" << part << "' in: " << e;
    };
    { auto d = back; d.pop_back(); refused(d, "has no row"); }
    { auto d = back; d.push_back(back[0]); refused(d, "has two rows"); }
    { auto d = back; d[0].member = 99999; refused(d, "is not a member"); }
    { auto d = back; d[0].b = 0.3; d[0].t = 0.03; refused(d, "more than 10 % apart"); }
    { auto d = back; d[0].net_area = 1.2; refused(d, "net area"); }
    { auto d = back; d[0].yield = -1.0; refused(d, "yield strength"); }
    auto bad_row = [&] (const std::string& row, const std::string& part) {
        { std::ofstream f(path); f << "# header\n" << row << "\n"; }
        std::vector<MemberDesign> x;
        const std::string e = read_member_designs(path, x);
        EXPECT_NE(e.find("line 2"), std::string::npos) << e;
        EXPECT_NE(e.find(part), std::string::npos) << "expected '" << part << "' in: " << e;
    };
    bad_row("1 strut 3.45e8 0.09 0.007 1 one both none", "the role must be");
    bad_row("1 leg 3.45e8 0.09 0.007 1 one both", "a row has 9 values");
    bad_row("1 leg 3.45e8 0.09 0.007 1 three both none", "Bolted must be");
    bad_row("1 leg 3.45e8 0.09 0.007 1 one twice none", "Ends must be");
    bad_row("1 leg 3.45e8 0.09 0.007 1 one both partly", "Restraint must be");
    bad_row("1.5 leg 3.45e8 0.09 0.007 1 one both none", "integer");
    std::filesystem::remove_all(dir);
}

TEST(MemberChecks, AColumnCarriesItsLoadAxially)
{
    FrameInputs in;
    in.file = "column";
    in.divisions = 3;
    in.joints = {FrameJoint{1, {{0.0, 0.0, 0.0}}}, FrameJoint{2, {{0.0, 0.0, 4.0}}}};
    const AngleProperties a = angle_properties(0.1, 0.01);
    FrameSection sec;
    sec.id = 1; sec.E = 2.0e11; sec.G = 7.7e10; sec.rho = 7850.0;
    sec.A = a.A; sec.Asx = a.A; sec.Asy = a.A; sec.Ixx = a.Iu; sec.Iyy = a.Iv; sec.J0 = a.Iu + a.Iv; sec.Jt = a.J;
    in.sections = {sec};
    FrameMember m;
    m.id = 1; m.joint_a = 1; m.joint_b = 2; m.section = 1;
    in.members = {m};
    in.supports = {FrameSupport{}};
    in.supports[0].joint = 1;
    const auto f = make(in);
    ASSERT_TRUE(f);
    MemberDesign d;
    d.member = 1; d.role = MemberRole::Leg; d.yield = 3.45e8; d.b = 0.1; d.t = 0.01; d.one_leg = false;
    for (const double P : {-3.0e4, 2.0e4}) {
        std::vector<double> load(f->num_dofs(), 0.0);
        load[6 * f->node_of_joint(2) + 2] = P;
        const FrameSolution s = f->solve(load, 0.0);
        for (const auto& e : s.element_force) { EXPECT_NEAR(e[8], P, 1e-9 * std::abs(P)); EXPECT_NEAR(-e[2], P, 1e-9 * std::abs(P)); }
        const auto c = check_members(*f, {d}, s.element_force, {});
        ASSERT_EQ(c.size(), 1u);
        EXPECT_NEAR(c[0].length, 4.0, 1e-12);
        EXPECT_NEAR(P < 0.0 ? c[0].compression : c[0].tension, std::abs(P), 1e-9 * std::abs(P));
        EXPECT_EQ(P < 0.0 ? c[0].tension : c[0].compression, 0.0);
        const double cap = (P < 0.0) ? c[0].compression_capacity : c[0].tension_capacity;
        EXPECT_NEAR(c[0].utilisation, std::abs(P) / cap, 1e-12);
    }
}

TEST(LatticeFrame, HasTheTowersGeometry)
{
    for (const bool crossed : {true, false}) {
        const LatticeSpec s = case_g(crossed);
        FrameInputs in;
        std::vector<MemberDesign> designs;
        std::vector<int> load_joints;
        const std::string err = lattice_frame(s, in, designs, &load_joints);
        ASSERT_TRUE(err.empty()) << err;
        // levels: 7 to the cross-arm's bottom, its top, the peak: 9 x 4 corners
        const int levels = s.panels + 3;
        const int faces = 4 * (levels - 1);
        const int arm_panels = static_cast<int>(std::lround(0.5 * (s.arm_length - s.top_width) / s.arm_depth));
        ASSERT_EQ(arm_panels, 4);
        const std::size_t joints = static_cast<std::size_t>(4 * levels + (crossed ? faces : 0) + 1 + 2 * (4 * (arm_panels - 1) + 1));
        EXPECT_EQ(in.joints.size(), joints);
        const std::size_t members = static_cast<std::size_t>(4 * (levels - 1) + 4 * (levels - 1) + (crossed ? 4 : 1) * faces + 8 +
                                                             2 * (4 * arm_panels + 5 * (arm_panels - 1) + 4 * (arm_panels - 1)));
        EXPECT_EQ(in.members.size(), members);
        EXPECT_EQ(designs.size(), members);
        EXPECT_EQ(in.supports.size(), 4u);
        // the base corners, the cross-arm's top corners and the peak's
        for (int c = 0; c < 4; ++c) {
            const auto& base = in.joints[static_cast<std::size_t>(c)].x;
            EXPECT_NEAR(std::abs(base[0]), 3.0, 1e-12);
            EXPECT_NEAR(std::abs(base[1]), 3.0, 1e-12);
            EXPECT_EQ(base[2], 0.0);
            const auto& top = in.joints[static_cast<std::size_t>(4 * (levels - 2) + c)].x;
            EXPECT_NEAR(std::abs(top[0]), 0.75, 1e-12);
            EXPECT_NEAR(top[2], 30.0, 1e-12);
            const auto& peak = in.joints[static_cast<std::size_t>(4 * (levels - 1) + c)].x;
            EXPECT_NEAR(std::abs(peak[1]), 0.75, 1e-12);
            EXPECT_NEAR(peak[2], 33.0, 1e-12);
        }
        // the bottom of the cross-arm at 28.5 m, on the taper
        const auto& bottom = in.joints[static_cast<std::size_t>(4 * (levels - 3))].x;
        EXPECT_NEAR(bottom[2], 28.5, 1e-12);
        EXPECT_NEAR(bottom[0], 0.5 * (6.0 + 28.5 / 30.0 * (1.5 - 6.0)), 1e-12);
        // the centre joint is the interface joint; the tips are the last joints of each cross-arm
        ASSERT_EQ(in.interface_joints.size(), 1u);
        const auto& centre = in.joints[static_cast<std::size_t>(in.interface_joints[0] - 1)].x;
        EXPECT_EQ(centre, (std::array<double,3>{{0.0, 0.0, 30.0}}));
        const auto& tip = in.joints.back().x;
        EXPECT_EQ(tip, (std::array<double,3>{{0.0, -6.0, 30.0}}));
        // loads may be tied to every joint but the crossings, which follow the corners in the joint list
        EXPECT_EQ(load_joints.size(), joints - static_cast<std::size_t>(crossed ? faces : 0));
        for (const int id : load_joints) {
            EXPECT_FALSE(crossed && id > 4 * levels && id <= 4 * levels + faces) << "joint " << id << " is a crossing";
        }
        // a link from the shaft's centre line at a crossing's height takes leg joints, not the crossings
        if (crossed) {
            const auto f = make(in);
            ASSERT_TRUE(f);
            std::vector<std::size_t> links;
            for (const int id : load_joints) { links.push_back(f->node_of_joint(id)); }
            const double zc = in.joints[static_cast<std::size_t>(4 * levels)].x[2];
            const RigidLink link(*f, {{0.0, 0.0, zc}}, 4, &links);
            for (const std::size_t n : link.nodes()) { EXPECT_LT(n, static_cast<std::size_t>(4 * levels)) << "node " << n; }
            const RigidLink any(*f, {{0.0, 0.0, zc}});
            EXPECT_GE(any.nodes()[0], static_cast<std::size_t>(4 * levels)) << "without the list the nearest node is a crossing";
        }
        // every crossing lies on both diagonals of its face
        if (crossed) {
            const auto& p = in.joints[static_cast<std::size_t>(4 * levels)].x;    // panel 0, face 0 -> 1
            const auto& a0 = in.joints[0].x;
            const auto& b1 = in.joints[5].x;
            const auto& a1 = in.joints[1].x;
            const auto& b0 = in.joints[4].x;
            auto off_line = [] (const std::array<double,3>& q, const std::array<double,3>& u, const std::array<double,3>& v) {
                const std::array<double,3> d{{v[0] - u[0], v[1] - u[1], v[2] - u[2]}}, w{{q[0] - u[0], q[1] - u[1], q[2] - u[2]}};
                const std::array<double,3> x{{d[1] * w[2] - d[2] * w[1], d[2] * w[0] - d[0] * w[2], d[0] * w[1] - d[1] * w[0]}};
                return std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
            };
            EXPECT_LT(off_line(p, a0, b1), 1e-12);
            EXPECT_LT(off_line(p, a1, b0), 1e-12);
        }
        // the legs and the cross-arms' chords are legs, the rest bracing
        int legs = 0;
        for (const auto& d : designs) { legs += (d.role == MemberRole::Leg) ? 1 : 0; }
        EXPECT_EQ(legs, 4 * (levels - 1) + 2 * 4 * arm_panels);
        EXPECT_TRUE(make(in));
    }
    // the dimensions are checked
    LatticeSpec bad = case_g();
    bad.arm_depth = 31.0;
    FrameInputs in;
    std::vector<MemberDesign> d;
    EXPECT_NE(lattice_frame(bad, in, d).find("cross-arm's depth"), std::string::npos);
    bad = case_g();
    bad.leg_t = 0.2;
    EXPECT_NE(lattice_frame(bad, in, d).find("legs' angle"), std::string::npos);
    bad = case_g();
    bad.arm_length = 1.0;
    EXPECT_NE(lattice_frame(bad, in, d).find("cross-arm's length"), std::string::npos);
}

TEST(LatticeFrame, ItsMembersBalanceTheLoadAtEveryLevel)
{
    FrameInputs in;
    std::vector<MemberDesign> designs;
    ASSERT_TRUE(lattice_frame(case_g(), in, designs).empty());
    const auto f = make(in);
    ASSERT_TRUE(f);
    // 10 kN along x and 4 kN down at the cross-arm's centre
    const std::size_t centre = f->node_of_joint(in.interface_joints[0]);
    std::vector<double> load(f->num_dofs(), 0.0);
    load[6 * centre] = 1.0e4;
    load[6 * centre + 2] = -4.0e3;
    const FrameSolution s = f->solve(load, 0.0);
    // cut at mid-height of each shaft panel below the cross-arm: the elements across it carry the load above
    for (const double zc : {1.0, 9.0, 21.0, 27.0}) {
        std::array<double,3> F{{0, 0, 0}}, M{{0, 0, 0}};
        double m_legs = 0.0;
        for (std::size_t e = 0; e < f->elements().size(); ++e) {
            const FrameElement& el = f->elements()[e];
            const auto& xa = f->node_position(el.node_a);
            const auto& xb = f->node_position(el.node_b);
            const bool a_below = xa[2] < zc, b_below = xb[2] < zc;
            if (a_below == b_below) { continue; }
            // the force the part above takes from the element: minus the force on the element at its upper node
            const std::size_t up = a_below ? 6 : 0;
            const auto& xu = a_below ? xb : xa;
            const auto& fe = s.element_force[e];
            std::array<double,3> fl{{fe[up], fe[up + 1], fe[up + 2]}}, ml{{fe[up + 3], fe[up + 4], fe[up + 5]}};
            std::array<double,3> fg{}, mg{};
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) { fg[i] += el.dc[3 * i + j] * fl[j]; mg[i] += el.dc[3 * i + j] * ml[j]; }
            }
            for (std::size_t i = 0; i < 3; ++i) { F[i] += fg[i]; }
            const std::array<double,3> r{{xu[0], xu[1], xu[2] - zc}};
            M[0] += r[1] * fg[2] - r[2] * fg[1] + mg[0];
            M[1] += r[2] * fg[0] - r[0] * fg[2] + mg[1];
            M[2] += r[0] * fg[1] - r[1] * fg[0] + mg[2];
            if (designs[el.member].role == MemberRole::Leg) { m_legs += r[2] * fg[0] - r[0] * fg[2]; }
        }
        // the elements' forces on the nodes above balance the load above: they equal the load
        const double moment = 1.0e4 * (30.0 - zc);
        EXPECT_NEAR(F[0], 1.0e4, 1e-6 * 1.0e4) << zc;
        EXPECT_NEAR(F[1], 0.0, 1e-6 * 1.0e4) << zc;
        EXPECT_NEAR(F[2], -4.0e3, 1e-6 * 1.0e4) << zc;
        EXPECT_NEAR(M[1], moment, 1e-6 * moment) << zc;
        EXPECT_NEAR(M[0], 0.0, 1e-6 * moment) << zc;
        // a truss: the legs carry most of the moment (81 to 96 % at these cuts), the diagonals the rest
        EXPECT_GT(m_legs, 0.75 * moment) << zc;
    }
}

TEST(LatticeFrame, AFileFramesLoadsAvoidTheCrossingsAsAGeneratedOnes)
{
    // read from a file, the generated frame has no roles: the joints a leg, chord or strut meets are
    // exactly those the generator lists, every one but the crossings of the diagonals
    FrameInputs in;
    std::vector<MemberDesign> designs;
    std::vector<int> listed;
    ASSERT_TRUE(lattice_frame(case_g(), in, designs, &listed).empty());
    std::vector<int> found = square_joints(in);
    std::sort(listed.begin(), listed.end());
    std::sort(found.begin(), found.end());
    ASSERT_LT(found.size(), in.joints.size()) << "the crossed bracing has crossings to leave out";
    EXPECT_EQ(found, listed);
}

TEST(LatticeFrame, AFrameTooLargeForTheDenseSolverIsRefused)
{
    // about 49 free degrees of freedom a panel: 150 panels are over 7000
    LatticeSpec s = case_g();
    s.panels = 150;
    FrameInputs in;
    std::vector<MemberDesign> designs;
    ASSERT_TRUE(lattice_frame(s, in, designs).empty());
    std::string err;
    EXPECT_FALSE(Frame::create(in, err));
    EXPECT_NE(err.find("more than the 6000 the dense frame solver takes"), std::string::npos) << err;
}

TEST(LatticeFrame, TheWrittenSubDynFileReadsBackTheSameFrame)
{
    FrameInputs in;
    std::vector<MemberDesign> designs;
    ASSERT_TRUE(lattice_frame(case_g(false), in, designs).empty());
    // a spring and a mass at one support, a concentrated mass and a spin, so that every table is written
    in.supports[1].fixed = {{false, false, false, false, false, false}};
    in.supports[1].stiffness[0] = 2.0e8; in.supports[1].stiffness[2] = 2.0e8; in.supports[1].stiffness[5] = 5.0e8;
    in.supports[1].stiffness[9] = 3.0e8; in.supports[1].stiffness[14] = 3.0e8; in.supports[1].stiffness[20] = 1.0e8;
    in.supports[1].mass[0] = 1.0e3;
    FrameMass cm;
    cm.joint = in.interface_joints[0]; cm.mass = 300.0; cm.inertia = {{10.0, 20.0, 30.0, 1.0, 0.5, 0.25}}; cm.offset = {{0.1, 0.0, -1.5}};
    in.masses.push_back(cm);
    in.members[3].spin = 0.3;
    const std::filesystem::path dir = erf_gtest_temp_path("lattice_frame_written");
    std::filesystem::create_directories(dir);
    const std::string path = (dir / "frame.dat").string();
    std::string err = write_subdyn(in, path, "test");
    ASSERT_TRUE(err.empty()) << err;
    FrameInputs back;
    err = read_subdyn(path, back);
    ASSERT_TRUE(err.empty()) << err;
    ASSERT_EQ(back.joints.size(), in.joints.size());
    for (std::size_t j = 0; j < in.joints.size(); ++j) {
        EXPECT_EQ(back.joints[j].id, in.joints[j].id);
        EXPECT_EQ(back.joints[j].x, in.joints[j].x);
    }
    ASSERT_EQ(back.members.size(), in.members.size());
    for (std::size_t m = 0; m < in.members.size(); ++m) {
        EXPECT_EQ(back.members[m].joint_a, in.members[m].joint_a);
        EXPECT_EQ(back.members[m].joint_b, in.members[m].joint_b);
        EXPECT_EQ(back.members[m].section, in.members[m].section);
        EXPECT_NEAR(back.members[m].spin, in.members[m].spin, 1e-15);
    }
    ASSERT_EQ(back.sections.size(), in.sections.size());
    for (std::size_t k = 0; k < in.sections.size(); ++k) {
        EXPECT_EQ(back.sections[k].A, in.sections[k].A);
        EXPECT_EQ(back.sections[k].Ixx, in.sections[k].Ixx);
        EXPECT_EQ(back.sections[k].Iyy, in.sections[k].Iyy);
        EXPECT_EQ(back.sections[k].Jt, in.sections[k].Jt);
    }
    EXPECT_EQ(back.supports[1].stiffness, in.supports[1].stiffness);
    EXPECT_EQ(back.supports[1].mass, in.supports[1].mass);
    EXPECT_EQ(back.supports[1].fixed, in.supports[1].fixed);
    ASSERT_EQ(back.masses.size(), 1u);
    EXPECT_EQ(back.masses[0].inertia, cm.inertia);
    EXPECT_EQ(back.masses[0].offset, cm.offset);
    EXPECT_EQ(back.interface_joints, in.interface_joints);
    std::filesystem::remove_all(dir);
}

TEST(LatticeFrame, TheAnglesLieOnTheirPrincipalAxes)
{
    for (const bool crossed : {true, false}) {
        SCOPED_TRACE(crossed ? "crossed bracing" : "single bracing");
        const LatticeSpec g = case_g(crossed);
        FrameInputs in;
        std::vector<MemberDesign> designs;
        ASSERT_TRUE(lattice_frame(g, in, designs).empty());
        auto at = [&] (int id) { return in.joints[static_cast<std::size_t>(in.joint_index(id))].x; };
        // the shaft's half width at a height, and a vector's part normal to a unit axis, made unit
        auto half = [&] (double z) { return 0.5 * (g.base_width + std::min(z, g.arm_height) / g.arm_height * (g.top_width - g.base_width)); };
        auto normal_to = [] (std::array<double,3> v, const std::array<double,3>& e) {
            const double d = v[0] * e[0] + v[1] * e[1] + v[2] * e[2];
            for (int i = 0; i < 3; ++i) { v[static_cast<std::size_t>(i)] -= d * e[static_cast<std::size_t>(i)]; }
            const double n = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
            for (auto& c : v) { c /= n; }
            return v;
        };
        int legs = 0, peak_legs = 0, braces = 0;
        for (std::size_t m = 0; m < in.members.size(); ++m) {
            const auto& mem = in.members[m];
            const auto a = at(mem.joint_a), b = at(mem.joint_b);
            const std::array<double,9> d = direction_cosines(a, b, mem.spin);
            const std::array<double,3> x{{d[0], d[3], d[6]}}, e{{d[2], d[5], d[8]}};
            // on the shaft or its peak: within its width at both ends (the peak keeps the top width)
            const bool on_shaft = std::abs(a[0]) <= half(a[2]) + 1.0e-9 && std::abs(a[1]) <= half(a[2]) + 1.0e-9 &&
                                  std::abs(b[0]) <= half(b[2]) + 1.0e-9 && std::abs(b[1]) <= half(b[2]) + 1.0e-9;
            if (designs[m].role == MemberRole::Leg && std::abs(e[2]) > 0.9 && on_shaft) {
                // a shaft leg: its axis of symmetry towards its corner
                const auto r = normal_to({{a[0] > 0.0 ? 1.0 : -1.0, a[1] > 0.0 ? 1.0 : -1.0, 0.0}}, e);
                EXPECT_NEAR(std::abs(x[0] * r[0] + x[1] * r[1] + x[2] * r[2]), 1.0, 1.0e-12) << "leg " << mem.id;
                ++((std::min(a[2], b[2]) < g.arm_height - 1.0e-9) ? legs : peak_legs);
                continue;
            }
            // a brace in a shaft face x = +-half(z) or y = +-half(z): its axis of symmetry at 45 degrees to the face
            for (int c = 0; c < 2; ++c) {
                const double sa = a[static_cast<std::size_t>(c)] / half(a[2]), sb = b[static_cast<std::size_t>(c)] / half(b[2]);
                if (!on_shaft || designs[m].role == MemberRole::Leg || std::abs(sa - sb) > 1.0e-9 || std::abs(std::abs(sa) - 1.0) > 1.0e-9) {
                    continue;
                }
                // the face's outward normal, leaning in with the taper below the cross-arm, upright on the peak; a
                // strut takes the face of the panel below it
                const bool strut = std::abs(a[2] - b[2]) < 1.0e-9;
                const bool leans = strut ? a[2] <= g.arm_height + 1.0e-9 : std::min(a[2], b[2]) < g.arm_height - 1.0e-9;
                std::array<double,3> n{{0.0, 0.0, leans ? 0.5 * (g.base_width - g.top_width) / g.arm_height : 0.0}};
                n[static_cast<std::size_t>(c)] = sa;
                const auto np = normal_to(n, e);
                EXPECT_NEAR(std::abs(n[0] * e[0] + n[1] * e[1] + n[2] * e[2]), 0.0, 1.0e-9) << "brace " << mem.id << " lies in its face";
                EXPECT_NEAR(std::abs(x[0] * np[0] + x[1] * np[1] + x[2] * np[2]), std::sqrt(0.5), 1.0e-9) << "brace " << mem.id;
                ++braces;
            }
        }
        EXPECT_EQ(legs, 4 * g.panels + 4) << "every panel's four legs and the cross-arm's";
        EXPECT_GE(peak_legs, 4) << "the peak's";
        EXPECT_GT(braces, (crossed ? 4 : 1) * 4 * g.panels) << "the diagonals (halves when crossed) and the struts";
        // the cross-arms, along y to their tips (0, +-arm_length/2, arm_height), from the shaft's corners at arm_height and
        // arm_depth below it (where the tapering shaft is wider): four plane faces per side through the tip, the top level,
        // the bottom rising to the tip, the sides closing in on it. Their outward normals from the spec: a chord's axis
        // of symmetry between its two faces' normals, a face member's at 45 degrees to its face
        const double H = g.arm_height, hw = 0.5 * g.top_width, hb = half(H - g.arm_depth);
        auto plane = [] (const std::array<double,3>& p, const std::array<double,3>& q, const std::array<double,3>& t,
                         const std::array<double,3>& out) {
            const std::array<double,3> u{{q[0] - p[0], q[1] - p[1], q[2] - p[2]}}, v{{t[0] - p[0], t[1] - p[1], t[2] - p[2]}};
            std::array<double,3> n{{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
            const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) * ((n[0] * out[0] + n[1] * out[1] + n[2] * out[2] < 0.0) ? -1.0 : 1.0);
            for (auto& c : n) { c /= l; }
            return n;
        };
        int chords = 0, arm_braces = 0;
        for (int sy = -1; sy <= 1; sy += 2) {
            const std::array<double,3> tip{{0.0, sy * 0.5 * g.arm_length, H}};
            const std::array<double,3> t0{{-hw, sy * hw, H}}, t1{{hw, sy * hw, H}}, b0{{-hb, sy * hb, H - g.arm_depth}}, b1{{hb, sy * hb, H - g.arm_depth}};
            // the top, the bottom, the x < 0 side and the x > 0 side
            const std::array<std::array<double,3>,4> fn{{plane(t0, t1, tip, {{0.0, 0.0, 1.0}}), plane(b0, b1, tip, {{0.0, 0.0, -1.0}}),
                                                         plane(t0, b0, tip, {{-1.0, 0.0, 0.0}}), plane(t1, b1, tip, {{1.0, 0.0, 0.0}})}};
            auto in_face = [&] (const std::array<double,3>& p, std::size_t f) {
                return std::abs(fn[f][0] * (p[0] - tip[0]) + fn[f][1] * (p[1] - tip[1]) + fn[f][2] * (p[2] - tip[2])) < 1.0e-9;
            };
            for (std::size_t m = 0; m < in.members.size(); ++m) {
                const auto& mem = in.members[m];
                const auto a = at(mem.joint_a), b = at(mem.joint_b);
                // on this side's cross-arm: from its root corners out, one end past the shaft
                if (std::min(a[2], b[2]) < H - g.arm_depth - 1.0e-9 || std::max(a[2], b[2]) > H + 1.0e-9) { continue; }
                if (sy * a[1] < hw - 1.0e-9 || sy * b[1] < hw - 1.0e-9 || std::max(sy * a[1], sy * b[1]) < hb + 1.0e-9) { continue; }
                const std::array<double,9> d = direction_cosines(a, b, mem.spin);
                const std::array<double,3> x{{d[0], d[3], d[6]}}, e{{d[2], d[5], d[8]}};
                std::vector<std::size_t> faces;
                for (std::size_t f = 0; f < 4; ++f) { if (in_face(a, f) && in_face(b, f)) { faces.push_back(f); } }
                if (faces.size() == 2) {
                    // a chord, on the edge of two faces
                    const auto& p = fn[faces[0]];
                    const auto& q = fn[faces[1]];
                    const auto rs = normal_to({{p[0] + q[0], p[1] + q[1], p[2] + q[2]}}, e);
                    EXPECT_NEAR(std::abs(x[0] * rs[0] + x[1] * rs[1] + x[2] * rs[2]), 1.0, 1.0e-9) << "chord " << mem.id;
                    ++chords;
                } else if (faces.size() == 1) {
                    const auto np = normal_to(fn[faces[0]], e);
                    EXPECT_NEAR(std::abs(x[0] * np[0] + x[1] * np[1] + x[2] * np[2]), std::sqrt(0.5), 1.0e-9) << "arm brace " << mem.id;
                    ++arm_braces;
                }
            }
        }
        EXPECT_GT(chords, 8) << "every cross-arm panel's four chords, both sides";
        EXPECT_GT(arm_braces, 8) << "the cross-arms' struts and face diagonals";
        // the turned angles keep the frame as symmetric as its members are: at the cross-arm's centre, a load along x
        // moves it along x, with y and a twist only as much as the frame with every MSpin 0 (the stiffness's
        // couplings, against its diagonal)
        {
            const auto f = make(in);
            ASSERT_TRUE(f);
            const auto k = condensed(*f, in.interface_joints[0]);
            auto coupling = [&] (std::size_t i, std::size_t j) { return std::abs(k[6 * i + j]) / std::sqrt(k[7 * i] * k[7 * j]); };
            RecordProperty(crossed ? "kxy_crossed" : "kxy_single", std::to_string(coupling(0, 1)));
            RecordProperty(crossed ? "kxrx_crossed" : "kxrx_single", std::to_string(coupling(0, 3)));
            RecordProperty(crossed ? "kxrz_crossed" : "kxrz_single", std::to_string(coupling(0, 5)));
            // the same couplings with every MSpin 0, for comparison
            LatticeSpec plain = g;
            plain.principal_axes = false;
            FrameInputs pin;
            ASSERT_TRUE(lattice_frame(plain, pin, designs).empty());
            const auto pf = make(pin);
            ASSERT_TRUE(pf);
            const auto k0 = condensed(*pf, pin.interface_joints[0]);
            auto coupling0 = [&] (std::size_t i, std::size_t j) { return std::abs(k0[6 * i + j]) / std::sqrt(k0[7 * i] * k0[7 * j]); };
            RecordProperty(crossed ? "kxy0_crossed" : "kxy0_single", std::to_string(coupling0(0, 1)));
            RecordProperty(crossed ? "kxrx0_crossed" : "kxrx0_single", std::to_string(coupling0(0, 3)));
            if (crossed) {
                // crossed bracing: within the generated frame's own asymmetry (the members without a mirror image,
                // each cross-arm station's one inner diagonal and the cross-arms' top and bottom diagonals), 2.2e-5
                // with every MSpin 0; angles turned without regard to the mirror images couple x with the twist
                // about x at 9e-5
                EXPECT_LT(coupling(0, 1), 1.0e-5) << "x with y";
                EXPECT_LT(coupling(0, 3), 1.0e-5) << "x with the twist about x";
                EXPECT_LT(coupling(0, 5), 1.0e-5) << "x with the twist about z";
                EXPECT_LT(coupling(1, 5), 1.0e-5) << "y with the twist about z";
                EXPECT_LT(coupling(1, 4), std::max(1.0e-5, coupling0(1, 4))) << "y with the twist about y";
            } else {
                // single bracing alternates its diagonals, which couples the frame by itself: the angles add little
                EXPECT_LT(coupling(0, 1), 1.5 * coupling0(0, 1)) << "x with y";
                EXPECT_LT(coupling(0, 3), 1.5 * coupling0(0, 3)) << "x with the twist about x";
            }
        }
    }
    const LatticeSpec g = case_g();
    std::vector<MemberDesign> designs;
    // without principal_axes every angle keeps MSpin 0
    LatticeSpec flat = g;
    flat.principal_axes = false;
    FrameInputs zero;
    ASSERT_TRUE(lattice_frame(flat, zero, designs).empty());
    for (const auto& mem : zero.members) { EXPECT_EQ(mem.spin, 0.0) << mem.id; }
}

TEST(LatticeFrame, TheStiffnessAtTheCrossArmIsSubDyns)
{
    // case G was written by write_subdyn from case_g() with every MSpin 0 and run through SubDyn's driver:
    // this checks the frame's stiffness against SubDyn's for the same members (the angles' turning is
    // LatticeFrame.TheAnglesLieOnTheirPrincipalAxes's)
    FrameInputs in;
    std::vector<MemberDesign> designs;
    LatticeSpec g = case_g();
    g.principal_axes = false;
    ASSERT_TRUE(lattice_frame(g, in, designs).empty());
    FrameInputs theirs;
    const std::string dir = std::string(ERF_FRAME_TEST_FILES) + "/caseG/";
    const std::string err = read_subdyn(dir + "towerG.dat", theirs);
    ASSERT_TRUE(err.empty()) << err;
    ASSERT_EQ(theirs.joints.size(), in.joints.size());
    // the file was written on another machine, where a coordinate can round one bit apart: the joints agree to 1e-12 m
    for (std::size_t j = 0; j < in.joints.size(); ++j) {
        for (std::size_t d = 0; d < 3; ++d) { EXPECT_NEAR(theirs.joints[j].x[d], in.joints[j].x[d], 1.0e-12) << "joint " << j + 1; }
    }
    ASSERT_EQ(theirs.members.size(), in.members.size());
    const auto f = make(in);
    ASSERT_TRUE(f);
    const auto k = condensed(*f, in.interface_joints[0]);
    const std::vector<double> kbbt = numbers_in(dir + "towerG_kbbt.txt");
    ASSERT_EQ(kbbt.size(), 36u);
    for (std::size_t i = 0; i < 6; ++i) {
        for (std::size_t j = 0; j < 6; ++j) {
            const double scale = std::sqrt(kbbt[7 * i] * kbbt[7 * j]);
            EXPECT_NEAR(k[6 * i + j], kbbt[6 * i + j], 2.0e-6 * scale) << "KBBt(" << i + 1 << "," << j + 1 << ")";
        }
    }
    const std::vector<double> freq = numbers_in(dir + "towerG_frequencies.txt");
    ASSERT_GE(freq.size(), 10u);
    FrameModes modes;
    const std::string merr = frame_modes(*f, 10, modes);
    ASSERT_TRUE(merr.empty()) << merr;
    for (std::size_t m = 0; m < 10; ++m) { EXPECT_NEAR(modes.frequency[m], freq[m], 2.0e-6 * freq[m]) << "mode " << m + 1; }
}

TEST(HeatedFrame, DeflectsByOneOverKE)
{
    FrameInputs in;
    std::vector<MemberDesign> designs;
    ASSERT_TRUE(lattice_frame(case_g(), in, designs).empty());
    const auto cold = make(in);
    FrameInputs hot_in = in;
    hot_in.temperature.assign(in.members.size(), 600.0);
    const auto hot = make(hot_in);
    ASSERT_TRUE(cold && hot);
    std::vector<double> load(cold->num_dofs(), 0.0);
    const std::size_t centre = cold->node_of_joint(in.interface_joints[0]);
    load[6 * centre] = 1.0e4;
    load[6 * centre + 1] = -3.0e3;
    const FrameSolution sc = cold->solve(load, 9.81), sh = hot->solve(load, 9.81);
    double umax = 0.0;
    for (const double u : sc.displacement) { umax = std::max(umax, std::abs(u)); }
    for (std::size_t i = 0; i < sc.displacement.size(); ++i) {
        EXPECT_NEAR(sh.displacement[i], sc.displacement[i] / 0.31, 1e-9 * umax / 0.31) << i;
    }
    // the member forces are the same (every member softened alike); the checks see the weaker steel
    const auto cc = check_members(*cold, designs, sc.element_force, {});
    const auto ch = check_members(*hot, designs, sh.element_force, hot_in.temperature);
    for (std::size_t m = 0; m < cc.size(); ++m) {
        EXPECT_NEAR(ch[m].compression, cc[m].compression, 1e-6 * (cc[m].compression + 1.0));
        EXPECT_NEAR(ch[m].tension_capacity, 0.47 * cc[m].tension_capacity, 1e-9 * cc[m].tension_capacity);
        EXPECT_GE(ch[m].utilisation, cc[m].utilisation);
    }
    // a member at 1200 C has no stiffness left: refused
    hot_in.temperature[0] = 1200.0;
    std::string err;
    EXPECT_FALSE(Frame::create(hot_in, err));
    EXPECT_NE(err.find("below 1200"), std::string::npos) << err;
    hot_in.temperature.pop_back();
    EXPECT_FALSE(Frame::create(hot_in, err));
    EXPECT_NE(err.find("member temperatures"), std::string::npos) << err;
}
