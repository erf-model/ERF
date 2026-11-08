// Contract of the lattice tower (shaft and cross-arm drag, foundation reactions) and of its type's inputs.
//
// - MemberDrag.OnlyTheFlowNormalToTheMemberAndRelativeToItLoadsIt: the drag on a drag node is
//   1/2 rho Cd w L |U_n| U_n, only the flow normal to the member and relative to it counting.
// - MemberDrag.TheLatticeForceCoefficientIsTheSquareTowerCurve: the lattice force coefficient is
//   the square-tower curve of ASCE 7 (American Society of Civil Engineers, minimum design loads).
// - TowerType.EveryValueOutsideItsRangeIsRefusedByName: validate() names the key of every value
//   outside its range, NaN and infinity included; the foundation stiffnesses and damping_ratio on a
//   tower that stands still, leg_spacing with a frame, a frame without damping, a cross-arm no longer than
//   the top width with arm_outside_shaft, and angle_principal_axes without frame_panels are refused; a type's own
//   drag width, epsilon, must be positive when given and is not checked (and is 0) when it is not.
// - Tower.TheBodyTapersUpToTheCrossArmWhichRunsAcrossTheLine: the shaft's drag nodes stand up the
//   tapering shaft with the direction factor's gain; the cross-arm's face is arm_depth deep out to the shaft's
//   corners at the arm's bottom and tapers linearly to the tips (arm_face_at), its nodes across the line outside the shaft at the
//   centroids of their pieces of the face with the face's mean depth; arm_outside_shaft = false and
//   diagonal_wind_factor = false give the old layout (the face's depth along the whole length).
// - Tower.InAUniformWindTheDragAndBaseMomentAreTheHandValues: the shaft's drag exactly (its width
//   is linear in height), the base moment to the midpoint rule's known error, the cross-arm's the closed
//   integrals of its tapering face outside the shaft, force and moment both exact; along a diagonal the
//   shaft's drag times ASCE 7's 1 + 0.75 phi, half that gain at 22.5 degrees.
// - Tower.InALogLawWindTheDragIsTheFineIntegral: within the segments' quadrature error.
// - Tower.TheLegsShareTheLoadAndResistTheOverturningMoment: the four legs share the downward load
//   equally and the overturning moment linearly, in equilibrium, whichever way the tower faces;
//   the corner leg takes sqrt(2) times more from a diagonal pull.
// - Tower.EachLegLoadIsFlaggedOverItsAllowable.
// - Tower.EachLinePullsWhereItHangsAndTheFootingsTakeThemAll: several lines on one tower.
// - Tower.ANonFiniteLinePullIsRefusedNamingTheTower and
//   Tower.AFlatOrBadlyAimedTowerIsRefusedNamingIt: the preconditions abort with the tower's name, and a cross-arm
//   no longer than the shaft's width at its face's middle, or whose face's middle is below the base, names its keys.

#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <gtest/gtest.h>

#include "ERF_GTestThrowOnAbort.H"
#include "ERF_MemberDrag.H"
#include "ERF_Tower.H"
#include "ERF_TowerInputs.H"

using amrex::Real;
using erf_towers::MemberDrag;
using erf_towers::MemberNode;
using erf_towers::Tower;
using erf_towers::TowerType;
using P3 = std::array<Real,3>;

namespace {

constexpr Real tol = (std::is_same<Real, float>::value) ? Real(1.0e-5) : Real(1.0e-12);
constexpr Real rho = 1.2;

TowerType lattice ()
{
    TowerType t;
    t.name = "lattice";
    t.base_width = 6.0; t.top_width = 1.5; t.solidity = 0.2; t.arm_length = 12.0; t.arm_depth = 1.2;
    return t;
}

// the drag on every node of a tower in a wind given as a function of the node's height above its base
template <typename F>
void blow (Tower& tw, F wind)
{
    std::vector<Real> u, v;
    for (const auto& n : tw.nodes()) {
        const P3 w = wind(n.pos[2] - tw.base()[2]);
        u.insert(u.end(), {w[0], w[1], w[2]});
        v.insert(v.end(), {0.0, 0.0, 0.0});
    }
    std::vector<Real> f;
    MemberDrag(rho).loads(tw.nodes(), u, v, f);
    tw.set_loads(f);
}

} // namespace

TEST(MemberDrag, OnlyTheFlowNormalToTheMemberAndRelativeToItLoadsIt)
{
    MemberNode n;
    n.axis = {{0.0, 0.0, 1.0}}; n.length = 2.0; n.drag_width = 0.5; n.drag_coefficient = 1.5;
    const P3 still{{0.0, 0.0, 0.0}};
    // normal: 1/2 rho Cd w L U^2 along the wind
    auto f = erf_towers::member_drag(n, {{10.0, 0.0, 0.0}}, still, rho);
    EXPECT_NEAR(f[0], 0.5 * 1.2 * 1.5 * 0.5 * 2.0 * 100.0, 1.0e3 * tol);
    EXPECT_EQ(f[1], 0.0); EXPECT_EQ(f[2], 0.0);
    // along the axis: nothing
    f = erf_towers::member_drag(n, {{0.0, 0.0, 10.0}}, still, rho);
    EXPECT_EQ(f[0], 0.0); EXPECT_EQ(f[1], 0.0); EXPECT_EQ(f[2], 0.0);
    // oblique, 30 degrees off the axis: the normal half of the speed, in the normal direction
    f = erf_towers::member_drag(n, {{Real(5.0), Real(0.0), Real(5.0 * std::sqrt(Real(3.0)))}}, still, rho);
    EXPECT_NEAR(f[0], 0.5 * 1.2 * 1.5 * 0.5 * 2.0 * 25.0, 1.0e3 * tol);
    EXPECT_NEAR(f[2], 0.0, 1.0e3 * tol);
    // a node moving with the wind feels none of it; one moving against it, the sum
    f = erf_towers::member_drag(n, {{10.0, 0.0, 0.0}}, {{10.0, 0.0, 0.0}}, rho);
    EXPECT_EQ(f[0], 0.0);
    f = erf_towers::member_drag(n, {{10.0, 0.0, 0.0}}, {{-10.0, 0.0, 0.0}}, rho);
    EXPECT_NEAR(f[0], 4.0 * 0.5 * 1.2 * 1.5 * 0.5 * 2.0 * 100.0, 1.0e4 * tol);
}

TEST(MemberDrag, TheLatticeForceCoefficientIsTheSquareTowerCurve)
{
    EXPECT_NEAR(erf_towers::lattice_force_coefficient(0.0), 4.0, tol);
    EXPECT_NEAR(erf_towers::lattice_force_coefficient(0.2), 4.0 * 0.04 - 5.9 * 0.2 + 4.0, 10 * tol);   // 2.98
    EXPECT_NEAR(erf_towers::lattice_force_coefficient(0.5), 2.05, 10 * tol);
    TowerType t = lattice();
    EXPECT_NEAR(t.force_coefficient(), 2.98, 10 * tol);
    t.drag_coefficient = 3.4;
    EXPECT_EQ(t.force_coefficient(), Real(3.4)) << "an explicit coefficient wins";
}

TEST(TowerType, EveryValueOutsideItsRangeIsRefusedByName)
{
    EXPECT_TRUE(lattice().validate().empty());
    auto bad = [](auto mutate, const std::string& key) {
        TowerType t = lattice();
        mutate(t);
        const std::string err = t.validate();
        EXPECT_NE(err.find("erf.conductors.lattice." + key), std::string::npos) << key << ": " << err;
    };
    bad([](TowerType& t) { t.base_width = 0.0; }, "base_width");
    bad([](TowerType& t) { t.top_width = 7.0; }, "top_width");
    bad([](TowerType& t) { t.top_width = 0.0; }, "top_width");
    bad([](TowerType& t) { t.solidity = 0.0; }, "solidity");
    bad([](TowerType& t) { t.solidity = 1.0; }, "solidity");
    bad([](TowerType& t) { t.arm_length = 0.0; }, "arm_length");
    bad([](TowerType& t) { t.arm_depth = -1.0; }, "arm_depth");
    bad([](TowerType& t) { t.peak = -1.0; }, "peak");
    bad([](TowerType& t) { t.drag_coefficient = -1.0; }, "drag_coefficient");
    bad([](TowerType& t) { t.segments = 0; }, "segments");
    bad([](TowerType& t) { t.weight = -1.0; }, "weight");
    bad([](TowerType& t) { t.leg_spacing = -1.0; }, "leg_spacing");
    bad([](TowerType& t) { t.allowable_uplift = -1.0; }, "allowable_uplift");
    bad([](TowerType& t) { t.allowable_compression = -1.0; }, "allowable_compression");
    bad([](TowerType& t) { t.frequency = -1.0; }, "frequency");
    bad([](TowerType& t) { t.frequency = 2.0; t.weight = 0.0; }, "frequency");
    bad([](TowerType& t) { t.damping_ratio = -0.1; }, "damping_ratio");
    bad([](TowerType& t) { t.damping_ratio = 1.0; }, "damping_ratio");
    bad([](TowerType& t) { t.foundation_rotational_stiffness = -1.0; }, "foundation_rotational_stiffness");
    bad([](TowerType& t) { t.foundation_lateral_stiffness = -1.0; }, "foundation_lateral_stiffness");
    bad([](TowerType& t) { t.epsilon_given = true; t.epsilon = 0.0; }, "epsilon must be positive");
    bad([](TowerType& t) { t.epsilon_given = true; t.epsilon = -1.0; }, "epsilon must be positive");
    {
        TowerType t = lattice();
        EXPECT_FALSE(t.epsilon_given);
        EXPECT_TRUE(t.validate().empty()) << "no epsilon of its own: the towers take erf.conductors.epsilon";
        t.epsilon_given = true;
        t.epsilon = 0.5;
        EXPECT_TRUE(t.validate().empty()) << t.validate();
    }
    // NaN passes a test written as x < 0 and infinity one written as x > 0: every key refuses both
    const Real nan = std::numeric_limits<Real>::quiet_NaN(), inf = std::numeric_limits<Real>::infinity();
    for (const Real v : {nan, inf}) {
        bad([v](TowerType& t) { t.base_width = v; }, "base_width");
        bad([v](TowerType& t) { t.top_width = v; }, "top_width");
        bad([v](TowerType& t) { t.solidity = v; }, "solidity");
        bad([v](TowerType& t) { t.arm_length = v; }, "arm_length");
        bad([v](TowerType& t) { t.arm_depth = v; }, "arm_depth");
        bad([v](TowerType& t) { t.peak = v; }, "peak");
        bad([v](TowerType& t) { t.drag_coefficient = v; }, "drag_coefficient");
        bad([v](TowerType& t) { t.weight = v; }, "weight");
        bad([v](TowerType& t) { t.leg_spacing = v; }, "leg_spacing");
        bad([v](TowerType& t) { t.allowable_uplift = v; }, "allowable_uplift");
        bad([v](TowerType& t) { t.allowable_compression = v; }, "allowable_compression");
        bad([v](TowerType& t) { t.frequency = v; }, "frequency");
        bad([v](TowerType& t) { t.damping_ratio = v; }, "damping_ratio");
        bad([v](TowerType& t) { t.foundation_rotational_stiffness = v; }, "foundation_rotational_stiffness");
        bad([v](TowerType& t) { t.foundation_lateral_stiffness = v; }, "foundation_lateral_stiffness");
        bad([v](TowerType& t) { t.epsilon_given = true; t.epsilon = v; }, "epsilon must be finite");
    }
    TowerType moving = lattice();
    moving.frequency = 2.0;
    moving.weight = 9.0e4;
    EXPECT_TRUE(moving.validate().empty()) << moving.validate();
    EXPECT_TRUE(moving.moves());
    EXPECT_FALSE(lattice().moves()) << "a tower stands still unless it has a frequency";
    TowerType d = lattice();
    d.arm_depth = 0.0;
    EXPECT_EQ(d.arm_face(), d.top_width) << "the arm's face defaults to the top width";
    // keys of a tower's motion on a tower that stands still, and leg_spacing with a frame, would do nothing
    bad([](TowerType& t) { t.foundation_rotational_stiffness = 1.0e9; },
        "foundation_rotational_stiffness needs a tower that moves in one mode");
    bad([](TowerType& t) { t.foundation_lateral_stiffness = 1.0e7; }, "foundation_lateral_stiffness needs a tower that moves in one mode");
    bad([](TowerType& t) { t.damping_given = true; }, "damping_ratio needs a tower that moves");
    {
        // with frequency = 0 given the motion is off on purpose: the keys are unused, a warning
        TowerType off = lattice();
        off.frequency_given = true;
        off.damping_given = true;
        off.foundation_lateral_stiffness = 1.0e7;
        EXPECT_TRUE(off.validate().empty()) << off.validate();
        const std::string note = off.unused_on_a_still_tower();
        EXPECT_NE(note.find("damping_ratio, foundation_lateral_stiffness not used"), std::string::npos) << note;
        EXPECT_TRUE(lattice().unused_on_a_still_tower().empty());
    }
    moving.foundation_rotational_stiffness = 1.0e9;
    moving.damping_given = true;
    EXPECT_TRUE(moving.validate().empty()) << moving.validate();
    bad([](TowerType& t) { t.arm_length = 1.5; }, "arm_length must exceed top_width for erf.conductors.lattice.arm_outside_shaft");
    bad([](TowerType& t) { t.arm_length = 1.5; }, "arm_outside_shaft = false for a cross-arm within");
    {
        TowerType narrow = lattice();
        narrow.arm_length = 1.5;
        narrow.arm_outside_shaft = false;
        EXPECT_TRUE(narrow.validate().empty()) << narrow.validate();
    }
    bad([](TowerType& t) { t.angle_axes_given = true; }, "angle_principal_axes needs erf.conductors.lattice.frame_panels");
    bad([](TowerType& t) { t.frame_panels = 4; t.leg_angle = {0.1, 0.01}; t.brace_angle = {0.08, 0.006}; t.leg_spacing = 5.0; },
        "leg_spacing is not given with a frame");
    // an undamped frame's higher modes grow under the members' drag
    bad([](TowerType& t) { t.frame_panels = 4; t.leg_angle = {0.1, 0.01}; t.brace_angle = {0.08, 0.006}; t.damping_ratio = 0.0; },
        "damping_ratio must be positive for a frame");
}

TEST(Tower, TheBodyTapersUpToTheCrossArmWhichRunsAcrossTheLine)
{
    TowerType t = lattice();
    t.peak = 4.5;
    const Tower tw("L1_t1", t, {{100.0, 200.0, 50.0}}, 30.0, {{0.0, 1.0, 0.0}});
    // ten body segments of 3 m, two peak segments of 2.25 m, and the cross-arm outside the shaft: two segments of
    // 2.6025 m on each side, from 0.795 m (half the shaft's width at the middle of the face's depth at the shaft,
    // 29.4 m: 6 - 4.5 x 29.4 / 30 = 1.59 m) to 6 m. The face is 1.2 m deep out to the shaft's corners at the arm's
    // bottom (28.8 m, half width 0.84 m) and tapers to 0 at the tips, d(s) = 1.2 (6 - s) / 5.16 beyond, so a piece
    // from a to b on the taper has the area (d(a) + d(b)) (b - a) / 2, its centroid at
    // a + (b - a) (d(a) + 2 d(b)) / (3 (d(a) + d(b))) along the arm and (d(a)^2 + d(a) d(b) + d(b)^2) / (3 (d(a) + d(b)))
    // below the top; a piece across 0.84 m adds its full-depth part. The inner piece [0.795, 3.3975]: 1.2 m deep to
    // 0.84 m, then from 1.2 to 0.60523 m, area 2.36244 m^2 (mean depth 0.90776 m), centroid 1.95178 m out and
    // 0.47066 m down; the tip piece [3.3975, 6]: d from 0.60523 to 0, area 0.78756 m^2 (0.30262 m), centroid
    // 4.265 m out, 0.20174 m down
    EXPECT_NEAR(tw.arm_taper_start(), 0.84, tol * 10) << "half the shaft's width at the arm's bottom";
    EXPECT_NEAR(tw.arm_face_at(0.0), 1.2, tol * 10);
    EXPECT_NEAR(tw.arm_face_at(-0.84), 1.2, tol * 10) << "the full depth out to the shaft's bottom corners";
    EXPECT_NEAR(tw.arm_face_at(3.42), 0.6, tol * 10) << "half way out along the taper, half the depth";
    EXPECT_NEAR(tw.arm_face_at(-3.42), 0.6, tol * 10);
    EXPECT_EQ(tw.arm_face_at(6.0), 0.0) << "nothing at the tip";
    EXPECT_EQ(tw.arm_face_at(7.0), 0.0) << "nothing beyond it";
    ASSERT_EQ(tw.num_body_nodes(), 12);
    ASSERT_EQ(tw.nodes().size(), 16u);
    const auto& n0 = tw.nodes()[0];
    EXPECT_NEAR(n0.pos[2], 51.5, tol * 100);
    EXPECT_NEAR(n0.length, 3.0, tol * 10);
    EXPECT_NEAR(n0.drag_width, 0.2 * (6.0 - 4.5 * 1.5 / 30.0), tol * 10) << "solidity x width at 1.5 m";
    EXPECT_NEAR(tw.nodes()[11].pos[2], 50.0 + 30.0 + 3.375, tol * 100);
    EXPECT_NEAR(tw.nodes()[11].drag_width, 0.2 * 1.5, tol * 10) << "the peak keeps the top width";
    Real arm = 0.0;
    for (std::size_t i = 12; i < 16; ++i) {
        const auto& n = tw.nodes()[i];
        EXPECT_NEAR(n.pos[0], 100.0, tol * 100);
        EXPECT_EQ(n.axis, (P3{{0.0, 1.0, 0.0}}));
        EXPECT_NEAR(n.length, 2.6025, tol * 10);
        EXPECT_GE(std::abs(n.pos[1] - 200.0), 0.795) << "outside the shaft";
        EXPECT_EQ(n.diagonal_gain, 0.0) << "the cross-arm takes no shaft direction factor";
        arm += n.pos[1] - 200.0;
    }
    EXPECT_NEAR(arm, 0.0, tol * 1000) << "the arm is centred on the body";
    // from one tip to the other: the tip piece, the inner piece, then the same on the other side
    EXPECT_NEAR(tw.nodes()[12].pos[1], 200.0 - 4.265, tol * 1000);
    EXPECT_NEAR(tw.nodes()[13].pos[1], 200.0 - 1.9517804569, tol * 1000);
    EXPECT_NEAR(tw.nodes()[14].pos[1], 200.0 + 1.9517804569, tol * 1000);
    EXPECT_NEAR(tw.nodes()[15].pos[1], 200.0 + 4.265, tol * 1000);
    for (const std::size_t i : {std::size_t(12), std::size_t(15)}) {
        EXPECT_NEAR(tw.nodes()[i].pos[2], 80.0 - 0.2017441860, tol * 1000) << i << ": the tip piece's centroid";
        EXPECT_NEAR(tw.nodes()[i].drag_width, 0.2 * 0.3026162791, tol * 100) << i << ": solidity x its mean depth";
    }
    for (const std::size_t i : {std::size_t(13), std::size_t(14)}) {
        EXPECT_NEAR(tw.nodes()[i].pos[2], 80.0 - 0.4706634005, tol * 1000) << i << ": the inner piece's centroid";
        EXPECT_NEAR(tw.nodes()[i].drag_width, 0.2 * 0.9077583607, tol * 100) << i;
    }
    // the four pieces' areas add up to the face outside the shaft: per side the full-depth strip from 0.795 to 0.84 m
    // plus the triangle from 0.84 to 6 m, d0 (0.045 + 5.16 / 2) = 3.15 m^2, 6.3 m^2 in all: 50.4 % of the 12.492 m^2
    // rectangle the model took before
    Real face = 0.0;
    for (std::size_t i = 12; i < 16; ++i) { face += tw.nodes()[i].drag_width * tw.nodes()[i].length / 0.2; }
    EXPECT_NEAR(face, 6.3, tol * 1000);
    // the shaft's and the peak's nodes take the direction factor, with the cross-arm's face normal
    EXPECT_NEAR(tw.nodes()[0].diagonal_gain, 0.15, tol) << "0.75 x solidity 0.2";
    EXPECT_EQ(tw.nodes()[0].face, (P3{{0.0, 1.0, 0.0}}));
    EXPECT_NEAR(tw.nodes()[11].diagonal_gain, 0.15, tol);
    // the old forms: the cross-arm along its whole length at its height, the face's depth all along, no direction factor
    t.arm_outside_shaft = false;
    t.diagonal_wind_factor = false;
    const Tower old("L1_t1", t, {{100.0, 200.0, 50.0}}, 30.0, {{0.0, 1.0, 0.0}});
    EXPECT_NEAR(old.nodes()[12].pos[1], 200.0 - 4.5, tol * 1000);
    EXPECT_NEAR(old.nodes()[13].pos[1], 200.0 - 1.5, tol * 1000);
    EXPECT_NEAR(old.nodes()[12].pos[2], 80.0, tol * 100);
    EXPECT_NEAR(old.nodes()[12].length, 3.0, tol * 10);
    for (std::size_t i = 12; i < 16; ++i) { EXPECT_NEAR(old.nodes()[i].drag_width, 0.2 * 1.2, tol * 10) << i; }
    EXPECT_EQ(old.nodes()[0].diagonal_gain, 0.0);
}

TEST(Tower, InAUniformWindTheDragAndBaseMomentAreTheHandValues)
{
    const TowerType t = lattice();
    const Real H = 30.0, U = 20.0, q = 0.5 * rho * U * U, cf = t.force_coefficient(), phi = t.solidity;
    Tower tw("T", t, {{0.0, 0.0, 10.0}}, H, {{0.0, 1.0, 0.0}});
    // wind along x, across the arm (along y): body and arm both loaded
    blow(tw, [&](Real) { return P3{{U, 0.0, 0.0}}; });
    const Real body = q * cf * phi * 0.5 * (t.base_width + t.top_width) * H;    // exact: the width is linear
    // the cross-arm's two parts outside the shaft, from the shaft's half width w(za) / 2 at the middle of the face's
    // depth at the shaft, za = H - d0 / 2, to the tips; the face is d0 deep out to the shaft's half width s0 at the
    // arm's bottom, H - d0, and tapers to the tip beyond, so each part's area is d0 (s0 - w(za)/2) + d0 (L/2 - s0) / 2,
    // 3.15 m^2, and the drag on it q cf phi times that, which the nodes' mean depths give exactly
    const Real za = H - 0.5 * t.arm_depth;
    const Real wza = t.base_width + (t.top_width - t.base_width) * za / H;
    const Real s0 = 0.5 * (t.base_width + (t.top_width - t.base_width) * (H - t.arm_depth) / H);
    const Real taper = 0.5 * t.arm_length - s0;
    const Real area1 = t.arm_depth * ((s0 - 0.5 * wza) + 0.5 * taper);
    const Real armF = q * cf * phi * 2.0 * area1;
    const auto F = tw.total_force();
    RecordProperty("body_drag_N", std::to_string(body));
    RecordProperty("arm_drag_N", std::to_string(armF));
    EXPECT_NEAR(F[0], body + armF, 1.0e-6 * (body + armF) + tol);
    EXPECT_NEAR(F[1], 0.0, tol); EXPECT_NEAR(F[2], 0.0, tol);
    EXPECT_NEAR(armF / (q * cf * phi * t.arm_depth * (t.arm_length - wza)), 0.5043227666, 100 * tol)
        << "50.4 % of the rectangle's drag the model took before";
    // the base moment: the body's q cf phi int w(z) z dz, with w linear, and the arm's q cf phi int d (H - d/2) ds
    // over both parts, int d^2 ds = d0^2 (s0 - w(za)/2) + d0^2 taper / 3 per part, which the nodes at their pieces'
    // centroids give exactly. The midpoint sum over n segments is exact on z and short by H^3 / (12 n^2) on z^2.
    const Real n = t.segments;
    const Real body_m = q * cf * phi * (t.base_width * H * H / 2.0 +
                                        (t.top_width - t.base_width) / H * (H * H * H / 3.0 - H * H * H / (12.0 * n * n)));
    const Real d2_1 = t.arm_depth * t.arm_depth * ((s0 - 0.5 * wza) + taper / 3.0);
    const Real arm_m = q * cf * phi * (2.0 * H * area1 - d2_1);
    EXPECT_NEAR(arm_m / armF, 29.5965714286, 1.0e4 * tol) << "the arm's drag acts 0.4034 m below the cross-arm";
    const auto M = tw.base_moment();
    EXPECT_NEAR(M[1], body_m + arm_m, 1.0e-6 * (body_m + arm_m) + tol);
    EXPECT_NEAR(M[0], 0.0, 1.0e-6 * M[1] + tol);
    // wind along the arm: only the body carries it
    blow(tw, [&](Real) { return P3{{0.0, U, 0.0}}; });
    EXPECT_NEAR(tw.total_force()[1], body, 1.0e-6 * body + tol);
    EXPECT_NEAR(tw.total_force()[0], 0.0, tol);
    // wind along the shaft's diagonal: ASCE 7's 1 + 0.75 phi on the shaft (1.15 here), half way at 22.5 degrees
    // (IEC 60826's sin^2 2 theta); the cross-arm takes the wind's part across it only
    const Real gain = 0.75 * phi;
    for (const Real deg : {Real(45.0), Real(22.5)}) {
        const Real th = deg * Real(3.14159265358979323846) / Real(180.0);
        blow(tw, [&](Real) { return P3{{U * std::cos(th), U * std::sin(th), 0.0}}; });
        Real fx = 0.0, fy = 0.0;
        for (int i = 0; i < tw.num_body_nodes(); ++i) {
            fx += tw.loads()[static_cast<std::size_t>(3 * i)];
            fy += tw.loads()[static_cast<std::size_t>(3 * i + 1)];
        }
        const Real s2 = std::sin(2.0 * th);
        EXPECT_NEAR(std::hypot(fx, fy) / (body * (1.0 + gain * s2 * s2)), 1.0, 10.0 * tol) << deg << " degrees";
        EXPECT_NEAR(std::atan2(fy, fx), th, 10.0 * tol) << "along the wind";
    }
    EXPECT_NEAR(erf_towers::lattice_diagonal_gain(0.5), 0.2, tol) << "ASCE 7's cap of 1.2";
}

TEST(Tower, InALogLawWindTheDragIsTheFineIntegral)
{
    const TowerType t = lattice();
    const Real H = 30.0, us = 1.26, z0 = 0.1, kappa = 0.4;
    auto u = [&](Real z) { return us / kappa * std::log((z + z0) / z0); };
    Tower tw("T", t, {{0.0, 0.0, 0.0}}, H, {{0.0, 1.0, 0.0}});
    blow(tw, [&](Real z) { return P3{{u(z), 0.0, 0.0}}; });
    // the body's drag and base moment by a 20000-point midpoint sum
    const Real k = 0.5 * rho * t.force_coefficient() * t.solidity;
    Real drag = 0.0, mom = 0.0;
    const int N = 20000;
    for (int i = 0; i < N; ++i) {
        const Real z = (i + 0.5) * H / N;
        const Real w = t.base_width + (t.top_width - t.base_width) * z / H;
        drag += k * w * u(z) * u(z) * H / N;
        mom += k * w * u(z) * u(z) * z * H / N;
    }
    // the cross-arm's two tapering parts outside the shaft, each point of the face at its own height, by a fine sum
    const Real za = H - 0.5 * t.arm_depth;
    const Real inner = 0.5 * (t.base_width + (t.top_width - t.base_width) * za / H);
    const Real s0 = 0.5 * (t.base_width + (t.top_width - t.base_width) * (H - t.arm_depth) / H);
    auto d = [&](Real s) { return (s <= s0) ? t.arm_depth : t.arm_depth * (0.5 * t.arm_length - s) / (0.5 * t.arm_length - s0); };
    Real arm = 0.0, arm_m = 0.0;
    for (int i = 0; i < N; ++i) {
        const Real s = inner + (i + 0.5) * (0.5 * t.arm_length - inner) / N, ds = (0.5 * t.arm_length - inner) / N;
        const Real zc = H - 0.5 * d(s);
        arm += 2.0 * k * d(s) * u(zc) * u(zc) * ds;
        arm_m += 2.0 * k * d(s) * u(zc) * u(zc) * zc * ds;
    }
    // ten segments resolve the log law's curvature near the ground to about 1 %
    EXPECT_NEAR(tw.total_force()[0] / (drag + arm), 1.0, 0.01);
    EXPECT_NEAR(tw.base_moment()[1] / (mom + arm_m), 1.0, 0.01);
}

namespace {
// the legs' reactions balance the tower's loads: their sum is the downward load and their moments
// about the base cancel the loads' overturning moment
void expect_equilibrium (const Tower& tw, const erf_towers::FoundationLoad& L)
{
    const Real a = 0.5 * tw.type().legs();
    const std::array<Real,2> across{{tw.across()[0], tw.across()[1]}}, along{{tw.across()[1], -tw.across()[0]}};
    Real sum = 0.0, mx = 0.0, my = 0.0;
    int i = 0;
    for (const Real su : {Real(1.0), Real(-1.0)}) {
        for (const Real sv : {Real(1.0), Real(-1.0)}) {
            const Real x = a * (su * along[0] + sv * across[0]), y = a * (su * along[1] + sv * across[1]);
            const Real R = L.legs[static_cast<std::size_t>(i++)];
            sum += R; mx += y * R; my -= x * R;
        }
    }
    const Real scale = std::abs(L.vertical) + L.overturning / a + 1.0;
    EXPECT_NEAR(sum, L.vertical, 1.0e2 * tol * scale);
    EXPECT_NEAR(mx, -L.moment[0], 1.0e2 * tol * scale * a);
    EXPECT_NEAR(my, -L.moment[1], 1.0e2 * tol * scale * a);
}
} // namespace

TEST(Tower, TheLegsShareTheLoadAndResistTheOverturningMoment)
{
    TowerType t = lattice();
    t.weight = 9.0e4;
    const Real H = 30.0, a = 3.0;   // legs 6 m apart, the base width
    Tower tw("T", t, {{50.0, 60.0, 10.0}}, H, {{0.0, 1.0, 0.0}});
    // the line pulls 10 kN along x and 8 kN down at the cross-arm, no wind
    const Real f = 1.0e4, down = 8.0e3;
    tw.set_line_load({{f, Real(0.0), Real(-down)}}, {{Real(50.0), Real(60.0), Real(10.0 + H)}});
    auto L = tw.foundation();
    const Real P = t.weight + down;
    EXPECT_NEAR(L.vertical, P, 1.0e2 * tol * P);
    EXPECT_NEAR(L.shear, f, 1.0e2 * tol * f);
    EXPECT_NEAR(L.overturning, H * f, 1.0e2 * tol * H * f);
    // face-on: two legs at P/4 + H f / (4 a), two at P/4 - H f / (4 a)
    EXPECT_NEAR(L.max_compression, P / 4 + H * f / (4 * a), 1.0e2 * tol * P);
    EXPECT_NEAR(L.max_uplift, std::max(Real(0.0), -(P / 4 - H * f / (4 * a))), 1.0e2 * tol * P);
    EXPECT_GT(L.max_uplift, 0.0) << "this pull lifts the windward legs: H f / 4a = 25 kN > P/4 = 24.5 kN";
    expect_equilibrium(tw, L);
    // the same pull along the diagonal: the corner leg takes sqrt(2) times the face-on moment share
    tw.set_line_load({{Real(f / std::sqrt(Real(2.0))), Real(f / std::sqrt(Real(2.0))), Real(-down)}}, {{Real(50.0), Real(60.0), Real(10.0 + H)}});
    L = tw.foundation();
    EXPECT_NEAR(L.max_compression - P / 4, std::sqrt(Real(2.0)) * H * f / (4 * a), 1.0e3 * tol * P);
    expect_equilibrium(tw, L);
    // a tower facing another way carries a pull across its line the same way
    Tower turned("T", t, {{50.0, 60.0, 10.0}}, H, {{1.0, 0.0, 0.0}});
    turned.set_line_load({{Real(0.0), f, Real(-down)}}, {{Real(50.0), Real(60.0), Real(10.0 + H)}});
    const auto Lt = turned.foundation();
    EXPECT_NEAR(Lt.max_compression, P / 4 + H * f / (4 * a), 1.0e2 * tol * P);
    expect_equilibrium(turned, Lt);
    // the wind's drag adds to the line's pull
    blow(tw, [&](Real) { return P3{{20.0, 0.0, 0.0}}; });
    tw.set_line_load({{f, Real(0.0), Real(-down)}}, {{Real(50.0), Real(60.0), Real(10.0 + H)}});
    L = tw.foundation();
    EXPECT_NEAR(L.force[0], f + tw.total_force()[0], 1.0e2 * tol * f);
    EXPECT_NEAR(L.moment[1], H * f + tw.base_moment()[1], 1.0e2 * tol * H * f);
    expect_equilibrium(tw, L);
}

TEST(Tower, EachLegLoadIsFlaggedOverItsAllowable)
{
    TowerType t = lattice();
    t.weight = 9.0e4;
    auto flagged = [&](Real uplift, Real compression) {
        TowerType u = t;
        u.allowable_uplift = uplift;
        u.allowable_compression = compression;
        Tower tw("T", u, {{0.0, 0.0, 0.0}}, 30.0, {{0.0, 1.0, 0.0}});
        tw.set_line_load({{1.0e4, 0.0, -8.0e3}}, {{0.0, 0.0, 30.0}});   // 49.5 kN compression, 0.5 kN uplift
        return tw.foundation().over_allowable;
    };
    EXPECT_FALSE(flagged(0.0, 0.0)) << "no allowable: not checked";
    EXPECT_TRUE(flagged(400.0, 0.0));
    EXPECT_FALSE(flagged(600.0, 0.0));
    EXPECT_TRUE(flagged(0.0, 4.9e4));
    EXPECT_FALSE(flagged(0.0, 5.0e4));
    EXPECT_TRUE(flagged(600.0, 4.9e4)) << "the compression alone is over";
    // legs default to the base width; a type gives its own spacing
    TowerType w = t;
    w.leg_spacing = 8.0;
    EXPECT_EQ(t.legs(), t.base_width);
    EXPECT_EQ(w.legs(), Real(8.0));
}

TEST(Tower, EachLinePullsWhereItHangsAndTheFootingsTakeThemAll)
{
    TowerType t = lattice();
    t.weight = 9.0e4;
    t.peak = 8.0;
    Tower tw("t", t, P3{{500.0, 500.0, 20.0}}, Real(30.0), P3{{0.0, 1.0, 0.0}});
    // three phases across the cross-arm and a shield wire on the peak
    const std::vector<P3> at{{{500.0, 500.0, 50.0}}, {{500.0, 505.5, 50.0}}, {{500.0, 494.5, 50.0}}, {{500.0, 500.0, 57.0}}};
    for (const auto& p : at) { tw.add_attachment(p); }
    EXPECT_EQ(tw.attachments().size(), 4u);
    const std::vector<P3> F{{{0.0, 1000.0, -5000.0}}, {{0.0, 1200.0, -5000.0}}, {{0.0, 800.0, -5000.0}}, {{0.0, 300.0, -1300.0}}};
    tw.set_line_loads(F, at);
    const P3 sum = tw.line_force();
    EXPECT_NEAR(sum[1], 3300.0, tol * 1.0e4);
    EXPECT_NEAR(sum[2], -16300.0, tol * 1.0e5);
    // each pull's own lever arm: M_x = sum (y_i - y_b) F_z - (z_i - z_b) F_y
    double Mx = 0.0;
    for (std::size_t i = 0; i < at.size(); ++i) { Mx += (at[i][1] - 500.0) * F[i][2] - (at[i][2] - 20.0) * F[i][1]; }
    const auto L = tw.foundation();
    EXPECT_NEAR(L.moment[0], Mx, 10.0 * tol * std::abs(Mx));
    EXPECT_NEAR(L.vertical, 9.0e4 + 16300.0, 1.0e-6 * L.vertical);
    expect_equilibrium(tw, L);
}

TEST(Tower, ANonFiniteLinePullIsRefusedNamingTheTower)
{
    TowerType t = lattice();
    t.weight = 9.0e4;
    Tower tw("L1_t3", t, P3{{500.0, 500.0, 20.0}}, Real(30.0), P3{{0.0, 1.0, 0.0}});
    tw.add_attachment(P3{{500.0, 500.0, 50.0}});
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    // a NaN pull, as a diverged MoorDyn line would give: refused, and the pulls already set are kept
    tw.set_line_loads({P3{{0.0, 1000.0, -5000.0}}}, {P3{{500.0, 500.0, 50.0}}});
    std::string msg = erf_gtest::abort_message([&] { tw.set_line_loads({P3{{nan, 0.0, 0.0}}}, {P3{{500.0, 500.0, 50.0}}}); });
    EXPECT_NE(msg.find("L1_t3"), std::string::npos) << msg;
    EXPECT_NE(msg.find("erf.conductors.moordyn_cfl"), std::string::npos) << msg;
    EXPECT_EQ(tw.line_forces()[0][1], Real(1000.0));
    // two pulls for one attachment
    msg = erf_gtest::abort_message([&] {
        tw.set_line_loads({P3{{0.0, 0.0, 0.0}}, P3{{0.0, 0.0, 0.0}}}, {P3{{500.0, 500.0, 50.0}}, P3{{500.0, 500.0, 50.0}}});
    });
    EXPECT_NE(msg.find("2 line pulls for 1 attachments"), std::string::npos) << msg;
}

TEST(Tower, AFlatOrBadlyAimedTowerIsRefusedNamingIt)
{
    const TowerType t = lattice();
    std::string msg = erf_gtest::abort_message([&] { Tower tw("T9", t, P3{{0.0, 0.0, 0.0}}, Real(0.0), P3{{0.0, 1.0, 0.0}}); });
    EXPECT_NE(msg.find("T9"), std::string::npos) << msg;
    EXPECT_NE(msg.find("cross-arm height"), std::string::npos) << msg;
    msg = erf_gtest::abort_message([&] { Tower tw("T9", t, P3{{0.0, 0.0, 0.0}}, Real(30.0), P3{{0.0, 0.6, 0.8}}); });
    EXPECT_NE(msg.find("horizontal unit vector"), std::string::npos) << msg;
    TowerType bad = t;
    bad.solidity = 0.0;
    msg = erf_gtest::abort_message([&] { Tower tw("T9", bad, P3{{0.0, 0.0, 0.0}}, Real(30.0), P3{{0.0, 1.0, 0.0}}); });
    EXPECT_NE(msg.find("erf.conductors.lattice.solidity"), std::string::npos) << msg;
    // a cross-arm no longer than the shaft is wide at the middle of its face (1.59 m at 29.4 m), though longer
    // than the top width, has no part outside the shaft
    TowerType stub = t;
    stub.arm_length = 1.55;
    msg = erf_gtest::abort_message([&] { Tower tw("T9", stub, P3{{0.0, 0.0, 0.0}}, Real(30.0), P3{{0.0, 1.0, 0.0}}); });
    EXPECT_NE(msg.find("erf.conductors.lattice.arm_length"), std::string::npos) << msg;
    EXPECT_NE(msg.find("arm_outside_shaft"), std::string::npos) << msg;
    // a cross-arm whose face's middle would lie below the base
    msg = erf_gtest::abort_message([&] { Tower tw("T9", t, P3{{0.0, 0.0, 0.0}}, Real(0.5), P3{{0.0, 1.0, 0.0}}); });
    EXPECT_NE(msg.find("erf.conductors.lattice.arm_depth"), std::string::npos) << msg;
    EXPECT_TRUE(erf_gtest::abort_message([&] { Tower tw("T9", t, P3{{0.0, 0.0, 0.0}}, Real(30.0), P3{{0.0, 1.0, 0.0}}); }).empty());
}
