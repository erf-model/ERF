// The frame model against the closed-form static response of beams. The frame works in double precision whatever
// amrex::Real is, so the tolerances here do not depend on ERF's precision.
//
// DirectionCosines.FollowSubDynsConvention: vertical, downward and horizontal elements get SubDyn's
//   local axes, every frame is orthonormal and right-handed, and the spin turns x towards y.
// BeamElement.IsSymmetricAndRigidMotionsCostNothing: an inclined Timoshenko element's 12 x 12
//   stiffness is symmetric and gives no force for the six rigid-body motions.
// Frame.ACantileverBendsStretchesAndTwistsAsBeamTheorySays: tip loads on a horizontal cantilever give
//   P L^3/(3 E I) (+ P L/(kappa G A) with shear deformation) in both bending planes, P L/(E A) and
//   T L/(G Jt), for one and for five elements per member.
// Frame.AVerticalCantileverCarriesItsWeightAndAMass: the base reaction equals the weight of the
//   member and of a mass at the top, the top sinks by rho g L^2/(2 E) + m g L/(E A), and the mass's
//   offset gives the base its moment.
// Frame.AHorizontalCantileverSagsUnderItsWeightAsBeamTheorySays: a cantilever under its own weight
//   sags w L^4/(8 E I) (+ w L^2/(2 kappa G A)) and turns w L^3/(6 E I) at its tip, and its base takes
//   the moment w L^2/2: the sign and size of SubDyn's gravity end moments.
// Frame.ASpringBaseAddsItsCompliance: on a 6 x 6 spring base with an x-translation, y-rotation
//   coupling, a tip load moves the tip by the beam's bending plus the base's translation and tilt
//   from the spring's own 2 x 2 solve.
// Frame.AFixedPortalSwaysAsSlopeDeflectionGives: a portal frame's sway under a lateral load equals
//   P h^2 (2a + 3b)/(12 a (a + 6b)), a = E Ic/h, b = E Ib/W.
// Frame.ATriangulatedTrussCarriesTheStaticallyDeterminateForces: slender members joined rigidly carry
//   the pin-jointed truss's axial forces and reactions.
// Frame.TheReactionsBalanceEveryLoadAndTheWeight: on the SubDyn lattice tower with a spring base,
//   the support reactions balance arbitrary joint loads and gravity in force and in moment.
// Frame.AMechanismIsRefusedNamingTheFreeDegreeOfFreedom: a tower free to spin about its axis, and a
//   frame without supports, are refused with the degree of freedom or the missing support named.

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "ERF_Frame.H"

using namespace erf_towers;

namespace {

constexpr double E_steel = 2.0e11;
constexpr double G_steel = 7.7e10;

FrameSection arbitrary (int id, double A, double Ixx, double Iyy, double Jt, double Asx, double Asy, double rho = 7850.0)
{
    FrameSection s;
    s.id = id;
    s.shape = SectionShape::Arbitrary;
    s.E = E_steel;
    s.G = G_steel;
    s.rho = rho;
    s.A = A;
    s.Asx = Asx;
    s.Asy = Asy;
    s.Ixx = Ixx;
    s.Iyy = Iyy;
    s.J0 = Ixx + Iyy;
    s.Jt = Jt;
    return s;
}

FrameJoint joint (int id, double x, double y, double z) { FrameJoint j; j.id = id; j.x = {{x, y, z}}; return j; }

FrameMember member (int id, int a, int b, int section, double spin = 0.0)
{
    FrameMember m;
    m.id = id;
    m.joint_a = a;
    m.joint_b = b;
    m.section = section;
    m.shape = SectionShape::Arbitrary;
    m.spin = spin;
    return m;
}

FrameSupport fixed_support (int j) { FrameSupport s; s.joint = j; return s; }

std::unique_ptr<Frame> build (const FrameInputs& in)
{
    std::string err;
    auto f = Frame::create(in, err);
    EXPECT_TRUE(f) << err;
    return f;
}

/** The 6 loads of joint id set to load, all others zero. */
std::vector<double> load_at (const Frame& f, int id, const std::array<double,6>& load)
{
    std::vector<double> l(f.num_dofs(), 0.0);
    const std::size_t n = f.node_of_joint(id);
    for (std::size_t d = 0; d < 6; ++d) { l[6 * n + d] = load[d]; }
    return l;
}

double dof (const Frame& f, const FrameSolution& s, int id, std::size_t d) { return s.displacement[6 * f.node_of_joint(id) + d]; }

/** Column j (local axis j) of a direction cosine matrix. */
std::array<double,3> col (const std::array<double,9>& dc, int j)
{
    const auto u = static_cast<std::size_t>(j);
    return {{dc[u], dc[3 + u], dc[6 + u]}};
}

} // namespace

TEST(DirectionCosines, FollowSubDynsConvention)
{
    const auto up = direction_cosines({{0, 0, 0}}, {{0, 0, 5}}, 0.0);
    const std::array<double,9> identity{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
    for (std::size_t i = 0; i < 9; ++i) { EXPECT_DOUBLE_EQ(up[i], identity[i]); }
    const auto down = direction_cosines({{0, 0, 5}}, {{0, 0, 0}}, 0.0);
    const std::array<double,9> flip{{1, 0, 0, 0, -1, 0, 0, 0, -1}};
    for (std::size_t i = 0; i < 9; ++i) { EXPECT_DOUBLE_EQ(down[i], flip[i]); }
    // along +x: local x = -y, local y = -z, local z = +x
    const auto ax = direction_cosines({{0, 0, 0}}, {{4, 0, 0}}, 0.0);
    const std::array<double,3> ex{{0, -1, 0}}, ey{{0, 0, -1}}, ez{{1, 0, 0}};
    for (int i = 0; i < 3; ++i) {
        EXPECT_NEAR(col(ax, 0)[static_cast<std::size_t>(i)], ex[static_cast<std::size_t>(i)], 1e-15);
        EXPECT_NEAR(col(ax, 1)[static_cast<std::size_t>(i)], ey[static_cast<std::size_t>(i)], 1e-15);
        EXPECT_NEAR(col(ax, 2)[static_cast<std::size_t>(i)], ez[static_cast<std::size_t>(i)], 1e-15);
    }
    // an inclined member: orthonormal, right-handed, z along the member; the spin turns x towards y
    const std::array<double,3> a{{1.0, -2.0, 0.5}}, b{{3.5, 1.0, 7.0}};
    const auto d0 = direction_cosines(a, b, 0.0);
    const double s = 0.7;
    const auto ds = direction_cosines(a, b, s);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            double dot = 0.0;
            for (int k = 0; k < 3; ++k) { dot += col(ds, i)[static_cast<std::size_t>(k)] * col(ds, j)[static_cast<std::size_t>(k)]; }
            EXPECT_NEAR(dot, i == j ? 1.0 : 0.0, 1e-14);
        }
    }
    const auto x = col(ds, 0), y = col(ds, 1), z = col(ds, 2);
    EXPECT_NEAR(x[1] * y[2] - x[2] * y[1], z[0], 1e-14);
    EXPECT_NEAR(x[2] * y[0] - x[0] * y[2], z[1], 1e-14);
    EXPECT_NEAR(x[0] * y[1] - x[1] * y[0], z[2], 1e-14);
    const double L = std::sqrt(2.5 * 2.5 + 3.0 * 3.0 + 6.5 * 6.5);
    EXPECT_NEAR(z[0], 2.5 / L, 1e-14);
    for (int k = 0; k < 3; ++k) {
        const auto u = static_cast<std::size_t>(k);
        EXPECT_NEAR(x[u], std::cos(s) * col(d0, 0)[u] + std::sin(s) * col(d0, 1)[u], 1e-14);
    }
}

TEST(BeamElement, IsSymmetricAndRigidMotionsCostNothing)
{
    const FrameSection sec = arbitrary(1, 4.0e-3, 1.2e-5, 0.6e-5, 4.0e-7, 2.0e-3, 2.5e-3);
    const BeamProperties p = beam_properties(sec, BeamTheory::Timoshenko);
    const std::array<double,3> a{{1.0, -2.0, 0.5}}, b{{3.5, 1.0, 7.0}};
    const double L = std::sqrt(2.5 * 2.5 + 3.0 * 3.0 + 6.5 * 6.5);
    const auto k = beam_stiffness(p, L, direction_cosines(a, b, 0.4));
    double kmax = 0.0;
    for (const double v : k) { kmax = std::max(kmax, std::abs(v)); }
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            EXPECT_NEAR(k[static_cast<std::size_t>(12 * i + j)], k[static_cast<std::size_t>(12 * j + i)], 1e-12 * kmax);
        }
    }
    // three translations and three rotations theta about the origin: u = theta x r at each node
    for (int m = 0; m < 6; ++m) {
        std::array<double,12> u{};
        for (int n = 0; n < 2; ++n) {
            const auto& r = (n == 0) ? a : b;
            if (m < 3) { u[static_cast<std::size_t>(6 * n + m)] = 1.0; continue; }
            std::array<double,3> th{};
            th[static_cast<std::size_t>(m - 3)] = 1.0;
            u[static_cast<std::size_t>(6 * n + 0)] = th[1] * r[2] - th[2] * r[1];
            u[static_cast<std::size_t>(6 * n + 1)] = th[2] * r[0] - th[0] * r[2];
            u[static_cast<std::size_t>(6 * n + 2)] = th[0] * r[1] - th[1] * r[0];
            for (int d = 0; d < 3; ++d) { u[static_cast<std::size_t>(6 * n + 3 + d)] = th[static_cast<std::size_t>(d)]; }
        }
        for (int i = 0; i < 12; ++i) {
            double f = 0.0;
            for (int j = 0; j < 12; ++j) { f += k[static_cast<std::size_t>(12 * i + j)] * u[static_cast<std::size_t>(j)]; }
            EXPECT_NEAR(f, 0.0, 1e-11 * kmax * 10.0) << "rigid motion " << m << ", row " << i;
        }
    }
}

TEST(Frame, ACantileverBendsStretchesAndTwistsAsBeamTheorySays)
{
    const double L = 6.0, A = 4.0e-3, Ixx = 1.2e-5, Iyy = 0.6e-5, Jt = 4.0e-7, Asx = 2.0e-3, Asy = 2.5e-3;
    for (const BeamTheory theory : {BeamTheory::EulerBernoulli, BeamTheory::Timoshenko}) {
        for (const int ndiv : {1, 5}) {
            FrameInputs in;
            in.theory = theory;
            in.divisions = ndiv;
            in.joints = {joint(1, 0, 0, 0), joint(2, L, 0, 0)};
            in.sections = {arbitrary(1, A, Ixx, Iyy, Jt, Asx, Asy)};
            in.members = {member(1, 1, 2, 1)};
            in.supports = {fixed_support(1)};
            auto f = build(in);
            ASSERT_TRUE(f);
            const bool shear = (theory == BeamTheory::Timoshenko);
            const double P = 1.0e3, T = 50.0;
            // along +x the local x axis is -y and the local y axis is -z: a z load bends about local x (Ixx,
            // shear along local y), a y load about local y (Iyy, shear along local x)
            const double dz = P * L * L * L / (3.0 * E_steel * Ixx) + (shear ? P * L / (Asy * G_steel) : 0.0);
            const double dy = P * L * L * L / (3.0 * E_steel * Iyy) + (shear ? P * L / (Asx * G_steel) : 0.0);
            const std::string tag = std::string(shear ? "Timoshenko" : "Euler-Bernoulli") + ", NDiv " + std::to_string(ndiv);
            auto s = f->solve(load_at(*f, 2, {{0, 0, P, 0, 0, 0}}), 0.0);
            EXPECT_NEAR(dof(*f, s, 2, 2), dz, 1e-10 * dz) << tag;
            EXPECT_NEAR(dof(*f, s, 2, 4), -P * L * L / (2.0 * E_steel * Ixx), 1e-10 * P * L * L / (E_steel * Ixx)) << tag;
            s = f->solve(load_at(*f, 2, {{0, P, 0, 0, 0, 0}}), 0.0);
            EXPECT_NEAR(dof(*f, s, 2, 1), dy, 1e-10 * dy) << tag;
            s = f->solve(load_at(*f, 2, {{P, 0, 0, 0, 0, 0}}), 0.0);
            EXPECT_NEAR(dof(*f, s, 2, 0), P * L / (E_steel * A), 1e-10 * P * L / (E_steel * A)) << tag;
            s = f->solve(load_at(*f, 2, {{0, 0, 0, T, 0, 0}}), 0.0);
            EXPECT_NEAR(dof(*f, s, 2, 3), T * L / (G_steel * Jt), 1e-10 * T * L / (G_steel * Jt)) << tag;
            // the base reaction holds the tip load: a force -P and the moment of P about the base
            s = f->solve(load_at(*f, 2, {{0, 0, P, 0, 0, 0}}), 0.0);
            EXPECT_NEAR(s.reaction[0][2], -P, 1e-9 * P) << tag;
            EXPECT_NEAR(s.reaction[0][4], P * L, 1e-9 * P * L) << tag;
        }
    }
}

TEST(Frame, AVerticalCantileverCarriesItsWeightAndAMass)
{
    const double L = 12.0, A = 4.0e-3, rho = 7850.0, g = 9.80665, m = 500.0, ox = 0.3;
    for (const int ndiv : {1, 4}) {
        FrameInputs in;
        in.divisions = ndiv;
        in.joints = {joint(1, 0, 0, 0), joint(2, 0, 0, L)};
        in.sections = {arbitrary(1, A, 1.2e-5, 0.6e-5, 4.0e-7, 2.0e-3, 2.0e-3, rho)};
        in.members = {member(1, 1, 2, 1)};
        in.supports = {fixed_support(1)};
        FrameMass cm;
        cm.joint = 2;
        cm.mass = m;
        cm.offset = {{ox, 0.0, 0.0}};
        in.masses = {cm};
        auto f = build(in);
        ASSERT_TRUE(f);
        EXPECT_NEAR(f->total_mass(), rho * A * L + m, 1e-9 * (rho * A * L + m));
        const auto s = f->solve({}, g);
        const double W = (rho * A * L + m) * g;
        EXPECT_NEAR(s.reaction[0][2], W, 1e-10 * W) << "NDiv " << ndiv;
        EXPECT_NEAR(s.reaction[0][0], 0.0, 1e-10 * W);
        // the mass's weight at x = ox: its moment about the base is (0, ox m g, 0), which the base balances
        EXPECT_NEAR(s.reaction[0][4], -ox * m * g, 1e-9 * W);
        // the top sinks by the bar's own weight and the mass's
        const double sink = rho * g * L * L / (2.0 * E_steel) + m * g * L / (E_steel * A);
        EXPECT_NEAR(dof(*f, s, 2, 2), -sink, 1e-9 * sink) << "NDiv " << ndiv;
    }
}

TEST(Frame, AHorizontalCantileverSagsUnderItsWeightAsBeamTheorySays)
{
    // consistent element loads make the nodal displacements of a uniformly loaded beam exact:
    // tip sag w L^4/(8 E I) (+ w L^2/(2 kappa G A) with shear deformation), tip rotation w L^3/(6 E I)
    const double L = 6.0, A = 4.0e-3, Ixx = 1.2e-5, Asy = 2.5e-3, rho = 7850.0, g = 9.80665, w = rho * A * g;
    for (const BeamTheory theory : {BeamTheory::EulerBernoulli, BeamTheory::Timoshenko}) {
        for (const int ndiv : {1, 5}) {
            FrameInputs in;
            in.theory = theory;
            in.divisions = ndiv;
            in.joints = {joint(1, 0, 0, 0), joint(2, L, 0, 0)};
            in.sections = {arbitrary(1, A, Ixx, 0.6e-5, 4.0e-7, 2.0e-3, Asy, rho)};
            in.members = {member(1, 1, 2, 1)};
            in.supports = {fixed_support(1)};
            auto f = build(in);
            ASSERT_TRUE(f);
            const auto s = f->solve({}, g);
            const bool shear = (theory == BeamTheory::Timoshenko);
            const double sag = w * std::pow(L, 4) / (8.0 * E_steel * Ixx) + (shear ? w * L * L / (2.0 * Asy * G_steel) : 0.0);
            const double tilt = w * std::pow(L, 3) / (6.0 * E_steel * Ixx);
            const std::string tag = std::string(shear ? "Timoshenko" : "Euler-Bernoulli") + ", NDiv " + std::to_string(ndiv);
            EXPECT_NEAR(dof(*f, s, 2, 2), -sag, 1e-10 * sag) << tag;
            // the tip slopes down towards +x: a positive rotation about y
            EXPECT_NEAR(dof(*f, s, 2, 4), tilt, 1e-10 * tilt) << tag;
            EXPECT_NEAR(s.reaction[0][4], -0.5 * w * L * L, 1e-9 * w * L * L) << tag;
        }
    }
}

TEST(Frame, ASpringBaseAddsItsCompliance)
{
    const double L = 10.0, Iyy = 0.6e-5, P = 2.0e3;
    const double kx = 2.0e7, kry = 3.0e7, kxty = 4.0e6;
    FrameInputs in;
    in.joints = {joint(1, 0, 0, 0), joint(2, 0, 0, L)};
    in.sections = {arbitrary(1, 4.0e-3, 1.2e-5, Iyy, 4.0e-7, 2.0e-3, 2.0e-3)};
    in.members = {member(1, 1, 2, 1)};
    FrameSupport s;
    s.joint = 1;
    s.fixed = {{false, false, false, false, false, false}};
    s.stiffness = {};
    s.stiffness[0] = kx;      // Kxx
    s.stiffness[2] = 5.0e7;   // Kyy
    s.stiffness[5] = 9.0e7;   // Kzz
    s.stiffness[9] = 6.0e7;   // Ktxtx
    s.stiffness[10] = kxty;   // Kxty: x translation with rotation about y
    s.stiffness[14] = kry;    // Ktyty
    s.stiffness[20] = 7.0e7;  // Ktztz
    s.ssi_file = "spring";
    in.supports = {s};
    auto f = build(in);
    ASSERT_TRUE(f);
    const auto sol = f->solve(load_at(*f, 2, {{P, 0, 0, 0, 0, 0}}), 0.0);
    // the base takes the force P and the moment P L about y: solve the spring's x / theta_y block
    const double fx = P, my = P * L, det = kx * kry - kxty * kxty;
    const double ux = (kry * fx - kxty * my) / det, ty = (kx * my - kxty * fx) / det;
    const double tip = P * L * L * L / (3.0 * E_steel * Iyy) + ux + ty * L;
    EXPECT_NEAR(dof(*f, sol, 1, 0), ux, 1e-9 * std::abs(ux));
    EXPECT_NEAR(dof(*f, sol, 1, 4), ty, 1e-9 * std::abs(ty));
    EXPECT_NEAR(dof(*f, sol, 2, 0), tip, 1e-9 * tip);
    // the reaction is the spring's force on the frame
    EXPECT_NEAR(sol.reaction[0][0], -P, 1e-9 * P);
    EXPECT_NEAR(sol.reaction[0][4], -P * L, 1e-9 * P * L);
}

TEST(Frame, AFixedPortalSwaysAsSlopeDeflectionGives)
{
    const double h = 4.0, W = 6.0, Ic = 8.0e-5, Ib = 2.0e-4, P = 1.0e4, A = 2.0;
    for (const int ndiv : {1, 4}) {
        FrameInputs in;
        in.divisions = ndiv;
        in.joints = {joint(1, 0, 0, 0), joint(2, W, 0, 0), joint(3, 0, 0, h), joint(4, W, 0, h)};
        // columns bend about global y (Iyy); the beam along x bends about its local x axis, global -y (Ixx)
        in.sections = {arbitrary(1, A, 1.0e-3, Ic, 1.0e-4, 1.0, 1.0), arbitrary(2, A, Ib, 1.0e-3, 1.0e-4, 1.0, 1.0)};
        in.members = {member(1, 1, 3, 1), member(2, 2, 4, 1), member(3, 3, 4, 2)};
        in.supports = {fixed_support(1), fixed_support(2)};
        auto f = build(in);
        ASSERT_TRUE(f);
        const auto s = f->solve(load_at(*f, 3, {{P, 0, 0, 0, 0, 0}}), 0.0);
        const double a = E_steel * Ic / h, b = E_steel * Ib / W;
        const double sway = P * h * h * (2.0 * a + 3.0 * b) / (12.0 * a * (a + 6.0 * b));
        // the members are 10^4 times stiffer axially than in bending: their stretch moves the sway by ~1e-5
        EXPECT_NEAR(dof(*f, s, 3, 0), sway, 1e-4 * sway) << "NDiv " << ndiv;
        EXPECT_NEAR(dof(*f, s, 4, 0), sway, 1e-4 * sway) << "NDiv " << ndiv;
        EXPECT_NEAR(s.reaction[0][0] + s.reaction[1][0], -P, 1e-9 * P);
    }
}

TEST(Frame, ATriangulatedTrussCarriesTheStaticallyDeterminateForces)
{
    // joints A (0,0,0) and B (4,0,0) on the ground, C (2,0,3) loaded; in-plane bending (about global y,
    // the members' local x axis) made 10^8 times softer than stretching, so the joints act as pins
    const double P = 1.0e4;
    FrameInputs in;
    in.joints = {joint(1, 0, 0, 0), joint(2, 4, 0, 0), joint(3, 2, 0, 3)};
    in.sections = {arbitrary(1, 1.0e-2, 1.0e-12, 1.0e-5, 1.0e-6, 5.0e-3, 5.0e-3)};
    in.members = {member(1, 1, 2, 1), member(2, 1, 3, 1), member(3, 2, 3, 1)};
    // out of the plane every joint is held (y translation, rotations about x and z); A is pinned, B on a roller
    auto planar = [] (int j, bool hold_x, bool hold_z) {
        FrameSupport s;
        s.joint = j;
        s.fixed = {{hold_x, true, hold_z, true, false, true}};
        return s;
    };
    in.supports = {planar(1, true, true), planar(2, false, true), planar(3, false, false)};
    auto f = build(in);
    ASSERT_TRUE(f);
    const auto s = f->solve(load_at(*f, 3, {{0, 0, -P, 0, 0, 0}}), 0.0);
    // member forces: N = the local z force on the element at its second node (tension positive)
    const double diag = -P * std::sqrt(13.0) / 6.0, bottom = P / 3.0;
    EXPECT_NEAR(s.element_force[0][8], bottom, 1e-6 * P);
    EXPECT_NEAR(s.element_force[1][8], diag, 1e-6 * P);
    EXPECT_NEAR(s.element_force[2][8], diag, 1e-6 * P);
    EXPECT_NEAR(s.reaction[0][2], 0.5 * P, 1e-6 * P);
    EXPECT_NEAR(s.reaction[1][2], 0.5 * P, 1e-6 * P);
    EXPECT_NEAR(s.reaction[0][0], 0.0, 1e-6 * P);
}

TEST(Frame, TheReactionsBalanceEveryLoadAndTheWeight)
{
    FrameInputs in;
    const std::string err = read_subdyn(std::string(ERF_FRAME_TEST_FILES) + "/caseB/towerB.dat", in);
    ASSERT_TRUE(err.empty()) << err;
    FrameMass cm;
    cm.joint = 18;
    cm.mass = 800.0;
    cm.offset = {{0.2, -0.1, -0.5}};
    in.masses = {cm};
    auto f = build(in);
    ASSERT_TRUE(f);
    const double g = 9.80665;
    std::vector<double> loads(f->num_dofs(), 0.0);
    for (std::size_t i = 0; i < loads.size(); ++i) { loads[i] = 100.0 * std::sin(1.7 * static_cast<double>(i) + 0.3); }
    const auto s = f->solve(loads, g);
    // every external load: the given ones, the members' weight as element loads, the mass's weight
    std::vector<double> ext = loads;
    for (const auto& e : f->elements()) {
        const auto fg = beam_gravity_load(e.prop, e.length, e.dc, g);
        for (int i = 0; i < 12; ++i) {
            ext[6 * (i < 6 ? e.node_a : e.node_b) + static_cast<std::size_t>(i % 6)] += fg[static_cast<std::size_t>(i)];
        }
    }
    const std::size_t nm = f->node_of_joint(18);
    ext[6 * nm + 2] -= cm.mass * g;
    ext[6 * nm + 3] += cm.offset[1] * (-cm.mass * g);
    ext[6 * nm + 4] -= cm.offset[0] * (-cm.mass * g);
    for (std::size_t k = 0; k < in.supports.size(); ++k) {
        const std::size_t n = f->node_of_joint(in.supports[k].joint);
        for (std::size_t d = 0; d < 6; ++d) { ext[6 * n + d] += s.reaction[k][d]; }
    }
    std::array<double,6> total{};
    double scale = 0.0;
    for (std::size_t n = 0; n < f->num_nodes(); ++n) {
        const auto& r = f->node_position(n);
        const double* q = &ext[6 * n];
        for (int d = 0; d < 3; ++d) { total[static_cast<std::size_t>(d)] += q[d]; total[static_cast<std::size_t>(d + 3)] += q[d + 3]; }
        total[3] += r[1] * q[2] - r[2] * q[1];
        total[4] += r[2] * q[0] - r[0] * q[2];
        total[5] += r[0] * q[1] - r[1] * q[0];
        for (int d = 0; d < 6; ++d) { scale = std::max(scale, std::abs(q[d])); }
    }
    for (std::size_t d = 0; d < 6; ++d) { EXPECT_NEAR(total[d], 0.0, 1e-9 * scale * 30.0) << "component " << d; }
}

TEST(Frame, AMechanismIsRefusedNamingTheFreeDegreeOfFreedom)
{
    FrameInputs in;
    in.file = "spinning";
    in.joints = {joint(1, 0, 0, 0), joint(2, 0, 0, 10)};
    in.sections = {arbitrary(1, 4.0e-3, 1.2e-5, 0.6e-5, 4.0e-7, 2.0e-3, 2.0e-3)};
    in.members = {member(1, 1, 2, 1)};
    FrameSupport s = fixed_support(1);
    s.fixed[5] = false;   // nothing holds the twist about z
    in.supports = {s};
    std::string err;
    EXPECT_FALSE(Frame::create(in, err));
    EXPECT_NE(err.find("mechanism"), std::string::npos) << err;
    EXPECT_NE(err.find("rotation about z"), std::string::npos) << err;
    in.supports.clear();
    EXPECT_FALSE(Frame::create(in, err));
    EXPECT_NE(err.find("base reaction"), std::string::npos) << err;
}
