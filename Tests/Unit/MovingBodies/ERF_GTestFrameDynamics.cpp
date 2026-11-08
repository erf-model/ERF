// The frame model's mass, natural modes and Newmark time response. The frame works in double
// precision whatever amrex::Real is, so the tolerances do not depend on ERF's precision.
//
// BeamMass.HoldsTheElementsMassAndTwistInertia: an inclined element's consistent mass moves rho A L
//   in every translation and rho J0 L in a twist about its axis, and is symmetric positive definite.
// RigidBodyMass.GivesTheKineticEnergyOfAnOffsetBody: for any joint velocity and rotation rate, the
//   6 x 6 concentrated mass gives the kinetic energy of the rigid body whose centre is offset.
// FrameModes.ACantileverRingsAtItsBeamTheoryFrequencies: the first two bending frequencies in both
//   planes, the axial and the torsional frequencies of a cantilever.
// FrameModes.ASpringBaseRingsWithItsOwnMass: a stiff, nearly massless post on a spring base with
//   masses rings at the base's frequencies, the coupled translation-rocking pair included.
// FrameModes.TheSturmCountIsTheNumberOfModesBelowTheShift: the negative pivots of K - sigma M count
//   the modes below sigma, the check frame_modes relies on.
// FrameModes.AreMassOrthonormalAndSolveTheEigenproblem: on the lattice tower with masses and a
//   spring base, phi_i^T M phi_j = delta_ij and K phi = omega^2 M phi.
// FrameModes.AreSubDynsFullFrequencies: the 20 lowest frequencies of three SubDyn lattice towers
//   (Euler-Bernoulli; Timoshenko with a spring base and a cluster of 8 close modes; and with
//   concentrated masses and a spring mass) equal SubDyn's to the 7 digits it prints.
// FrameDynamics.AModeSwingsAtNewmarksFrequencyAndKeepsItsEnergy: released in its first mode, a
//   cantilever's tip follows cos(n 2 atan(omega h/2)) and the energy stays constant.
// FrameDynamics.RayleighDampingDecaysAModeAsTheOscillatorItIs: with Rayleigh damping, the mode follows
//   Newmark's recursion for the damped one-degree-of-freedom oscillator of its frequency and ratio.
// FrameDynamics.AChangingStepFollowsTheOscillatorOfEachStep: one dynamics stepping h1, h2, h1, ... and
//   another on the same frame stepping h2 throughout (sharing its factorisation) each follow the modal
//   oscillator's Newmark recursion with their own steps, and share the factorisation of a step they both take.
// FrameDynamics.AHeldLoadDoesTheWorkTheEnergyGains: under loads held over each step (step_held), changing
//   from step to step, the loads' work f (u' - u) is the gain of kinetic and strain energy plus the
//   Rayleigh damping's work h vm^T C vm (vm = (u' - u)/h) to round-off, and the end acceleration
//   satisfies M a = f - C v - K u; Newmark's step to the same loads misses the balance.
// FrameDynamics.AStepLoadSettlesAtTheStaticResponse: a damped tower under a sudden load and gravity
//   comes to rest at the static solution.
// FrameDynamics.TheStateRestoresTheSameMotion: a state saved and restored continues bit for bit.
// FrameDynamics.RefusesBadStepsStatesAndAMasslessNode: a step <= 0 aborts, a state of the wrong
//   size or with a non-finite value is refused, a node without mass is refused at the start.

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "ERF_FrameDynamics.H"
#include "ERF_GTestThrowOnAbort.H"

using namespace erf_towers;

namespace {

constexpr double pi = 3.14159265358979323846;
constexpr double E_steel = 2.0e11, G_steel = 7.7e10, rho_steel = 7850.0;

FrameSection arbitrary (double A, double Ixx, double Iyy, double J0, double Jt, double rho = rho_steel)
{
    FrameSection s;
    s.id = 1;
    s.shape = SectionShape::Arbitrary;
    s.E = E_steel;
    s.G = G_steel;
    s.rho = rho;
    s.A = A;
    s.Asx = 0.5 * A;
    s.Asy = 0.5 * A;
    s.Ixx = Ixx;
    s.Iyy = Iyy;
    s.J0 = J0;
    s.Jt = Jt;
    return s;
}

/** A vertical cantilever of length L, fixed at z = 0, in ndiv elements. */
FrameInputs cantilever (double L, int ndiv, const FrameSection& s)
{
    FrameInputs in;
    in.divisions = ndiv;
    FrameJoint a, b;
    a.id = 1;
    b.id = 2;
    b.x = {{0.0, 0.0, L}};
    in.joints = {a, b};
    in.sections = {s};
    FrameMember m;
    m.id = 1;
    m.joint_a = 1;
    m.joint_b = 2;
    m.section = 1;
    in.members = {m};
    FrameSupport base;
    base.joint = 1;
    in.supports = {base};
    return in;
}

std::unique_ptr<Frame> build (const FrameInputs& in)
{
    std::string err;
    auto f = Frame::create(in, err);
    EXPECT_TRUE(f) << err;
    return f;
}

std::string case_file (const std::string& c, const std::string& suffix)
{
    return std::string(ERF_FRAME_TEST_FILES) + "/case" + c + "/tower" + c + suffix;
}

std::vector<double> subdyn_frequencies (const std::string& c)
{
    std::ifstream f(case_file(c, "_frequencies.txt"));
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

double dot (const std::vector<double>& a, const std::vector<double>& b)
{
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) { s += a[i] * b[i]; }
    return s;
}

/** The tip's x displacement (degree of freedom 0 of joint 2). */
double tip_x (const Frame& f, const std::vector<double>& u) { return u[6 * f.node_of_joint(2)]; }

} // namespace

TEST(BeamMass, HoldsTheElementsMassAndTwistInertia)
{
    const FrameSection s = arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7);
    const BeamProperties p = beam_properties(s, BeamTheory::Timoshenko);
    const std::array<double,3> a{{1.0, -2.0, 0.5}}, b{{3.5, 1.0, 7.0}};
    const double L = std::sqrt(2.5 * 2.5 + 3.0 * 3.0 + 6.5 * 6.5);
    const auto dc = direction_cosines(a, b, 0.4);
    const auto m = beam_mass(p, L, dc);
    auto quad = [&m] (const std::array<double,12>& x) {
        double q = 0.0;
        for (std::size_t i = 0; i < 12; ++i) {
            for (std::size_t j = 0; j < 12; ++j) { q += x[i] * m[12 * i + j] * x[j]; }
        }
        return q;
    };
    for (int d = 0; d < 3; ++d) {
        std::array<double,12> t{};
        t[static_cast<std::size_t>(d)] = 1.0;
        t[static_cast<std::size_t>(6 + d)] = 1.0;
        EXPECT_NEAR(quad(t), p.rho * p.A * L, 1e-12 * p.rho * p.A * L) << "translation " << d;
    }
    std::array<double,12> twist{};
    for (int d = 0; d < 3; ++d) {
        twist[static_cast<std::size_t>(3 + d)] = dc[static_cast<std::size_t>(3 * d + 2)];
        twist[static_cast<std::size_t>(9 + d)] = dc[static_cast<std::size_t>(3 * d + 2)];
    }
    EXPECT_NEAR(quad(twist), p.rho * p.J0 * L, 1e-12 * p.rho * p.J0 * L);
    double mmax = 0.0;
    for (const double v : m) { mmax = std::max(mmax, std::abs(v)); }
    for (int i = 0; i < 12; ++i) {
        for (int j = 0; j < 12; ++j) {
            EXPECT_NEAR(m[static_cast<std::size_t>(12 * i + j)], m[static_cast<std::size_t>(12 * j + i)], 1e-14 * mmax);
        }
    }
    DenseCholesky c;
    EXPECT_LT(c.factor(std::vector<double>(m.begin(), m.end()), 12), 0) << "the element mass must be positive definite";
}

TEST(RigidBodyMass, GivesTheKineticEnergyOfAnOffsetBody)
{
    FrameMass c;
    c.mass = 400.0;
    c.inertia = {{20.0, 30.0, 25.0, 2.0, -1.0, 1.5}};   // Jxx Jyy Jzz Jxy Jxz Jyz, tensor entries about the centre
    c.offset = {{0.1, -0.2, -1.5}};
    const auto m66 = rigid_body_mass(c);
    const std::array<std::array<double,6>,3> motions{{{{0.3, -0.7, 0.2, 0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0, 1.1, -0.4, 0.6}},
                                                       {{-0.5, 0.25, 0.9, 0.3, 0.8, -1.2}}}};
    for (const auto& x : motions) {
        // the centre moves with u + w x r; the body spins at w
        const std::array<double,3> u{{x[0], x[1], x[2]}}, w{{x[3], x[4], x[5]}};
        const auto& r = c.offset;
        const std::array<double,3> vg{{u[0] + w[1] * r[2] - w[2] * r[1],
                                       u[1] + w[2] * r[0] - w[0] * r[2],
                                       u[2] + w[0] * r[1] - w[1] * r[0]}};
        const double jx = c.inertia[0] * w[0] + c.inertia[3] * w[1] + c.inertia[4] * w[2];
        const double jy = c.inertia[3] * w[0] + c.inertia[1] * w[1] + c.inertia[5] * w[2];
        const double jz = c.inertia[4] * w[0] + c.inertia[5] * w[1] + c.inertia[2] * w[2];
        const double twice_ke = c.mass * (vg[0] * vg[0] + vg[1] * vg[1] + vg[2] * vg[2]) + w[0] * jx + w[1] * jy + w[2] * jz;
        double q = 0.0;
        for (std::size_t i = 0; i < 6; ++i) { for (std::size_t j = 0; j < 6; ++j) { q += x[i] * m66[6 * i + j] * x[j]; } }
        EXPECT_NEAR(q, twice_ke, 1e-12 * twice_ke);
    }
}

TEST(FrameModes, ACantileverRingsAtItsBeamTheoryFrequencies)
{
    const double L = 10.0, A = 4.0e-3, Ixx = 1.2e-5, Iyy = 0.6e-5, J0 = 1.8e-5, Jt = 4.0e-7;
    auto f = build(cantilever(L, 20, arbitrary(A, Ixx, Iyy, J0, Jt)));
    ASSERT_TRUE(f);
    FrameModes modes;
    // 20 modes: the torsion's odd harmonics (3, 5, 7, 9 times its fundamental) come before the axial mode
    const std::string err = frame_modes(*f, 20, modes);
    ASSERT_TRUE(err.empty()) << err;
    // Euler-Bernoulli cantilever: (beta L)^2 / (2 pi L^2) sqrt(E I/(rho A)); the consistent mass's rotary
    // inertia rho I lowers them by ~(beta r)^2/2 ~ 1e-4 here. A bar and a shaft fixed at one end: c/(4 L).
    const double b1 = 1.875104069, b2 = 4.694091133;
    auto bend = [&] (double beta, double I) { return beta * beta / (2.0 * pi * L * L) * std::sqrt(E_steel * I / (rho_steel * A)); };
    const std::vector<std::pair<double,double>> expected = {
        {bend(b1, Iyy), 2e-4}, {bend(b1, Ixx), 2e-4}, {bend(b2, Iyy), 5e-4}, {bend(b2, Ixx), 5e-4},
        {std::sqrt(E_steel / rho_steel) / (4.0 * L), 1e-3}, {std::sqrt(G_steel * Jt / (rho_steel * J0)) / (4.0 * L), 1e-3}};
    for (const auto& e : expected) {
        double best = 1e30;
        for (const double fr : modes.frequency) { best = std::min(best, std::abs(fr - e.first) / e.first); }
        std::string all;
        for (const double fr : modes.frequency) { all += " " + std::to_string(fr); }
        EXPECT_LT(best, e.second) << "no mode near " << e.first << " Hz; the modes:" << all;
    }
}

TEST(FrameModes, ASpringBaseRingsWithItsOwnMass)
{
    // a short, stiff and nearly massless post on a spring base with masses: its lowest modes are the
    // base's, the generalized eigenvalues of the spring's and the base mass's blocks
    FrameInputs in = cantilever(1.0, 1, arbitrary(1.0, 1.0, 1.0, 2.0, 1.0, 1.0e-6));
    FrameSupport& b = in.supports[0];
    b.fixed = {{false, false, false, false, false, false}};
    b.ssi_file = "spring";
    const double kx = 2.0e8, kty = 3.0e8, kxty = 1.0e7, kz = 5.0e8, ktz = 1.0e8;
    const double mx = 2.0e3, mty = 5.0e2, mxty = 1.0e2, mz = 2.0e3, mtz = 4.0e2;
    b.stiffness = {};
    b.stiffness[0] = kx; b.stiffness[2] = kx; b.stiffness[5] = kz; b.stiffness[9] = kty; b.stiffness[14] = kty;
    b.stiffness[20] = ktz; b.stiffness[10] = kxty; b.stiffness[7] = -kxty;    // Kxty: x with theta_y; Kytx: y with theta_x
    b.mass = {};
    b.mass[0] = mx; b.mass[2] = mx; b.mass[5] = mz; b.mass[9] = mty; b.mass[14] = mty; b.mass[20] = mtz;
    b.mass[10] = mxty; b.mass[7] = -mxty;
    auto f = build(in);
    ASSERT_TRUE(f);
    FrameModes modes;
    const std::string err = frame_modes(*f, 6, modes);
    ASSERT_TRUE(err.empty()) << err;
    // the x / theta_y block (and its mirror y / theta_x): det(K - w^2 M) = 0
    const double a = mx * mty - mxty * mxty, bb = -(kx * mty + kty * mx - 2.0 * kxty * mxty), c = kx * kty - kxty * kxty;
    const double disc = std::sqrt(bb * bb - 4.0 * a * c);
    std::vector<double> expected = {std::sqrt((-bb - disc) / (2.0 * a)), std::sqrt((-bb - disc) / (2.0 * a)),
                                    std::sqrt((-bb + disc) / (2.0 * a)), std::sqrt((-bb + disc) / (2.0 * a)),
                                    std::sqrt(kz / mz), std::sqrt(ktz / mtz)};
    for (double& e : expected) { e /= 2.0 * pi; }
    std::sort(expected.begin(), expected.end());
    for (std::size_t k = 0; k < 6; ++k) { EXPECT_NEAR(modes.frequency[k], expected[k], 1e-6 * expected[k]) << "mode " << k + 1; }
}

TEST(FrameModes, TheSturmCountIsTheNumberOfModesBelowTheShift)
{
    auto f = build(cantilever(10.0, 20, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7)));
    ASSERT_TRUE(f);
    FrameModes modes;
    ASSERT_TRUE(frame_modes(*f, 14, modes).empty());
    for (const std::size_t below : {std::size_t(1), std::size_t(5), std::size_t(12)}) {
        const double fs = 0.5 * (modes.frequency[below - 1] + modes.frequency[below]);
        const double sigma = std::pow(2.0 * pi * fs, 2);
        EXPECT_EQ(negative_pivots(f->assemble_free(1.0, -sigma), f->num_free_dofs()), below) << "shift at " << fs << " Hz";
    }
}

TEST(FrameModes, AreMassOrthonormalAndSolveTheEigenproblem)
{
    FrameInputs in;
    ASSERT_TRUE(read_subdyn(case_file("C", ".dat"), in).empty());
    auto f = build(in);
    ASSERT_TRUE(f);
    FrameModes modes;
    const std::string err = frame_modes(*f, 10, modes);
    ASSERT_TRUE(err.empty()) << err;
    std::vector<std::vector<double>> mphi;
    for (const auto& phi : modes.shape) { mphi.push_back(f->apply_mass(phi)); }
    for (std::size_t i = 0; i < 10; ++i) {
        for (std::size_t j = 0; j < 10; ++j) { EXPECT_NEAR(dot(modes.shape[i], mphi[j]), i == j ? 1.0 : 0.0, 1e-9) << i << "," << j; }
        const double w2 = std::pow(2.0 * pi * modes.frequency[i], 2);
        const std::vector<double> kphi = f->apply_stiffness(modes.shape[i]);
        double r2 = 0.0, k2 = 0.0;
        for (const std::size_t d : f->free_dofs()) { r2 += std::pow(kphi[d] - w2 * mphi[i][d], 2); k2 += kphi[d] * kphi[d]; }
        // the eigenvalues converge to 1e-12, the vectors to about its square root
        EXPECT_LT(std::sqrt(r2 / k2), 1e-5) << "mode " << i;
    }
}

TEST(FrameModes, AreSubDynsFullFrequencies)
{
    for (const std::string c : {"A", "B", "C"}) {
        FrameInputs in;
        const std::string rerr = read_subdyn(case_file(c, ".dat"), in);
        ASSERT_TRUE(rerr.empty()) << rerr;
        auto f = build(in);
        ASSERT_TRUE(f);
        const std::vector<double> theirs = subdyn_frequencies(c);
        ASSERT_GE(theirs.size(), 20u);
        FrameModes modes;
        const std::string err = frame_modes(*f, 20, modes);
        ASSERT_TRUE(err.empty()) << "case " << c << ": " << err;
        for (std::size_t k = 0; k < 20; ++k) {
            EXPECT_NEAR(modes.frequency[k], theirs[k], 1e-6 * theirs[k]) << "case " << c << ", mode " << k + 1;
        }
    }
}

TEST(FrameDynamics, AModeSwingsAtNewmarksFrequencyAndKeepsItsEnergy)
{
    auto f = build(cantilever(10.0, 10, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7)));
    ASSERT_TRUE(f);
    FrameModes modes;
    ASSERT_TRUE(frame_modes(*f, 1, modes).empty());
    const double w = 2.0 * pi * modes.frequency[0], h = (2.0 * pi / w) / 40.0;
    std::vector<double> u0 = modes.shape[0];
    for (double& x : u0) { x *= 0.01; }
    FrameDynamics dyn(*f, 0.0, 0.0);
    ASSERT_TRUE(dyn.start(u0, std::vector<double>(u0.size(), 0.0), {}, 0.0).empty());
    const double amp = tip_x(*f, u0), e0 = dyn.kinetic_energy() + dyn.strain_energy();
    const double wh = 2.0 * std::atan(0.5 * w * h);     // Newmark's phase advance per step
    for (int n = 1; n <= 120; ++n) {
        dyn.step(h, {}, 0.0);
        EXPECT_NEAR(tip_x(*f, dyn.displacement()), amp * std::cos(n * wh), 1e-5 * std::abs(amp)) << "step " << n;
        EXPECT_NEAR(dyn.kinetic_energy() + dyn.strain_energy(), e0, 1e-11 * e0) << "step " << n;
    }
    // Newmark's phase lags the exact one by 0.04 rad after 120 steps, 4000 times the tolerance near a zero
    // crossing: the check is not vacuous
    EXPECT_GT(std::abs(120 * w * h - 120 * wh), 0.01);
}

TEST(FrameDynamics, RayleighDampingDecaysAModeAsTheOscillatorItIs)
{
    auto f = build(cantilever(10.0, 10, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7)));
    ASSERT_TRUE(f);
    FrameModes modes;
    ASSERT_TRUE(frame_modes(*f, 3, modes).empty());
    double a0 = 0.0, a1 = 0.0;
    rayleigh_coefficients(modes.frequency[0], 0.02, modes.frequency[2], 0.05, a0, a1);
    const double w = 2.0 * pi * modes.frequency[0];
    const double zeta = 0.5 * (a0 / w + a1 * w);
    EXPECT_NEAR(zeta, 0.02, 1e-14);
    const double w3 = 2.0 * pi * modes.frequency[2];
    EXPECT_NEAR(0.5 * (a0 / w3 + a1 * w3), 0.05, 1e-14);
    const double h = (2.0 * pi / w) / 25.0;
    std::vector<double> u0 = modes.shape[0];
    for (double& x : u0) { x *= 0.01; }
    FrameDynamics dyn(*f, a0, a1);
    ASSERT_TRUE(dyn.start(u0, std::vector<double>(u0.size(), 0.0), {}, 0.0).empty());
    // the same Newmark step on the modal oscillator q'' + 2 zeta w q' + w^2 q = 0, q(0) = 1
    double q = 1.0, qd = 0.0, qdd = -w * w * q;
    const double amp = tip_x(*f, u0);
    for (int n = 1; n <= 250; ++n) {
        dyn.step(h, {}, 0.0);
        const double c0 = 4.0 / (h * h), c1 = 2.0 / h, c2 = 4.0 / h;
        const double keff = w * w + c1 * 2.0 * zeta * w + c0;
        const double rhs = c0 * q + c2 * qd + qdd + 2.0 * zeta * w * (c1 * q + qd);
        const double qn = rhs / keff;
        const double qddn = c0 * (qn - q) - c2 * qd - qdd;
        qd += 0.5 * h * (qdd + qddn);
        qdd = qddn;
        q = qn;
        EXPECT_NEAR(tip_x(*f, dyn.displacement()), amp * q, 1e-5 * std::abs(amp)) << "step " << n;
    }
    EXPECT_LT(std::abs(q), 0.4) << "the mode must have decayed over its 10 periods: the check is not vacuous";
}

TEST(FrameDynamics, AChangingStepFollowsTheOscillatorOfEachStep)
{
    auto f = build(cantilever(10.0, 10, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7)));
    ASSERT_TRUE(f);
    FrameModes modes;
    ASSERT_TRUE(frame_modes(*f, 3, modes).empty());
    double a0 = 0.0, a1 = 0.0;
    rayleigh_coefficients(modes.frequency[0], 0.02, modes.frequency[2], 0.05, a0, a1);
    const double w = 2.0 * pi * modes.frequency[0];
    const double zeta = 0.5 * (a0 / w + a1 * w);
    const double h1 = (2.0 * pi / w) / 25.0, h2 = 0.7 * h1;
    std::vector<double> u0 = modes.shape[0];
    for (double& x : u0) { x *= 0.01; }
    const double amp = tip_x(*f, u0);
    FrameDynamics changing(*f, a0, a1), steady(*f, a0, a1);
    ASSERT_TRUE(changing.start(u0, std::vector<double>(u0.size(), 0.0), {}, 0.0).empty());
    ASSERT_TRUE(steady.start(u0, std::vector<double>(u0.size(), 0.0), {}, 0.0).empty());
    // Newmark's step of the modal oscillator q'' + 2 zeta w q' + w^2 q = 0
    struct Osc { double q{1.0}, qd{0.0}, qdd{0.0}; };
    auto advance = [&] (Osc& o, double h) {
        const double c0 = 4.0 / (h * h), c1 = 2.0 / h, c2 = 4.0 / h;
        const double qn = (c0 * o.q + c2 * o.qd + o.qdd + 2.0 * zeta * w * (c1 * o.q + o.qd)) / (w * w + c1 * 2.0 * zeta * w + c0);
        const double qddn = c0 * (qn - o.q) - c2 * o.qd - o.qdd;
        o.qd += 0.5 * h * (o.qdd + qddn);
        o.qdd = qddn;
        o.q = qn;
    };
    Osc a, b;
    a.qdd = b.qdd = -w * w;
    for (int n = 1; n <= 200; ++n) {
        const double h = (n % 2) ? h1 : h2;
        changing.step(h, {}, 0.0);
        steady.step(h2, {}, 0.0);
        advance(a, h);
        advance(b, h2);
        EXPECT_NEAR(tip_x(*f, changing.displacement()), amp * a.q, 1e-5 * std::abs(amp)) << "step " << n;
        EXPECT_NEAR(tip_x(*f, steady.displacement()), amp * b.q, 1e-5 * std::abs(amp)) << "step " << n;
    }
    EXPECT_GT(std::abs(a.q - b.q), 1e-3) << "the two step sequences must part: the check is not vacuous";
    // both took h2 last: one factorisation serves the two
    ASSERT_NE(steady.effective_stiffness_factor(), nullptr);
    EXPECT_EQ(changing.effective_stiffness_factor(), steady.effective_stiffness_factor());
}

TEST(FrameDynamics, AHeldLoadDoesTheWorkTheEnergyGains)
{
    FrameInputs in;
    ASSERT_TRUE(read_subdyn(case_file("C", ".dat"), in).empty());
    auto f = build(in);
    ASSERT_TRUE(f);
    FrameModes modes;
    ASSERT_TRUE(frame_modes(*f, 2, modes).empty());
    double a0 = 0.0, a1 = 0.0;
    rayleigh_coefficients(modes.frequency[0], 0.02, 10.0 * modes.frequency[0], 0.02, a0, a1);
    const double h = 1.0 / (20.0 * modes.frequency[0]);
    // swinging in its first two modes, under a pull on an arm tip and on the peak that changes every step
    std::vector<double> u0(modes.shape[0].size()), v0(u0.size());
    for (std::size_t i = 0; i < u0.size(); ++i) {
        u0[i] = 0.01 * modes.shape[0][i];
        v0[i] = 0.3 * modes.shape[1][i];
    }
    auto load = [&] (int n) {
        std::vector<double> l(f->num_dofs(), 0.0);
        l[6 * f->node_of_joint(18)] = 2.0e4 * (1.0 + 0.5 * std::sin(0.7 * n));
        l[6 * f->node_of_joint(17) + 2] = -5.0e3 * (1.0 + 0.5 * std::cos(1.3 * n));
        return l;
    };
    // the loads' work less the energy gain and the damping's work over one step
    auto imbalance = [&] (FrameDynamics& d, int n, bool held, double& scale) {
        const std::vector<double> u = d.displacement();
        const double e = d.kinetic_energy() + d.strain_energy();
        const std::vector<double> l = load(n);
        if (held) { d.step_held(h, l, 0.0); } else { d.step(h, l, 0.0); }
        std::vector<double> du(u.size()), vm(u.size());
        for (std::size_t i = 0; i < u.size(); ++i) { du[i] = d.displacement()[i] - u[i]; vm[i] = du[i] / h; }
        const std::vector<double> mv = f->apply_mass(vm), kv = f->apply_stiffness(vm);
        double damping = 0.0;
        for (std::size_t i = 0; i < vm.size(); ++i) { damping += h * vm[i] * (a0 * mv[i] + a1 * kv[i]); }
        const double work = dot(f->load_vector(l, 0.0), du);
        const double gain = d.kinetic_energy() + d.strain_energy() - e;
        scale = std::abs(work) + std::abs(gain) + std::abs(damping);
        return work - gain - damping;
    };
    FrameDynamics held(*f, a0, a1), newmark(*f, a0, a1);
    ASSERT_TRUE(held.start(u0, v0, load(0), 0.0).empty());
    ASSERT_TRUE(newmark.start(u0, v0, load(0), 0.0).empty());
    double worst_newmark = 0.0;
    for (int n = 1; n <= 40; ++n) {
        // imbalance() sets scale: call it before the tolerance reads scale (the order in which a call's
        // arguments are evaluated is unspecified; GCC read scale first, as 0)
        double scale = 0.0;
        const double r_held = imbalance(held, n, true, scale);
        EXPECT_NEAR(r_held, 0.0, 1e-10 * scale) << "step " << n;
        // the end acceleration is the equation of motion's under the held loads (on the free degrees of
        // freedom: the fixed ones carry the support reactions)
        const std::vector<double> fl = f->load_vector(load(n), 0.0);
        const std::vector<double> ma = f->apply_mass(held.acceleration());
        const std::vector<double> mv = f->apply_mass(held.velocity()), kv = f->apply_stiffness(held.velocity());
        const std::vector<double> ku = f->apply_stiffness(held.displacement());
        double rmax = 0.0, fmax = 0.0;
        for (const std::size_t i : f->free_dofs()) {
            rmax = std::max(rmax, std::abs(ma[i] + a0 * mv[i] + a1 * kv[i] + ku[i] - fl[i]));
            fmax = std::max(fmax, std::abs(ku[i]) + std::abs(fl[i]));
        }
        EXPECT_LT(rmax, 1e-9 * fmax) << "step " << n;
        double s = 0.0;
        const double r_newmark = imbalance(newmark, n, false, s);
        worst_newmark = std::max(worst_newmark, std::abs(r_newmark) / s);
    }
    // Newmark's step takes the mean of the last step's loads and these: the balance is not met by construction
    EXPECT_GT(worst_newmark, 1e-3);
}

TEST(FrameDynamics, AStepLoadSettlesAtTheStaticResponse)
{
    FrameInputs in;
    ASSERT_TRUE(read_subdyn(case_file("C", ".dat"), in).empty());
    auto f = build(in);
    ASSERT_TRUE(f);
    FrameModes modes;
    const std::string err = frame_modes(*f, 1, modes);
    ASSERT_TRUE(err.empty()) << err;
    double a0 = 0.0, a1 = 0.0;
    rayleigh_coefficients(modes.frequency[0], 0.3, 10.0 * modes.frequency[0], 0.3, a0, a1);
    std::vector<double> load(f->num_dofs(), 0.0);
    load[6 * f->node_of_joint(18)] = 2.0e4;       // a lateral pull on an arm tip
    load[6 * f->node_of_joint(17) + 1] = -5.0e3;  // and on the peak
    const double g = 9.80665;
    FrameDynamics dyn(*f, a0, a1);
    dyn.start_static({}, 0.0);
    const double h = 1.0 / (20.0 * modes.frequency[0]);
    for (int n = 0; n < 2000; ++n) { dyn.step(h, load, g); }
    const FrameSolution s = f->solve(load, g);
    double umax = 0.0;
    for (const double x : s.displacement) { umax = std::max(umax, std::abs(x)); }
    for (std::size_t d = 0; d < s.displacement.size(); ++d) {
        EXPECT_NEAR(dyn.displacement()[d], s.displacement[d], 1e-6 * umax) << "dof " << d;
    }
}

TEST(FrameDynamics, TheStateRestoresTheSameMotion)
{
    FrameInputs in;
    ASSERT_TRUE(read_subdyn(case_file("B", ".dat"), in).empty());
    auto f = build(in);
    ASSERT_TRUE(f);
    std::vector<double> load(f->num_dofs(), 0.0);
    load[6 * f->node_of_joint(17)] = 1.0e4;
    FrameDynamics a(*f, 0.1, 1.0e-3);
    a.start_static({}, 9.80665);
    for (int n = 0; n < 10; ++n) { a.step(0.01, load, 9.80665); }
    const std::vector<double> saved = a.state();
    for (int n = 0; n < 10; ++n) { a.step(0.01, load, 9.80665); }
    FrameDynamics b(*f, 0.1, 1.0e-3);
    ASSERT_TRUE(b.set_state(saved));
    for (int n = 0; n < 10; ++n) { b.step(0.01, load, 9.80665); }
    EXPECT_EQ(a.state(), b.state());
    EXPECT_DOUBLE_EQ(b.time(), saved[0] + 0.1);
}

TEST(FrameDynamics, RefusesBadStepsStatesAndAMasslessNode)
{
    auto f = build(cantilever(10.0, 2, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7)));
    ASSERT_TRUE(f);
    FrameDynamics dyn(*f, 0.0, 0.0);
    dyn.start_static({}, 0.0);
    EXPECT_NE(erf_gtest::abort_message([&] { dyn.step(0.0, {}, 0.0); }).find("positive"), std::string::npos);
    std::vector<double> s = dyn.state();
    EXPECT_FALSE(dyn.set_state(std::vector<double>(s.begin(), s.end() - 1)));
    s[3] = std::nan("");
    EXPECT_FALSE(dyn.set_state(s));
    EXPECT_FALSE(dyn.start(std::vector<double>(3, 0.0), std::vector<double>(3, 0.0), {}, 0.0).empty());
    // a member of zero density leaves its nodes without mass
    auto massless = build(cantilever(10.0, 2, arbitrary(4.0e-3, 1.2e-5, 0.6e-5, 1.8e-5, 4.0e-7, 0.0)));
    ASSERT_TRUE(massless);
    FrameDynamics md(*massless, 0.0, 0.0);
    const std::vector<double> zero(massless->num_dofs(), 0.0);
    EXPECT_NE(md.start(zero, zero, {}, 0.0).find("mass"), std::string::npos);
    FrameModes modes;
    EXPECT_FALSE(frame_modes(*f, 0, modes).empty());
}
