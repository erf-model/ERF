// A lattice tower's frame model generated from its dimensions.

#include "ERF_LatticeFrame.H"
#include "ERF_Frame.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <sstream>

namespace erf_towers {

namespace {

using Point = std::array<double,3>;

/** The corners of a level in order: (+,+), (-,+), (-,-), (+,-). */
constexpr int sx[4] = {1, -1, -1, 1};
constexpr int sy[4] = {1, 1, -1, -1};

/** The point where segments a0-a1 and b0-b1 come closest: the crossing of a face's two diagonals. */
Point crossing (const Point& a0, const Point& a1, const Point& b0, const Point& b1)
{
    Point d1{}, d2{}, w{};
    for (int i = 0; i < 3; ++i) { d1[i] = a1[i] - a0[i]; d2[i] = b1[i] - b0[i]; w[i] = a0[i] - b0[i]; }
    const double a = d1[0] * d1[0] + d1[1] * d1[1] + d1[2] * d1[2];
    const double b = d1[0] * d2[0] + d1[1] * d2[1] + d1[2] * d2[2];
    const double c = d2[0] * d2[0] + d2[1] * d2[1] + d2[2] * d2[2];
    const double d = d1[0] * w[0] + d1[1] * w[1] + d1[2] * w[2];
    const double e = d2[0] * w[0] + d2[1] * w[1] + d2[2] * w[2];
    const double s = (b * e - c * d) / (a * c - b * b);
    return {{a0[0] + s * d1[0], a0[1] + s * d1[1], a0[2] + s * d1[2]}};
}

} // namespace

std::string LatticeSpec::validate () const
{
    std::ostringstream m;
    auto finite = [] (double v) { return std::isfinite(v); };
    if (!(finite(base_width) && base_width > 0.0)) { return "the base width must be > 0 (m)"; }
    if (!(finite(top_width) && top_width > 0.0 && top_width <= base_width)) { return "the top width must be in (0, base width] (m)"; }
    if (!(finite(arm_height) && arm_height > 0.0)) { return "the cross-arm's height must be > 0 (m)"; }
    if (!(finite(arm_depth) && arm_depth > 0.0 && arm_depth < arm_height)) {
        m << "the cross-arm's depth (" << arm_depth << " m) must be in (0, its height " << arm_height << " m)";
        return m.str();
    }
    if (!(finite(arm_length) && arm_length > top_width)) {
        m << "the cross-arm's length (" << arm_length << " m) must exceed the top width (" << top_width << " m)";
        return m.str();
    }
    if (!(finite(peak) && peak >= 0.0)) { return "the peak must be >= 0 (m)"; }
    if (panels < 1 || panels > 200) { return "the number of panels must be in [1, 200]"; }
    if (!(finite(leg_b) && finite(leg_t) && leg_t > 0.0 && leg_t < leg_b)) { return "the legs' angle needs 0 < t < b (m)"; }
    if (!(finite(brace_b) && finite(brace_t) && brace_t > 0.0 && brace_t < brace_b)) { return "the bracing's angle needs 0 < t < b (m)"; }
    if (!(finite(E) && E > 0.0 && finite(G) && G > 0.0 && finite(rho) && rho >= 0.0 && finite(yield) && yield > 0.0)) {
        return "the steel needs E > 0, G > 0, rho >= 0 and a yield strength > 0";
    }
    return std::string();
}

std::string lattice_frame (const LatticeSpec& s, FrameInputs& in, std::vector<MemberDesign>& designs, std::vector<int>* load_joints)
{
    const std::string err = s.validate();
    if (!err.empty()) { return err; }
    in = FrameInputs();
    designs.clear();
    in.theory = BeamTheory::EulerBernoulli;
    in.divisions = 1;

    // the sections: 1 the legs' angle, 2 the bracing's; principal axes along the members' local axes
    for (int id = 1; id <= 2; ++id) {
        const double b = (id == 1) ? s.leg_b : s.brace_b, t = (id == 1) ? s.leg_t : s.brace_t;
        const AngleProperties a = angle_properties(b, t);
        FrameSection sec;
        sec.id = id;
        sec.shape = SectionShape::Arbitrary;
        sec.E = s.E; sec.G = s.G; sec.rho = s.rho;
        sec.A = a.A; sec.Asx = 0.5 * a.A; sec.Asy = 0.5 * a.A;
        sec.Ixx = a.Iu; sec.Iyy = a.Iv; sec.J0 = a.Iu + a.Iv; sec.Jt = a.J;
        in.sections.push_back(sec);
    }

    std::vector<Point> x;
    auto joint = [&] (const Point& p) { x.push_back(p); return static_cast<int>(x.size()); };
    // the spin that turns a member's local x (the axis of Ixx, an angle's major axis: its axis of symmetry)
    // towards sym, projected normal to the member; 0 without principal_axes or a direction
    auto spin_towards = [&] (int a, int b, const Point& sym) {
        if (!s.principal_axes || (sym[0] == 0.0 && sym[1] == 0.0 && sym[2] == 0.0)) { return 0.0; }
        const std::array<double,9> d = direction_cosines(x[static_cast<std::size_t>(a - 1)], x[static_cast<std::size_t>(b - 1)], 0.0);
        const double px = sym[0] * d[0] + sym[1] * d[3] + sym[2] * d[6];
        const double py = sym[0] * d[1] + sym[1] * d[4] + sym[2] * d[7];
        return std::atan2(py, px);
    };
    // a face member's axis of symmetry: half way between the face's outward normal and the face's direction
    // normal to the member (one leg flat on the face, the other standing out of it), both normal to the member
    auto face_sym = [&] (int a, int b, const Point& normal, std::size_t h) {
        const Point& pa = x[static_cast<std::size_t>(a - 1)];
        const Point& pb = x[static_cast<std::size_t>(b - 1)];
        Point e{{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]}};
        const double le = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
        for (auto& v : e) { v /= le; }
        const double ne = normal[0] * e[0] + normal[1] * e[1] + normal[2] * e[2];
        Point n{{normal[0] - ne * e[0], normal[1] - ne * e[1], normal[2] - ne * e[2]}};
        const double ln = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
        for (auto& v : n) { v /= ln; }
        Point p{{n[1] * e[2] - n[2] * e[1], n[2] * e[0] - n[0] * e[2], n[0] * e[1] - n[1] * e[0]}};
        // the leg flat on the face points away from the member, along p or -p: the same way for a member and its mirror
        // images in the tower's planes of symmetry x = 0 and y = 0, so the frame keeps the tower's symmetry. Towards
        // the face's line of symmetry, across the horizontal axis h (0: x, 1: y) the face's own layout gives (y on a
        // shaft face facing x and on a cross-arm's side, x on a shaft face facing y and on a cross-arm's top and
        // bottom); for a member centred on it, down; for a level one, out along the other horizontal axis g
        const std::size_t g = 1 - h;
        const double mh = 0.5 * (pa[h] + pb[h]), mg = 0.5 * (pa[g] + pb[g]);
        constexpr double tiny = 1.0e-9;
        double sense = 1.0;
        if (std::abs(mh) > tiny && std::abs(p[h]) > tiny) { sense = (p[h] * mh < 0.0) ? 1.0 : -1.0; }
        else if (std::abs(p[2]) > tiny) { sense = (p[2] < 0.0) ? 1.0 : -1.0; }
        else if (std::abs(mg) > tiny) { sense = (p[g] * mg > 0.0) ? 1.0 : -1.0; }
        for (auto& v : p) { v *= sense; }
        return Point{{n[0] + p[0], n[1] + p[1], n[2] + p[2]}};
    };
    auto member = [&] (int a, int b, bool leg, MemberRole role = MemberRole::Bracing, const Point& sym = Point{{0.0, 0.0, 0.0}}) {
        FrameMember m;
        m.id = static_cast<int>(in.members.size()) + 1;
        m.joint_a = a;
        m.joint_b = b;
        m.section = leg ? 1 : 2;
        m.shape = SectionShape::Arbitrary;
        m.spin = spin_towards(a, b, sym);
        in.members.push_back(m);
        MemberDesign d;
        d.member = m.id;
        d.yield = s.yield;
        d.b = leg ? s.leg_b : s.brace_b;
        d.t = leg ? s.leg_t : s.brace_t;
        d.role = leg ? MemberRole::Leg : role;
        d.one_leg = !leg;
        d.ends = leg ? EndLoading::Concentric : EndLoading::BothEccentric;
        d.restraint = EndRestraint::None;
        designs.push_back(d);
    };

    // the shaft's levels and their corners
    const double H = s.arm_height, bottom = H - s.arm_depth;
    std::vector<double> z;
    for (int k = 0; k <= s.panels; ++k) { z.push_back(bottom * k / s.panels); }
    z.push_back(H);
    // the peak in panels no taller than the shaft's
    const int peak_panels = (s.peak > 0.0) ? static_cast<int>(std::ceil(s.peak / (bottom / s.panels) - 1.0e-9)) : 0;
    for (int k = 1; k <= peak_panels; ++k) { z.push_back(H + s.peak * k / peak_panels); }
    const int levels = static_cast<int>(z.size());
    auto half = [&] (double zz) {
        const double f = std::min(zz, H) / H;
        return 0.5 * (s.base_width + f * (s.top_width - s.base_width));
    };
    std::vector<std::array<int,4>> corner(static_cast<std::size_t>(levels));
    for (int k = 0; k < levels; ++k) {
        const double h = half(z[static_cast<std::size_t>(k)]);
        const double zk = z[static_cast<std::size_t>(k)];
        for (int c = 0; c < 4; ++c) { corner[static_cast<std::size_t>(k)][c] = joint({{sx[c] * h, sy[c] * h, zk}}); }
    }
    auto at = [&] (int id) -> const Point& { return x[static_cast<std::size_t>(id - 1)]; };
    // the crossings of each panel face's diagonals
    std::vector<std::array<int,4>> cross(static_cast<std::size_t>(levels - 1), std::array<int,4>{{0, 0, 0, 0}});
    if (s.crossed) {
        for (int k = 0; k + 1 < levels; ++k) {
            const auto& lo = corner[static_cast<std::size_t>(k)];
            const auto& hi = corner[static_cast<std::size_t>(k + 1)];
            for (int c = 0; c < 4; ++c) {
                const int n = (c + 1) % 4;
                cross[static_cast<std::size_t>(k)][c] = joint(crossing(at(lo[c]), at(hi[n]), at(lo[n]), at(hi[c])));
            }
        }
    }
    const int centre = joint({{0.0, 0.0, H}});
    const int top = levels - 1 - peak_panels;                 // the cross-arm's top level
    const int arm_bottom = top - 1;                           // the cross-arm's bottom level

    // the cross-arms: four chords from the shaft's corners at the cross-arm to the tip, in panels about
    // as long as the cross-arm is deep
    const double reach = 0.5 * (s.arm_length - s.top_width);
    const int arm_panels = std::max(1, static_cast<int>(std::lround(reach / s.arm_depth)));
    // per side, per station: the chords' points (top 0, top 1, bottom 0, bottom 1)
    std::array<std::vector<std::array<int,4>>,2> station;
    std::array<int,2> tip{{0, 0}};
    for (int side = 0; side < 2; ++side) {
        const int c0 = (side == 0) ? 0 : 2, c1 = (side == 0) ? 1 : 3;
        std::array<int,4> root{{corner[static_cast<std::size_t>(top)][c0], corner[static_cast<std::size_t>(top)][c1],
                                corner[static_cast<std::size_t>(arm_bottom)][c0], corner[static_cast<std::size_t>(arm_bottom)][c1]}};
        station[static_cast<std::size_t>(side)].push_back(root);
        const Point tp{{0.0, (side == 0 ? 0.5 : -0.5) * s.arm_length, H}};
        for (int j = 1; j < arm_panels; ++j) {
            const double f = static_cast<double>(j) / arm_panels;
            std::array<int,4> st{};
            for (int q = 0; q < 4; ++q) {
                const Point r = at(root[static_cast<std::size_t>(q)]);
                st[static_cast<std::size_t>(q)] =
                    joint({{r[0] + f * (tp[0] - r[0]), r[1] + f * (tp[1] - r[1]), r[2] + f * (tp[2] - r[2])}});
            }
            station[static_cast<std::size_t>(side)].push_back(st);
        }
        tip[static_cast<std::size_t>(side)] = joint(tp);
    }
    for (const auto& p : x) { FrameJoint j; j.id = static_cast<int>(in.joints.size()) + 1; j.x = p; in.joints.push_back(j); }

    // the outward normal of the face between corners c and c + 1 of the panel above level k: leaning in with the
    // shaft's taper below the cross-arm, upright on the peak
    const double lean = 0.5 * (s.base_width - s.top_width) / H;
    auto face_normal = [&] (int c, int k) {
        const int n = (c + 1) % 4;
        const double nx = 0.5 * (sx[c] + sx[n]), ny = 0.5 * (sy[c] + sy[n]);
        const double nz = (z[static_cast<std::size_t>(k)] < H - 1.0e-9) ? lean : 0.0;
        const double l = std::sqrt(nx * nx + ny * ny + nz * nz);
        return Point{{nx / l, ny / l, nz / l}};
    };
    // the horizontal axis along the face between corners c and c + 1: y on a face facing x, x on one facing y
    auto face_axis = [&] (int c) { return (sx[c] == sx[(c + 1) % 4]) ? std::size_t(1) : std::size_t(0); };
    // legs, their axis of symmetry towards the corner they stand on
    for (int k = 0; k + 1 < levels; ++k) {
        for (int c = 0; c < 4; ++c) {
            member(corner[static_cast<std::size_t>(k)][c], corner[static_cast<std::size_t>(k + 1)][c], true, MemberRole::Leg,
                   Point{{static_cast<double>(sx[c]), static_cast<double>(sy[c]), 0.0}});
        }
    }
    // struts at every level above the ground: redundant between crossed diagonals, which carry the shear without them
    const MemberRole strut = s.crossed ? MemberRole::Redundant : MemberRole::Bracing;
    for (int k = 1; k < levels; ++k) {
        for (int c = 0; c < 4; ++c) {
            const int a = corner[static_cast<std::size_t>(k)][c], b = corner[static_cast<std::size_t>(k)][(c + 1) % 4];
            member(a, b, false, strut, face_sym(a, b, face_normal(c, k - 1), face_axis(c)));   // in the face of the panel below
        }
    }
    // face bracing: crossed diagonals in halves, or one diagonal per face, alternating panel by panel
    for (int k = 0; k + 1 < levels; ++k) {
        const auto& lo = corner[static_cast<std::size_t>(k)];
        const auto& hi = corner[static_cast<std::size_t>(k + 1)];
        for (int c = 0; c < 4; ++c) {
            const int n = (c + 1) % 4;
            const Point fn = face_normal(c, k);
            auto brace = [&] (int a, int b) { member(a, b, false, MemberRole::Bracing, face_sym(a, b, fn, face_axis(c))); };
            if (s.crossed) {
                const int xc = cross[static_cast<std::size_t>(k)][c];
                brace(lo[c], xc);
                brace(xc, hi[n]);
                brace(lo[n], xc);
                brace(xc, hi[c]);
            } else if (k % 2 == 0) {
                brace(lo[c], hi[n]);
            } else {
                brace(lo[n], hi[c]);
            }
        }
    }
    // the hanger: the centre of the cross-arm to its four top corners and its four bottom corners
    for (int c = 0; c < 4; ++c) { member(centre, corner[static_cast<std::size_t>(top)][c], false); }
    for (int c = 0; c < 4; ++c) { member(centre, corner[static_cast<std::size_t>(arm_bottom)][c], false); }
    // the cross-arms: chords, then at each inner station its frame (four struts and a diagonal), then the
    // faces' diagonals of every panel but the last (whose faces the chords and the struts triangulate)
    for (int side = 0; side < 2; ++side) {
        const auto& st = station[static_cast<std::size_t>(side)];
        // the corners the side's chords start from: (top c0, top c1, bottom c0, bottom c1), and the outward normals
        // of its faces (top, bottom, the c0 side, the c1 side). Each face is a plane through two chords that meet at
        // the tip: the top is level, the bottom rises to the tip and the sides close in on it
        const int c0 = (side == 0) ? 0 : 2, c1 = (side == 0) ? 1 : 3;
        const std::array<Point,4> outward{{Point{{0.0, 0.0, 1.0}}, Point{{0.0, 0.0, -1.0}},
                                           Point{{static_cast<double>(sx[c0]), 0.0, 0.0}}, Point{{static_cast<double>(sx[c1]), 0.0, 0.0}}}};
        static constexpr int face_corners[4][2] = {{0, 1}, {2, 3}, {0, 2}, {1, 3}};
        // the faces' horizontal axes for face_sym: across the cross-arm (x) on its top and bottom, along it (y) on its sides
        static constexpr std::array<std::size_t,4> arm_axis{{0, 0, 1, 1}};
        std::array<Point,4> normal;
        for (std::size_t f = 0; f < 4; ++f) {
            const Point& pa = x[static_cast<std::size_t>(st[0][static_cast<std::size_t>(face_corners[f][0])] - 1)];
            const Point& pb = x[static_cast<std::size_t>(st[0][static_cast<std::size_t>(face_corners[f][1])] - 1)];
            const Point& pt = x[static_cast<std::size_t>(tip[static_cast<std::size_t>(side)] - 1)];
            const Point u{{pb[0] - pa[0], pb[1] - pa[1], pb[2] - pa[2]}}, v{{pt[0] - pa[0], pt[1] - pa[1], pt[2] - pa[2]}};
            Point n{{u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]}};
            const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            const double sense = (n[0] * outward[f][0] + n[1] * outward[f][1] + n[2] * outward[f][2] < 0.0) ? -1.0 : 1.0;
            for (auto& c : n) { c *= sense / l; }
            normal[f] = n;
        }
        // a chord's axis of symmetry towards the cross-arm's corner it runs along, between its two faces' normals
        auto between = [] (const Point& p, const Point& q) { return Point{{p[0] + q[0], p[1] + q[1], p[2] + q[2]}}; };
        const std::array<Point,4> corner_sym{{between(normal[0], normal[2]), between(normal[0], normal[3]),
                                              between(normal[1], normal[2]), between(normal[1], normal[3])}};
        for (int j = 1; j <= arm_panels; ++j) {
            for (std::size_t q = 0; q < 4; ++q) {
                const int to = (j < arm_panels) ? st[static_cast<std::size_t>(j)][q] : tip[static_cast<std::size_t>(side)];
                member(st[static_cast<std::size_t>(j - 1)][q], to, true, MemberRole::Leg, corner_sym[q]);
            }
        }
        for (int j = 1; j < arm_panels; ++j) {
            const auto& p = st[static_cast<std::size_t>(j)];
            member(p[0], p[1], false, MemberRole::Bracing, face_sym(p[0], p[1], normal[0], arm_axis[0]));
            member(p[2], p[3], false, MemberRole::Bracing, face_sym(p[2], p[3], normal[1], arm_axis[1]));
            member(p[0], p[2], false, MemberRole::Bracing, face_sym(p[0], p[2], normal[2], arm_axis[2]));
            member(p[1], p[3], false, MemberRole::Bracing, face_sym(p[1], p[3], normal[3], arm_axis[3]));
            member(p[0], p[3], false);
        }
        // the four faces of the cross-arm: (top 0, top 1), (bottom 0, bottom 1), (top 0, bottom 0), (top 1, bottom 1)
        const auto& faces = face_corners;
        for (int j = 1; j < arm_panels; ++j) {
            const auto& a = st[static_cast<std::size_t>(j - 1)];
            const auto& b = st[static_cast<std::size_t>(j)];
            for (std::size_t f = 0; f < 4; ++f) {
                const int u = (j % 2 == 1) ? a[faces[f][0]] : a[faces[f][1]];
                const int v = (j % 2 == 1) ? b[faces[f][1]] : b[faces[f][0]];
                member(u, v, false, MemberRole::Bracing, face_sym(u, v, normal[f], arm_axis[f]));
            }
        }
    }

    for (int c = 0; c < 4; ++c) {
        FrameSupport sup;
        sup.joint = corner[0][c];
        in.supports.push_back(sup);
    }
    in.interface_joints.push_back(centre);
    if (load_joints != nullptr) {
        load_joints->clear();
        std::vector<bool> crossing(x.size() + 1, false);
        for (const auto& k : cross) { for (const int id : k) { if (id > 0) { crossing[static_cast<std::size_t>(id)] = true; } } }
        for (const auto& j : in.joints) { if (!crossing[static_cast<std::size_t>(j.id)]) { load_joints->push_back(j.id); } }
    }
    return in.validate();
}

std::vector<int> square_joints (const FrameInputs& in)
{
    const double deg = 3.14159265358979323846 / 180.0;
    const double cmax = std::cos(20.0 * deg), cmin = std::sin(20.0 * deg);
    std::vector<bool> square(in.joints.size(), false);
    for (const auto& m : in.members) {
        const int a = in.joint_index(m.joint_a), b = in.joint_index(m.joint_b);
        if (a < 0 || b < 0) { continue; }
        const auto& xa = in.joints[static_cast<std::size_t>(a)].x;
        const auto& xb = in.joints[static_cast<std::size_t>(b)].x;
        const double len = std::sqrt((xb[0] - xa[0]) * (xb[0] - xa[0]) + (xb[1] - xa[1]) * (xb[1] - xa[1]) +
                                     (xb[2] - xa[2]) * (xb[2] - xa[2]));
        if (!(len > 0.0)) { continue; }
        const double ez = std::abs(xb[2] - xa[2]) / len;
        if (ez >= cmax || ez <= cmin) { square[static_cast<std::size_t>(a)] = square[static_cast<std::size_t>(b)] = true; }
    }
    std::vector<int> ids;
    for (std::size_t j = 0; j < in.joints.size(); ++j) { if (square[j]) { ids.push_back(in.joints[j].id); } }
    return ids;
}

} // namespace erf_towers
