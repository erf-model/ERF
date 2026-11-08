// The frame model: SubDyn's beam element, assembly, supports, the static solve and its recovery.

#include "ERF_Frame.H"
#include "ERF_MemberChecks.H"

#include <algorithm>
#include <cmath>
#include <utility>

#include <AMReX.H>
#include <AMReX_BLassert.H>

namespace erf_towers {

namespace {

constexpr double pi = 3.14159265358979323846;

/** The (row, column) of each of SubDyn's 21 SSI entries in the 6 x 6 matrix: the upper triangle column by column. */
std::array<std::pair<int,int>,21> ssi_entries ()
{
    std::array<std::pair<int,int>,21> e{};
    int k = 0;
    for (int j = 0; j < 6; ++j) { for (int i = 0; i <= j; ++i) { e[static_cast<std::size_t>(k++)] = {i, j}; } }
    return e;
}

/** The frame degree of freedom of an element's local index i (0..11: node_a's six, then node_b's). */
std::size_t element_dof (const FrameElement& e, int i)
{
    return 6 * (i < 6 ? e.node_a : e.node_b) + static_cast<std::size_t>(i % 6);
}

/** Rotate 12 values block by block: out = T^T v (to local axes) or out = T v (to frame axes). */
std::array<double,12> rotate12 (const std::array<double,9>& dc, const std::array<double,12>& v, bool to_local)
{
    std::array<double,12> out{};
    for (int b = 0; b < 4; ++b) {
        for (int i = 0; i < 3; ++i) {
            double s = 0.0;
            for (int j = 0; j < 3; ++j) {
                s += (to_local ? dc[static_cast<std::size_t>(3 * j + i)] : dc[static_cast<std::size_t>(3 * i + j)]) *
                     v[static_cast<std::size_t>(3 * b + j)];
            }
            out[static_cast<std::size_t>(3 * b + i)] = s;
        }
    }
    return out;
}

/** T A T^T for a 12 x 12 element matrix A in local axes, T block-diagonal with four dc. */
std::array<double,144> to_frame (const std::array<double,144>& local, const std::array<double,9>& dc)
{
    std::array<double,144> tmp{}, out{};
    for (int j = 0; j < 12; ++j) {          // rotate every column to the frame
        std::array<double,12> col{};
        for (int i = 0; i < 12; ++i) { col[static_cast<std::size_t>(i)] = local[static_cast<std::size_t>(12 * i + j)]; }
        col = rotate12(dc, col, false);
        for (int i = 0; i < 12; ++i) { tmp[static_cast<std::size_t>(12 * i + j)] = col[static_cast<std::size_t>(i)]; }
    }
    for (int i = 0; i < 12; ++i) {          // then every row
        std::array<double,12> rowv{};
        for (int j = 0; j < 12; ++j) { rowv[static_cast<std::size_t>(j)] = tmp[static_cast<std::size_t>(12 * i + j)]; }
        rowv = rotate12(dc, rowv, false);
        for (int j = 0; j < 12; ++j) { out[static_cast<std::size_t>(12 * i + j)] = rowv[static_cast<std::size_t>(j)]; }
    }
    return out;
}

} // namespace

BeamProperties beam_properties (const FrameSection& s, BeamTheory theory)
{
    BeamProperties p;
    p.E = s.E;
    p.G = s.G;
    p.rho = s.rho;
    p.shear = (theory == BeamTheory::Timoshenko);
    const double nu = s.E / (2.0 * s.G) - 1.0;
    if (s.shape == SectionShape::Circular) {
        const double r1 = 0.5 * s.D;
        const double r2 = (s.t == 0.0) ? 0.0 : r1 - s.t;
        p.A = pi * (r1 * r1 - r2 * r2);
        p.Ixx = 0.25 * pi * (std::pow(r1, 4) - std::pow(r2, 4));
        p.Iyy = p.Ixx;
        p.J0 = 2.0 * p.Ixx;
        p.Jt = p.J0;
        if (p.shear) {
            // Steinboeck et al., eq. 13 of the SubDyn theory manual, with D_inner = D - 2t as SubDyn takes it
            const double ratio_sq = std::pow((s.D - 2.0 * s.t) / s.D, 2);
            const double k = 6.0 * std::pow(1.0 + nu, 2) * std::pow(1.0 + ratio_sq, 2) /
                             (std::pow(1.0 + ratio_sq, 2) * (7.0 + 14.0 * nu + 8.0 * nu * nu) +
                              4.0 * ratio_sq * (5.0 + 10.0 * nu + 4.0 * nu * nu));
            p.kappa_x = k;
            p.kappa_y = k;
        }
    } else if (s.shape == SectionShape::Rectangular) {
        const double sa1 = s.Sa, sb1 = s.Sb;
        const double sa2 = (s.t == 0.0) ? 0.0 : sa1 - 2.0 * s.t;
        const double sb2 = (s.t == 0.0) ? 0.0 : sb1 - 2.0 * s.t;
        p.A = sa1 * sb1 - sa2 * sb2;
        p.Ixx = (sa1 * std::pow(sb1, 3) - sa2 * std::pow(sb2, 3)) / 12.0;
        p.Iyy = (std::pow(sa1, 3) * sb1 - std::pow(sa2, 3) * sb2) / 12.0;
        p.J0 = (sa1 * sb1 * (sa1 * sa1 + sb1 * sb1) - sa2 * sb2 * (sa2 * sa2 + sb2 * sb2)) / 12.0;
        const bool solid = (sa2 == 0.0 || sb2 == 0.0);
        if (solid) {
            const double a = std::max(sa1, sb1), b = std::min(sa1, sb1);
            p.Jt = a * std::pow(b, 3) / 16.0 * (16.0 / 3.0 - 3.36 * b / a * (1.0 - std::pow(b, 4) / 12.0 / std::pow(a, 4)));
        } else {
            const double t = s.t;
            p.Jt = 2.0 * t * std::pow(sa1 - t, 2) * std::pow(sb1 - t, 2) / (sa1 + sb1 - 2.0 * t);
        }
        if (p.shear) {
            if (solid) {
                p.kappa_x = 10.0 * (1.0 + nu) / (12.0 + 11.0 * nu);
                p.kappa_y = p.kappa_x;
            } else {
                auto kappa = [nu] (double ratio) {
                    return 10.0 * (1.0 + nu) * std::pow(1.0 + 3.0 * ratio, 2) /
                           ((12.0 + 72.0 * ratio + 150.0 * ratio * ratio + 90.0 * std::pow(ratio, 3)) +
                            nu * (11.0 + 66.0 * ratio + 135.0 * ratio * ratio + 90.0 * std::pow(ratio, 3)) +
                            10.0 * ratio * ratio * ((3.0 + nu) * ratio + 3.0 * ratio * ratio));
                };
                p.kappa_x = kappa(sb2 / sa1);
                p.kappa_y = kappa(sa2 / sb1);
            }
        }
    } else {
        p.A = s.A;
        p.Ixx = s.Ixx;
        p.Iyy = s.Iyy;
        p.J0 = s.J0;
        p.Jt = s.Jt;
        // Euler-Bernoulli beams take no shear area (SubDyn ignores it): 1, so that no shear
        // ratio divides by zero, even where the compiler evaluates both sides of p.shear's test
        p.kappa_x = p.shear ? s.Asx / s.A : 1.0;
        p.kappa_y = p.shear ? s.Asy / s.A : 1.0;
    }
    return p;
}

std::array<double,9> direction_cosines (const std::array<double,3>& a, const std::array<double,3>& b, double spin)
{
    const double dx = b[0] - a[0], dy = b[1] - a[1], dz = b[2] - a[2];
    const double dxy = std::sqrt(dx * dx + dy * dy);
    const double L = std::sqrt(dx * dx + dy * dy + dz * dz);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(L > 0.0, "direction_cosines: an element needs two distinct points");
    std::array<double,9> d{};
    if (dxy <= 1.0e-14 * L) {
        // a vertical element: x kept along the frame's x; pointing down flips y and z
        d = {{1.0, 0.0, 0.0, 0.0, (dz < 0.0) ? -1.0 : 1.0, 0.0, 0.0, 0.0, (dz < 0.0) ? -1.0 : 1.0}};
    } else {
        d = {{ dy / dxy,  dx * dz / (L * dxy), dx / L,
              -dx / dxy,  dz * dy / (L * dxy), dy / L,
               0.0,      -dxy / L,             dz / L}};
    }
    // the spin turns the local x and y axes about z: dc <- dc R(spin)
    const double c = std::cos(spin), s = std::sin(spin);
    std::array<double,9> r = d;
    for (int i = 0; i < 3; ++i) {
        r[static_cast<std::size_t>(3 * i)]     =  c * d[static_cast<std::size_t>(3 * i)] + s * d[static_cast<std::size_t>(3 * i + 1)];
        r[static_cast<std::size_t>(3 * i + 1)] = -s * d[static_cast<std::size_t>(3 * i)] + c * d[static_cast<std::size_t>(3 * i + 1)];
    }
    return r;
}

std::array<double,144> beam_stiffness_local (const BeamProperties& p, double L)
{
    std::array<double,144> k{};
    auto K = [&k] (int i, int j) -> double& { return k[static_cast<std::size_t>(12 * (i - 1) + (j - 1))]; };   // SubDyn's 1-based K(i,j)
    const double kx = p.shear ? 12.0 * p.E * p.Iyy / (p.G * p.kappa_x * p.A * L * L) : 0.0;
    const double ky = p.shear ? 12.0 * p.E * p.Ixx / (p.G * p.kappa_y * p.A * L * L) : 0.0;
    const double E = p.E, G = p.G;
    K(9, 9) = E * p.A / L;
    K(7, 7) = 12.0 * E * p.Iyy / (L * L * L * (1.0 + kx));
    K(8, 8) = 12.0 * E * p.Ixx / (L * L * L * (1.0 + ky));
    K(12, 12) = G * p.Jt / L;
    K(10, 10) = (4.0 + ky) * E * p.Ixx / (L * (1.0 + ky));
    K(11, 11) = (4.0 + kx) * E * p.Iyy / (L * (1.0 + kx));
    K(2, 4) = -6.0 * E * p.Ixx / (L * L * (1.0 + ky));
    K(1, 5) = 6.0 * E * p.Iyy / (L * L * (1.0 + kx));
    K(4, 10) = (2.0 - ky) * E * p.Ixx / (L * (1.0 + ky));
    K(5, 11) = (2.0 - kx) * E * p.Iyy / (L * (1.0 + kx));
    K(3, 3) = K(9, 9);
    K(1, 1) = K(7, 7);
    K(2, 2) = K(8, 8);
    K(6, 6) = K(12, 12);
    K(4, 4) = K(10, 10);
    K(5, 5) = K(11, 11);
    K(4, 2) = K(2, 4);
    K(5, 1) = K(1, 5);
    K(10, 4) = K(4, 10);
    K(11, 5) = K(5, 11);
    K(12, 6) = -K(6, 6);
    K(10, 2) = K(4, 2);
    K(11, 1) = K(5, 1);
    K(9, 3) = -K(3, 3);
    K(7, 1) = -K(1, 1);
    K(8, 2) = -K(2, 2);
    K(6, 12) = -K(6, 6);
    K(2, 10) = K(4, 2);
    K(1, 11) = K(5, 1);
    K(3, 9) = -K(3, 3);
    K(1, 7) = -K(1, 1);
    K(2, 8) = -K(2, 2);
    K(11, 7) = -K(5, 1);
    K(10, 8) = -K(4, 2);
    K(7, 11) = -K(5, 1);
    K(8, 10) = -K(4, 2);
    K(7, 5) = -K(5, 1);
    K(5, 7) = -K(5, 1);
    K(8, 4) = -K(4, 2);
    K(4, 8) = -K(4, 2);
    return k;
}

std::array<double,144> beam_stiffness (const BeamProperties& p, double L, const std::array<double,9>& dc)
{
    return to_frame(beam_stiffness_local(p, L), dc);
}

std::array<double,144> beam_mass_local (const BeamProperties& p, double L)
{
    std::array<double,144> m{};
    auto M = [&m] (int i, int j) -> double& { return m[static_cast<std::size_t>(12 * (i - 1) + (j - 1))]; };   // SubDyn's 1-based M(i,j)
    const double t = p.rho * p.A * L, rx = p.rho * p.Ixx, ry = p.rho * p.Iyy, po = p.rho * p.J0 * L;
    M(9, 9) = t / 3.0;
    M(7, 7) = 13.0 * t / 35.0 + 6.0 * ry / (5.0 * L);
    M(8, 8) = 13.0 * t / 35.0 + 6.0 * rx / (5.0 * L);
    M(12, 12) = po / 3.0;
    M(10, 10) = t * L * L / 105.0 + 2.0 * L * rx / 15.0;
    M(11, 11) = t * L * L / 105.0 + 2.0 * L * ry / 15.0;
    M(2, 4) = -11.0 * t * L / 210.0 - rx / 10.0;
    M(1, 5) = 11.0 * t * L / 210.0 + ry / 10.0;
    M(3, 9) = t / 6.0;
    M(5, 7) = 13.0 * t * L / 420.0 - ry / 10.0;
    M(4, 8) = -13.0 * t * L / 420.0 + rx / 10.0;
    M(6, 12) = po / 6.0;
    M(2, 10) = 13.0 * t * L / 420.0 - rx / 10.0;
    M(1, 11) = -13.0 * t * L / 420.0 + ry / 10.0;
    M(8, 10) = 11.0 * t * L / 210.0 + rx / 10.0;
    M(7, 11) = -11.0 * t * L / 210.0 - ry / 10.0;
    M(1, 7) = 9.0 * t / 70.0 - 6.0 * ry / (5.0 * L);
    M(2, 8) = 9.0 * t / 70.0 - 6.0 * rx / (5.0 * L);
    M(4, 10) = -L * L * t / 140.0 - rx * L / 30.0;
    M(5, 11) = -L * L * t / 140.0 - ry * L / 30.0;
    M(3, 3) = M(9, 9);
    M(1, 1) = M(7, 7);
    M(2, 2) = M(8, 8);
    M(6, 6) = M(12, 12);
    M(4, 4) = M(10, 10);
    M(5, 5) = M(11, 11);
    // the lower triangle mirrors the upper
    for (int i = 1; i <= 12; ++i) { for (int j = i + 1; j <= 12; ++j) { M(j, i) = M(i, j); } }
    return m;
}

std::array<double,144> beam_mass (const BeamProperties& p, double L, const std::array<double,9>& dc)
{
    return to_frame(beam_mass_local(p, L), dc);
}

std::array<double,36> rigid_body_mass (const FrameMass& c)
{
    const double m = c.mass, x = c.offset[0], y = c.offset[1], z = c.offset[2];
    const double jxx = c.inertia[0], jyy = c.inertia[1], jzz = c.inertia[2], jxy = c.inertia[3], jxz = c.inertia[4], jyz = c.inertia[5];
    return {{ m,      0.0,    0.0,    0.0,                       z * m,                     -y * m,
              0.0,    m,      0.0,   -z * m,                     0.0,                        x * m,
              0.0,    0.0,    m,      y * m,                    -x * m,                      0.0,
              0.0,   -z * m,  y * m,  jxx + m * (y * y + z * z), jxy - m * x * y,            jxz - m * x * z,
              z * m,  0.0,   -x * m,  jxy - m * x * y,           jyy + m * (x * x + z * z),  jyz - m * y * z,
             -y * m,  x * m,  0.0,    jxz - m * x * z,           jyz - m * y * z,            jzz + m * (x * x + y * y)}};
}

std::array<double,12> beam_gravity_load (const BeamProperties& p, double L, const std::array<double,9>& dc, double gravity)
{
    std::array<double,12> f{};
    const double w = p.rho * p.A * gravity;      // weight per length (N/m)
    const double m = L * L * w / 12.0;
    f[2] = -0.5 * L * w;
    f[8] = f[2];
    f[3] = -m * dc[5];                           // dc[5] = dy/L
    f[4] = m * dc[2];                            // dc[2] = dx/L
    f[9] = -f[3];
    f[10] = -f[4];
    return f;
}

std::unique_ptr<Frame> Frame::create (const FrameInputs& in, std::string& err)
{
    err = in.validate();
    if (!err.empty()) { return nullptr; }
    std::unique_ptr<Frame> f(new Frame());
    f->m_in = in;
    for (const auto& j : in.joints) {
        f->m_x.push_back(j.x);
        f->m_node_name.push_back("joint " + std::to_string(j.id));
    }
    for (std::size_t m = 0; m < in.members.size(); ++m) {
        const FrameMember& mem = in.members[m];
        const std::size_t a = static_cast<std::size_t>(in.joint_index(mem.joint_a));
        const std::size_t b = static_cast<std::size_t>(in.joint_index(mem.joint_b));
        const auto& xa = in.joints[a].x;
        const auto& xb = in.joints[b].x;
        const std::array<double,9> dc = direction_cosines(xa, xb, mem.spin);
        BeamProperties prop = beam_properties(*in.section(mem.section, mem.shape), in.theory);
        if (!in.temperature.empty()) {
            // the steel's stiffness at the member's temperature (EN 1993-1-2); G keeps its ratio to E
            double k_y = 1.0, k_E = 1.0;
            steel_reduction(in.temperature[m], k_y, k_E);
            prop.E *= k_E;
            prop.G *= k_E;
        }
        const int n = in.divisions;
        std::size_t prev = a;
        for (int k = 1; k <= n; ++k) {
            std::size_t node = b;
            if (k < n) {
                const double s = static_cast<double>(k) / n;
                f->m_x.push_back({{xa[0] + s * (xb[0] - xa[0]), xa[1] + s * (xb[1] - xa[1]), xa[2] + s * (xb[2] - xa[2])}});
                f->m_node_name.push_back("member " + std::to_string(mem.id) + " node " + std::to_string(k));
                node = f->m_x.size() - 1;
            }
            FrameElement e;
            e.node_a = prev;
            e.node_b = node;
            e.member = m;
            const auto& pa = f->m_x[prev];
            const auto& pb = f->m_x[node];
            e.length = std::sqrt((pb[0] - pa[0]) * (pb[0] - pa[0]) + (pb[1] - pa[1]) * (pb[1] - pa[1]) + (pb[2] - pa[2]) * (pb[2] - pa[2]));
            e.dc = dc;
            e.prop = prop;
            f->m_elems.push_back(e);
            prev = node;
        }
    }
    // the free degrees of freedom: every one a support does not fix
    const std::size_t ndof = f->num_dofs();
    std::vector<bool> fixed(ndof, false);
    for (const auto& s : in.supports) {
        const std::size_t node = static_cast<std::size_t>(in.joint_index(s.joint));
        for (std::size_t d = 0; d < 6; ++d) { if (s.fixed[d]) { fixed[6 * node + d] = true; } }
    }
    f->m_free_index.assign(ndof, -1);
    for (std::size_t d = 0; d < ndof; ++d) {
        if (!fixed[d]) { f->m_free_index[d] = static_cast<long>(f->m_free.size()); f->m_free.push_back(d); }
    }
    if (f->m_free.empty()) { err = in.file + ": every degree of freedom is fixed; nothing is left to solve"; return nullptr; }
    if (f->m_free.size() > max_free_dofs) {
        const double mb = 8.0 * static_cast<double>(f->m_free.size()) * static_cast<double>(f->m_free.size()) / 1.0e6;
        err = in.file + ": the frame has " + std::to_string(f->m_free.size()) + " free degrees of freedom, more than the " +
              std::to_string(max_free_dofs) + " the dense frame solver takes (each of its matrices would hold " +
              std::to_string(static_cast<long>(mb)) + " MB on every rank); use fewer panels, members or NDiv";
        return nullptr;
    }
    // the stiffness of the free degrees of freedom: the elements, then the support springs
    const std::size_t nf = f->m_free.size();
    const std::vector<double> kf = f->assemble_free(1.0, 0.0);
    const long bad = f->m_chol.factor(kf, nf);
    if (bad >= 0) {
        err = in.file + ": the frame is a mechanism: nothing holds the " + f->describe_dof(f->m_free[static_cast<std::size_t>(bad)]) +
              " (add a support, a member or a bracing there)";
        return nullptr;
    }
    return f;
}

std::vector<double> Frame::assemble_free (double cK, double cM) const
{
    const std::size_t nf = m_free.size();
    std::vector<double> a(nf * nf, 0.0);
    auto add = [&] (std::size_t gi, std::size_t gj, double v) {
        const long fi = m_free_index[gi], fj = m_free_index[gj];
        if (fi >= 0 && fj >= 0) { a[static_cast<std::size_t>(fi) * nf + static_cast<std::size_t>(fj)] += v; }
    };
    for (const auto& e : m_elems) {
        const std::array<double,144> ke = beam_stiffness(e.prop, e.length, e.dc);
        const std::array<double,144> me = (cM != 0.0) ? beam_mass(e.prop, e.length, e.dc) : std::array<double,144>{};
        for (int i = 0; i < 12; ++i) {
            for (int j = 0; j < 12; ++j) {
                const std::size_t k = static_cast<std::size_t>(12 * i + j);
                add(element_dof(e, i), element_dof(e, j), cK * ke[k] + cM * me[k]);
            }
        }
    }
    const auto entries = ssi_entries();
    for (const auto& s : m_in.supports) {
        const std::size_t node = static_cast<std::size_t>(m_in.joint_index(s.joint));
        for (std::size_t k = 0; k < 21; ++k) {
            const double v = cK * s.stiffness[k] + cM * s.mass[k];
            if (v == 0.0) { continue; }
            const std::size_t gi = 6 * node + static_cast<std::size_t>(entries[k].first);
            const std::size_t gj = 6 * node + static_cast<std::size_t>(entries[k].second);
            add(gi, gj, v);
            if (gi != gj) { add(gj, gi, v); }
        }
    }
    if (cM != 0.0) {
        for (const auto& c : m_in.masses) {
            const std::size_t node = static_cast<std::size_t>(m_in.joint_index(c.joint));
            const std::array<double,36> m66 = rigid_body_mass(c);
            for (std::size_t i = 0; i < 6; ++i) {
                for (std::size_t j = 0; j < 6; ++j) { add(6 * node + i, 6 * node + j, cM * m66[6 * i + j]); }
            }
        }
    }
    return a;
}

std::vector<double> Frame::apply_stiffness (const std::vector<double>& u, bool with_springs) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(u.size() == num_dofs(), "Frame::apply_stiffness: one value per degree of freedom is needed");
    std::vector<double> r(u.size(), 0.0);
    for (const auto& e : m_elems) {
        const std::array<double,144> ke = beam_stiffness(e.prop, e.length, e.dc);
        for (int i = 0; i < 12; ++i) {
            double s = 0.0;
            for (int j = 0; j < 12; ++j) { s += ke[static_cast<std::size_t>(12 * i + j)] * u[element_dof(e, j)]; }
            r[element_dof(e, i)] += s;
        }
    }
    if (!with_springs) { return r; }
    const auto entries = ssi_entries();
    for (const auto& sp : m_in.supports) {
        const std::size_t node = static_cast<std::size_t>(m_in.joint_index(sp.joint));
        for (std::size_t k = 0; k < 21; ++k) {
            const double v = sp.stiffness[k];
            if (v == 0.0) { continue; }
            const std::size_t gi = 6 * node + static_cast<std::size_t>(entries[k].first);
            const std::size_t gj = 6 * node + static_cast<std::size_t>(entries[k].second);
            r[gi] += v * u[gj];
            if (gi != gj) { r[gj] += v * u[gi]; }
        }
    }
    return r;
}

std::vector<double> Frame::apply_mass (const std::vector<double>& a, bool with_springs) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a.size() == num_dofs(), "Frame::apply_mass: one value per degree of freedom is needed");
    std::vector<double> r(a.size(), 0.0);
    for (const auto& e : m_elems) {
        const std::array<double,144> me = beam_mass(e.prop, e.length, e.dc);
        for (int i = 0; i < 12; ++i) {
            double s = 0.0;
            for (int j = 0; j < 12; ++j) { s += me[static_cast<std::size_t>(12 * i + j)] * a[element_dof(e, j)]; }
            r[element_dof(e, i)] += s;
        }
    }
    if (with_springs) {
        const auto entries = ssi_entries();
        for (const auto& sp : m_in.supports) {
            const std::size_t node = static_cast<std::size_t>(m_in.joint_index(sp.joint));
            for (std::size_t k = 0; k < 21; ++k) {
                const double v = sp.mass[k];
                if (v == 0.0) { continue; }
                const std::size_t gi = 6 * node + static_cast<std::size_t>(entries[k].first);
                const std::size_t gj = 6 * node + static_cast<std::size_t>(entries[k].second);
                r[gi] += v * a[gj];
                if (gi != gj) { r[gj] += v * a[gi]; }
            }
        }
    }
    for (const auto& c : m_in.masses) {
        const std::size_t node = static_cast<std::size_t>(m_in.joint_index(c.joint));
        const std::array<double,36> m66 = rigid_body_mass(c);
        for (std::size_t i = 0; i < 6; ++i) {
            for (std::size_t j = 0; j < 6; ++j) { r[6 * node + i] += m66[6 * i + j] * a[6 * node + j]; }
        }
    }
    return r;
}

std::vector<double> Frame::point_loads (const std::vector<double>& node_loads, double gravity) const
{
    const std::size_t ndof = num_dofs();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node_loads.empty() || node_loads.size() == ndof, "Frame: 6 loads per node are needed");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(gravity) && gravity >= 0.0, "Frame: gravity must be finite and >= 0 (m/s^2)");
    std::vector<double> point(ndof, 0.0);
    for (std::size_t d = 0; d < node_loads.size(); ++d) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::isfinite(node_loads[d]), "Frame: a node load is not finite");
        point[d] = node_loads[d];
    }
    if (gravity > 0.0) {
        for (const auto& c : m_in.masses) {
            const std::size_t node = static_cast<std::size_t>(m_in.joint_index(c.joint));
            const double fz = -c.mass * gravity;
            point[6 * node + 2] += fz;
            point[6 * node + 3] += c.offset[1] * fz;     // r x F with F = (0, 0, fz)
            point[6 * node + 4] -= c.offset[0] * fz;
        }
    }
    return point;
}

std::vector<double> Frame::load_vector (const std::vector<double>& node_loads, double gravity) const
{
    std::vector<double> f = point_loads(node_loads, gravity);
    if (gravity > 0.0) {
        for (const auto& e : m_elems) {
            const std::array<double,12> fg = beam_gravity_load(e.prop, e.length, e.dc, gravity);
            for (int i = 0; i < 12; ++i) { f[element_dof(e, i)] += fg[static_cast<std::size_t>(i)]; }
        }
    }
    return f;
}

std::string Frame::describe_dof (std::size_t dof) const
{
    static const char* names[6] = {"x translation", "y translation", "z translation",
                                   "rotation about x", "rotation about y", "rotation about z"};
    return names[dof % 6] + std::string(" of ") + m_node_name[dof / 6];
}

const std::array<double,3>& Frame::node_position (std::size_t n) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n < m_x.size(), "Frame::node_position: no such node");
    return m_x[n];
}

std::size_t Frame::node_of_joint (int id) const
{
    const int i = m_in.joint_index(id);
    if (i < 0) { amrex::Abort("Frame: " + m_in.file + " has no joint " + std::to_string(id)); }
    return static_cast<std::size_t>(i);
}

double Frame::total_mass () const
{
    double m = 0.0;
    for (const auto& e : m_elems) { m += e.prop.rho * e.prop.A * e.length; }
    for (const auto& c : m_in.masses) { m += c.mass; }
    return m;
}

std::vector<double> Frame::stiffness () const
{
    const std::size_t n = num_dofs();
    std::vector<double> k(n * n, 0.0);
    for (const auto& e : m_elems) {
        const std::array<double,144> ke = beam_stiffness(e.prop, e.length, e.dc);
        for (int i = 0; i < 12; ++i) {
            const std::size_t gi = element_dof(e, i);
            for (int j = 0; j < 12; ++j) {
                const std::size_t gj = element_dof(e, j);
                k[gi * n + gj] += ke[static_cast<std::size_t>(12 * i + j)];
            }
        }
    }
    return k;
}

std::array<double,12> Frame::end_forces (const FrameElement& el, const std::vector<double>& u, double gravity) const
{
    const std::array<double,144> ke = beam_stiffness(el.prop, el.length, el.dc);
    const std::array<double,12> fg = (gravity > 0.0) ? beam_gravity_load(el.prop, el.length, el.dc, gravity) : std::array<double,12>{};
    std::array<double,12> fe{};
    for (int i = 0; i < 12; ++i) {
        double s = -fg[static_cast<std::size_t>(i)];
        for (int j = 0; j < 12; ++j) { s += ke[static_cast<std::size_t>(12 * i + j)] * u[element_dof(el, j)]; }
        fe[static_cast<std::size_t>(i)] = s;
    }
    return fe;
}

std::vector<std::array<double,12>> Frame::element_forces (const std::vector<double>& u, double gravity) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(u.size() == num_dofs(), "Frame::element_forces: one displacement per degree of freedom is needed");
    std::vector<std::array<double,12>> out(m_elems.size());
    for (std::size_t e = 0; e < m_elems.size(); ++e) { out[e] = rotate12(m_elems[e].dc, end_forces(m_elems[e], u, gravity), true); }
    return out;
}

FrameSolution Frame::solve (const std::vector<double>& node_loads, double gravity) const
{
    const std::size_t ndof = num_dofs();
    // the loads applied at the nodes, and with them the members' weight as element loads
    const std::vector<double> point = point_loads(node_loads, gravity);
    const std::vector<double> f = load_vector(node_loads, gravity);
    std::vector<double> b(m_free.size());
    for (std::size_t i = 0; i < m_free.size(); ++i) { b[i] = f[m_free[i]]; }
    m_chol.solve(b);
    FrameSolution sol;
    sol.displacement.assign(ndof, 0.0);
    for (std::size_t i = 0; i < m_free.size(); ++i) { sol.displacement[m_free[i]] = b[i]; }
    // each element's end forces, and their sum at every node (the force the elements take from it)
    std::vector<double> taken(ndof, 0.0);
    sol.element_force.resize(m_elems.size());
    for (std::size_t e = 0; e < m_elems.size(); ++e) {
        const FrameElement& el = m_elems[e];
        const std::array<double,12> fe = end_forces(el, sol.displacement, gravity);
        for (int i = 0; i < 12; ++i) { taken[element_dof(el, i)] += fe[static_cast<std::size_t>(i)]; }
        sol.element_force[e] = rotate12(el.dc, fe, true);
    }
    // a support's reaction: what its node's elements take from it less what is applied there (springs included)
    for (const auto& s : m_in.supports) {
        const std::size_t node = static_cast<std::size_t>(m_in.joint_index(s.joint));
        std::array<double,6> r{};
        for (std::size_t d = 0; d < 6; ++d) { r[d] = taken[6 * node + d] - point[6 * node + d]; }
        sol.reaction.push_back(r);
    }
    return sol;
}

} // namespace erf_towers
