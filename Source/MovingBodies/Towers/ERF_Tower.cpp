// Tower: drag nodes, line loads and the foundation's reactions.

#include "ERF_Tower.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <utility>

#include <AMReX.H>
#include <AMReX_BLassert.H>

using amrex::Real;

namespace erf_towers {

namespace {

/** The face area of a piece of the cross-arm and the centroid of that area, from arm_segment(). */
struct ArmPiece
{
    Real area{0.0};   //!< m^2
    Real s{0.0};      //!< the centroid's distance from the shaft's axis along the cross-arm (m)
    Real drop{0.0};   //!< the centroid's depth below the cross-arm's top (m)
};

/**
 * The face area between s = a and s = b (a < b, both within the cross-arm) and its centroid, for the depth d(s)
 * of Tower::arm_face_at(): linear in s between the breaks at the tips, where the taper starts on either side
 * (Tower::arm_taper_start()) and 0, so on each piece between breaks d runs from d1 to d2 and the integrals of d,
 * s d and d^2 are closed forms. The centroid lies at (int s d ds / A, int d^2 / 2 ds / A) below the top,
 * A = int d ds. The tips are breaks too, so a b an ulp past the tip adds no length.
 */
ArmPiece arm_segment (const Tower& tower, Real a, Real b)
{
    const Real s0 = tower.arm_taper_start(), half_len = Real(0.5) * tower.type().arm_length;
    std::array<Real,7> at{{a, -half_len, -s0, Real(0.0), s0, half_len, b}};
    std::sort(at.begin(), at.end());
    Real area = 0.0, moment_s = 0.0, moment_d2 = 0.0;
    for (std::size_t i = 0; i + 1 < at.size(); ++i) {
        const Real p = std::max(a, at[i]), q = std::min(b, at[i + 1]);
        if (!(q > p)) { continue; }
        const Real h = q - p, d1 = tower.arm_face_at(p), d2 = tower.arm_face_at(q);
        area += h * (d1 + d2) / Real(2.0);
        moment_s += p * h * (d1 + d2) / Real(2.0) + h * h * (d1 + Real(2.0) * d2) / Real(6.0);
        moment_d2 += h * (d1 * d1 + d1 * d2 + d2 * d2) / Real(3.0);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(area > Real(0.0), "Tower: a cross-arm segment has no face area");
    return ArmPiece{area, moment_s / area, Real(0.5) * moment_d2 / area};
}

} // namespace

Real
Tower::arm_taper_start () const
{
    // half the shaft's width at the cross-arm's bottom, H - d0 (the base when the face is deeper than the arm is
    // high, which arm_outside_shaft refuses): where the generated frame's bottom chords leave the shaft's corners
    const Real zb = std::max(Real(0.0), m_arm_height - m_type.arm_face());
    return Real(0.5) * (m_type.base_width + (m_type.top_width - m_type.base_width) * zb / m_arm_height);
}

Real
Tower::arm_face_at (Real s) const
{
    const Real half_len = Real(0.5) * m_type.arm_length, s0 = arm_taper_start();
    const Real d0 = m_type.arm_face();
    const Real r = std::abs(s);
    if (r >= half_len) { return Real(0.0); }
    if (r <= s0) { return d0; }
    // here half_len > s0, so the taper's length is positive
    return d0 * (half_len - r) / (half_len - s0);
}

Tower::Tower (std::string name, const TowerType& type, const std::array<Real,3>& base, Real arm_height,
              const std::array<Real,3>& across)
    : m_name(std::move(name)), m_type(type), m_base(base), m_arm_height(arm_height), m_across(across)
{
    // preconditions the callers establish: a validated type, a cross-arm above the base, a horizontal unit cross-arm
    const std::string terr = m_type.validate();
    if (!terr.empty()) { amrex::Abort("Tower " + m_name + ": " + terr); }
    if (!(std::isfinite(m_arm_height) && m_arm_height > 0.0)) {
        amrex::Abort("Tower " + m_name + ": the cross-arm height above the base (" + std::to_string(m_arm_height) + " m) must be positive");
    }
    const Real across_norm = std::sqrt(m_across[0] * m_across[0] + m_across[1] * m_across[1] + m_across[2] * m_across[2]);
    if (!(std::abs(across_norm - Real(1.0)) < Real(1.0e-4) && std::abs(m_across[2]) < Real(1.0e-4))) {
        amrex::Abort("Tower " + m_name + ": the cross-arm direction must be a horizontal unit vector");
    }
    const Real phi = m_type.solidity;
    const Real cf = m_type.force_coefficient();
    // the shaft is a square lattice whose faces lie along the cross-arm and across it: wind along a diagonal
    // meets more of its members than wind normal to a face
    const Real gain = m_type.diagonal_wind_factor ? lattice_diagonal_gain(phi) : Real(0.0);
    // the shaft up to the cross-arm, tapering, then the peak at the top width
    const int n = m_type.segments;
    const Real dz = m_arm_height / n;
    for (int i = 0; i < n; ++i) {
        const Real z = (Real(i) + Real(0.5)) * dz;
        const Real w = m_type.base_width + (m_type.top_width - m_type.base_width) * z / m_arm_height;
        m_nodes.push_back(MemberNode{{{m_base[0], m_base[1], m_base[2] + z}}, {{0.0, 0.0, 1.0}}, dz, phi * w, cf, m_across, gain});
    }
    if (m_type.peak > 0.0) {
        const int np = std::max(1, static_cast<int>(std::ceil(m_type.peak / dz)));
        const Real dp = m_type.peak / np;
        for (int i = 0; i < np; ++i) {
            const Real z = m_arm_height + (Real(i) + Real(0.5)) * dp;
            m_nodes.push_back(MemberNode{{{m_base[0], m_base[1], m_base[2] + z}}, {{0.0, 0.0, 1.0}}, dp, phi * m_type.top_width, cf,
                                         m_across, gain});
        }
    }
    m_nbody = static_cast<int>(m_nodes.size());
    // the cross-arm, centred on the shaft: with arm_outside_shaft its two parts outside the shaft (the shaft's
    // nodes already carry the drag where the arm passes through it, the shaft as wide there as at the middle
    // of the cross-arm's face at the shaft), each node at the centroid of its segment's face, which is arm_face()
    // deep out to the shaft's half width at the arm's bottom and tapers to the tip (arm_face_at), with the face's
    // mean depth over the segment; else the earlier layout, along its whole length at the cross-arm's height,
    // arm_face() deep all along
    const int na = TowerType::arm_segments;
    const bool outside = m_type.arm_outside_shaft;
    const std::string key = "erf.conductors." + m_type.name + ".";
    const Real za = m_arm_height - (outside ? Real(0.5) * m_type.arm_face() : Real(0.0));
    if (!(za > Real(0.0))) {
        amrex::Abort("Tower " + m_name + ": the middle of the cross-arm's face lies " + std::to_string(-za) + " m below the base: " + key +
                     "arm_depth (or top_width, when arm_depth is 0) must be under twice the cross-arm's height, or set " + key +
                     "arm_outside_shaft = false");
    }
    const Real inner = outside ? Real(0.5) * (m_type.base_width + (m_type.top_width - m_type.base_width) * za / m_arm_height) : Real(0.0);
    if (!(Real(0.5) * m_type.arm_length > inner)) {
        amrex::Abort("Tower " + m_name + ": " + key + "arm_length (" + std::to_string(m_type.arm_length) + " m) must exceed the shaft's width " +
                     "at the middle of the cross-arm's face (" + std::to_string(2.0 * inner) + " m) for " + key + "arm_outside_shaft");
    }
    auto arm_node = [&] (Real s, Real z, Real length, Real depth) {
        m_nodes.push_back(MemberNode{{{m_base[0] + s * m_across[0], m_base[1] + s * m_across[1], m_base[2] + z}},
                                     m_across, length, phi * depth, cf});
    };
    if (outside) {
        // half the nodes on each side, from the tip inwards on the first, outwards on the second: the nodes run
        // from one tip to the other; each segment's node at its face's centroid with the face's mean depth
        const int nside = na / 2;
        const Real da = (Real(0.5) * m_type.arm_length - inner) / nside;
        auto piece_node = [&] (Real a, Real b) {
            const ArmPiece p = arm_segment(*this, a, b);
            arm_node(p.s, m_arm_height - p.drop, b - a, p.area / (b - a));
        };
        for (int i = 0; i < nside; ++i) {
            const Real b = -(inner + Real(nside - 1 - i) * da);
            piece_node(b - da, b);
        }
        for (int i = 0; i < nside; ++i) {
            const Real a = inner + Real(i) * da;
            piece_node(a, a + da);
        }
    } else {
        const Real da = m_type.arm_length / na;
        for (int i = 0; i < na; ++i) {
            arm_node(-Real(0.5) * m_type.arm_length + (Real(i) + Real(0.5)) * da, za, da, m_type.arm_face());
        }
    }
    m_force.assign(3 * m_nodes.size(), 0.0);
    m_disp.assign(3 * m_nodes.size(), 0.0);
    m_vel.assign(3 * m_nodes.size(), 0.0);
    m_inertia.assign(3 * m_nodes.size(), 0.0);
}

void Tower::set_motion (const std::vector<Real>& displacement, const std::vector<Real>& velocity, const std::vector<Real>& inertia)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(displacement.size() == 3 * m_nodes.size() && velocity.size() == 3 * m_nodes.size() &&
                                     inertia.size() == 3 * m_nodes.size(), "Tower::set_motion: 3 values per drag node are needed");
    m_disp = displacement;
    m_vel = velocity;
    m_inertia = inertia;
}

std::size_t Tower::add_attachment (const std::array<Real,3>& at)
{
    m_attach.push_back(at);
    return m_attach.size() - 1;
}

void Tower::set_loads (const std::vector<Real>& force)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(force.size() == 3 * m_nodes.size(), "Tower::set_loads: 3 forces per drag node are needed");
    m_force = force;
}

void Tower::set_line_loads (const std::vector<std::array<Real,3>>& force, const std::vector<std::array<Real,3>>& at)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(force.size() == at.size(), "Tower::set_line_loads: one point per force is needed");
    const bool fits = m_attach.empty() ? force.size() <= 1 : force.size() == m_attach.size();
    if (!fits) {
        amrex::Abort("Tower " + m_name + ": " + std::to_string(force.size()) + " line pulls for " + std::to_string(m_attach.size()) +
                     " attachments; one pull per attachment is needed");
    }
    for (std::size_t a = 0; a < force.size(); ++a) {
        for (int d = 0; d < 3; ++d) {
            if (!std::isfinite(force[a][d]) || !std::isfinite(at[a][d])) {
                amrex::Abort("Tower " + m_name + ": the pull of the line at attachment " + std::to_string(a) +
                             " (0-based) or its point is not finite; MoorDyn's line integration may have diverged: "
                             "reduce erf.conductors.moordyn_cfl or erf.conductors.moordyn_dt");
            }
        }
    }
    m_line_force = force;
    m_line_at = at;
}

std::array<Real,3> Tower::line_force () const
{
    std::array<Real,3> F{{0.0, 0.0, 0.0}};
    for (const auto& f : m_line_force) { for (int d = 0; d < 3; ++d) { F[d] += f[d]; } }
    return F;
}

std::vector<MemberNode> Tower::current_nodes () const
{
    std::vector<MemberNode> nodes = m_nodes;
    for (std::size_t i = 0; i < nodes.size(); ++i) { for (int d = 0; d < 3; ++d) { nodes[i].pos[d] += m_disp[3*i+d]; } }
    return nodes;
}

std::array<Real,3> Tower::arm_displacement () const
{
    // the cross-arm's drag nodes follow the shaft's
    const std::size_t i = static_cast<std::size_t>(m_nbody);
    return {{m_disp[3*i], m_disp[3*i+1], m_disp[3*i+2]}};
}

std::array<Real,3> Tower::total_force () const
{
    std::array<Real,3> F{{0.0, 0.0, 0.0}};
    for (std::size_t i = 0; i < m_nodes.size(); ++i) { for (int d = 0; d < 3; ++d) { F[d] += m_force[3*i+d]; } }
    return F;
}

std::array<Real,3> Tower::base_moment () const
{
    std::array<Real,3> M{{0.0, 0.0, 0.0}};
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {
        const Real rx = m_nodes[i].pos[0] - m_base[0], ry = m_nodes[i].pos[1] - m_base[1], rz = m_nodes[i].pos[2] - m_base[2];
        const Real fx = m_force[3*i], fy = m_force[3*i+1], fz = m_force[3*i+2];
        M[0] += ry * fz - rz * fy;
        M[1] += rz * fx - rx * fz;
        M[2] += rx * fy - ry * fx;
    }
    return M;
}

FoundationLoad Tower::foundation () const
{
    FoundationLoad L;
    auto Fd = total_force();
    auto Md = base_moment();
    // a moving tower's nodes also carry their inertial forces to the foundation, the cross-arm's from its height,
    // where its mass is (OneModeTower), wherever its drag nodes stand
    for (std::size_t i = 0; i < m_nodes.size(); ++i) {
        const Real rx = m_nodes[i].pos[0] - m_base[0], ry = m_nodes[i].pos[1] - m_base[1];
        const Real rz = ((static_cast<int>(i) < m_nbody) ? m_nodes[i].pos[2] : m_base[2] + m_arm_height) - m_base[2];
        const Real fx = m_inertia[3*i], fy = m_inertia[3*i+1], fz = m_inertia[3*i+2];
        Fd[0] += fx; Fd[1] += fy; Fd[2] += fz;
        Md[0] += ry * fz - rz * fy;
        Md[1] += rz * fx - rx * fz;
        Md[2] += rx * fy - ry * fx;
    }
    L.force = Fd;
    L.moment = Md;
    // each line's pull where it hangs from the tower
    for (std::size_t a = 0; a < m_line_force.size(); ++a) {
        const auto& F = m_line_force[a];
        const Real rx = m_line_at[a][0] - m_base[0], ry = m_line_at[a][1] - m_base[1], rz = m_line_at[a][2] - m_base[2];
        for (int d = 0; d < 3; ++d) { L.force[d] += F[d]; }
        L.moment[0] += ry * F[2] - rz * F[1];
        L.moment[1] += rz * F[0] - rx * F[2];
        L.moment[2] += rx * F[1] - ry * F[0];
    }
    L.shear = std::hypot(L.force[0], L.force[1]);
    L.overturning = std::hypot(L.moment[0], L.moment[1]);
    L.vertical = m_type.weight - L.force[2];
    // the legs in the frame of the line (along) and the cross-arm (across), a = s / 2 from the centre
    const Real a = Real(0.5) * m_type.legs();
    const std::array<Real,2> across{{m_across[0], m_across[1]}};
    const std::array<Real,2> along{{m_across[1], -m_across[0]}};
    int i = 0;
    for (const Real su : {Real(1.0), Real(-1.0)}) {
        for (const Real sv : {Real(1.0), Real(-1.0)}) {
            const Real x = a * (su * along[0] + sv * across[0]);
            const Real y = a * (su * along[1] + sv * across[1]);
            // the leg compressions R_i (N) balance the load: sum R_i = P, sum x_i R_i = M_y, sum y_i R_i = -M_x,
            // x_i and y_i the leg offsets in the ERF frame
            L.legs[static_cast<std::size_t>(i++)] = Real(0.25) * L.vertical + (L.moment[1] * x - L.moment[0] * y) / (Real(4.0) * a * a);
        }
    }
    L.max_compression = std::max({L.legs[0], L.legs[1], L.legs[2], L.legs[3], Real(0.0)});
    L.max_uplift = std::max(Real(0.0), -std::min({L.legs[0], L.legs[1], L.legs[2], L.legs[3]}));
    L.over_allowable = (m_type.allowable_uplift > 0.0 && L.max_uplift > m_type.allowable_uplift) ||
                       (m_type.allowable_compression > 0.0 && L.max_compression > m_type.allowable_compression);
    return L;
}

} // namespace erf_towers
