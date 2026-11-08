// A lattice tower that bends as its frame model: rigid links between the tower's points and the
// frame's nodes, the Newmark step, and the footings' loads from the support reactions.

#include "ERF_FrameTower.H"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>

#include <AMReX.H>
#include <AMReX_BLassert.H>

using amrex::Real;

namespace erf_towers {

namespace {

std::array<double,3> cross (const std::array<double,3>& a, const std::array<double,3>& b)
{
    return {{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}};
}

/** The pseudo-inverse of a symmetric positive semi-definite 3 x 3 matrix (row-major), by Jacobi rotations; eigenvalues below 1e-9 of the largest are dropped. */
std::array<double,9> pseudo_inverse (std::array<double,9> a)
{
    std::array<double,9> v{{1, 0, 0, 0, 1, 0, 0, 0, 1}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        const double off = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
        if (off <= 1.0e-30 * (a[0] * a[0] + a[4] * a[4] + a[8] * a[8] + 1.0e-300)) { break; }
        for (int p = 0; p < 3; ++p) {
            for (int q = p + 1; q < 3; ++q) {
                const double apq = a[static_cast<std::size_t>(3 * p + q)];
                if (apq == 0.0) { continue; }
                const double theta = (a[static_cast<std::size_t>(4 * q)] - a[static_cast<std::size_t>(4 * p)]) / (2.0 * apq);
                const double t = (theta >= 0.0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
                const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const double akp = a[static_cast<std::size_t>(3 * k + p)], akq = a[static_cast<std::size_t>(3 * k + q)];
                    a[static_cast<std::size_t>(3 * k + p)] = c * akp - s * akq;
                    a[static_cast<std::size_t>(3 * k + q)] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const double apk = a[static_cast<std::size_t>(3 * p + k)], aqk = a[static_cast<std::size_t>(3 * q + k)];
                    a[static_cast<std::size_t>(3 * p + k)] = c * apk - s * aqk;
                    a[static_cast<std::size_t>(3 * q + k)] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const double vkp = v[static_cast<std::size_t>(3 * k + p)], vkq = v[static_cast<std::size_t>(3 * k + q)];
                    v[static_cast<std::size_t>(3 * k + p)] = c * vkp - s * vkq;
                    v[static_cast<std::size_t>(3 * k + q)] = s * vkp + c * vkq;
                }
            }
        }
    }
    const double lmax = std::max({a[0], a[4], a[8], 0.0});
    std::array<double,9> r{};
    for (int k = 0; k < 3; ++k) {
        const double l = a[static_cast<std::size_t>(4 * k)];
        if (!(l > 1.0e-9 * lmax)) { continue; }
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                r[static_cast<std::size_t>(3 * i + j)] += v[static_cast<std::size_t>(3 * i + k)] * v[static_cast<std::size_t>(3 * j + k)] / l;
            }
        }
    }
    return r;
}

std::array<double,3> times (const std::array<double,9>& m, const std::array<double,3>& x)
{
    return {{m[0] * x[0] + m[1] * x[1] + m[2] * x[2], m[3] * x[0] + m[4] * x[1] + m[5] * x[2], m[6] * x[0] + m[7] * x[1] + m[8] * x[2]}};
}

} // namespace

RigidLink::RigidLink (const Frame& frame, const std::array<double,3>& p, std::size_t count, const std::vector<std::size_t>* candidates)
{
    std::vector<std::size_t> order;
    if (candidates != nullptr && !candidates->empty()) {
        order = *candidates;
    } else {
        order.resize(frame.num_nodes());
        std::iota(order.begin(), order.end(), std::size_t(0));
    }
    const std::size_t n = order.size();
    std::vector<double> dist(frame.num_nodes(), 0.0);
    for (const std::size_t i : order) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(i < frame.num_nodes(), "RigidLink: a candidate node is not a node of the frame");
        const auto& x = frame.node_position(i);
        dist[i] = std::sqrt((x[0] - p[0]) * (x[0] - p[0]) + (x[1] - p[1]) * (x[1] - p[1]) + (x[2] - p[2]) * (x[2] - p[2]));
    }
    std::stable_sort(order.begin(), order.end(), [&] (std::size_t a, std::size_t b) { return dist[a] < dist[b]; });
    m_distance = dist[order[0]];
    const std::size_t k = (m_distance <= 1.0e-6) ? 1 : std::min(count, n);
    m_nodes.assign(order.begin(), order.begin() + static_cast<long>(k));
    std::array<double,3> c{{0.0, 0.0, 0.0}};
    for (const std::size_t i : m_nodes) { for (int d = 0; d < 3; ++d) { c[static_cast<std::size_t>(d)] += frame.node_position(i)[static_cast<std::size_t>(d)] / static_cast<double>(k); } }
    std::array<double,9> inertia{};
    for (const std::size_t i : m_nodes) {
        const auto& x = frame.node_position(i);
        const std::array<double,3> d{{x[0] - c[0], x[1] - c[1], x[2] - c[2]}};
        m_d.push_back(d);
        const double d2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                inertia[static_cast<std::size_t>(3 * a + b)] += (a == b ? d2 : 0.0) - d[static_cast<std::size_t>(a)] * d[static_cast<std::size_t>(b)];
            }
        }
    }
    m_r = {{p[0] - c[0], p[1] - c[1], p[2] - c[2]}};
    m_iinv = pseudo_inverse(inertia);
}

void RigidLink::add_load (const std::array<double,3>& force, std::vector<double>& loads) const
{
    const std::size_t k = m_nodes.size();
    const std::array<double,3> w = times(m_iinv, cross(m_r, force));
    for (std::size_t j = 0; j < k; ++j) {
        const std::array<double,3> extra = cross(w, m_d[j]);
        for (std::size_t d = 0; d < 3; ++d) { loads[6 * m_nodes[j] + d] += force[d] / static_cast<double>(k) + extra[d]; }
    }
}

std::array<double,3> RigidLink::motion (const std::vector<double>& u) const
{
    const std::size_t k = m_nodes.size();
    std::array<double,3> mean{{0.0, 0.0, 0.0}};
    for (const std::size_t i : m_nodes) { for (std::size_t d = 0; d < 3; ++d) { mean[d] += u[6 * i + d] / static_cast<double>(k); } }
    std::array<double,3> s{{0.0, 0.0, 0.0}};
    for (std::size_t j = 0; j < k; ++j) {
        const std::array<double,3> du{{u[6 * m_nodes[j]] - mean[0], u[6 * m_nodes[j] + 1] - mean[1], u[6 * m_nodes[j] + 2] - mean[2]}};
        const std::array<double,3> t = cross(m_d[j], du);
        for (std::size_t d = 0; d < 3; ++d) { s[d] += t[d]; }
    }
    const std::array<double,3> theta = times(m_iinv, s);
    const std::array<double,3> turn = cross(theta, m_r);
    return {{mean[0] + turn[0], mean[1] + turn[1], mean[2] + turn[2]}};
}

FrameTower::FrameTower (const Tower& tower, std::shared_ptr<const Frame> frame, double gravity,
                        std::shared_ptr<const std::vector<MemberDesign>> designs, const std::string& source,
                        const std::vector<std::size_t>& link_nodes, double first_frequency)
    : m_frame(std::move(frame)), m_name(tower.name()), m_designs(std::move(designs)), m_gravity(gravity)
{
    const TowerType& type = tower.type();
    const std::string key = source.empty() ? "erf.conductors." + type.name + ".frame_file = " + type.frame_file : source;
    m_source = key;
    const auto& ac = tower.across();
    m_across = {{static_cast<double>(ac[0]), static_cast<double>(ac[1]), 0.0}};
    m_along = {{m_across[1], -m_across[0], 0.0}};
    const auto& base = tower.base();
    auto local_point = [&] (const std::array<Real,3>& p) {
        return to_local({{p[0] - base[0], p[1] - base[1], p[2] - base[2]}});
    };
    const double reach = static_cast<double>(type.base_width);
    for (std::size_t i = 0; i < tower.nodes().size(); ++i) {
        m_drag.emplace_back(*m_frame, local_point(tower.nodes()[i].pos), 4, &link_nodes);
        if (!(m_drag.back().distance() <= reach)) {
            amrex::Abort("tower " + m_name + ": drag node " + std::to_string(i) + " lies " + std::to_string(m_drag.back().distance()) +
                         " m from the nearest node of " + key + ", more than the base width; the frame's axes must be "
                         "tower-local (origin at the base centre, x along the line, y along the cross-arm, z up)");
        }
    }
    for (std::size_t a = 0; a < tower.attachments().size(); ++a) {
        m_attach.emplace_back(*m_frame, local_point(tower.attachments()[a]), 4, &link_nodes);
        if (!(m_attach.back().distance() <= reach)) {
            amrex::Abort("tower " + m_name + ": the line attachment " + std::to_string(a) + " lies " +
                         std::to_string(m_attach.back().distance()) + " m from the nearest node of " + key +
                         ", more than the base width; the cross-arm of the frame must be where the lines hang");
        }
    }
    {
        // the frame's cross-arm (its interface joint, the cross-arm's centre) at the tower's cross-arm height:
        // a frame built for another height would tie the lines to the shaft
        const auto& iface = m_frame->inputs().interface_joints;
        if (iface.size() != 1) {
            amrex::Abort("tower " + m_name + ": " + key + " has " + std::to_string(iface.size()) + " interface joints; a tower's "
                         "frame has one, the cross-arm's centre, where the lines' pulls are tied");
        }
        const double zi = m_frame->node_position(m_frame->node_of_joint(iface.front()))[2];
        const double za = static_cast<double>(tower.arm_height());
        if (!(std::abs(zi - za) <= 0.5 * static_cast<double>(type.arm_face()))) {
            amrex::Abort("tower " + m_name + ": its cross-arm stands " + std::to_string(za) + " m above the base, and the "
                         "cross-arm of " + key + " (its interface joint " + std::to_string(iface.front()) + ") " + std::to_string(zi) +
                         " m, more than half the cross-arm's face depth apart: the frame's interface joint must be its "
                         "cross-arm's centre (not the peak, as SubDyn's own cases have it), and a frame is built for one "
                         "cross-arm height; give a frame file per tower height, or frame_panels");
        }
    }
    // the cross-arm's centre, on the tower's axis, for the cross-arm's displacement in the logs
    m_arm.emplace_back(*m_frame, std::array<double,3>{{0.0, 0.0, static_cast<double>(tower.arm_height())}}, 4, &link_nodes);
    // the four supports, one per quadrant of the base, at z = 0
    const auto& sup = m_frame->inputs().supports;
    if (sup.size() != 4) { amrex::Abort("tower " + m_name + ": " + key + " has " + std::to_string(sup.size()) + " supports; a lattice tower stands on four legs"); }
    std::array<int,4> seen{{-1, -1, -1, -1}};
    for (std::size_t s = 0; s < 4; ++s) {
        const auto& x = m_frame->node_position(m_frame->node_of_joint(sup[s].joint));
        if (!(std::abs(x[2]) <= 1.0e-6 && x[0] != 0.0 && x[1] != 0.0)) {
            amrex::Abort("tower " + m_name + ": support joint " + std::to_string(sup[s].joint) + " of " + key +
                         " must stand at z = 0 off both axes (a leg at a corner of the base)");
        }
        // Tower::foundation()'s order: (+along, +across), (+along, -across), (-along, +across), (-along, -across)
        const std::size_t leg = (x[0] > 0.0 ? 0 : 2) + (x[1] > 0.0 ? 0 : 1);
        if (seen[leg] >= 0) { amrex::Abort("tower " + m_name + ": " + key + " has two supports in one quadrant of the base"); }
        seen[leg] = static_cast<int>(s);
        m_leg_support[leg] = s;
    }
    const std::size_t nm = m_frame->inputs().members.size();
    if (m_designs && m_designs->size() != nm) {
        amrex::Abort("tower " + m_name + ": " + key + " has " + std::to_string(nm) + " members but " +
                     std::to_string(m_designs->size()) + " rows of design data");
    }
    m_theta = m_frame->inputs().temperature;
    if (m_theta.empty()) { m_theta.assign(nm, 20.0); }
    // the static equilibrium under the frame's weight, about which it moves
    const FrameSolution s0 = m_frame->solve({}, gravity);
    m_static = s0.reaction;
    m_static_force = s0.element_force;
    m_sag = s0.displacement;
    m_sag0 = m_sag;
    if (first_frequency > 0.0) {
        m_f1 = first_frequency;
    } else {
        FrameModes modes;
        const std::string err = frame_modes(*m_frame, 1, modes);
        if (!err.empty()) { amrex::Abort("tower " + m_name + ": " + key + ": " + err); }
        m_f1 = modes.frequency[0];
    }
    double a0 = 0.0, a1 = 0.0;
    const double zeta = static_cast<double>(type.damping_ratio);
    if (zeta > 0.0) { rayleigh_coefficients(m_f1, zeta, 10.0 * m_f1, zeta, a0, a1); }
    m_a0 = a0;
    m_a1 = a1;
    m_dyn = std::make_unique<FrameDynamics>(*m_frame, a0, a1);
    m_dyn->start_static({}, 0.0);
    m_load.assign(m_frame->num_dofs(), 0.0);
}

std::array<double,3> FrameTower::to_local (const std::array<Real,3>& v) const
{
    const double x = v[0], y = v[1], z = v[2];
    return {{x * m_along[0] + y * m_along[1], x * m_across[0] + y * m_across[1], z}};
}

std::array<Real,3> FrameTower::to_erf (const std::array<double,3>& v) const
{
    return {{static_cast<Real>(v[0] * m_along[0] + v[1] * m_across[0]), static_cast<Real>(v[0] * m_along[1] + v[1] * m_across[1]),
             static_cast<Real>(v[2])}};
}

void FrameTower::step (Real dt, const std::vector<Real>& node_force, const std::vector<std::array<Real,3>>& line_force)
{
    advance(dt, node_force, line_force, false);
}

void FrameTower::advance (Real dt, const std::vector<Real>& node_force, const std::vector<std::array<Real,3>>& line_force, bool held)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node_force.size() == 3 * m_drag.size(), "FrameTower::step: 3 forces per drag node are needed");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(line_force.size() == m_attach.size(), "FrameTower::step: one line pull per attachment is needed");
    std::vector<double> load(m_frame->num_dofs(), 0.0);
    for (std::size_t i = 0; i < m_drag.size(); ++i) {
        m_drag[i].add_load(to_local({{node_force[3 * i], node_force[3 * i + 1], node_force[3 * i + 2]}}), load);
    }
    for (std::size_t a = 0; a < m_attach.size(); ++a) { m_attach[a].add_load(to_local(line_force[a]), load); }
    for (const double f : load) {
        if (!std::isfinite(f)) { amrex::Abort("tower " + m_name + ": a load on its frame is not finite (the wind's drag or a line's pull)"); }
    }
    if (held) { m_dyn->step_held(static_cast<double>(dt), load, 0.0); }
    else { m_dyn->step(static_cast<double>(dt), load, 0.0); }
    m_load = load;
    m_stepped = true;
}

void FrameTower::step_between (Real dt, const std::vector<Real>& node_force, const std::vector<std::array<Real,3>>& line_start,
                               const std::vector<std::array<Real,3>>& line_end)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(line_start.size() == m_attach.size() && line_end.size() == m_attach.size(),
                                     "FrameTower::step_between: one start and one end pull per attachment are needed");
    std::vector<std::array<Real,3>> mean(line_end.size());
    for (std::size_t a = 0; a < mean.size(); ++a) {
        for (std::size_t d = 0; d < 3; ++d) { mean[a][d] = Real(0.5) * (line_start[a][d] + line_end[a][d]); }
    }
    advance(dt, node_force, mean, true);
}

std::array<Real,3> FrameTower::displacement (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_drag.size(), "FrameTower: no such drag node");
    return to_erf(moved(m_drag[node], m_dyn->displacement()));
}

std::array<Real,3> FrameTower::velocity (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_drag.size(), "FrameTower: no such drag node");
    return to_erf(m_drag[node].motion(m_dyn->velocity()));
}

std::array<Real,3> FrameTower::inertial_force (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_drag.size(), "FrameTower: no such drag node");
    return {{0.0, 0.0, 0.0}};
}

std::array<Real,3> FrameTower::attachment_displacement (std::size_t a) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a < m_attach.size(), "FrameTower: no such attachment");
    return to_erf(moved(m_attach[a], m_dyn->displacement()));
}

std::array<Real,3> FrameTower::attachment_velocity (std::size_t a) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a < m_attach.size(), "FrameTower: no such attachment");
    return to_erf(m_attach[a].motion(m_dyn->velocity()));
}

std::array<double,3> FrameTower::moved (const RigidLink& link, const std::vector<double>& u) const
{
    std::array<double,3> p = link.motion(u);
    if (m_sag != m_sag0) {
        const std::array<double,3> now = link.motion(m_sag), first = link.motion(m_sag0);
        for (std::size_t d = 0; d < 3; ++d) { p[d] += now[d] - first[d]; }
    }
    return p;
}

bool FrameTower::arm_centre_displacement (std::array<Real,3>& x) const
{
    x = to_erf(moved(m_arm.front(), m_dyn->displacement()));
    return true;
}

std::vector<double> FrameTower::state () const
{
    std::vector<double> s = m_dyn->state();
    s.insert(s.end(), m_load.begin(), m_load.end());
    s.insert(s.end(), m_theta.begin(), m_theta.end());
    return s;
}

bool FrameTower::set_state (const std::vector<double>& s)
{
    const std::size_t n = m_frame->num_dofs();
    const std::size_t nm = m_theta.size();
    // a state without the members' temperatures (1 + 4 num_dofs() values) keeps the present ones
    if (s.size() != 1 + 4 * n + nm && s.size() != 1 + 4 * n) { return false; }
    for (const double x : s) { if (!std::isfinite(x)) { return false; } }
    if (s.size() == 1 + 4 * n + nm) {
        const std::vector<double> theta(s.begin() + static_cast<long>(1 + 4 * n), s.end());
        if (theta != m_theta && !set_temperature(theta).empty()) { return false; }
    }
    if (!m_dyn->set_state(std::vector<double>(s.begin(), s.begin() + static_cast<long>(1 + 3 * n)))) { return false; }
    m_load.assign(s.begin() + static_cast<long>(1 + 3 * n), s.begin() + static_cast<long>(1 + 4 * n));
    // a state from before the first step (time 0, the frame at rest) keeps the static footings
    m_stepped = (s[0] > 0.0);
    return true;
}

std::string FrameTower::restart_mismatch (const std::vector<double>& saved) const
{
    const std::size_t n = m_frame->num_dofs();
    if (saved.size() != 1 + 4 * n + m_theta.size()) { return std::string(); }
    for (std::size_t m = 0; m < m_theta.size(); ++m) {
        const double was = saved[1 + 4 * n + m];
        if (was != m_theta[m]) {
            return "its member " + std::to_string(m + 1) + " is at " + std::to_string(was) + " C in the checkpoint and " +
                   std::to_string(m_theta[m]) + " C from the inputs; a restart cannot change the steel's temperature "
                   "(steel_temperature of its type)";
        }
    }
    return std::string();
}

std::vector<double> FrameTower::present_loads (const Tower& tower) const
{
    std::vector<double> load(m_frame->num_dofs(), 0.0);
    const auto& drag = tower.loads();
    for (std::size_t i = 0; i < m_drag.size() && 3 * i + 2 < drag.size(); ++i) {
        m_drag[i].add_load(to_local({{drag[3 * i], drag[3 * i + 1], drag[3 * i + 2]}}), load);
    }
    const auto& pulls = tower.line_forces();
    for (std::size_t a = 0; a < m_attach.size() && a < pulls.size(); ++a) { m_attach[a].add_load(to_local(pulls[a]), load); }
    return load;
}

std::vector<std::array<double,12>> FrameTower::element_forces (const Tower& tower) const
{
    if (!m_stepped) { return m_frame->solve(present_loads(tower), m_gravity).element_force; }
    const auto& u = m_dyn->displacement();
    const auto& v = m_dyn->velocity();
    std::vector<double> w(u.size());
    for (std::size_t i = 0; i < u.size(); ++i) { w[i] = u[i] + m_a1 * v[i]; }
    std::vector<std::array<double,12>> f = m_frame->element_forces(w, 0.0);
    for (std::size_t e = 0; e < f.size(); ++e) {
        for (std::size_t i = 0; i < 12; ++i) { f[e][i] += m_static_force[e][i]; }
    }
    return f;
}

std::vector<MemberCheck> FrameTower::member_checks (const Tower& tower) const
{
    if (!m_designs) { return {}; }
    return check_members(*m_frame, *m_designs, element_forces(tower), m_theta);
}

std::vector<double> FrameTower::node_displacements () const
{
    std::vector<double> u = m_dyn->displacement();
    for (std::size_t i = 0; i < u.size(); ++i) { u[i] += m_sag[i] - m_sag0[i]; }
    return u;
}

std::string FrameTower::set_temperature (const std::vector<double>& theta)
{
    const std::size_t n = m_frame->num_dofs();
    if (theta.size() != m_theta.size()) {
        return "tower " + m_name + ": " + std::to_string(theta.size()) + " temperatures for " + std::to_string(m_theta.size()) + " members";
    }
    FrameInputs in = m_frame->inputs();
    in.temperature = theta;
    std::string err;
    std::shared_ptr<const Frame> frame(Frame::create(in, err));
    if (!frame) { return "tower " + m_name + ": " + m_source + ": " + err; }
    FrameModes modes;
    err = frame_modes(*frame, 1, modes);
    if (!err.empty()) { return "tower " + m_name + ": " + m_source + ": " + err; }
    const FrameSolution s0 = frame->solve({}, m_gravity);
    // the same position and velocity: the motion about the heated frame's static equilibrium makes up the change in sag
    std::vector<double> u = m_dyn->displacement();
    for (std::size_t i = 0; i < n; ++i) { u[i] += m_sag[i] - s0.displacement[i]; }
    auto dyn = std::make_unique<FrameDynamics>(*frame, m_a0, m_a1);
    err = dyn->start(u, m_dyn->velocity(), m_load, 0.0, m_dyn->time());
    if (!err.empty()) { return "tower " + m_name + ": " + m_source + ": " + err; }
    m_dyn = std::move(dyn);
    m_frame = frame;
    m_static = s0.reaction;
    m_static_force = s0.element_force;
    m_sag = s0.displacement;
    m_f1 = modes.frequency[0];
    m_theta = theta;
    return std::string();
}

bool FrameTower::foundation (const Tower& tower, FoundationLoad& L) const
{
    const std::size_t n = m_frame->num_dofs();
    const auto& sup = m_frame->inputs().supports;
    // per support, the reaction beyond the one under the frame's weight (frame axes)
    std::vector<std::array<double,6>> reaction(sup.size());
    if (m_stepped) {
        // the dynamic reactions at the supports: what the members (and concentrated masses) take
        // from each support node, K u + a1 K v + M a, less the load applied there; a spring's force
        // reaches the ground through the support, so the support springs and their masses stay out
        const auto& u = m_dyn->displacement();
        const auto& v = m_dyn->velocity();
        std::vector<double> w(n);
        for (std::size_t i = 0; i < n; ++i) { w[i] = u[i] + m_a1 * v[i]; }
        const std::vector<double> kw = m_frame->apply_stiffness(w, false);
        const std::vector<double> ma = m_frame->apply_mass(m_dyn->acceleration(), false);
        for (std::size_t s = 0; s < sup.size(); ++s) {
            const std::size_t node = m_frame->node_of_joint(sup[s].joint);
            for (std::size_t d = 0; d < 6; ++d) { reaction[s][d] = kw[6 * node + d] + ma[6 * node + d] - m_load[6 * node + d]; }
        }
    } else {
        // at rest before the first step: the static reactions under the tower's present loads
        reaction = m_frame->solve(present_loads(tower), 0.0).reaction;
    }
    // their resultant about the base centre
    std::array<double,3> force{{0.0, 0.0, 0.0}}, moment{{0.0, 0.0, 0.0}};
    std::array<double,4> up{};
    for (std::size_t s = 0; s < sup.size(); ++s) {
        const std::size_t node = m_frame->node_of_joint(sup[s].joint);
        const std::array<double,6>& r = reaction[s];
        const auto& x = m_frame->node_position(node);
        const std::array<double,3> f{{r[0], r[1], r[2]}};
        const std::array<double,3> m = cross(x, f);
        for (std::size_t d = 0; d < 3; ++d) { force[d] += f[d]; moment[d] += m[d] + r[3 + d]; }
        for (std::size_t leg = 0; leg < 4; ++leg) { if (m_leg_support[leg] == s) { up[leg] = r[2] + m_static[s][2]; } }
    }
    // the loads on the tower but its weight balance the dynamic reactions
    const auto F = to_erf({{-force[0], -force[1], -force[2]}});
    const auto M = to_erf({{-moment[0], -moment[1], -moment[2]}});
    L = FoundationLoad();
    L.force = F;
    L.moment = M;
    L.shear = std::hypot(L.force[0], L.force[1]);
    L.overturning = std::hypot(L.moment[0], L.moment[1]);
    double vertical = 0.0;
    for (std::size_t leg = 0; leg < 4; ++leg) { L.legs[leg] = static_cast<Real>(up[leg]); vertical += up[leg]; }
    L.vertical = static_cast<Real>(vertical);
    L.max_compression = std::max({L.legs[0], L.legs[1], L.legs[2], L.legs[3], Real(0.0)});
    L.max_uplift = std::max(Real(0.0), -std::min({L.legs[0], L.legs[1], L.legs[2], L.legs[3]}));
    const TowerType& type = tower.type();
    L.over_allowable = (type.allowable_uplift > 0.0 && L.max_uplift > type.allowable_uplift) ||
                       (type.allowable_compression > 0.0 && L.max_compression > type.allowable_compression);
    return true;
}

} // namespace erf_towers
