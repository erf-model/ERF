// OneModeTower: assembly of the mode and its exact step.

#include "ERF_TowerDynamics.H"

#include <algorithm>
#include <cmath>

#include <AMReX.H>
#include <AMReX_BLassert.H>

using amrex::Real;

namespace erf_towers {

namespace {
constexpr double two_pi = 2.0 * 3.14159265358979323846;
}

OneModeTower::OneModeTower (const Tower& tower, Real gravity)
    : m_name(tower.name()), m_type(tower.type().name)
{
    const TowerType& type = tower.type();
    const std::string key = "erf.conductors." + m_type + ".";
    const std::string who = "OneModeTower " + m_name + ": ";
    if (!type.moves()) { amrex::Abort(who + key + "frequency must be positive (Hz) for a tower that moves"); }
    if (!(type.weight > 0.0)) { amrex::Abort(who + key + "frequency > 0 needs " + key + "weight > 0 (N), which sets the mass"); }
    if (!(gravity > 0.0)) { amrex::Abort(who + "the gravitational acceleration must be positive (m/s^2)"); }
    if (!(type.damping_ratio >= 0.0 && type.damping_ratio < 1.0)) { amrex::Abort(who + key + "damping_ratio must be in [0, 1)"); }
    const auto& nodes = tower.nodes();
    const double H = tower.arm_height();
    if (!(H > 0.0)) { amrex::Abort(who + "the cross-arm height above the base must be positive (m)"); }
    // the mass is spread along the members by length: the shaft's at its drag nodes, and the cross-arm's whole length
    // (the part within the shaft too) at the cross-arm's height, shared by its drag nodes, wherever arm_outside_shaft
    // puts them; so the drag nodes' layout leaves the modal mass, and with it the stiffness, as they are
    const auto nbody = static_cast<std::size_t>(tower.num_body_nodes());
    const double arm_share = static_cast<double>(type.arm_length) / static_cast<double>(TowerType::arm_segments);
    double length = 0.0;
    for (std::size_t i = 0; i < nodes.size(); ++i) { length += (i < nbody) ? static_cast<double>(nodes[i].length) : arm_share; }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(length > 0.0, "OneModeTower: the drag nodes stand for no length of member");
    // the bending shape's generalized mass on a rigid foundation sets the bending stiffness
    const double mass = static_cast<double>(type.weight) / static_cast<double>(gravity);
    const Real arm_z = tower.base()[2] + tower.arm_height();
    double Mb = 0.0;
    std::vector<double> zh(nodes.size()), zm(nodes.size());
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        m_mass.push_back(mass * ((i < nbody) ? static_cast<double>(nodes[i].length) : arm_share) / length);
        zh[i] = (nodes[i].pos[2] - tower.base()[2]) / H;
        zm[i] = (i < nbody) ? zh[i] : (arm_z - tower.base()[2]) / H;
        Mb += m_mass[i] * std::pow(zm[i], 4);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(Mb > 0.0, "OneModeTower: no drag node stands above the base");
    m_Kb = Mb * std::pow(two_pi * static_cast<double>(type.frequency), 2);
    // the compliances in series under a load at the cross-arm; a stiffness of 0 is a rigid foundation. Each
    // division runs on the stiffness floored away from zero and a factor selects the rigid case, rather than
    // a guarded division: an optimising compiler may evaluate a guarded arm speculatively, and H^2 / 0 would
    // then trip amrex.fpe_trap_zero on a rigid footing although the result is discarded. The floor is far below
    // any physical stiffness, so a positive stiffness gives the same compliance as the plain division.
    constexpr double k_floor = 1.0e-30;
    const double kr = static_cast<double>(type.foundation_rotational_stiffness);
    const double kl = static_cast<double>(type.foundation_lateral_stiffness);
    const double cb = 1.0 / m_Kb;
    const double cr = ((kr > 0.0) ? 1.0 : 0.0) * (H * H / std::max(kr, k_floor));
    const double cl = ((kl > 0.0) ? 1.0 : 0.0) * (1.0 / std::max(kl, k_floor));
    const double c = cb + cr + cl;
    auto shape = [&] (double z) { return (cb * z * z + cr * z + cl) / c; };
    for (std::size_t i = 0; i < nodes.size(); ++i) {
        m_phi.push_back(shape(zh[i]));
        m_phi_mass.push_back(shape(zm[i]));
        m_M += m_mass[i] * m_phi_mass[i] * m_phi_mass[i];
    }
    // the lines hang from the tower's attachments, or from the centre of its cross-arm
    for (const auto& at : tower.attachments()) { m_phi_att.push_back(shape((at[2] - tower.base()[2]) / H)); }
    if (m_phi_att.empty()) { m_phi_att.push_back(1.0); }
    m_K = 1.0 / c;
    m_omega = std::sqrt(m_K / m_M);
    m_zeta = type.damping_ratio;
}

std::array<double,2> OneModeTower::generalized_force (const std::vector<Real>& node_force,
                                                      const std::vector<std::array<Real,3>>& line_force) const
{
    // each node's horizontal load and each line's pull by the shape where it acts
    std::array<double,2> Q{{0.0, 0.0}};
    for (int d = 0; d < 2; ++d) {
        for (std::size_t a = 0; a < m_phi_att.size(); ++a) { Q[d] += m_phi_att[a] * line_force[a][d]; }
        for (std::size_t i = 0; i < m_phi.size(); ++i) { Q[d] += m_phi[i] * node_force[3*i+d]; }
    }
    if (!std::isfinite(Q[0]) || !std::isfinite(Q[1])) {
        amrex::Abort("tower " + m_name + ": the generalized load on its mode is not finite (the lines' pull or the members' drag); "
                     "MoorDyn's line integration may have diverged: reduce erf.conductors.moordyn_cfl or erf.conductors.moordyn_dt");
    }
    return Q;
}

void OneModeTower::step_between (Real dt, const std::vector<Real>& node_force, const std::vector<std::array<Real,3>>& line_start,
                                 const std::vector<std::array<Real,3>>& line_end)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(line_start.size() == m_phi_att.size() && line_end.size() == m_phi_att.size(),
                                     "OneModeTower::step_between: one start and one end pull per attachment are needed");
    TowerModel::step_between(dt, node_force, line_start, line_end);
    m_Q_end = generalized_force(node_force, line_end);
}

void OneModeTower::step (Real dt, const std::vector<Real>& node_force, const std::vector<std::array<Real,3>>& line_force)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node_force.size() == 3 * m_phi.size(), "OneModeTower::step: 3 forces per drag node are needed");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(line_force.size() == m_phi_att.size(),
                                     "OneModeTower::step: one line pull per attachment is needed");
    if (!(std::isfinite(dt) && dt > 0.0)) {
        amrex::Abort("OneModeTower " + m_name + ": the step must be finite and positive (s), " + std::to_string(dt) + " given");
    }
    m_Q = generalized_force(node_force, line_force);
    m_Q_end = m_Q;
    // the damped oscillation about the static displacement Q / K, exact for a load held over the step:
    // a = zeta omega (1/s), wd = omega sqrt(1 - zeta^2) (rad/s), x0 the offset from Q/K (m)
    const double h = dt;
    const double a = m_zeta * m_omega;
    const double wd = m_omega * std::sqrt(1.0 - m_zeta * m_zeta);
    const double e = std::exp(-a * h), cs = std::cos(wd * h), sn = std::sin(wd * h);
    for (int d = 0; d < 2; ++d) {
        const double x0 = m_q[d] - m_Q[d] / m_K;
        const double v0 = m_v[d];
        m_q[d] = m_Q[d] / m_K + e * (x0 * cs + (v0 + a * x0) / wd * sn);
        m_v[d] = e * (v0 * cs - (a * v0 + m_omega * m_omega * x0) / wd * sn);
    }
}

std::array<Real,3> OneModeTower::displacement (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_phi.size(), "OneModeTower::displacement: no such drag node");
    return {{static_cast<Real>(m_phi[node] * m_q[0]), static_cast<Real>(m_phi[node] * m_q[1]), Real(0.0)}};
}

std::array<Real,3> OneModeTower::velocity (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_phi.size(), "OneModeTower::velocity: no such drag node");
    return {{static_cast<Real>(m_phi[node] * m_v[0]), static_cast<Real>(m_phi[node] * m_v[1]), Real(0.0)}};
}

std::array<Real,3> OneModeTower::inertial_force (std::size_t node) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(node < m_phi.size(), "OneModeTower::inertial_force: no such drag node");
    const double m = m_mass[node] * m_phi_mass[node];
    return {{static_cast<Real>(-m * acceleration(0)), static_cast<Real>(-m * acceleration(1)), Real(0.0)}};
}

std::array<Real,3> OneModeTower::attachment_displacement (std::size_t a) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a < m_phi_att.size(), "OneModeTower::attachment_displacement: no such attachment");
    return {{static_cast<Real>(m_phi_att[a] * m_q[0]), static_cast<Real>(m_phi_att[a] * m_q[1]), Real(0.0)}};
}

std::array<Real,3> OneModeTower::attachment_velocity (std::size_t a) const
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a < m_phi_att.size(), "OneModeTower::attachment_velocity: no such attachment");
    return {{static_cast<Real>(m_phi_att[a] * m_v[0]), static_cast<Real>(m_phi_att[a] * m_v[1]), Real(0.0)}};
}

bool OneModeTower::arm_centre_displacement (std::array<Real,3>& x) const
{
    x = {{static_cast<Real>(m_q[0]), static_cast<Real>(m_q[1]), Real(0.0)}};
    return true;
}

Real OneModeTower::frequency () const { return static_cast<Real>(m_omega / two_pi); }

std::vector<double> OneModeTower::state () const
{
    return {m_q[0], m_q[1], m_v[0], m_v[1], m_Q[0], m_Q[1], m_Q_end[0], m_Q_end[1]};
}

bool OneModeTower::set_state (const std::vector<double>& s)
{
    // a state without the end load (6 values, an older checkpoint) takes the held one
    if (s.size() != 8 && s.size() != 6) { return false; }
    for (const double x : s) { if (!std::isfinite(x)) { return false; } }
    m_q = {{s[0], s[1]}};
    m_v = {{s[2], s[3]}};
    m_Q = {{s[4], s[5]}};
    m_Q_end = (s.size() == 8) ? std::array<double,2>{{s[6], s[7]}} : m_Q;
    return true;
}

} // namespace erf_towers
