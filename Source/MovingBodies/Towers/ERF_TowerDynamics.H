#ifndef ERF_TOWER_DYNAMICS_H
#define ERF_TOWER_DYNAMICS_H

/** The structural model of a tower (the TowerModel interface) and its one-mode implementation, OneModeTower. */

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include <AMReX.H>
#include <AMReX_REAL.H>

#include "ERF_Tower.H"

namespace erf_towers {

/**
 * One tower's structural model behind a narrow interface: the loads on its drag nodes (the
 * members' drag, 3 per node) and each line's pull at its attachment (one per attachment) go in,
 * held over a step; the motion of every drag node and of every attachment comes out, with each
 * node's inertial force for the foundation's reactions. The one-mode tower implements it; a frame
 * model of the lattice can stand in its place without touching the callers.
 *
 * Indices: drag nodes 0-based in the order of Tower::nodes(); attachments 0-based in the order
 * of Tower::add_attachment(). Frame: ERF, z up.
 */
class TowerModel
{
public:
    virtual ~TowerModel () = default;
    /** The number of attachments the lines hang from. */
    virtual std::size_t num_attachments () const = 0;
    /**
     * Advance by dt under loads held over the step.
     *
     * @param dt          the step (s, > 0)
     * @param node_force  3 values per drag node: the force ON the tower (N, ERF frame)
     * @param line_force  one per attachment: the pull ON the tower BY the line there (N, ERF frame)
     *
     * Only the horizontal components drive OneModeTower.
     */
    virtual void step (amrex::Real dt, const std::vector<amrex::Real>& node_force,
                       const std::vector<std::array<amrex::Real,3>>& line_force) = 0;
    /** step() for a tower with a single attachment. */
    void step (amrex::Real dt, const std::vector<amrex::Real>& node_force, const std::array<amrex::Real,3>& line_force)
    {
        step(dt, node_force, std::vector<std::array<amrex::Real,3>>{line_force});
    }
    /**
     * Advance by dt with the lines' pulls at the start and at the end of the step (the drag
     * held over it). The tower takes the mean of the two held over the step: OneModeTower
     * integrates a held load exactly, FrameTower by the implicit midpoint rule. Either way the
     * work the tower takes from the pulls is their mean times the attachments' displacement,
     * which the coupling with MoorDyn makes the work the lines give up, by handing the lines'
     * mean pull over the step as both the start and the end pull.
     */
    virtual void step_between (amrex::Real dt, const std::vector<amrex::Real>& node_force,
                               const std::vector<std::array<amrex::Real,3>>& line_start,
                               const std::vector<std::array<amrex::Real,3>>& line_end)
    {
        std::vector<std::array<amrex::Real,3>> mean(line_end.size());
        for (std::size_t a = 0; a < mean.size(); ++a) {
            for (std::size_t d = 0; d < 3; ++d) { mean[a][d] = amrex::Real(0.5) * (line_start[a][d] + line_end[a][d]); }
        }
        step(dt, node_force, mean);
    }
    /** A drag node's displacement from where it stands unloaded (m, ERF frame; horizontal only, z = 0, in OneModeTower). */
    virtual std::array<amrex::Real,3> displacement (std::size_t node) const = 0;
    /** A drag node's velocity (m/s, ERF frame; horizontal only in OneModeTower). */
    virtual std::array<amrex::Real,3> velocity (std::size_t node) const = 0;
    /** A drag node's inertial force, minus its mass times its acceleration (N, ERF frame), under the loads of the last step. */
    virtual std::array<amrex::Real,3> inertial_force (std::size_t node) const = 0;
    /** The displacement of attachment a (m, ERF frame), the first by default. */
    virtual std::array<amrex::Real,3> attachment_displacement (std::size_t a = 0) const = 0;
    /** The velocity of attachment a (m/s, ERF frame), the first by default. */
    virtual std::array<amrex::Real,3> attachment_velocity (std::size_t a = 0) const = 0;
    /**
     * The highest natural frequency (Hz, undamped) that the coupling must resolve: the caller
     * takes coupling_steps_per_period coupling steps per period of the fastest tower in a group.
     * For OneModeTower, the frequency of its one mode on its compliant foundation, sqrt(K/M) / (2 pi).
     */
    virtual amrex::Real frequency () const = 0;
    /** The state to checkpoint. OneModeTower: {q_x, q_y, v_x, v_y, Q_x, Q_y, Q_end_x, Q_end_y} (m, m/s, N). */
    virtual std::vector<double> state () const = 0;
    /**
     * The displacement of the cross-arm's centre at the cross-arm's height (m, ERF frame) into x: false when the
     * model leaves it to Tower::arm_displacement(). OneModeTower: its modal displacement, the mode shape being 1
     * at the cross-arm; FrameTower: the frame's cross-arm centre.
     */
    virtual bool arm_centre_displacement (std::array<amrex::Real,3>& /*x*/) const { return false; }
    /**
     * Why a checkpointed state cannot continue this run although it fits the model (state() of the run
     * that wrote it): an empty string when it can. A frame tower refuses members heated differently from
     * the inputs'; a one-mode tower takes any state that fits.
     */
    virtual std::string restart_mismatch (const std::vector<double>& /*saved*/) const { return std::string(); }
    /** Restore a state from state(): returns false, changing nothing, when the values do not fit the model or are not finite. */
    virtual bool set_state (const std::vector<double>& s) = 0;
    /**
     * The load on the tower's footings, when the model knows its support reactions: fills load and
     * returns true. Otherwise returns false and Tower::foundation() estimates it from the drag nodes'
     * loads and inertial forces (d'Alembert). Frame: ERF, z up.
     */
    virtual bool foundation (const Tower& tower, FoundationLoad& load) const
    {
        amrex::ignore_unused(tower, load);
        return false;
    }
};

/**
 * One tower bending in its first mode, the same way in x and y (a square tower), on a foundation
 * that can tilt and slide. Its mass, the type's weight over gravity, is spread by length of member:
 * the shaft's over its drag nodes, and the whole cross-arm's (arm_length, the part within the shaft
 * too) at the cross-arm's height, shared by the cross-arm's drag nodes wherever they stand, so the
 * drag nodes' layout (arm_outside_shaft) leaves the mode's mass as it is. A drag node at height z above the base moves
 * horizontally by phi(z) times the cross-arm's displacement q, where
 *
 *   phi(z) = a_b (z/H)^2 + a_r z/H + a_l,   H the cross-arm's height,
 *
 * combines an assumed quadratic bending shape (zero slope at the footing; the static shape of a
 * uniform cantilever under a tip load is (z/H)^2 (3 - z/H)/2, so K_b is defined from the
 * frequency, not from EI), the tilt of the footing and its slide, each in proportion to its
 * compliance under a load at the cross-arm: a_b = (1/K_b)/c, a_r = (H^2/k_r)/c, a_l = (1/k_l)/c,
 * c = 1/K_b + H^2/k_r + 1/k_l, a rigid foundation (stiffness 0) adding nothing. The generalized
 * mass is M = sum m_i phi_i^2 and the stiffness K = 1/c, the strain energy of that shape exactly;
 * K_b is the bending stiffness that gives the type's frequency on a rigid foundation,
 * K_b = sum m_i (z_i/H)^4 (2 pi f)^2. In each horizontal direction
 *
 *   M q'' + 2 zeta omega M q' + K q = Q,   Q = sum_i phi_i F_i + sum_a phi_a F_a,   omega^2 = K / M,
 *
 * with F_i the drag on drag node i and F_a the pull of the line at attachment a (horizontal
 * components, N), and phi_a the shape at the attachment's height. The load Q is held over a step
 * and integrated exactly, so any step is stable. Each line pulls at its attachment (the tower's
 * attachments, or the centre of the cross-arm when it has none), which moves by phi_a q:
 * phi = 1 at the cross-arm, more on a peak above it. The tower does not twist, rise or sink.
 *
 * Symbols: H (m), m_i the drag nodes' masses (kg), K_b and k_l (N/m), k_r (N m/rad), f the type's
 * frequency (Hz), zeta its damping_ratio (-), q the cross-arm's horizontal displacement (m).
 */
class OneModeTower final : public TowerModel
{
public:
    /**
     * Build the mode of a tower whose type moves (frequency > 0, weight > 0), with every
     * attachment already added to the tower. Aborts naming the tower and its type's keys otherwise.
     *
     * @param tower    the tower
     * @param gravity  the gravitational acceleration (m/s^2, > 0)
     */
    OneModeTower (const Tower& tower, amrex::Real gravity);

    using TowerModel::step;
    std::size_t num_attachments () const override { return m_phi_att.size(); }
    /** Aborts naming the tower when dt is not positive or the generalized load is not finite. */
    void step (amrex::Real dt, const std::vector<amrex::Real>& node_force,
               const std::vector<std::array<amrex::Real,3>>& line_force) override;
    /**
     * The mean pull held over the step, integrated exactly; the inertia the footings take at the
     * step's end (inertial_force()) is that under line_end. The coupling with MoorDyn hands the
     * lines' mean pull over the step as line_end too, so the footings carry that mean.
     */
    void step_between (amrex::Real dt, const std::vector<amrex::Real>& node_force,
                       const std::vector<std::array<amrex::Real,3>>& line_start,
                       const std::vector<std::array<amrex::Real,3>>& line_end) override;
    std::array<amrex::Real,3> displacement (std::size_t node) const override;
    std::array<amrex::Real,3> velocity (std::size_t node) const override;
    std::array<amrex::Real,3> inertial_force (std::size_t node) const override;
    std::array<amrex::Real,3> attachment_displacement (std::size_t a = 0) const override;
    std::array<amrex::Real,3> attachment_velocity (std::size_t a = 0) const override;
    bool arm_centre_displacement (std::array<amrex::Real,3>& x) const override;
    amrex::Real frequency () const override;
    std::vector<double> state () const override;
    bool set_state (const std::vector<double>& s) override;

    /** The drag node's mass (kg). */
    double node_mass (std::size_t node) const { return m_mass[node]; }
    /** The drag node's value of the mode shape (-). */
    double mode_shape (std::size_t node) const { return m_phi[node]; }
    /** The mode shape at attachment a (-). */
    double attachment_shape (std::size_t a) const { return m_phi_att[a]; }
    /** The generalized mass of the mode (kg), the cross-arm's displacement its coordinate. */
    double generalized_mass () const { return m_M; }
    /** The generalized stiffness of the mode (N/m). */
    double stiffness () const { return m_K; }
    /** The stiffness of the bending alone, on a rigid foundation (N/m). */
    double bending_stiffness () const { return m_Kb; }
    /** The cross-arm's horizontal displacement (m), x and y. */
    const std::array<double,2>& q () const { return m_q; }
    /** The cross-arm's horizontal velocity (m/s), x and y. */
    const std::array<double,2>& v () const { return m_v; }

private:
    /** The generalized force of the drag on the nodes and the lines' pull at the attachments (N); aborts when not finite. */
    std::array<double,2> generalized_force (const std::vector<amrex::Real>& node_force,
                                            const std::vector<std::array<amrex::Real,3>>& line_force) const;
    /** The acceleration of the mode at the end of the last step, under the load then (m/s^2). */
    double acceleration (int d) const { return (m_Q_end[d] - 2.0 * m_zeta * m_omega * m_M * m_v[d] - m_K * m_q[d]) / m_M; }

    std::string m_name;                    //!< the tower's name, for messages
    std::string m_type;                    //!< its type's name, for messages
    std::vector<double> m_mass;            //!< per drag node (kg)
    std::vector<double> m_phi;             //!< the mode shape at each drag node (-)
    std::vector<double> m_phi_mass;        //!< the mode shape where each drag node's mass is (the cross-arm's at its height) (-)
    std::vector<double> m_phi_att;         //!< the mode shape at each attachment (-)
    double m_M{0.0};                       //!< generalized mass (kg)
    double m_K{0.0};                       //!< generalized stiffness (N/m)
    double m_Kb{0.0};                      //!< bending stiffness on a rigid foundation (N/m)
    double m_omega{0.0};                   //!< natural angular frequency (rad/s)
    double m_zeta{0.0};                    //!< damping ratio (-)
    std::array<double,2> m_q{{0.0, 0.0}};  //!< cross-arm displacement (m)
    std::array<double,2> m_v{{0.0, 0.0}};  //!< cross-arm velocity (m/s)
    std::array<double,2> m_Q{{0.0, 0.0}};  //!< generalized force held over the last step (N)
    std::array<double,2> m_Q_end{{0.0, 0.0}};  //!< generalized force at the end of the last step (N)
};

} // namespace erf_towers

#endif
