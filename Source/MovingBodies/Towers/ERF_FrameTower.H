#ifndef ERF_FRAME_TOWER_H
#define ERF_FRAME_TOWER_H

/**
 * A lattice tower that bends as its frame model (ERF_Frame.H, read from a SubDyn input file): the
 * TowerModel the conductor coupling drives when a tower type gives erf.conductors.<type>.frame_file.
 * The wind's drag on the tower's drag nodes and the lines' pull at its attachments are carried onto
 * the frame's nodes by rigid links, the frame advances by Newmark's method (ERF_FrameDynamics.H), and
 * the nodes' motion comes back to the drag nodes and attachments through the same links. The frame
 * moves about its static equilibrium under its own weight; its support reactions give the footings'
 * loads and its element forces the members' utilisation (ERF_MemberChecks.H), when the members'
 * design data is given. The steel's temperature softens the frame (EN 1993-1-2): set_temperature()
 * rebuilds and re-factors it.
 *
 * Axes: the frame file's are tower-local, origin at the base centre, x along the line (along =
 * (across_y, -across_x, 0)), y along the cross-arm (Tower::across()), z up. Everything this class takes
 * and returns is in ERF's axes.
 */

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "ERF_Frame.H"
#include "ERF_FrameDynamics.H"
#include "ERF_MemberChecks.H"
#include "ERF_Tower.H"
#include "ERF_TowerDynamics.H"

namespace erf_towers {

/**
 * A point tied rigidly to a few frame nodes. A force at the point is shared among the nodes so that
 * the force and its moment are kept: F_i = F/n + w x d_i, with d_i the node's offset from the nodes'
 * centroid c and w = I^+ ((p - c) x F), I = sum (|d_i|^2 1 - d_i d_i^T) and ^+ the pseudo-inverse. The
 * point's motion is the adjoint, so the links do no work: u_p = mean(u_i) + theta x (p - c), theta =
 * I^+ sum d_i x (u_i - mean(u)). A moment about a line through collinear nodes is not carried. A node
 * within 1e-6 m of the point takes the force alone and gives its motion.
 */
class RigidLink
{
public:
    RigidLink () = default;
    /**
     * Tie point p (m, frame axes) to the count frame nodes nearest to it (all of them when there are fewer),
     * among candidates (0-based node indices; null or empty: every node).
     */
    RigidLink (const Frame& frame, const std::array<double,3>& p, std::size_t count = 4,
               const std::vector<std::size_t>* candidates = nullptr);
    /** The distance (m) from the point to its nearest node. */
    double distance () const { return m_distance; }
    /** The frame nodes the point is tied to (0-based). */
    const std::vector<std::size_t>& nodes () const { return m_nodes; }
    /** Add a force at the point (N, frame axes) to loads, 6 values per frame node (the translations only). */
    void add_load (const std::array<double,3>& force, std::vector<double>& loads) const;
    /** The point's translation (m, frame axes) from u, 6 values per frame node. */
    std::array<double,3> motion (const std::vector<double>& u) const;

private:
    std::vector<std::size_t> m_nodes;
    std::vector<std::array<double,3>> m_d;      //!< each node's offset from the centroid (m)
    std::array<double,3> m_r{{0.0, 0.0, 0.0}};  //!< the point's offset from the centroid (m)
    std::array<double,9> m_iinv{};              //!< I^+ (1/m^2), row-major
    double m_distance{0.0};
};

/** One tower bending as its frame model. */
class FrameTower final : public TowerModel
{
public:
    /**
     * Stand the frame on the tower: tie the tower's drag nodes and attachments to the frame's nodes
     * (each must lie within the type's base_width of a frame node), match the frame's four supports to
     * the tower's legs (one per quadrant of the base, all at z = 0), solve the static equilibrium under
     * the frame's weight (gravity in m/s^2), find the first natural frequency and set the Rayleigh
     * damping of the type's damping_ratio at that frequency and ten times it. designs, when not null,
     * holds the design data of every frame member in the frame's member order (match_designs()), for
     * member_checks(). link_nodes are the frame nodes (0-based) the drag nodes and attachments may be tied
     * to (empty: every node). first_frequency (Hz) is the frame's first natural frequency when its
     * caller has found it for the towers that share the frame (0: found here). Aborts, naming the tower
     * and its frame (source, for messages), when the frame does not fit the tower: a drag node or an
     * attachment far from the frame, a frame with other than one interface joint, or its interface joint
     * (its cross-arm's centre) more than half the cross-arm's face depth from the tower's cross-arm height
     * (a frame built for another height).
     */
    FrameTower (const Tower& tower, std::shared_ptr<const Frame> frame, double gravity,
                std::shared_ptr<const std::vector<MemberDesign>> designs = nullptr, const std::string& source = std::string(),
                const std::vector<std::size_t>& link_nodes = {}, double first_frequency = 0.0);

    using TowerModel::step;
    std::size_t num_attachments () const override { return m_attach.size(); }
    /** Advance by dt to the loads given, those at the end of the step (Newmark's method). */
    void step (amrex::Real dt, const std::vector<amrex::Real>& node_force,
               const std::vector<std::array<amrex::Real,3>>& line_force) override;
    /**
     * The mean of the pulls at the start and the end held over the step with the drag, by the
     * implicit midpoint rule (FrameDynamics::step_held): the work the frame takes from the pulls is
     * their mean times the attachments' displacement, the work MoorDyn's lines give up when they
     * move the attachments at a constant velocity and the mean is their mean pull over the step.
     */
    void step_between (amrex::Real dt, const std::vector<amrex::Real>& node_force,
                       const std::vector<std::array<amrex::Real,3>>& line_start,
                       const std::vector<std::array<amrex::Real,3>>& line_end) override;
    /** A drag node's displacement from where it stands under the frame's weight (m, ERF axes). */
    std::array<amrex::Real,3> displacement (std::size_t node) const override;
    std::array<amrex::Real,3> velocity (std::size_t node) const override;
    /** Zero: the footings' loads come from the frame's support reactions (foundation()), not from the drag nodes' inertia. */
    std::array<amrex::Real,3> inertial_force (std::size_t node) const override;
    std::array<amrex::Real,3> attachment_displacement (std::size_t a = 0) const override;
    std::array<amrex::Real,3> attachment_velocity (std::size_t a = 0) const override;
    /** The frame's first natural frequency (Hz), which the coupling resolves; Newmark's method integrates the higher modes stably. */
    amrex::Real frequency () const override { return static_cast<amrex::Real>(m_f1); }
    /**
     * The Newmark state (time, u, v, a of every frame degree of freedom), the frame loads of the last step and the
     * members' temperatures (C).
     */
    std::vector<double> state () const override;
    /** The motion of the frame where the cross-arm's centre stands (the tower's axis at its arm height). */
    bool arm_centre_displacement (std::array<amrex::Real,3>& x) const override;
    /** Members heated differently from the inputs' (the steel_temperature of the type) cannot continue a restart. */
    std::string restart_mismatch (const std::vector<double>& saved) const override;
    /**
     * Restore a state from state(), reheating the frame first when its members' temperatures differ from the
     * present ones (heating within a run; a restart refuses different temperatures through restart_mismatch()).
     */
    bool set_state (const std::vector<double>& s) override;
    /**
     * The footings' load from the frame's support reactions: the static reactions under its weight plus
     * the dynamic ones, K u + a1 K v + M a - f at the supports (the mass-proportional damping excluded).
     * Before the frame's first step, the static reactions under the tower's present drag (Tower::loads())
     * and line pulls (Tower::line_forces()), as Tower::foundation() reckons a tower at rest.
     * force and moment (about the base centre) are the loads on the tower but its weight, as
     * Tower::foundation() defines them; legs[i] is the upward reaction at the support of leg i
     * (compression positive), in Tower::foundation()'s leg order.
     */
    bool foundation (const Tower& tower, FoundationLoad& load) const override;

    const Frame& frame () const { return *m_frame; }
    /** Where the frame came from (its key and file, or how it was generated), for messages. */
    const std::string& source () const { return m_source; }
    double mass () const { return m_frame->total_mass(); }

    /** Whether the members' design data was given, so member_checks() has something to check. */
    bool has_member_checks () const { return static_cast<bool>(m_designs); }
    /**
     * Every member's element forces (FrameSolution::element_force: per element, the forces on it at its
     * nodes in its local axes, N and N m): the static ones under the frame's weight plus K_e (u + a1 v), the
     * elastic and the stiffness-proportional damping forces of the motion. Before the frame's first step, the
     * static forces under the frame's weight and the tower's present drag and line pulls, as foundation().
     */
    std::vector<std::array<double,12>> element_forces (const Tower& tower) const;
    /** The members' strength checks under element_forces() at the members' temperatures; empty without design data. */
    std::vector<MemberCheck> member_checks (const Tower& tower) const;
    /**
     * Every frame node's displacement from where it stands under the frame's weight at its first temperature
     * (6 per node, m and rad; frame axes).
     */
    std::vector<double> node_displacements () const;

    /** The members' temperatures (C, one per member in the frame's order). */
    const std::vector<double>& temperature () const { return m_theta; }
    /**
     * Heat (or cool) the members to theta (C, one per member, each below 1200): rebuild the frame with the steel's
     * stiffness at those temperatures (EN 1993-1-2), re-factor it, and carry on from the same position and
     * velocity: the static equilibrium under the weight, the reactions under it, the first natural frequency
     * and the Newmark step's factorisation follow the heated frame's stiffness; the Rayleigh coefficients stay those
     * of the frame at its first temperature. Returns an empty string, or the problem (and changes nothing).
     */
    std::string set_temperature (const std::vector<double>& theta);

private:
    std::array<double,3> to_local (const std::array<amrex::Real,3>& v) const;
    std::array<amrex::Real,3> to_erf (const std::array<double,3>& v) const;
    /** The frame loads of the tower's present drag and line pulls (6 per node, frame axes). */
    std::vector<double> present_loads (const Tower& tower) const;
    /** Advance by dt under the drag and line pulls given: at the end of the step (Newmark), or held over it (held, the midpoint rule). */
    void advance (amrex::Real dt, const std::vector<amrex::Real>& node_force,
                  const std::vector<std::array<amrex::Real,3>>& line_force, bool held);
    /** A point's motion through its link: of the dynamic displacement u, plus the weight's extra sag since the first temperature. */
    std::array<double,3> moved (const RigidLink& link, const std::vector<double>& u) const;

    std::shared_ptr<const Frame> m_frame;
    std::string m_name;                          //!< the tower's name, for messages
    std::array<double,3> m_along{}, m_across{};  //!< the tower-local x and y axes in ERF's axes
    std::vector<RigidLink> m_drag;               //!< per drag node of the tower
    std::vector<RigidLink> m_attach;             //!< per attachment of the tower
    std::vector<RigidLink> m_arm;                //!< the cross-arm's centre (one link)
    std::array<std::size_t,4> m_leg_support{};   //!< per leg (Tower::foundation()'s order): the index into the frame's supports
    std::vector<std::array<double,6>> m_static;  //!< per support: the reaction under the frame's weight (N, N m; frame axes)
    std::vector<std::array<double,12>> m_static_force; //!< per element: its end forces under the frame's weight (local axes)
    std::vector<double> m_sag;                   //!< the displacement under the frame's weight (6 per node, frame axes)
    std::vector<double> m_sag0;                  //!< the same at the first temperature, which the tower's motion is measured from
    std::shared_ptr<const std::vector<MemberDesign>> m_designs;
    std::vector<double> m_theta;                 //!< the members' temperatures (C)
    std::string m_source;                        //!< where the frame came from, for messages
    double m_gravity{0.0};                       //!< m/s^2
    double m_f1{0.0};                            //!< first natural frequency (Hz)
    double m_a0{0.0};                            //!< the mass-proportional Rayleigh coefficient (1/s)
    double m_a1{0.0};                            //!< the stiffness-proportional Rayleigh coefficient (s)
    std::unique_ptr<FrameDynamics> m_dyn;
    std::vector<double> m_load;                  //!< the frame loads of the last step (6 per node, frame axes)
    bool m_stepped{false};                       //!< whether the frame has taken a step (or continues a saved state)
};

} // namespace erf_towers

#endif
