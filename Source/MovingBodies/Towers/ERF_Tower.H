#ifndef ERF_TOWER_H
#define ERF_TOWER_H

/**
 * A lattice tower standing on the terrain: its drag nodes, the loads on them, the lines' pull at
 * its attachments, and the reactions at its four footings.
 */

#include <array>
#include <string>
#include <vector>

#include <AMReX_REAL.H>

#include "ERF_MemberDrag.H"
#include "ERF_TowerInputs.H"

namespace erf_towers {

/**
 * One tower's foundation load: the resultant on the tower about the centre of its base, and the
 * axial force in each of its four legs, a rigid square of footings sharing the vertical load and
 * resisting the overturning moment (compression positive). Frame: ERF, z up.
 */
struct FoundationLoad
{
    /**
     * the members' drag, the lines' pull and, for a moving tower, the nodes' inertial forces -m a (d'Alembert);
     * the weight excluded (N)
     */
    std::array<amrex::Real,3> force{{0.0, 0.0, 0.0}};
    std::array<amrex::Real,3> moment{{0.0, 0.0, 0.0}};  //!< the moment of those forces about the base centre (N m)
    amrex::Real shear{0.0};                             //!< the magnitude of the horizontal force, the base shear (N)
    amrex::Real overturning{0.0};                       //!< the magnitude of the horizontal moment (N m)
    amrex::Real vertical{0.0};                          //!< the downward load on the footings, the weight included (N)
    /** each leg's axial force, compression positive, in the order given at Tower::foundation() (N) */
    std::array<amrex::Real,4> legs{{0.0, 0.0, 0.0, 0.0}};
    amrex::Real max_compression{0.0};                   //!< the largest leg compression (N)
    amrex::Real max_uplift{0.0};                        //!< the largest leg tension, 0 when no leg is in tension (N)
    bool over_allowable{false};                         //!< the largest uplift or compression exceeds its allowable (where one is given)
};

/**
 * One lattice tower of a TowerType standing on the terrain: a vertical shaft (the tower's "body")
 * from its base up to the cross-arm, and the peak above it, tapering linearly from base_width to
 * top_width; and a horizontal cross-arm of arm_length centred on the shaft at arm_height above
 * the base, along the horizontal unit vector `across`. Each is cut into drag nodes of equal
 * segments: a shaft node stands at the middle of its segment, with the shaft's width there times
 * the solidity as its drag width; the cross-arm's face is arm_face() deep out to the shaft's
 * corners at the arm's bottom and tapers linearly to nothing at the tips, as the generated frame's
 * truss of chords converging on its tip (arm_face_at()), so a cross-arm node stands at the centroid of its segment's face, with
 * the face's mean depth over the segment times the solidity as its drag width. The shaft's nodes
 * come first (the peak's included), then the cross-arm's. Frame: ERF domain coordinates, z up.
 */
class Tower
{
public:
    /**
     * Build the tower and its drag nodes.
     *
     * @param name        the tower's name, used in messages
     * @param type        its type (validated: TowerType::validate() is empty)
     * @param base        centre of the base on the terrain (m, ERF frame)
     * @param arm_height  height of the cross-arm above the base (m, > 0)
     * @param across      horizontal unit vector along the cross-arm (ERF frame)
     */
    Tower (std::string name, const TowerType& type, const std::array<amrex::Real,3>& base,
           amrex::Real arm_height, const std::array<amrex::Real,3>& across);

    const std::string& name () const { return m_name; }
    const TowerType& type () const { return m_type; }
    /** The centre of the base, on the terrain (m, ERF frame). */
    const std::array<amrex::Real,3>& base () const { return m_base; }
    /** The height of the cross-arm above the base (m). */
    amrex::Real arm_height () const { return m_arm_height; }
    /**
     * The depth of the cross-arm's face (m) at the distance s (m) from the shaft's axis along the cross-arm:
     * the type's arm_face() out to arm_taper_start(), then tapering linearly to 0 at the tips (|s| = half
     * the arm length), the shape of the generated frame's cross-arm (ERF_LatticeFrame.H), whose four chords
     * run from the shaft's corners at the cross-arm's top and bottom to the tip, so the side view keeps the
     * full depth out to the bottom corners and tapers from there. 0 at and beyond the tips, whatever the
     * layout. The earlier layout (arm_outside_shaft = false) does not use it: its nodes carry arm_face()
     * along the whole length.
     */
    amrex::Real arm_face_at (amrex::Real s) const;
    /**
     * Where the cross-arm's face starts to taper (m from the shaft's axis): half the shaft's width at the
     * cross-arm's bottom, arm_height - arm_face() (at the base when the face is deeper than the arm is high).
     */
    amrex::Real arm_taper_start () const;
    /** The horizontal unit vector along the cross-arm (ERF frame). */
    const std::array<amrex::Real,3>& across () const { return m_across; }
    /** The drag nodes where the tower stands unloaded: the shaft's (the peak's included), then the cross-arm's. */
    const std::vector<MemberNode>& nodes () const { return m_nodes; }
    /**
     * The number of drag nodes of the shaft, the peak's included; the cross-arm's nodes have the
     * 0-based indices [num_body_nodes(), nodes().size()).
     */
    int num_body_nodes () const { return m_nbody; }

    /**
     * Hand over how the tower moves: every drag node's displacement from where it stands unloaded
     * (m), its velocity (m/s) and its inertial force, minus its mass times its acceleration (N).
     * Each vector holds 3 values per node, node i at [3i, 3i+3) (0-based, the order of nodes()).
     * Frame: ERF. A tower that stands still keeps them zero.
     */
    void set_motion (const std::vector<amrex::Real>& displacement, const std::vector<amrex::Real>& velocity,
                     const std::vector<amrex::Real>& inertia);
    /** The drag nodes at their current positions: each unloaded position plus its displacement (m, ERF frame). */
    std::vector<MemberNode> current_nodes () const;
    /** The drag nodes' velocities (3 per node, m/s, ERF frame), zero for a tower that stands still. */
    const std::vector<amrex::Real>& node_velocities () const { return m_vel; }
    /**
     * The displacement of the cross-arm's first drag node (m, ERF frame), which stands at the centroid of the
     * tip piece of the cross-arm's face on one side (arm_outside_shaft) or at the cross-arm's height; the tower models give the cross-arm
     * centre's displacement at its height (TowerModel::arm_centre_displacement()), which the logs report. A line
     * hanging from an attachment moves by TowerModel::attachment_displacement().
     */
    std::array<amrex::Real,3> arm_displacement () const;

    /** Set the force ON every drag node BY the air (3 per node, N, ERF frame): the members' drag. */
    void set_loads (const std::vector<amrex::Real>& force);
    /** The force on every drag node from set_loads() (3 per node, N, ERF frame). */
    const std::vector<amrex::Real>& loads () const { return m_force; }
    /** The sum of the members' drag (N, ERF frame). */
    std::array<amrex::Real,3> total_force () const;
    /**
     * The moment of the members' drag about the centre of the base (N m, ERF frame), with lever
     * arms from the unloaded drag-node positions (the sway's P-delta is neglected).
     */
    std::array<amrex::Real,3> base_moment () const;

    /**
     * Add an attachment, the point one line hangs from, where it stands on the unloaded tower
     * (m, ERF frame), and return its 0-based index: a circuit's phases and shield wire each hang
     * from their own attachment of one tower. OneModeTower treats a tower with no attachments as
     * having one attachment at the centre of the cross-arm (shape 1); Tower itself applies each
     * pull where set_line_loads() places it.
     */
    std::size_t add_attachment (const std::array<amrex::Real,3>& at);
    /** The attachments, where they stand on the unloaded tower (m, ERF frame), in the order of add_attachment(). */
    const std::vector<std::array<amrex::Real,3>>& attachments () const { return m_attach; }
    /**
     * Set the force ON the tower BY each line (N, ERF frame) and the point where that line hangs
     * at present (m, ERF frame), one per attachment in the order of add_attachment() (at most one
     * for a tower without attachments). Aborts, naming the tower, when a force or point is not finite.
     */
    void set_line_loads (const std::vector<std::array<amrex::Real,3>>& force, const std::vector<std::array<amrex::Real,3>>& at);
    /** set_line_loads() for a tower with a single line. */
    void set_line_load (const std::array<amrex::Real,3>& force, const std::array<amrex::Real,3>& at)
    {
        set_line_loads({force}, {at});
    }
    /** The sum of the lines' pull on the tower (N, ERF frame). */
    std::array<amrex::Real,3> line_force () const;
    /** The pull of each line on the tower, one per attachment (N, ERF frame). */
    const std::vector<std::array<amrex::Real,3>>& line_forces () const { return m_line_force; }
    /**
     * Return the foundation's load from the members' drag and the lines' pull together, with the
     * tower's weight and, for a tower that moves, the inertial forces of its drag nodes (the
     * reactions balance the loads less the mass times the acceleration of every node).
     *
     * The legs stand at +-s/2 along the line and +-s/2 along the cross-arm (s the leg spacing)
     * and share the downward load P equally and the overturning moment linearly:
     * N_i = P/4 + (M_y x_i - M_x y_i) / (4 (s/2)^2), with x_i, y_i the leg's offset from the base
     * centre and M_x, M_y the moment about the base centre, all in the ERF frame (the formula is
     * the same in any horizontal frame because the four legs form a square). Legs in order:
     * (+along,+across), (+along,-across), (-along,+across), (-along,-across), with
     * along = (across_y, -across_x).
     */
    FoundationLoad foundation () const;

private:
    std::string m_name;
    TowerType m_type;
    std::array<amrex::Real,3> m_base;                      //!< m
    amrex::Real m_arm_height{0.0};                         //!< m
    std::array<amrex::Real,3> m_across;
    std::vector<MemberNode> m_nodes;
    int m_nbody{0};
    std::vector<amrex::Real> m_force;     //!< 3 per drag node: the air's force on it (N)
    std::vector<amrex::Real> m_disp;      //!< 3 per drag node: the displacement from where it stands unloaded (m)
    std::vector<amrex::Real> m_vel;       //!< 3 per drag node (m/s)
    std::vector<amrex::Real> m_inertia;   //!< 3 per drag node: minus mass times acceleration (N)
    std::vector<std::array<amrex::Real,3>> m_attach;       //!< per attachment: where it stands on the unloaded tower (m)
    std::vector<std::array<amrex::Real,3>> m_line_force;   //!< per attachment: the line's pull ON the tower (N)
    std::vector<std::array<amrex::Real,3>> m_line_at;      //!< per attachment: where the line pulls at present (m)
};

} // namespace erf_towers

#endif
