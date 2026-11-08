#ifndef ERF_MEMBER_DRAG_H
#define ERF_MEMBER_DRAG_H

/**
 * Aerodynamic drag on slender structural members (tower shafts and cross-arms) from the wind
 * relative to the member, and the interface through which a tower gets its member loads.
 */

#include <array>
#include <vector>

#include <AMReX_REAL.H>

namespace erf_towers {

/**
 * One drag node of a structural member: where it is, the member's axis there, the length of
 * member it stands for, and the member's drag per unit length, given by the width the wind acts
 * on (a diameter, or a lattice face's width times its solidity) and the drag coefficient on that
 * width.
 */
struct MemberNode
{
    std::array<amrex::Real,3> pos{{0.0, 0.0, 0.0}};    //!< position (m, ERF domain coordinates, z up)
    std::array<amrex::Real,3> axis{{0.0, 0.0, 1.0}};   //!< unit vector along the member (ERF frame)
    amrex::Real length{0.0};             //!< length of member the node stands for (m)
    amrex::Real drag_width{0.0};         //!< width the wind acts on, the area per unit length (m)
    amrex::Real drag_coefficient{1.0};   //!< drag coefficient on drag_width (-)
    std::array<amrex::Real,3> face{{0.0, 0.0, 0.0}};   //!< a square lattice's face normal (horizontal unit vector), for diagonal_gain
    amrex::Real diagonal_gain{0.0};      //!< g of the direction factor 1 + g sin^2(2 theta) (lattice_diagonal_gain); 0: none
};

/**
 * Return the drag on one drag node: 1/2 rho Cd w L |U_n| U_n, where U_n is the wind relative to
 * the node with its component along the member's axis removed (only the normal flow loads a
 * slender member), w the drag width and L the node's length of member. A node of a square lattice
 * shaft (diagonal_gain > 0) takes the factor 1 + g sin^2(2 theta), theta the angle of U_n's
 * horizontal part from the node's face normal: 1 for wind normal to a face, 1 + g along a diagonal.
 *
 * @param node      the drag node
 * @param wind      fluid velocity at the node (m/s, ERF frame)
 * @param velocity  the node's own velocity (m/s, ERF frame)
 * @param rho       air density (kg/m^3)
 * @return the force ON the member BY the air (N, ERF frame)
 */
std::array<amrex::Real,3> member_drag (const MemberNode& node, const std::array<amrex::Real,3>& wind,
                                       const std::array<amrex::Real,3>& velocity, amrex::Real rho);

/**
 * Return the force coefficient (-) of a square lattice tower of flat-sided members on the
 * members' projected area of one face, from the face's solidity phi (the members' area over the
 * face's outline): Cf = 4 phi^2 - 5.9 phi + 4 (ASCE 7, the American Society of Civil Engineers'
 * minimum design loads, for trussed towers of square cross-section). It counts the windward and
 * leeward faces together, so the drag of a lattice shaft is Cf q phi A, with q = rho U^2 / 2 the
 * dynamic pressure (Pa) and A the outline area of one face (m^2).
 */
amrex::Real lattice_force_coefficient (amrex::Real solidity);

/**
 * Return g of a square lattice shaft's wind-direction factor 1 + g sin^2(2 theta) from its faces'
 * solidity phi: g = min(0.75 phi, 0.2). Along a diagonal the factor is ASCE 7's 1 + 0.75 phi, which
 * ASCE 7 caps at 1.2; between a face and a diagonal it varies as IEC 60826's 1 + 0.2 sin^2(2 theta),
 * whose 0.2 is that cap.
 */
amrex::Real lattice_diagonal_gain (amrex::Real solidity);

/**
 * The aerodynamics of a structure's members behind a narrow interface: the drag nodes, the wind
 * at the nodes and the nodes' own velocities go in, the force on each node comes out. The ERF
 * member drag implements it; another aerodynamic model, such as an external tower-drag library,
 * can stand in its place without touching the callers.
 */
class MemberAero
{
public:
    virtual ~MemberAero () = default;
    /**
     * Compute the force on every drag node.
     *
     * @param nodes     the drag nodes
     * @param wind      fluid velocity, 3 values per node, node i at [3i, 3i+3) (0-based, the order of nodes) (m/s)
     * @param velocity  the nodes' own velocities, laid out as wind (m/s)
     * @param force     resized to 3 values per node: the force ON node i BY the air (N, ERF frame)
     */
    virtual void loads (const std::vector<MemberNode>& nodes, const std::vector<amrex::Real>& wind,
                        const std::vector<amrex::Real>& velocity, std::vector<amrex::Real>& force) const = 0;
};

/** ERF's member drag: member_drag() at every drag node with the air density it is given. */
class MemberDrag final : public MemberAero
{
public:
    /** Precondition: rho finite and > 0 (kg/m^3). */
    explicit MemberDrag (amrex::Real rho);
    void loads (const std::vector<MemberNode>& nodes, const std::vector<amrex::Real>& wind,
                const std::vector<amrex::Real>& velocity, std::vector<amrex::Real>& force) const override;
    /** The air density (kg/m^3). */
    amrex::Real air_density () const { return m_rho; }
private:
    amrex::Real m_rho{1.225};
};

} // namespace erf_towers

#endif
