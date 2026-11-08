#ifndef ERF_FRAME_H
#define ERF_FRAME_H

/**
 * A linear 3D frame model of a lattice tower: two-node Euler-Bernoulli or Timoshenko beam
 * elements with six degrees of freedom per node, assembled in double precision, fixed or
 * spring-supported at its base joints, and solved for the static response to nodal loads and
 * gravity. The element formulation, the local axes and the gravity loads follow OpenFAST's SubDyn
 * (FEM.f90: GetDirCos, ElemK_Beam, ElemG; SD_FEM.f90: SetElementProperties), so that a SubDyn
 * input file gives the same stiffness here and there.
 *
 * Conventions used throughout:
 *   Frame: the coordinates of the input file (m), z up; gravity acts along -z.
 *   Degrees of freedom: 6 per node, node n at [6n, 6n+6): x, y, z translations (m), then rotations
 *     about x, y, z (rad). Loads in the same order: forces (N), then moments (N m).
 *   Nodes: the joints first, in input order, then each member's interior nodes when it is divided
 *     into more than one element, member by member, from its joint_a towards its joint_b.
 *   Local axes of an element: z from its first node to its second; x and y as SubDyn's GetDirCos
 *     builds them, then turned about z by the member's spin.
 */

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "ERF_DenseCholesky.H"
#include "ERF_FrameInputs.H"

namespace erf_towers {

/** The properties of one beam element, derived from its section as SubDyn derives them. Units: SI. */
struct BeamProperties
{
    double E{0.0};            //!< Young's modulus (Pa)
    double G{0.0};            //!< shear modulus (Pa)
    double rho{0.0};          //!< density (kg/m^3)
    double A{0.0};            //!< area (m^2)
    double Ixx{0.0};          //!< second moment about the local x axis (m^4): bending in the local y-z plane
    double Iyy{0.0};          //!< second moment about the local y axis (m^4): bending in the local x-z plane
    double J0{0.0};           //!< polar moment (m^4), for the rotary inertia
    double Jt{0.0};           //!< torsion constant (m^4), for the torsional stiffness
    double kappa_x{0.0};      //!< shear coefficient for shear along local x (-): shear area kappa_x A
    double kappa_y{0.0};      //!< shear coefficient for shear along local y (-)
    bool shear{false};        //!< Timoshenko (true) or Euler-Bernoulli (false)
};

/**
 * A section's beam properties under a beam theory, as SubDyn's SetElementProperties computes them:
 *   Circular:    A, Ixx = Iyy, J0 = Jt = 2 Ixx from D and t; Timoshenko kappa from Steinboeck et al.
 *                with the ratio of inner to outer diameter, D_inner = D - 2t (so 1 for a solid section, t = 0, as SubDyn takes it).
 *   Rectangular: A, Ixx = Sa Sb^3/12, Iyy = Sa^3 Sb/12 (minus the hollow), J0 polar; Jt and kappa
 *                for a solid or a thin-walled section.
 *   Arbitrary:   as given; kappa_x = Asx/A, kappa_y = Asy/A.
 */
BeamProperties beam_properties (const FrameSection& s, BeamTheory theory);

/**
 * The direction cosine matrix of an element from point a to point b (m) with the given spin (rad),
 * row-major: dc[3i+j] is the i-th frame component of the element's local axis j (x, y, z), so a
 * local vector v maps to the frame as dc v. As SubDyn's GetDirCos: a vertical element pointing up
 * has the identity, one pointing down flips y and z; otherwise local x is horizontal.
 * Precondition: a != b.
 */
std::array<double,9> direction_cosines (const std::array<double,3>& a, const std::array<double,3>& b, double spin);

/** The 12 x 12 stiffness of a beam element of length L (m) in its local axes, row-major (SubDyn's ElemK_Beam before rotation). */
std::array<double,144> beam_stiffness_local (const BeamProperties& p, double L);

/** The 12 x 12 stiffness in the frame's axes: T K_local T^T with T the block-diagonal of four dc. */
std::array<double,144> beam_stiffness (const BeamProperties& p, double L, const std::array<double,9>& dc);

/**
 * The nodal loads equivalent to an element's own weight (frame axes, 12 values, N and N m), as
 * SubDyn's ElemG: half the weight rho A g L at each node and the end moments of a uniformly loaded
 * fixed-fixed beam, w L^2/12 with w = rho A g, about the horizontal axes.
 */
std::array<double,12> beam_gravity_load (const BeamProperties& p, double L, const std::array<double,9>& dc, double gravity);

/**
 * The 12 x 12 consistent mass of a beam element of length L (m) in its local axes, row-major, as
 * SubDyn's ElemM_Beam: the translational mass rho A L with cubic shape functions, the rotary
 * inertia rho Ixx and rho Iyy of the cross-section, and the torsional inertia rho J0 L (kg, kg m, kg m^2).
 */
std::array<double,144> beam_mass_local (const BeamProperties& p, double L);

/** The element mass in the frame's axes: T M_local T^T. */
std::array<double,144> beam_mass (const BeamProperties& p, double L, const std::array<double,9>& dc);

/**
 * The 6 x 6 mass matrix (row-major) a concentrated mass adds at its joint, as SubDyn's
 * rigidBodyMassMatrix: the rigid body of mass m whose centre sits at the offset (x, y, z) from
 * the joint, with the inertia tensor entries Jxx Jyy Jzz Jxy Jxz Jyz about its centre, moved to
 * the joint (parallel axes; the translation-rotation coupling m [r]x). Units: kg, kg m, kg m^2.
 */
std::array<double,36> rigid_body_mass (const FrameMass& c);

/** One element: its nodes (0-based), the member it belongs to (0-based, input order), length (m), axes and properties. */
struct FrameElement
{
    std::size_t node_a{0};
    std::size_t node_b{0};
    std::size_t member{0};
    double length{0.0};
    std::array<double,9> dc{};
    BeamProperties prop;
};

/** The static response of a frame. */
struct FrameSolution
{
    std::vector<double> displacement;                    //!< 6 per node (m, rad), frame axes
    /** Per support (input order): the force and moment ON the frame BY the ground at its joint (N, N m), frame axes. */
    std::vector<std::array<double,6>> reaction;
    /** Per element: the forces and moments ON the element at node_a, then at node_b (N, N m), in the element's local axes. */
    std::vector<std::array<double,12>> element_force;
};

/**
 * A frame built from FrameInputs: nodes, elements, the supports applied, and the stiffness of the
 * free degrees of freedom factored once.
 */
class Frame
{
public:
    /**
     * Build and factor the frame, each member's E and G reduced for its temperature (FrameInputs::temperature,
     * EN 1993-1-2's k_E). Returns null with err set when the inputs are invalid
     * (FrameInputs::validate), the structure is a mechanism, naming the joint (or member node)
     * and the degree of freedom that nothing holds, or it has more than max_free_dofs free degrees of
     * freedom.
     */
    static std::unique_ptr<Frame> create (const FrameInputs& in, std::string& err);
    /**
     * The most free degrees of freedom a frame may have: its stiffness and mass are dense matrices, a
     * few of them held at once on every rank (288 MB each at this size) and factored in O(n^3).
     */
    static constexpr std::size_t max_free_dofs = 6000;

    const FrameInputs& inputs () const { return m_in; }
    std::size_t num_nodes () const { return m_x.size(); }
    std::size_t num_dofs () const { return 6 * m_x.size(); }
    /** The number of degrees of freedom not fixed by a support. */
    std::size_t num_free_dofs () const { return m_free.size(); }
    /** Node n's position (m). */
    const std::array<double,3>& node_position (std::size_t n) const;
    /** The node of the joint with this id (its index in the input). Aborts when there is none. */
    std::size_t node_of_joint (int id) const;
    const std::vector<FrameElement>& elements () const { return m_elems; }
    /** The total mass of the members and the concentrated masses (kg). */
    double total_mass () const;

    /** The assembled stiffness of every degree of freedom (num_dofs()^2, row-major), supports and springs left out. */
    std::vector<double> stiffness () const;

    /** The free degrees of freedom (indices into the num_dofs() vectors), ascending. */
    const std::vector<std::size_t>& free_dofs () const { return m_free; }
    /**
     * cK K + cM M on the free degrees of freedom, row-major (num_free_dofs()^2): K the elements'
     * stiffness and the support springs, M the elements' consistent mass, the concentrated masses
     * and the support springs' masses.
     */
    std::vector<double> assemble_free (double cK, double cM) const;
    /**
     * K u for a displacement u of every degree of freedom (num_dofs()): N, N m. The support springs
     * are included unless with_springs is false.
     */
    std::vector<double> apply_stiffness (const std::vector<double>& u, bool with_springs = true) const;
    /**
     * M a for an acceleration a of every degree of freedom (num_dofs(), concentrated masses
     * included): N, N m. The support springs' masses are included unless with_springs is false.
     */
    std::vector<double> apply_mass (const std::vector<double>& a, bool with_springs = true) const;
    /**
     * The load on every degree of freedom (num_dofs(), N and N m): node_loads (empty: none) and,
     * when gravity > 0 (m/s^2), the members' weight as element loads and the concentrated masses'
     * weight with its moment about the joint.
     */
    std::vector<double> load_vector (const std::vector<double>& node_loads, double gravity) const;
    /** Solve K x = b on the free degrees of freedom with the factored stiffness: b holds num_free_dofs() values, x on return. */
    void solve_free (std::vector<double>& b) const { m_chol.solve(b); }

    /**
     * The static response to node_loads (6 per node, frame axes, N and N m; empty: none) plus, when
     * gravity > 0 (m/s^2), the weight of the members and the concentrated masses along -z.
     */
    FrameSolution solve (const std::vector<double>& node_loads, double gravity) const;

    /**
     * The end forces of every element at the displacement u (num_dofs() values, m and rad): K_e u_e less the
     * element's weight as nodal loads when gravity > 0 (m/s^2), in the element's local axes, as
     * FrameSolution::element_force (the forces on the element at node_a, then at node_b; N, N m).
     */
    std::vector<std::array<double,12>> element_forces (const std::vector<double>& u, double gravity) const;

private:
    Frame () = default;
    std::string describe_dof (std::size_t dof) const;
    /** One element's end forces at the displacement u, in the frame's axes (12 values, N and N m). */
    std::array<double,12> end_forces (const FrameElement& el, const std::vector<double>& u, double gravity) const;
    /** node_loads and the concentrated masses' weight, the loads applied at the nodes themselves. */
    std::vector<double> point_loads (const std::vector<double>& node_loads, double gravity) const;

    FrameInputs m_in;
    std::vector<std::array<double,3>> m_x;          //!< node positions (m)
    std::vector<std::string> m_node_name;           //!< "joint <id>" or "member <id> node <k>"
    std::vector<FrameElement> m_elems;
    std::vector<std::size_t> m_free;                //!< the free degrees of freedom, ascending
    std::vector<long> m_free_index;                 //!< per degree of freedom: its index in m_free, or -1 when fixed
    DenseCholesky m_chol;                           //!< the factored stiffness of the free degrees of freedom, springs included
};

} // namespace erf_towers

#endif
