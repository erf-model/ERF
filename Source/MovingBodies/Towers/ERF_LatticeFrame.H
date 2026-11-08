#ifndef ERF_LATTICE_FRAME_H
#define ERF_LATTICE_FRAME_H

/**
 * A lattice tower's frame model generated from a few dimensions: a square shaft of four angle legs
 * tapering from the base to the cross-arm, horizontal struts at every panel level and braced faces;
 * a cross-arm on each side, a truss of four chords converging on its tip; a hanger at the centre of
 * the cross-arm, where a line attached at the tower's centre hangs; and an optional peak above the
 * cross-arm. Every member is an equal-leg angle; the result is the FrameInputs a SubDyn file would
 * give (ERF_FrameInputs.H) and the design data of every member (ERF_MemberChecks.H).
 *
 * Axes: tower-local, as a frame file for a tower (ERF_FrameTower.H): origin at the base centre, x
 * along the line, y along the cross-arm, z up. SI units.
 *
 * Joints (ids from 1, in this order):
 *   the shaft's levels from the ground up, four corners each, corner c at (sx h, sy h) with
 *     (sx, sy) = (+,+), (-,+), (-,-), (+,-) and h the half-width at the level's height: the panels'
 *     levels, equally spaced from the ground to the bottom of the cross-arm (arm_height - arm_depth),
 *     then the cross-arm's top (arm_height), then the peak's levels when peak > 0, in equal panels no
 *     taller than the shaft's, up to arm_height + peak;
 *   with crossed bracing, the point where each panel face's two diagonals cross;
 *   the centre of the cross-arm (0, 0, arm_height), the interface joint SubDyn needs;
 *   the cross-arms' panel points, +y side first, and their tips (0, +-arm_length/2, arm_height).
 * Members (ids from 1, in this order): the legs; the struts; the face bracing; the hanger; the
 * cross-arms. The four base corners are the supports, fixed in all six degrees of freedom.
 */

#include <string>
#include <vector>

#include "ERF_FrameInputs.H"
#include "ERF_MemberChecks.H"

namespace erf_towers {

/** The dimensions, the angles and the steel of a generated lattice tower. */
struct LatticeSpec
{
    double base_width{0.0};   //!< the shaft's width at the ground (m, > 0)
    double top_width{0.0};    //!< the shaft's width at the cross-arm (m, in (0, base_width]); the peak keeps it
    double arm_height{0.0};   //!< the height of the cross-arm's top, where the lines hang (m, > 0)
    double arm_length{0.0};   //!< tip to tip (m, > top_width)
    double arm_depth{0.0};    //!< the cross-arm's depth at the shaft (m, in (0, arm_height))
    double peak{0.0};         //!< the shaft above the cross-arm (m, >= 0)
    int panels{0};            //!< the shaft's panels below the cross-arm (>= 1)
    bool crossed{true};       //!< crossed diagonals on every panel face, bolted where they cross; else one diagonal, alternating
    bool principal_axes{true};//!< turn each angle (MSpin) so that its axis of symmetry lies as on a real tower; else MSpin 0
    double leg_b{0.0};        //!< the legs' and the cross-arms' chords' angle: leg width (m)
    double leg_t{0.0};        //!< and thickness (m)
    double brace_b{0.0};      //!< every other member's angle: leg width (m)
    double brace_t{0.0};      //!< and thickness (m)
    double E{2.0e11};         //!< Young's modulus (Pa)
    double G{7.7e10};         //!< shear modulus (Pa)
    double rho{7850.0};       //!< density (kg/m^3)
    double yield{3.45e8};     //!< yield strength (Pa)

    /** An empty string when the dimensions describe a tower, else the first problem. */
    std::string validate () const;
};

/**
 * Generate the tower's frame (Euler-Bernoulli beams, one element per member) and the design data of its
 * members: the legs and the cross-arms' chords as legs, bolted in both faces; every other member an
 * angle bolted by one leg with a normal framing eccentricity at both ends and no rotational restraint
 * (ASCE 10-15 curves 3 and 4), as bracing, but for the struts between crossed diagonals, which are
 * redundant. load_joints, when not null, receives the ids of the joints a load on the tower may be tied
 * to: every joint but the crossings of diagonals, which only the diagonals' bending holds out of their
 * face; a panel's wind then acts at its leg joints, as lattice tower analyses apply it. Returns an empty
 * string, or the problem with the dimensions.
 */
std::string lattice_frame (const LatticeSpec& spec, FrameInputs& in, std::vector<MemberDesign>& designs,
                           std::vector<int>* load_joints = nullptr);

/**
 * The ids of the joints of a frame read from a file that a load on the tower may be tied to: those a
 * member within 20 degrees of the vertical or of the horizontal meets (a leg, a chord or a strut),
 * not the crossings of diagonals, which only the diagonals' bending holds out of their face.
 */
std::vector<int> square_joints (const FrameInputs& in);

} // namespace erf_towers

#endif
