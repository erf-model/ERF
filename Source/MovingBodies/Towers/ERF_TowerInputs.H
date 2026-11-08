#ifndef ERF_TOWER_INPUTS_H
#define ERF_TOWER_INPUTS_H

/** The inputs of a lattice tower type (erf.conductors.<type>.*) and their range checks. */

#include <string>
#include <vector>

#include <AMReX_REAL.H>

namespace erf_towers {

/**
 * One lattice tower type: the geometry, weight, damping and foundation shared by every tower of
 * that type, read from the block erf.conductors.<type>.* for each name in erf.conductors.tower_types.
 * The ParmParse reader that fills this struct belongs to the body model that owns the towers (the
 * conductor lines); here the struct is filled by its caller and checked by validate().
 *
 * The tower's shaft (its "body") is a square lattice tapering linearly from base_width at the
 * ground to top_width at the cross-arm, with an optional peak above the cross-arm. The cross-arm
 * is a lattice of arm_length across the line, at the height the conductor hangs from, whose face
 * is arm_depth deep out to the shaft's corners at the arm's bottom and tapers to the tips (a truss
 * of chords converging on each tip, as the generated frame's). The wind acts on the members'
 * projected area: the solidity times the outline of a face.
 *
 *   erf.conductors.tower_types       = lattice
 *   erf.conductors.lattice.base_width = 6.     (m, the shaft's width at the ground; required, > 0)
 *   erf.conductors.lattice.top_width  = 1.5    (m, the shaft's width at the cross-arm; required, in (0, base_width])
 *   erf.conductors.lattice.solidity   = 0.2    (-, the members' area over a face's outline, shaft and cross-arm;
 *                                               required, in (0, 1))
 *   erf.conductors.lattice.arm_length = 12.    (m, the cross-arm's length, centred on the shaft; required, > 0,
 *                                               and with arm_outside_shaft longer than the shaft is wide at the
 *                                               middle of the arm's face, at least top_width)
 *   erf.conductors.lattice.arm_depth  = 1.5    (m, the depth of the cross-arm's face at the shaft, kept out to
 *                                               the shaft's corners at the arm's bottom, then tapering linearly
 *                                               to 0 at the tips; >= 0, default 0: top_width;
 *                                               with arm_outside_shaft under twice the cross-arm's height)
 *   erf.conductors.lattice.peak       = 0.     (m, the shaft above the cross-arm, at top_width; >= 0, default 0)
 *   erf.conductors.lattice.drag_coefficient = 0. (-, the force coefficient on the members' area; >= 0,
 *                                                 default 0: 4 phi^2 - 5.9 phi + 4 from the solidity phi)
 *   erf.conductors.lattice.segments   = 10     (the number of equal segments of the shaft up to the cross-arm,
 *                                               one drag node at the middle of each; >= 1, default 10)
 *   erf.conductors.lattice.diagonal_wind_factor = true (the shaft's and the peak's drag times 1 + g sin^2(2 theta), theta the
 *                                               wind's angle from a face's normal and g = min(0.75 solidity, 0.2):
 *                                               ASCE 7's factor along a diagonal; default true, false: none)
 *   erf.conductors.lattice.arm_outside_shaft = true (the cross-arm's drag nodes on its two outer parts only, from
 *                                               the shaft's half width at the middle of the face's depth at the
 *                                               shaft to the tips, each at the centroid of its segment's tapering
 *                                               face; default true, and then arm_length must exceed the shaft's
 *                                               width there; false, the earlier layout: along its whole length,
 *                                               through the shaft, at the cross-arm's height, arm_depth deep all
 *                                               along)
 *   erf.conductors.lattice.weight     = 9.0e4  (N, the tower's own weight on its foundation; >= 0, default 0;
 *                                               required > 0 when frequency > 0, since it sets the mass)
 *   erf.conductors.lattice.leg_spacing = 6.    (m, between the four legs' footings, a square; >= 0,
 *                                               default 0: base_width; not with a frame, whose supports are its footings)
 *   erf.conductors.lattice.allowable_uplift      = 1.0e5 (N, a footing's pull-out capacity; >= 0, default 0: not checked)
 *   erf.conductors.lattice.allowable_compression = 3.0e5 (N, a footing's bearing capacity; >= 0, default 0: not checked)
 *   erf.conductors.lattice.frequency  = 2.0    (Hz, the first bending frequency on a rigid foundation; >= 0,
 *                                               default 0: the tower stands still and the lines hang from
 *                                               fixed MoorDyn points; > 0 needs weight > 0)
 *   erf.conductors.lattice.damping_ratio = 0.02 (-, the fraction of critical damping of that mode; in [0, 1),
 *                                               default 0.02; needs a tower that moves, and > 0 with a frame;
 *                                               on a still tower with frequency = 0 given: a warning)
 *   erf.conductors.lattice.foundation_rotational_stiffness = 1.0e9 (N m/rad, the footing's resistance to tilting;
 *                                                                   >= 0, default 0: rigid; > 0 needs frequency,
 *                                                                   or with frequency = 0 given: a warning)
 *   erf.conductors.lattice.foundation_lateral_stiffness    = 1.0e8 (N/m, the footing's resistance to sliding;
 *                                                                   >= 0, default 0: rigid; > 0 needs frequency,
 *                                                                   or with frequency = 0 given: a warning)
 *   erf.conductors.lattice.frame_file = tower.dat (a SubDyn input file of the tower's frame in tower-local axes:
 *                                       origin at the base centre, x along the line, y along the cross-arm,
 *                                       z up; default none. The tower then bends as that frame model, which
 *                                       gives its stiffness, mass and footings: frequency, weight and the
 *                                       foundation stiffnesses are not given with it. damping_ratio sets the
 *                                       frame's Rayleigh damping at its first natural frequency and ten times it)
 *   erf.conductors.lattice.member_file = members.dat (the design data of every member of frame_file's frame, for
 *                                       the members' strength checks (ERF_MemberChecks.H: role, yield strength,
 *                                       angle, net area, end conditions); default none: no checks; needs frame_file)
 *   erf.conductors.lattice.frame_panels = 10 (generate the tower's frame (ERF_LatticeFrame.H) from this type's
 *                                       dimensions, with this many panels of the shaft below the cross-arm;
 *                                       >= 0, at most 200, default 0: none. Not given with frame_file; like it,
 *                                       not with frequency, weight or the foundation stiffnesses)
 *   erf.conductors.lattice.leg_angle   = 0.15 0.012 (m, the leg width and thickness of the equal-leg angle of the
 *                                       generated frame's legs and cross-arm chords; required with frame_panels,
 *                                       0 < thickness < width)
 *   erf.conductors.lattice.brace_angle = 0.09 0.007 (m, the same for every other member; required with frame_panels)
 *   erf.conductors.lattice.bracing     = crossed (the generated frame's face bracing: crossed, two diagonals per
 *                                       panel face bolted where they cross, or single, one diagonal alternating
 *                                       from panel to panel; default crossed; only with frame_panels)
 *   erf.conductors.lattice.yield_strength = 3.45e8 (Pa, the generated frame's steel; > 0, default 0: 3.45e8;
 *                                          only with frame_panels)
 *   erf.conductors.lattice.angle_principal_axes = true (turn every generated angle (SubDyn's MSpin) so that its
 *                                          principal axes lie as on a real tower: a leg's or a chord's axis of
 *                                          symmetry towards the corner it stands on, a brace's at 45 degrees to
 *                                          its face; default true; false: MSpin 0, the hangers always; only with
 *                                          frame_panels)
 *   erf.conductors.lattice.steel_temperature = 20. (C, the steel's temperature, which lowers its stiffness and
 *                                          yield strength by EN 1993-1-2; above -273.15 and below 1200, default 20;
 *                                          only with a frame model, frame_file or frame_panels)
 *   erf.conductors.lattice.epsilon = 3. (> 0; the width of the kernel that spreads the drag of this type's towers'
 *                                          drag nodes into the flow, in units of the anchor level's x spacing dx;
 *                                          default erf.conductors.epsilon; needs erf.conductors.drag_on_flow = true
 *                                          as that key does, a warning with drag_on_flow = false given)
 *
 * Every value must also be finite. The members mirror the keys, with the same units.
 */
struct TowerType
{
    std::string name;
    amrex::Real base_width{0.0};                       //!< m
    amrex::Real top_width{0.0};                        //!< m
    amrex::Real solidity{0.0};                         //!< -
    amrex::Real arm_length{0.0};                       //!< m
    amrex::Real arm_depth{0.0};                        //!< m; 0: top_width
    amrex::Real peak{0.0};                             //!< m
    amrex::Real drag_coefficient{0.0};                 //!< -; 0: the lattice curve of the solidity
    int segments{10};                                  //!< equal segments of the shaft up to the cross-arm
    amrex::Real weight{0.0};                           //!< N
    amrex::Real leg_spacing{0.0};                      //!< m; 0: base_width
    amrex::Real allowable_uplift{0.0};                 //!< N; 0: not checked
    amrex::Real allowable_compression{0.0};            //!< N; 0: not checked
    amrex::Real frequency{0.0};                        //!< Hz; 0: the tower stands still
    bool frequency_given{false};                       //!< frequency was given (0 then turns the motion off on purpose)
    amrex::Real damping_ratio{0.02};                   //!< - (fraction of critical)
    bool damping_given{false};                         //!< damping_ratio was given (it needs a tower that moves)
    amrex::Real foundation_rotational_stiffness{0.0};  //!< N m/rad; 0: rigid
    amrex::Real foundation_lateral_stiffness{0.0};     //!< N/m; 0: rigid
    std::string frame_file;                            //!< SubDyn input file of the tower's frame; empty: none
    std::string member_file;                           //!< the design data of frame_file's members; empty: none
    int frame_panels{0};                               //!< the generated frame's shaft panels; 0: no generated frame
    std::vector<amrex::Real> leg_angle;                //!< m: width, thickness
    std::vector<amrex::Real> brace_angle;              //!< m: width, thickness
    std::string bracing;                               //!< crossed or single; empty: crossed
    amrex::Real yield_strength{0.0};                   //!< Pa; 0: 3.45e8
    amrex::Real steel_temperature{20.0};               //!< C
    bool diagonal_wind_factor{true};                   //!< the shaft's and peak's wind-direction factor (lattice_diagonal_gain)
    bool arm_outside_shaft{true};                      //!< the cross-arm's drag outside the shaft, at its pieces' centroids
    bool angle_principal_axes{true};                   //!< the generated angles turned to their principal axes
    bool angle_axes_given{false};                      //!< angle_principal_axes was given (it needs frame_panels)
    amrex::Real epsilon{0.0};                          //!< the drag kernel's width (anchor level's dx); used when epsilon_given
    bool epsilon_given{false};                         //!< epsilon was given; else the towers take erf.conductors.epsilon

    /** The number of drag nodes on the cross-arm, one per equal segment at its face's centroid (half on each side). */
    static constexpr int arm_segments = 4;
    static_assert(arm_segments % 2 == 0, "the cross-arm's drag nodes come in pairs, one per side");

    /** The depth of the cross-arm's face (m): arm_depth, or top_width when arm_depth is 0. */
    amrex::Real arm_face () const { return arm_depth > 0.0 ? arm_depth : top_width; }
    /** The spacing of the legs' footings (m): leg_spacing, or base_width when leg_spacing is 0. */
    amrex::Real legs () const { return leg_spacing > 0.0 ? leg_spacing : base_width; }
    /** Whether the tower bends as a frame model, read from frame_file or generated (frame_panels > 0). */
    bool has_frame () const { return !frame_file.empty() || frame_panels > 0; }
    /** Whether the tower bends under its loads (a frequency or a frame model is given) or stands still. */
    bool moves () const { return frequency > 0.0 || has_frame(); }
    /** The generated frame's yield strength (Pa): yield_strength, or 3.45e8 when it is 0. */
    amrex::Real yield () const { return yield_strength > 0.0 ? yield_strength : amrex::Real(3.45e8); }
    /** The force coefficient on the members' area (-): drag_coefficient, or the lattice value of the solidity. */
    amrex::Real force_coefficient () const;
    /**
     * Check every range: returns an empty string when the type is acceptable, else the message, naming the key.
     * Keys of a tower's motion on a tower that stands still are refused when frequency was never given.
     */
    std::string validate () const;
    /** With frequency = 0 given and no frame, the keys of a tower's motion that are given but unused (a warning). */
    std::string unused_on_a_still_tower () const;
};

} // namespace erf_towers

#endif
