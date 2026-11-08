#ifndef ERF_MEMBER_CHECKS_H
#define ERF_MEMBER_CHECKS_H

/**
 * The strength checks of a lattice tower's members: the axial force each member of a frame model
 * carries against its design strength in tension and in compression, as ASCE 10-15 (Design of
 * Latticed Steel Transmission Structures) sets them for angle members, and the reduction of the
 * steel's yield strength and stiffness with its temperature, as EN 1993-1-2 (Eurocode 3, structural
 * fire design) Table 3.1 gives them.
 *
 * Conventions: SI units (m, N, Pa); axial forces are tension positive; temperatures in degrees
 * Celsius. A member's utilisation is its axial force over its design strength (1: at the strength).
 * The design strengths are nominal: the loads are the frame's as computed, without load factors.
 */

#include <array>
#include <string>
#include <vector>

#include "ERF_Frame.H"
#include "ERF_FrameInputs.H"

namespace erf_towers {

/** The section properties of an equal-leg angle of leg width b and thickness t (m), without its root fillet. */
struct AngleProperties
{
    double A{0.0};        //!< area t (2b - t) (m^2)
    double Ig{0.0};       //!< second moment about a centroidal axis parallel to a leg (m^4)
    double Iu{0.0};       //!< about the major principal axis (m^4)
    double Iv{0.0};       //!< about the minor principal axis (m^4), the smallest
    double J{0.0};        //!< Saint-Venant torsion constant of the two thin legs, (2b - t) t^3 / 3 (m^4)
    double rv{0.0};       //!< the least radius of gyration, sqrt(Iv/A) (m)
    double centroid{0.0}; //!< the centroid's distance from the outer face of each leg (m)
};

/** The properties of the equal-leg angle b x t (m; 0 < t < b), the two legs as rectangles meeting at the heel. */
AngleProperties angle_properties (double b, double t);

/** What a member does in the tower, which sets its effective slenderness and its slenderness limit (ASCE 10-15). */
enum class MemberRole {
    Leg,        //!< a main member bolted in both faces at its connections (the legs, the cross-arm's chords): KL/r = L/r, L/r <= 150
    Bracing,    //!< a member carrying load between the main members: KL/r from its end conditions, KL/r <= 200
    Redundant   //!< a member bracing others, with no calculated load of its own: KL/r from its end conditions, KL/r <= 250
};

/** How a bracing or redundant member's ends are loaded, for L/r <= 120 (ASCE 10-15 curves 1 to 3). */
enum class EndLoading {
    Concentric,      //!< concentric at both ends:                  KL/r = L/r            (curve 1)
    OneEccentric,    //!< concentric at one end, normal framing eccentricity at the other: KL/r = 30 + 0.75 L/r (curve 2)
    BothEccentric    //!< normal framing eccentricity at both ends:  KL/r = 60 + 0.5 L/r   (curve 3)
};

/** How a bracing or redundant member's ends are restrained against rotation, for L/r > 120 (ASCE 10-15 curves 4 to 6). */
enum class EndRestraint {
    None,       //!< unrestrained at both ends:            KL/r = L/r              (curve 4)
    OneEnd,     //!< partially restrained at one end:      KL/r = 28.6 + 0.762 L/r (curve 5)
    BothEnds    //!< partially restrained at both ends:    KL/r = 46.2 + 0.615 L/r (curve 6)
};

/**
 * The design data of one member: its role and end conditions, the steel's yield strength, and the
 * equal-leg angle it is made of. The frame's section gives the member's stiffness and mass; this
 * gives its strength. A member with b = t = 0 is not an angle: its radius of gyration is the least
 * of its frame section's, its area the section's, and it has no local buckling check.
 */
struct MemberDesign
{
    int member{0};                                  //!< the frame member's id
    MemberRole role{MemberRole::Bracing};
    double yield{0.0};                              //!< the yield strength Fy at 20 C (Pa, > 0)
    double b{0.0};                                  //!< the angle's leg width (m; 0: not an angle)
    double t{0.0};                                  //!< the angle's thickness (m; 0: not an angle)
    double net_area{1.0};                           //!< the net area over the gross area in tension, the bolt holes taken out (in (0, 1])
    bool one_leg{true};                             //!< an angle bolted by one leg only: its tension stress is 0.9 Fy on the net area
    EndLoading ends{EndLoading::BothEccentric};
    EndRestraint restraint{EndRestraint::None};
};

/** The result of checking one member under one set of loads. */
struct MemberCheck
{
    int member{0};                  //!< the frame member's id
    MemberRole role{MemberRole::Bracing};
    double length{0.0};             //!< joint to joint, the unbraced length L (m)
    double r{0.0};                  //!< the radius of gyration for buckling (m)
    double slenderness{0.0};        //!< L/r
    double effective{0.0};          //!< KL/r
    double limit{0.0};              //!< the largest slenderness ASCE 10-15 allows the role (L/r for legs, KL/r otherwise)
    double w_t{0.0};                //!< the angle leg's flat width over its thickness, (b - t)/t (0: not an angle)
    double yield{0.0};              //!< Fy at the member's temperature (Pa)
    double E{0.0};                  //!< Young's modulus at the member's temperature (Pa)
    double compression_stress{0.0}; //!< the design compressive stress Fa (Pa)
    double tension_capacity{0.0};   //!< the design tensile strength (N)
    double compression_capacity{0.0}; //!< the design compressive strength, Fa times the gross area (N)
    double tension{0.0};            //!< the largest tension along the member (N, >= 0)
    double compression{0.0};        //!< the largest compression along the member (N, >= 0)
    double utilisation{0.0};        //!< max(tension / tension_capacity, compression / compression_capacity)
    bool slender{false};            //!< the slenderness exceeds the role's limit
    bool thin{false};               //!< w/t exceeds 25, ASCE 10-15's largest for angles
};

/** The effective slenderness KL/r of a member of the role, end loading and restraint, at slenderness L/r (ASCE 10-15). */
double effective_slenderness (MemberRole role, EndLoading ends, EndRestraint restraint, double slenderness);

/** The slenderness limit of the role (ASCE 10-15): 150 (legs, on L/r), 200 (bracing) or 250 (redundant), on KL/r. */
double slenderness_limit (MemberRole role);

/**
 * The yield strength reduced for local buckling of an angle leg of flat width over thickness w_t
 * (ASCE 10-15): Fy up to (w/t)_lim = 80 psi/sqrt(Fy), then [1.677 - 0.677 (w/t)/(w/t)_lim] Fy up to
 * 144 psi/sqrt(Fy), then 0.0332 pi^2 E/(w/t)^2. Written with sqrt(E/Fy) in place of psi/sqrt(Fy), so
 * the limits are 80 sqrt(E/Fy)/sqrt(29000) and 144 sqrt(E/Fy)/sqrt(29000): the standard's values for
 * E = 29000 ksi (200 GPa), and the same in any units and at any temperature's E. Pa in and out.
 */
double local_buckling_stress (double w_t, double yield, double E);

/**
 * The design compressive stress Fa at effective slenderness KL/r (ASCE 10-15): [1 - (KL/r / Cc)^2 / 2] Fy
 * up to Cc = pi sqrt(2 E / Fy), then pi^2 E / (KL/r)^2. yield is the local-buckling stress where that
 * governs (local_buckling_stress()). Pa in and out.
 */
double compression_stress (double effective, double yield, double E);

/**
 * The reduction factors of carbon steel at temperature theta (C), EN 1993-1-2 Table 3.1, linear
 * between its rows: k_y of the effective yield strength and k_E of the slope of the linear elastic
 * range. Both are 1 up to 100 C (and below 20 C); both reach 0 at 1200 C.
 */
void steel_reduction (double theta, double& k_y, double& k_E);

/**
 * Check one member. length (m) is its unbraced length; A (m^2) and r_section (m) the area and the
 * least radius of gyration of its frame section (used when the design is not an angle); E (Pa) its
 * Young's modulus at 20 C; theta (C) its temperature; tension and compression (N, >= 0) the largest
 * of each along it.
 */
MemberCheck check_member (const MemberDesign& d, double length, double A, double r_section, double E, double theta,
                          double tension, double compression);

/**
 * The checks of every member of a frame under the element forces f (one per element, as
 * FrameSolution::element_force: the forces on the element at its nodes in its local axes, N).
 * designs holds one row per frame member in the frame's member order (match_designs()); theta the
 * members' temperatures (C, one per member; empty: 20 C).
 */
std::vector<MemberCheck> check_members (const Frame& frame, const std::vector<MemberDesign>& designs,
                                        const std::vector<std::array<double,12>>& f, const std::vector<double>& theta);

/**
 * Put designs in the order of in.members and check them: one row per member, every row a member of
 * the frame, the yield strength, the angle and the net area in range, and an angle's area within 10 %
 * of its frame section's (else the strength and the stiffness describe different members). Returns an
 * empty string, or the first problem, naming the file (where) and the member.
 */
std::string match_designs (const FrameInputs& in, std::vector<MemberDesign>& designs, const std::string& where);

/**
 * Read a member design file: one row per member,
 *   MemberID  Role  Fy(Pa)  b(m)  t(m)  NetArea(-)  Bolted  Ends  Restraint
 * with Role leg | bracing | redundant; Bolted one | both (the legs of the angle bolted at its ends);
 * Ends concentric | one | both (the ends with a normal framing eccentricity, ASCE 10-15 curves 1 to 3);
 * Restraint none | one | both (the ends partially restrained against rotation, curves 4 to 6). Blank
 * lines and lines starting with # or ! are skipped. Returns an empty string, or the problem with its line.
 */
std::string read_member_designs (const std::string& path, std::vector<MemberDesign>& designs);

/**
 * Write designs in the layout read_member_designs() reads, with a header naming the columns. Returns false
 * when the file cannot be written. Not collective: it writes on the calling rank, so call it on the I/O rank only.
 */
bool write_member_designs (const std::string& path, const std::vector<MemberDesign>& designs, const std::string& title);

/** The word a member file uses for a role: leg, bracing or redundant. */
const char* role_name (MemberRole role);

} // namespace erf_towers

#endif
