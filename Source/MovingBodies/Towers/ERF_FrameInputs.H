#ifndef ERF_FRAME_INPUTS_H
#define ERF_FRAME_INPUTS_H

/**
 * The inputs of a lattice tower's frame model: joints, beam members, cross-sections, supports and
 * concentrated masses, read from a SubDyn input file (OpenFAST's SubDyn module, the subset a
 * free-standing lattice tower needs), and their checks. All values in SI units, double precision.
 */

#include <array>
#include <string>
#include <vector>

namespace erf_towers {

/** The beam theory of every element (SubDyn's FEMMod): Euler-Bernoulli (1) or Timoshenko with shear deformation (3). */
enum class BeamTheory { EulerBernoulli, Timoshenko };

/** The shape a cross-section is given by (SubDyn's member type MType: 1c circular, 1r rectangular, 4 arbitrary). */
enum class SectionShape { Circular, Rectangular, Arbitrary };

/** One joint (SubDyn's STRUCTURE JOINTS row): an id and its position (m, the frame's coordinates, z up). */
struct FrameJoint
{
    int id{0};
    std::array<double,3> x{{0.0, 0.0, 0.0}};
};

/**
 * One cross-section property set (a row of SubDyn's circular, rectangular or arbitrary table).
 * Material: E (Pa), G (Pa), rho (kg/m^3). Geometry by shape:
 *   Circular:    D outer diameter (m), t wall thickness (m; 0: solid)
 *   Rectangular: Sa along the local x axis (m), Sb along the local y axis (m), t wall (m; 0: solid)
 *   Arbitrary:   A (m^2), shear areas Asx, Asy (m^2), Ixx, Iyy second moments about the local x and y
 *                axes (m^4), J0 polar moment (m^4), Jt torsion constant (m^4)
 */
struct FrameSection
{
    int id{0};
    SectionShape shape{SectionShape::Arbitrary};
    double E{0.0}, G{0.0}, rho{0.0};
    double D{0.0}, t{0.0}, Sa{0.0}, Sb{0.0};
    double A{0.0}, Asx{0.0}, Asy{0.0}, Ixx{0.0}, Iyy{0.0}, J0{0.0}, Jt{0.0};
};

/**
 * One beam member (a SubDyn MEMBERS row): from joint joint_a to joint_b (ids), with section id
 * section of the given shape, rotated about its own axis by spin (rad, SubDyn's MSpin in degrees
 * on input). The member's local z axis runs from joint_a to joint_b.
 */
struct FrameMember
{
    int id{0};
    int joint_a{0};
    int joint_b{0};
    int section{0};
    SectionShape shape{SectionShape::Arbitrary};
    double spin{0.0};
};

/**
 * One support (a SubDyn BASE REACTION JOINTS row): the joint, which of its six degrees of freedom
 * (x, y, z translations, then rotations about x, y, z) are fixed, and a soil-structure spring on
 * the free ones from the row's SSI file: the upper triangle of a symmetric 6 x 6 stiffness in
 * SubDyn's order Kxx Kxy Kyy Kxz Kyz Kzz Kxtx Kytx Kztx Ktxtx Kxty Kyty Kzty Ktxty Ktyty Kxtz Kytz
 * Kztz Ktxtz Ktytz Ktztz (N/m, N/rad, N m/rad), and the same for its mass (kg, kg m, kg m^2).
 */
struct FrameSupport
{
    int joint{0};
    std::array<bool,6> fixed{{true, true, true, true, true, true}};
    std::string ssi_file;
    std::array<double,21> stiffness{};
    std::array<double,21> mass{};
};

/** A concentrated mass on a joint (a SubDyn row): mass (kg), its centre's offset from the joint (m). Inertia is read for dynamics. */
struct FrameMass
{
    int joint{0};
    double mass{0.0};
    std::array<double,6> inertia{};             //!< Jxx Jyy Jzz Jxy Jxz Jyz about the mass centre (kg m^2)
    std::array<double,3> offset{{0.0, 0.0, 0.0}};
};

/** Everything the frame model is built from. */
struct FrameInputs
{
    std::string file;                            //!< the SubDyn file read (for messages)
    BeamTheory theory{BeamTheory::EulerBernoulli};
    int divisions{1};                            //!< elements per member (SubDyn's NDiv), >= 1
    std::vector<FrameJoint> joints;
    std::vector<FrameSection> sections;
    std::vector<FrameMember> members;
    std::vector<FrameSupport> supports;
    std::vector<int> interface_joints;           //!< SubDyn's INTERFACE JOINTS ids (a tower's frame: one, its cross-arm's centre)
    std::vector<FrameMass> masses;
    /**
     * The steel's temperature of each member (C, input order; empty: 20 C everywhere), which reduces its
     * Young's and shear moduli by EN 1993-1-2's k_E (steel_reduction() in ERF_MemberChecks.H). Not part
     * of a SubDyn file.
     */
    std::vector<double> temperature;

    /** The index (0-based, input order) of the joint with this id, or -1. */
    int joint_index (int id) const;
    /** The section with this id and shape, or null. */
    const FrameSection* section (int id, SectionShape shape) const;
    /** An empty string when the inputs describe a frame the model can build, else the first problem found, naming the file and item. */
    std::string validate () const;
};

/**
 * Read a whole text file into text: the I/O rank reads it and broadcasts it, so every rank parses the
 * same contents and the file system serves one reader. Collective: every rank calls it. False, on
 * every rank, when the file cannot be read.
 */
bool read_text_file (const std::string& path, std::string& text);

/**
 * Read a SubDyn input file (OpenFAST 5.0 layout). Returns an empty string on success, else a
 * message naming the file, the section and the problem. Read: FEMMod (1 or 3) and NDiv; joints
 * (cantilever joints only); base reaction joints with fixed or free DOFs and an optional SSI file
 * (relative to the SubDyn file's directory); interface joints; beam members of type 1/1c, 1r and 4
 * with equal property sets at both ends; the circular, rectangular and arbitrary property tables;
 * concentrated masses. Cable, rigid-link and spring members are refused; their tables, the
 * direction cosine matrices, the damping inputs and the output sections are skipped. Collective:
 * the files are read through read_text_file().
 */
std::string read_subdyn (const std::string& path, FrameInputs& in);

/**
 * Write the frame as a SubDyn input file (OpenFAST 5.0 layout) that read_subdyn() and SubDyn read
 * back to the same frame: the joints, the supports (a support with a spring writes its SSI file next
 * to path, named <path without .dat>_ssi_<joint>.dat), the interface joints, the members and their
 * sections, the concentrated masses; every value at 17 significant digits. The members' temperatures
 * are not written. title is the file's second line. Returns an empty string, or the problem.
 * Not collective: it writes on the calling rank, so call it on the I/O rank only.
 */
std::string write_subdyn (const FrameInputs& in, const std::string& path, const std::string& title);

} // namespace erf_towers

#endif
