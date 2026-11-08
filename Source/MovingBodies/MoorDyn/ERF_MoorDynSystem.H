#ifndef ERF_MOORDYN_SYSTEM_H
#define ERF_MOORDYN_SYSTEM_H

/**
 * ERF's wrapper of one MoorDyn-C system (the v2 C API, from the installed library or the bundled
 * stub), with ERF's error handling, and helpers for the library version and the FPE traps.
 */

#include <array>
#include <cfenv>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include <AMReX.H>

#include <moordyn/MoorDyn2.h>

namespace erf_moordyn {

/** The MoorDyn-C version ERF was configured against (from its CMake package), "stub" for the bundled stand-in. */
std::string library_version ();
/** Whether the bundled stub library stands in for MoorDyn-C. */
bool is_stub ();
/** The name of a MoorDyn error code (MOORDYN_SUCCESS, MOORDYN_INVALID_INPUT_FILE, ...). */
std::string error_name (int code);
/**
 * Masks the floating-point traps (amrex.fpe_trap_invalid, amrex.fpe_trap_zero and
 * amrex.fpe_trap_overflow, as set on the calling thread) for its lifetime, with
 * amrex::disableFPExcept; on destruction it discards the exception flags raised meanwhile, puts
 * back the flags raised before it and restores the traps with amrex::setFPExcept, so that a
 * masked exception is not reported later from unrelated code.
 *
 * MoorDyn-C is not trap clean: in 2.7.1 the stationary initial-condition solver
 * (StationaryScheme::Step, under MoorDyn_Init) bounds its pseudo-step by cfl * length / v for
 * every point, whose length is DBL_MAX (no length), with cfl <= 1 and v a residual speed: the
 * quotient overflows to +inf when v < cfl (and divides by zero when v = 0), and the std::min it
 * goes into discards the infinity, so the masked result is the untrapped one. With
 * amrex.fpe_trap_overflow on, the run stopped there (SIGFPE on Linux, SIGILL on Apple arm64).
 * MoorDynSystem::init and MoorDynSystem::step run MoorDyn under this mask, so a run with lines
 * may keep the traps on for ERF's own code. MoorDyn's results are checked instead: a line
 * integration that diverges (its internal step too long) overflows and makes NaNs inside
 * MoorDyn_Step, which then reports MOORDYN_NAN_ERROR and ERF aborts naming the internal step,
 * where the traps stopped the run with a SIGFPE inside MoorDyn; step() also aborts on a
 * non-finite coupled force. The traps are a per-thread setting: the mask applies to the calling
 * thread, which is the thread MoorDyn runs on (MoorDyn-C built without OpenMP, as
 * Build/setup_moordyn.sh builds it, starts no threads). Linux (glibc) and macOS only, as
 * amrex::disableFPExcept; elsewhere it masks nothing.
 */
class FPTrapMask
{
public:
    FPTrapMask ();
    ~FPTrapMask ();
    FPTrapMask (const FPTrapMask&) = delete;
    FPTrapMask& operator= (const FPTrapMask&) = delete;

private:
    amrex::FPExcept m_traps;
    std::fexcept_t m_flags{};
};

/**
 * Check the OPTIONS block of a MoorDyn input file for `1 WaveKin`, the external fluid kinematics
 * ERF's wrapper relies on. Returns an empty string when the option is there, else the message,
 * naming the file. MoorDyn-C 2.7.1 reports kinematics points whatever WaveKin is and then
 * ignores the velocities ERF sets unless WaveKin = 1, so the file itself is checked. Every rank
 * that calls it reads the file (as MoorDyn_Create does on every rank that creates a system).
 */
std::string check_wave_kinematics_option (const std::string& input_file);

/**
 * One MoorDyn-C system (one input file) behind the v2 C API, with ERF's error handling:
 * creation, initialisation and external_kinematics_init() report their failure in a string, so a
 * body can name its input; every other call aborts with the MoorDyn function and error code,
 * since nothing sensible can follow it.
 *
 * A "line" here is a MoorDyn line (for overhead conductor lines, one span or one insulator string;
 * the conductor line is the whole system). MoorDyn lines and points are numbered from 1, as in
 * MoorDyn. Node indices are 0-based: 0 at a MoorDyn line's first end (end A) to
 * line_num_nodes(line) - 1 at its second end (end B).
 *
 * The fluid velocity the MoorDyn lines feel comes through MoorDyn's external wave-kinematics
 * calls: external_kinematics_init() returns the number of kinematics points MoorDyn wants the
 * velocity at (the nodes of MoorDyn line 1, nodes 0..N, then of line 2, and so on; then the nodes
 * of any rods; then every point, fixed, coupled or free, in input order; then every body, the
 * ground body at MoorDyn's origin first. For ERF's inputs, which have no rods or bodies, that is
 * the line nodes, the points and one entry at the origin), kinematics_points() their current
 * coordinates, and set_kinematics() hands the velocity and acceleration there before each step.
 * The input file must set `WaveKin = 1` for this (create() checks it), and the nodes must lie
 * below z = 0 in MoorDyn's frame, which is where MoorDyn applies its fluid loads.
 *
 * MoorDyn's frame is the coordinates of the input file. A body model for lines in air writes ERF
 * coordinates lowered by an offset (erf.conductors.surface_offset for the conductor lines), so that
 * the lines lie below MoorDyn's free surface z = 0.
 */
class MoorDynSystem
{
public:
    /**
     * Create a system from a MoorDyn input file. Returns null and sets err when MoorDyn cannot
     * read the file or when its OPTIONS do not set `WaveKin = 1`. log_level is a
     * MOORDYN_*_LEVEL value; log_file, when not empty, receives MoorDyn's log, written by every
     * rank that passes it (and MoorDyn-C writes its own outputs next to the input file when the
     * input asks for them): give a per-rank name, or pass it on one rank only.
     */
    static std::unique_ptr<MoorDynSystem> create (const std::string& input_file,
                                                  const std::string& log_file,
                                                  int log_level,
                                                  std::string& err);
    ~MoorDynSystem ();
    MoorDynSystem (const MoorDynSystem&) = delete;
    MoorDynSystem& operator= (const MoorDynSystem&) = delete;

    const std::string& input_file () const { return m_file; }

    /** The number of coupled degrees of freedom: 3 per coupled point (ERF's inputs couple only points); 0 when no point is coupled. */
    unsigned num_coupled_dof () const;

    /**
     * Initialise the system with the coupled positions (m) and velocities (m/s), num_coupled_dof()
     * finite values each in MoorDyn's frame (empty when there are none). With compute_ic MoorDyn
     * solves for the static initial shape; without it the state is expected from load() or
     * deserialize(). MoorDyn runs with the floating-point traps masked (FPTrapMask). Returns an
     * error message, empty on success.
     */
    std::string init (const std::vector<double>& x, const std::vector<double>& xd, bool compute_ic = true);

    /**
     * Start the external kinematics; returns the number of kinematics points MoorDyn takes the
     * fluid velocity at. The stub returns 0 (and err is set) when the input file does not set
     * WaveKin = 1; MoorDyn-C 2.7.1 returns the full count whatever WaveKin is (create() has
     * checked the option). Also 0 with err set when MoorDyn reports an error.
     */
    unsigned external_kinematics_init (std::string& err);
    unsigned num_kinematics_points () const { return m_nkin; }
    /** 3*n coordinates of the kinematics points (m, MoorDyn's frame). */
    std::vector<double> kinematics_points () const;
    /**
     * Set the fluid velocity (m/s) and acceleration (m/s^2), 3*n finite values each, at the
     * kinematics points, valid at time t (s). MoorDyn-C 2.7.1 and the stub ignore t; the values
     * hold until the next call. Aborts, naming the input file, on a non-finite value.
     */
    void set_kinematics (const std::vector<double>& U, const std::vector<double>& Ud, double t);

    /**
     * Advance the system by dt (MoorDyn sub-steps internally with its own dtM): the coupled points
     * start the step at x and move at xd over it (MoorDyn places them at x + xd s, s the time into
     * the step), and f receives the coupled forces at its end. t receives MoorDyn's clock, the
     * time since it was initialised.
     *
     * Units: x m, xd m/s, f N, t and dt s; MoorDyn's frame. f is the force ON each coupled point
     * BY the lines ending there (plus the point's own weight, buoyancy and drag, zero for the
     * massless points ERF writes). MoorDyn runs with the floating-point traps masked (FPTrapMask).
     * Aborts when dt is not positive or a returned force is not finite.
     */
    void step (const std::vector<double>& x, const std::vector<double>& xd, std::vector<double>& f, double& t, double dt);

    /** MoorDyn's internal step dtM (s), the bound on the sub-steps step() takes. */
    double dt () const;
    /** Set MoorDyn's internal step dtM (s, finite and > 0). */
    void set_dt (double dt);

    /** The number of MoorDyn lines. */
    unsigned num_lines () const;
    /** The number of nodes of MoorDyn line `line` (1-based): its segments + 1. */
    unsigned line_num_nodes (unsigned line) const;
    /** The unstretched length of MoorDyn line `line` (m). */
    double line_unstretched_length (unsigned line) const;
    /** The position of node `node` (0-based) of MoorDyn line `line` (1-based) (m, MoorDyn's frame). */
    std::array<double,3> line_node_position (unsigned line, unsigned node) const;
    /** The velocity of node `node` of MoorDyn line `line` (m/s). */
    std::array<double,3> line_node_velocity (unsigned line, unsigned node) const;
    /**
     * The tension at node `node` (N, MoorDyn's frame, internal damping included): at an inner
     * node the mean of the two adjacent segments, at an end node the end segment's. It points
     * from node i towards node i+1.
     */
    std::array<double,3> line_node_tension (unsigned line, unsigned node) const;
    /** The fluid's drag on a node (N): the force the air exerts on the MoorDyn line there. */
    std::array<double,3> line_node_drag (unsigned line, unsigned node) const;
    /**
     * The net force on a node (N): the tensions of the segments either side, the weight, the
     * fluid's drag and the rest MoorDyn adds up to move the node. At an end node on a fixed or
     * coupled point it is the force the MoorDyn line puts on that point, since the point holds the node.
     */
    std::array<double,3> line_node_force (unsigned line, unsigned node) const;
    /** The magnitude of the tension at the MoorDyn line's last node (MoorDyn's fairlead tension, N). */
    double line_end_tension (unsigned line) const;
    /** The largest node tension magnitude (N): over nodes 0..N-1 in MoorDyn-C 2.7.1, over 0..N in the stub. */
    double line_max_tension (unsigned line) const;

    /** The number of MoorDyn points. */
    unsigned num_points () const;
    /** The type of point `point` (1-based): -1 coupled, 0 free, 1 fixed. */
    int point_type (unsigned point) const;
    /** The position of point `point` (m, MoorDyn's frame). */
    std::array<double,3> point_position (unsigned point) const;
    /**
     * The net force on a point (N). MoorDyn computes it for the points it integrates (free) and
     * for the coupled ones; a fixed point reports zero. The pull of a MoorDyn line on a fixed
     * point is line_node_tension(line, 0) where the line starts at the point and
     * -line_node_tension(line, N) where it ends there (N = line_num_nodes(line) - 1): MoorDyn's
     * node tension points from node i towards node i+1.
     */
    std::array<double,3> point_force (unsigned point) const;

    /** Write the whole system state to a file (MoorDyn_Save). */
    void save (const std::string& path);
    /** Restore a state written by save() into a system created from the same input and initialised without an IC solve. */
    void load (const std::string& path);
    /**
     * The whole system state in memory (MoorDyn_Serialize), to redo a step: MoorDyn's integrator
     * (with its clock), points and lines. The fluid kinematics last set by set_kinematics() are
     * not part of it, so set them again before the next step.
     */
    std::vector<std::uint64_t> serialize () const;
    /** Restore a state from serialize() of this system (MoorDyn_Deserialize); aborts when its size differs. */
    void deserialize (const std::vector<std::uint64_t>& data);

private:
    MoorDynSystem (MoorDyn sys, std::string file);
    void check (int rc, const std::string& what) const;
    MoorDynLine line (unsigned l) const;
    MoorDynPoint point (unsigned p) const;

    MoorDyn m_sys;
    std::string m_file;
    unsigned m_nkin{0};
    mutable std::size_t m_serial_words{0};   //!< size of the last serialize(), in 64-bit words; 0 before the first
};

} // namespace erf_moordyn

#endif
