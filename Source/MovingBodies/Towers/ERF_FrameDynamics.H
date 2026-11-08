#ifndef ERF_FRAME_DYNAMICS_H
#define ERF_FRAME_DYNAMICS_H

/**
 * The dynamics of a frame model: its natural modes by subspace iteration, checked by a Sturm
 * sequence count, and its time response by Newmark's average-acceleration method with Rayleigh
 * damping. Double precision; the frame's conventions (ERF_Frame.H) hold: 6 degrees of freedom per
 * node, frame axes, SI units.
 */

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

#include "ERF_DenseCholesky.H"
#include "ERF_Frame.H"

namespace erf_towers {

/** The lowest natural modes of a frame. */
struct FrameModes
{
    std::vector<double> frequency;              //!< undamped natural frequencies (Hz), ascending
    /** Per mode: one value per degree of freedom (num_dofs(); fixed ones 0), mass-normalised: phi^T M phi = 1. */
    std::vector<std::vector<double>> shape;
};

/**
 * The count lowest modes of K phi = omega^2 M phi on the frame's free degrees of freedom, by
 * subspace iteration on the factored stiffness, converged to a relative change of 1e-12 in every
 * eigenvalue. A Sturm sequence count of K - sigma M, sigma between the last mode returned and the
 * next, proves that no lower mode was missed. Returns an empty string on success, else the reason
 * (no convergence, a missed mode, a singular mass).
 */
std::string frame_modes (const Frame& f, std::size_t count, FrameModes& modes);

/**
 * The Rayleigh coefficients a0 (1/s) and a1 (s) of the damping C = a0 M + a1 K that give the
 * damping ratios zeta1 at frequency f1 and zeta2 at frequency f2 (Hz, f1 != f2). A mode of
 * angular frequency omega then has the ratio (a0/omega + a1 omega)/2.
 */
void rayleigh_coefficients (double f1, double zeta1, double f2, double zeta2, double& a0, double& a1);

/**
 * The time response of a frame by Newmark's average-acceleration method (beta = 1/4, gamma = 1/2):
 * unconditionally stable, second-order accurate, and without numerical damping, so an undamped free
 * vibration keeps its energy exactly; each mode advances at the angular frequency
 * (2/h) atan(omega h/2) per step h. The damping is Rayleigh's, C = a0 M + a1 K. The loads (forces on
 * the nodes, gravity) are those at the end of each step. The effective stiffness
 * K + (2/h) C + (4/h^2) M is factored once per step size, and shared by every FrameDynamics of the same
 * frame, Rayleigh coefficients and step (the towers of a type share one frame).
 *
 * State vectors hold one value per degree of freedom (num_dofs(); fixed ones 0), frame axes:
 * displacement (m, rad), velocity (m/s, rad/s), acceleration (m/s^2, rad/s^2).
 */
class FrameDynamics
{
public:
    /** The frame must outlive this object. a0 (1/s) and a1 (s) are the Rayleigh coefficients, both >= 0. */
    FrameDynamics (const Frame& frame, double a0, double a1);

    /** Start at rest in static equilibrium under node_loads (empty: none) and gravity (m/s^2), at time t (s). */
    void start_static (const std::vector<double>& node_loads, double gravity, double t = 0.0);
    /**
     * Start from a given displacement and velocity under node_loads and gravity at time t (s); the
     * acceleration follows from M a = f - C v - K u. Returns an empty string, or the reason (sizes,
     * non-finite values, a singular mass).
     */
    std::string start (const std::vector<double>& u, const std::vector<double>& v, const std::vector<double>& node_loads,
                       double gravity, double t = 0.0);
    /** Advance by h (s, > 0) to the loads node_loads and gravity at the end of the step. */
    void step (double h, const std::vector<double>& node_loads, double gravity);
    /**
     * Advance by h (s, > 0) under node_loads and gravity held over the step: the implicit midpoint
     * rule, (u' - u)/h = (v + v')/2 and M (v' - v)/h + C (v + v')/2 + K (u + u')/2 = f, which is
     * Newmark's method started from the acceleration M a = f - C v - K u under the held loads. The
     * work f (u' - u) of the held loads then equals the change of the kinetic and strain energy plus
     * the energy the damping takes, exactly; a partner that moves the frame's loaded points at the
     * constant velocity (u' - u)/h and returns its mean force over the step exchanges energy with the
     * frame without loss or gain. Uses the step's K_eff, and the mass matrix factored once per frame
     * for the acceleration at the end of the step.
     */
    void step_held (double h, const std::vector<double>& node_loads, double gravity);

    double time () const { return m_t; }
    const std::vector<double>& displacement () const { return m_u; }
    const std::vector<double>& velocity () const { return m_v; }
    const std::vector<double>& acceleration () const { return m_a; }
    /** The factorisation of K_eff the last step used (null before the first): shared by the dynamics of one frame and step. */
    const DenseCholesky* effective_stiffness_factor () const { return m_keff.get(); }
    /** The kinetic energy v^T M v / 2 (J). */
    double kinetic_energy () const;
    /** The strain energy u^T K u / 2 (J), the support springs included. */
    double strain_energy () const;

    /** The state to checkpoint: the time, then u, v and a (1 + 3 num_dofs() values). */
    std::vector<double> state () const;
    /** Restore a state from state(): returns false, changing nothing, when its size does not fit the frame or a value is not finite. */
    bool set_state (const std::vector<double>& s);

private:
    std::vector<double> to_free (const std::vector<double>& full) const;
    std::vector<double> to_full (const std::vector<double>& free) const;

    const Frame& m_frame;
    double m_a0{0.0};
    double m_a1{0.0};
    double m_t{0.0};
    std::vector<double> m_u, m_v, m_a;
    double m_h{0.0};                  //!< the step the effective stiffness is factored for (s); 0: none yet
    std::shared_ptr<const DenseCholesky> m_keff;   //!< K_eff factored for m_h, shared with the dynamics of the same frame
    std::shared_ptr<const DenseCholesky> m_mass;   //!< the free mass matrix factored (step_held), shared with the dynamics of the same frame
};

} // namespace erf_towers

#endif
