#ifndef ERF_ACTUATOR_SPREADING_H
#define ERF_ACTUATOR_SPREADING_H

/** Spreading of point forces onto the face-centred momentum sources with a normalised Gaussian kernel, and the integral check. */

#include <vector>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace erf_actuator {

/**
 * Spread point forces onto the face-centred momentum sources with a 3-D Gaussian kernel.
 *
 * Each point p carries a force F_p on the fluid (N) and has its own kernel width epsilon_p
 * (the conductor lines and the towers of a run may differ). Its share of the source at a face is
 * F_p w_p(face) / S_p, where w_p = exp(-|x_face - x_p|^2 / epsilon_p^2) is cut off beyond
 * 3 epsilon_p and S_p = sum over the faces of that component of w_p dV is the discrete
 * normalisation with the face volumes (dx dy dz detJ on a terrain-following mesh). The sum
 * of src dV over the mesh therefore equals F_p exactly for every point and component,
 * whatever the resolution and however much of the kernel is cut off by the ground or the
 * domain top; the kernel's shape is Gaussian only where it is resolved. The faces sit at
 * prob_lo + i dx with the volumes dx dy dz detJ, with no map factors: a body model must refuse a
 * run whose map factors are not 1 until the core carries them. The normalisation sums are
 * accumulated with atomic adds, so on a GPU the last bits of the source can differ between runs;
 * CPU runs are bitwise reproducible.
 *
 * Every rank accumulates the faces its boxes own (each physical face once, however the boxes
 * share their boundary faces) and the normalisations are summed across ranks, so the result
 * does not depend on the domain decomposition. Every periodic image of a point within
 * 3 epsilon_p of the domain (in x and y; z is not imaged) adds its kernel, so a kernel wraps across
 * a periodic boundary, and overlaps itself when the period is under 6 epsilon_p. The sources are
 * set (not added to) on every valid face, both copies of a shared face alike; the faces on a
 * non-periodic domain boundary get no source.
 *
 * Aborts when a position or force is not finite. z_phys_nd and detJ_cc must share the boxes and
 * distribution of src_x; src_x, src_y and src_z are face-centred in x, y and z.
 *
 * @param[in] pos 3*n point coordinates (m, ERF domain coordinates)
 * @param[in] force 3*n forces ON the fluid BY the body (N, ERF frame): minus the fluid's force on the body
 * @param[in] epsilon n kernel widths, one per point (m, > 0 and finite)
 * @param[in] z_phys_nd nodal heights (m), or nullptr on a uniform-dz mesh
 * @param[in] detJ_cc cell Jacobians (1 ghost cell), or nullptr on a uniform-dz mesh
 * @param[in] geom the level's geometry
 * @param[out] src_x,src_y,src_z face-centred momentum sources (N/m^3), one component each
 */
void spread_forces (const std::vector<amrex::Real>& pos,
                    const std::vector<amrex::Real>& force,
                    const std::vector<amrex::Real>& epsilon,
                    const amrex::MultiFab* z_phys_nd,
                    const amrex::MultiFab* detJ_cc,
                    const amrex::Geometry& geom,
                    amrex::MultiFab& src_x,
                    amrex::MultiFab& src_y,
                    amrex::MultiFab& src_z);

/** The same with one kernel width epsilon (m, > 0 and finite) for every point. */
void spread_forces (const std::vector<amrex::Real>& pos,
                    const std::vector<amrex::Real>& force,
                    amrex::Real epsilon,
                    const amrex::MultiFab* z_phys_nd,
                    const amrex::MultiFab* detJ_cc,
                    const amrex::Geometry& geom,
                    amrex::MultiFab& src_x,
                    amrex::MultiFab& src_y,
                    amrex::MultiFab& src_z);

/**
 * Return the integral of a face-centred source over the mesh, sum of src dV (the total force, N).
 *
 * @param dir 0, 1 or 2: src holds the x-, y- or z-face component
 */
amrex::Real integrate_source (int dir,
                              const amrex::MultiFab& src,
                              const amrex::MultiFab* detJ_cc,
                              const amrex::Geometry& geom);

} // namespace erf_actuator

#endif
