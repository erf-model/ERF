#ifndef ERF_LOG_SHEAR_H_
#define ERF_LOG_SHEAR_H_

#include <cmath>

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include "ERF_NumericalConstants.H"

/**
 * erf.rans_shear_stencil = log_law (k-eqn RANS, issue 5261): the vertical gradients of u and v in
 * the strain invariant S_mn S_mn, exact for a logarithmic profile over the ground.
 *
 * ComputeSmnSmn averages the strains of the faces below and above a cell, which for du/dz is the
 * central difference (u_{k+1} - u_{k-1}) / (z_{k+1} - z_{k-1}) on uniform layers. For the profile
 * u = (u_* / kappa) ln((d + z0) / z0) that difference is larger than du/dz at the cell center by
 * the factor (d_k + z0) ln((d_{k+1} + z0) / (d_{k-1} + z0)) / (d_{k+1} - d_{k-1}): 1.21 in the
 * second cell of uniform layers, 1.06 in the third and 1.03 in the fourth (z0 << dz), and
 * 1 + dz^2 / (3 d^2) far from the ground. The shear production of k, 2 nu_t S_mn S_mn, is then
 * 46 %, 12 % and 6 % high there (45 %, 12 % and 6 % with z0 = 0.1 m and dz = 15.625 m).
 *
 * Here, in a cell k above the wall cell and below the top cell,
 *
 *   du/dz = (u_{k+1} - u_{k-1}) / [ (d_k + z0) ln((d_{k+1} + z0) / (d_{k-1} + z0)) ]
 *
 * and the same for v, with u and v averaged to the cell centers, d_k the closure's wall distance,
 * d_{k-1} and d_{k+1} that distance less and plus the height differences of the cell centers, and
 * z0 the surface layer's roughness length of the column. S_mn S_mn changes by the difference of
 * its (1,3) and (2,3) terms; dw/dx, dw/dy and the other components stay as ComputeStrain_* formed
 * them. A cell whose d_{k-1} + z0 would not be positive is left as it is.
 *
 * The stencil is exact for the log profile only. For a profile that is linear in d it returns the
 * gradient divided by the factor above (17 %, 6 % and 3 % low in the cells 2, 3 and 4), so it is
 * less accurate than the central difference where the profile is closer to linear than to
 * logarithmic: in a stable layer above about L / 5 (L the Obukhov length), and well above the
 * ground, where both errors are below 1 %. On sloping ground d is the distance to the wall and
 * the height differences are vertical, so the form is exact to leading order in the slope only.
 */
namespace log_shear {

/**
 * Inverse of the vertical distance over which ComputeStrain_{N,S,T} difference u on the x-z edge
 * (i,j,k), the edge below u(i,j,k): the height difference of the u points k and k-1.
 * It is dz on uniform layers, (dz_k + dz_{k-1}) / 2 on stretched layers and
 * dz / h_zeta (Compute_h_zeta_AtEdgeCenterJ) on a terrain-fitted mesh.
 *
 * @param[in] i x index of the edge
 * @param[in] j y index of the edge
 * @param[in] k z index of the edge
 * @param[in] z_nd nodal physical heights
 * @return inverse spacing [1/m]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
inv_dz_xedge (int i, int j, int k, const amrex::Array4<amrex::Real const>& z_nd) noexcept
{
    return amrex::Real(4) / ( z_nd(i,j,k+1) + z_nd(i,j+1,k+1) - z_nd(i,j,k-1) - z_nd(i,j+1,k-1) );
}

/**
 * As inv_dz_xedge for v on the y-z edge (i,j,k) (Compute_h_zeta_AtEdgeCenterI).
 *
 * @param[in] i x index of the edge
 * @param[in] j y index of the edge
 * @param[in] k z index of the edge
 * @param[in] z_nd nodal physical heights
 * @return inverse spacing [1/m]
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
inv_dz_yedge (int i, int j, int k, const amrex::Array4<amrex::Real const>& z_nd) noexcept
{
    return amrex::Real(4) / ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) - z_nd(i,j,k-1) - z_nd(i+1,j,k-1) );
}

/**
 * S_mn S_mn of cell (i,j,k) with the log-consistent vertical gradients of u and v.
 *
 * The caller keeps k above the wall cell and below the top cell of the domain: the stencil reads
 * the velocities of the cells k-1 and k+1, the edges k and k+1 (none of them a boundary edge)
 * and the nodes k-1 to k+2.
 *
 * @param[in] i cell index in x
 * @param[in] j cell index in y
 * @param[in] k cell index in z
 * @param[in] SmnSmn S_mn S_mn of the cell from ComputeSmnSmn
 * @param[in] u x velocity (x faces)
 * @param[in] v y velocity (y faces)
 * @param[in] s13 (1,3) strain on the x-z edges, as ComputeSmnSmn read it
 * @param[in] s23 (2,3) strain on the y-z edges, as ComputeSmnSmn read it
 * @param[in] d wall distance of the cells [m], read in the cell (i,j,k) only
 * @param[in] z0 roughness length of the column [m]
 * @param[in] z_nd nodal physical heights
 * @return the corrected S_mn S_mn, not below zero
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
corrected_SmnSmn (int i, int j, int k, amrex::Real SmnSmn,
                  const amrex::Array4<amrex::Real const>& u,
                  const amrex::Array4<amrex::Real const>& v,
                  const amrex::Array4<amrex::Real const>& s13,
                  const amrex::Array4<amrex::Real const>& s23,
                  const amrex::Array4<amrex::Real const>& d,
                  amrex::Real z0,
                  const amrex::Array4<amrex::Real const>& z_nd) noexcept
{
    // Heights of the cell centers k-1 and k+1 relative to the center k, from the nodes: the wall
    // distance is read in the cell itself only, so the result does not depend on how the boxes
    // are cut or on ghost cells of the distance.
    const amrex::Real dz_lo = amrex::Real(0.125) * ( z_nd(i,j,k+1) + z_nd(i+1,j,k+1) + z_nd(i,j+1,k+1) + z_nd(i+1,j+1,k+1)
                                                   - z_nd(i,j,k-1) - z_nd(i+1,j,k-1) - z_nd(i,j+1,k-1) - z_nd(i+1,j+1,k-1) );
    const amrex::Real dz_hi = amrex::Real(0.125) * ( z_nd(i,j,k+2) + z_nd(i+1,j,k+2) + z_nd(i,j+1,k+2) + z_nd(i+1,j+1,k+2)
                                                   - z_nd(i,j,k  ) - z_nd(i+1,j,k  ) - z_nd(i,j+1,k  ) - z_nd(i+1,j+1,k  ) );
    const amrex::Real dc = d(i,j,k) + z0;
    const amrex::Real dm = dc - dz_lo;
    const amrex::Real dp = dc + dz_hi;
    // Not where the cell below would lie at or under the roughness height (a wall distance
    // shorter than the height above the ground). The negated form also leaves a cell with a NaN
    // distance as it is.
    if (!(dm > zero && dc > dm && dp > dc)) { return SmnSmn; }
    const amrex::Real inv_lg = one / (dc * std::log(dp / dm));

    // The log-consistent gradients, from the velocities averaged to the cell centers
    const amrex::Real gu_new = myhalf * ( u(i,j,k+1) + u(i+1,j,k+1) - u(i,j,k-1) - u(i+1,j,k-1) ) * inv_lg;
    const amrex::Real gv_new = myhalf * ( v(i,j,k+1) + v(i,j+1,k+1) - v(i,j,k-1) - v(i,j+1,k-1) ) * inv_lg;

    // The gradients the four edges of ComputeSmnSmn hold: du/dz of s13, dv/dz of s23
    const amrex::Real gu_old = fourth * ( (u(i  ,j,k  ) - u(i  ,j,k-1)) * inv_dz_xedge(i  ,j,k  ,z_nd)
                                        + (u(i  ,j,k+1) - u(i  ,j,k  )) * inv_dz_xedge(i  ,j,k+1,z_nd)
                                        + (u(i+1,j,k  ) - u(i+1,j,k-1)) * inv_dz_xedge(i+1,j,k  ,z_nd)
                                        + (u(i+1,j,k+1) - u(i+1,j,k  )) * inv_dz_xedge(i+1,j,k+1,z_nd) );
    const amrex::Real gv_old = fourth * ( (v(i,j  ,k  ) - v(i,j  ,k-1)) * inv_dz_yedge(i,j  ,k  ,z_nd)
                                        + (v(i,j  ,k+1) - v(i,j  ,k  )) * inv_dz_yedge(i,j  ,k+1,z_nd)
                                        + (v(i,j+1,k  ) - v(i,j+1,k-1)) * inv_dz_yedge(i,j+1,k  ,z_nd)
                                        + (v(i,j+1,k+1) - v(i,j+1,k  )) * inv_dz_yedge(i,j+1,k+1,z_nd) );

    // The cell-centered strains of ComputeSmnSmn, and the same with the new gradients
    const amrex::Real s13_old = fourth * ( s13(i,j,k) + s13(i,j,k+1) + s13(i+1,j,k) + s13(i+1,j,k+1) );
    const amrex::Real s23_old = fourth * ( s23(i,j,k) + s23(i,j,k+1) + s23(i,j+1,k) + s23(i,j+1,k+1) );
    const amrex::Real s13_new = s13_old + myhalf * (gu_new - gu_old);
    const amrex::Real s23_new = s23_old + myhalf * (gv_new - gv_old);

    // Not below zero: SmnSmn holds the two old terms to round-off only
    return amrex::max(amrex::Real(zero), SmnSmn + two * (s13_new*s13_new - s13_old*s13_old)
                                   + two * (s23_new*s23_new - s23_old*s23_old));
}

/**
 * Replace S_mn S_mn by corrected_SmnSmn in the cells of bx above the wall cell and below the top
 * cell of the domain. The wall cell keeps the strain of its faces (the lower one is the one-sided
 * difference to the ghost value under the surface layer) and the top cell that of the top boundary.
 *
 * @param[in] bx cells ComputeSmnSmn filled
 * @param[in] domain domain of the level
 * @param[in,out] SmnSmn S_mn S_mn of the cells
 * @param[in] u x velocity the strains were formed from
 * @param[in] v y velocity the strains were formed from
 * @param[in] s13 (1,3) strain on the x-z edges
 * @param[in] s23 (2,3) strain on the y-z edges
 * @param[in] d wall distance of the cells [m]
 * @param[in] z0 roughness length of the surface layer (one layer, k = 0) [m]
 * @param[in] z_nd nodal physical heights
 */
inline void
correct_box (const amrex::Box& bx, const amrex::Box& domain,
             const amrex::Array4<amrex::Real>& SmnSmn,
             const amrex::Array4<amrex::Real const>& u,
             const amrex::Array4<amrex::Real const>& v,
             const amrex::Array4<amrex::Real const>& s13,
             const amrex::Array4<amrex::Real const>& s23,
             const amrex::Array4<amrex::Real const>& d,
             const amrex::Array4<amrex::Real const>& z0,
             const amrex::Array4<amrex::Real const>& z_nd)
{
    amrex::Box cbx = bx;
    cbx.setSmall(2, amrex::max(bx.smallEnd(2), domain.smallEnd(2) + 1));
    cbx.setBig  (2, amrex::min(bx.bigEnd(2)  , domain.bigEnd(2)   - 1));
    if (!cbx.ok()) { return; }
    amrex::ParallelFor(cbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        SmnSmn(i,j,k) = corrected_SmnSmn(i, j, k, SmnSmn(i,j,k), u, v, s13, s23, d, z0(i,j,0), z_nd);
    });
}

} // namespace log_shear

#endif
