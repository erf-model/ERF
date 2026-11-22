#ifndef ERF_TIMESTEP_UTILS_H_
#define ERF_TIMESTEP_UTILS_H_

#include <ERF_Constants.H>
#include <ERF_TerrainMetrics.H>
#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>
#include <cmath>

/**
 * Pointwise inverse time step (1/dt, before the CFL is applied) for the
 * compressible solver.
 *
 * The terrain metric terms are computed here from z_nd, so the same routine
 * returns either the terrain aware or the terrain unaware estimate:
 *
 *   terrain aware   : h_xi and h_eta are retained and the vertical spacing is
 *                     the local dxinv[2]/h_zeta -- w should be the contravariant
 *                     (Omega) vertical velocity
 *   terrain unaware : h_xi = h_eta = 0 and the vertical spacing is dxinv[2]
 *                     -- w should be the Cartesian vertical velocity
 *
 * Note that dxinv[2]*(h_xi/h_zeta) == h_xi*idz, which is why h_zeta only ever
 * enters through idz.
 *
 * The physical horizontal cell widths are dx/mf_x and dy/mf_y, so every
 * horizontal term (advective and acoustic, including the h_xi and h_eta slope
 * contributions) is scaled by the map factor in that direction.  The map
 * factors do not touch the vertical.
 *
 * @param[in] i               x-index
 * @param[in] j               y-index
 * @param[in] k               z-index
 * @param[in] c               sound speed
 * @param[in] u               cell-centered x-velocity
 * @param[in] v               cell-centered y-velocity
 * @param[in] w               cell-centered vertical velocity (Omega if terrain aware)
 * @param[in] z_nd            nodal physical height field (not accessed if terrain unaware)
 * @param[in] mf_x            x map factor at the cell center
 * @param[in] mf_y            y map factor at the cell center
 * @param[in] dxinv           inverse cell spacing of the computational mesh
 * @param[in] l_terrain_aware true to include the metric terms
 * @param[in] l_substepping   true if implicit acoustic substepping is used
 * @param[in] nxc             number of cells in x on the level 0 domain
 * @param[in] nyc             number of cells in y on the level 0 domain
 * @return                    1/dt for this cell
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_InvDt_Compressible (const int i, const int j, const int k,
                            const amrex::Real c,
                            const amrex::Real u,
                            const amrex::Real v,
                            const amrex::Real w,
                            const amrex::Array4<const amrex::Real>& z_nd,
                            const amrex::Real mf_x,
                            const amrex::Real mf_y,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                            const bool l_terrain_aware,
                            const bool l_substepping,
                            const int  nxc=2,
                            const int  nyc=2)
{
    amrex::Real h_xi  = zero;
    amrex::Real h_eta = zero;
    amrex::Real idz   = dxinv[2];

    // Inverse physical horizontal spacing
    const amrex::Real idx = mf_x * dxinv[0];
    const amrex::Real idy = mf_y * dxinv[1];

    if (l_terrain_aware) {
        h_xi  =   Compute_h_xi_AtCellCenter(i,j,k,dxinv,z_nd);
        h_eta =  Compute_h_eta_AtCellCenter(i,j,k,dxinv,z_nd);
        idz  /= Compute_h_zeta_AtCellCenter(i,j,k,dxinv,z_nd);
    }

    const amrex::Real abs_u = amrex::Math::abs(u);
    const amrex::Real abs_v = amrex::Math::abs(v);
    const amrex::Real abs_w = amrex::Math::abs(w);

    amrex::Real inv_dt = -bogus_large_value;

    // If we are doing implicit acoustic substepping, then the z-direction is not constrained
    //    by the speed of sound for the computation of the time step
    if (l_substepping) {
        // Implicit substepping removes the g33 metric
        // Acoustic contribution limits dtau, whose max is in the first RK stage (1/3 dt)
        const amrex::Real Uacoustic = third * c * idx;
        const amrex::Real Vacoustic = third * c * idy;
        if ((nxc > 1) && (nyc == 1)) {
            // 2-D in x-z
            inv_dt = abs_u*idx + Uacoustic
                   + abs_w*idz;
        } else if ((nyc > 1) && (nxc == 1)) {
            // 2-D in y-z
            inv_dt = abs_v*idy + Vacoustic
                   + abs_w*idz;
        } else {
            // 3-D
            inv_dt = abs_u*idx + Uacoustic
                   + abs_v*idy + Vacoustic
                   + abs_w*idz;
        }

    // If we are not doing implicit acoustic substepping, then the z-direction is constrained
    //    by the speed of sound for the computation of the time step
    } else {
        // All the metric terms
        const amrex::Real Uacoustic = c * amrex::Math::abs(mf_x * (dxinv[0] + idz * h_xi ));
        const amrex::Real Vacoustic = c * amrex::Math::abs(mf_y * (dxinv[1] + idz * h_eta));
        const amrex::Real Oacoustic = c * amrex::Math::abs(idz);
        if (nxc > 1 && nyc > 1) {
            // 3-D
            inv_dt = abs_u*idx + Uacoustic
                   + abs_v*idy + Vacoustic
                   + abs_w*idz      + Oacoustic;
        } else if (nxc > 1) {
            // 2-D in x-z
            inv_dt = abs_u*idx + Uacoustic
                   + abs_w*idz      + Oacoustic;
        } else if (nyc > 1) {
            // 2-D in y-z
            inv_dt = abs_v*idy + Vacoustic
                   + abs_w*idz      + Oacoustic;
        } else {
            // 1-D in z
            inv_dt = abs_w*idz + Oacoustic;
        }
    }

    return inv_dt;
}

/**
 * Pointwise inverse time step (1/dt, before the CFL is applied) for the
 * anelastic (low Mach) solver -- this is a purely advective constraint.
 *
 * As above, the metric terms are computed here from z_nd:
 *
 *   terrain aware   : the vertical spacing is the local dxinv[2]/h_zeta and w
 *                     should be the contravariant (Omega) vertical velocity
 *   terrain unaware : the vertical spacing is dxinv[2] and w should be the
 *                     Cartesian vertical velocity
 *
 * The slope metrics h_xi and h_eta do not appear explicitly: for the advective
 * constraint their effect is already carried by Omega, which is the vertical
 * velocity relative to the terrain-following surfaces.
 *
 * As in the compressible estimate, the horizontal terms use the physical cell
 * widths dx/mf_x and dy/mf_y, and a lateral direction with only one cell on the
 * level 0 domain does not contribute, so that dx or dy there sets nothing.
 *
 * @param[in] i               x-index
 * @param[in] j               y-index
 * @param[in] k               z-index
 * @param[in] u               cell-centered x-velocity
 * @param[in] v               cell-centered y-velocity
 * @param[in] w               cell-centered vertical velocity (Omega if terrain aware)
 * @param[in] z_nd            nodal physical height field (not accessed if terrain unaware)
 * @param[in] mf_x            x map factor at the cell center
 * @param[in] mf_y            y map factor at the cell center
 * @param[in] dxinv           inverse cell spacing of the computational mesh
 * @param[in] l_terrain_aware true to include the metric terms
 * @param[in] nxc             number of cells in x on the level 0 domain
 * @param[in] nyc             number of cells in y on the level 0 domain
 * @return                    1/dt for this cell
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_InvDt_Anelastic (const int i, const int j, const int k,
                         const amrex::Real u,
                         const amrex::Real v,
                         const amrex::Real w,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::Real mf_x,
                         const amrex::Real mf_y,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                         const bool l_terrain_aware,
                         const int  nxc=2,
                         const int  nyc=2)
{
    amrex::Real idz = dxinv[2];
    if (l_terrain_aware) {
        idz /= Compute_h_zeta_AtCellCenter(i,j,k,dxinv,z_nd);
    }

    const amrex::Real idx = (nxc > 1) ? mf_x*dxinv[0] : zero;
    const amrex::Real idy = (nyc > 1) ? mf_y*dxinv[1] : zero;

    return ( amrex::Math::abs(u)*idx
           + amrex::Math::abs(v)*idy
           + amrex::Math::abs(w)*idz );
}

/**
 * Horizontal acoustic part of the inverse time step of a cell: the sound speed over the
 * physical horizontal spacing, summed over the horizontal directions the run has.  The
 * product of this with the acoustic substep is the horizontal acoustic Courant number of
 * the substep.  It holds no velocity, so it depends on the grid and the sound speed only.
 *
 * @param[in] c      sound speed (m/s)
 * @param[in] mf_x   map factor in x
 * @param[in] mf_y   map factor in y
 * @param[in] dxinv  inverse computational cell size (1/m)
 * @param[in] nxc    number of cells in x on this level
 * @param[in] nyc    number of cells in y on this level
 * @return           c (m_x/dx + m_y/dy), with only its own direction in a two-dimensional run
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_InvDt_HorizAcoustic (const amrex::Real c,
                             const amrex::Real mf_x,
                             const amrex::Real mf_y,
                             const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                             const int nxc,
                             const int nyc)
{
    const amrex::Real idx = mf_x * dxinv[0];
    const amrex::Real idy = mf_y * dxinv[1];
    if ((nxc > 1) && (nyc == 1)) {
        return c * idx;
    } else if ((nyc > 1) && (nxc == 1)) {
        return c * idy;
    }
    return c * (idx + idy);
}

/**
 * Divergence damping coefficient of an acoustic substep.
 *
 * The substep projects (rho theta) forward by beta times its last change before taking its
 * horizontal gradient.  That damps the horizontal divergence with a diffusivity
 * beta c^2 dtau, which falls with the substep.  So that a run with more substeps is not
 * damped less per unit time, erf.beta_d is the coefficient at a Courant number of one and the
 * coefficient used is erf.beta_d divided by the horizontal acoustic Courant number of the
 * substep,
 *
 *     Cr = dtau inv_dt,
 *
 * with inv_dt the largest value on the level of Compute_InvDt_HorizAcoustic.  The
 * diffusivity is then beta_d c^2 / inv_dt, which does not depend on the number of substeps
 * or on the flow.
 *
 * @param[in] beta_d  erf.beta_d
 * @param[in] dtau    acoustic substep (s)
 * @param[in] inv_dt  largest horizontal acoustic inverse time step on the level (1/s);
 *                    zero where there is none, and beta_d is then returned
 * @return            the coefficient the substep uses
 */
AMREX_FORCE_INLINE
amrex::Real
DivergenceDampingCoeff (const amrex::Real beta_d, const amrex::Real dtau, const amrex::Real inv_dt)
{
    const amrex::Real Cr = dtau * inv_dt;
    return (Cr > zero) ? beta_d / Cr : beta_d;
}

#endif
