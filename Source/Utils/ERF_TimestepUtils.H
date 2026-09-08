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
 * @param[in] i               x-index
 * @param[in] j               y-index
 * @param[in] k               z-index
 * @param[in] c               sound speed
 * @param[in] u               cell-centered x-velocity
 * @param[in] v               cell-centered y-velocity
 * @param[in] w               cell-centered vertical velocity (Omega if terrain aware)
 * @param[in] z_nd            nodal physical height field (not accessed if terrain unaware)
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
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                            const bool l_terrain_aware,
                            const bool l_substepping,
                            const int  nxc=2,
                            const int  nyc=2)
{
    amrex::Real h_xi  = zero;
    amrex::Real h_eta = zero;
    amrex::Real idz   = dxinv[2];
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
        const amrex::Real Uacoustic = third * c * std::sqrt( dxinv[0] * (dxinv[0] + two * idz * amrex::Math::abs(h_xi )) );
        const amrex::Real Vacoustic = third * c * std::sqrt( dxinv[1] * (dxinv[1] + two * idz * amrex::Math::abs(h_eta)) );
        if ((nxc > 1) && (nyc == 1)) {
            // 2-D in x-z
            inv_dt = amrex::max(abs_u*dxinv[0] + Uacoustic,
                                abs_w*idz                 );
        } else if ((nyc > 1) && (nxc == 1)) {
            // 2-D in y-z
            inv_dt = amrex::max(abs_v*dxinv[1] + Vacoustic,
                                abs_w*idz                 );
        } else {
            // 3-D
            inv_dt = amrex::max(abs_u*dxinv[0] + Uacoustic,
                                abs_v*dxinv[1] + Vacoustic,
                                abs_w*idz                 );
        }

    // If we are not doing implicit acoustic substepping, then the z-direction is constrained
    //    by the speed of sound for the computation of the time step
    } else {
        // All the metric terms
        const amrex::Real Uacoustic = c * amrex::Math::abs(dxinv[0] + idz * h_xi );
        const amrex::Real Vacoustic = c * amrex::Math::abs(dxinv[1] + idz * h_eta);
        const amrex::Real Oacoustic = c * amrex::Math::abs(idz);
        if (nxc > 1 && nyc > 1) {
            // 3-D
            inv_dt = amrex::max(abs_u*dxinv[0] + Uacoustic,
                                abs_v*dxinv[1] + Vacoustic,
                                abs_w*idz      + Oacoustic);
        } else if (nxc > 1) {
            // 2-D in x-z
            inv_dt = amrex::max(abs_u*dxinv[0] + Uacoustic,
                                abs_w*idz      + Oacoustic);
        } else if (nyc > 1) {
            // 2-D in y-z
            inv_dt = amrex::max(abs_v*dxinv[1] + Vacoustic,
                                abs_w*idz      + Oacoustic);
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
 * @param[in] i               x-index
 * @param[in] j               y-index
 * @param[in] k               z-index
 * @param[in] u               cell-centered x-velocity
 * @param[in] v               cell-centered y-velocity
 * @param[in] w               cell-centered vertical velocity (Omega if terrain aware)
 * @param[in] z_nd            nodal physical height field (not accessed if terrain unaware)
 * @param[in] dxinv           inverse cell spacing of the computational mesh
 * @param[in] l_terrain_aware true to include the metric terms
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
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv,
                         const bool l_terrain_aware)
{
    amrex::Real idz = dxinv[2];
    if (l_terrain_aware) {
        idz /= Compute_h_zeta_AtCellCenter(i,j,k,dxinv,z_nd);
    }

    return amrex::max(amrex::Math::abs(u)*dxinv[0],
                      amrex::Math::abs(v)*dxinv[1],
                      amrex::Math::abs(w)*idz     );
}

#endif
