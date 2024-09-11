#include <ERF_Advection.H>

using namespace amrex;

/** Compute advection tendencies for momentum normal to BC */
void
AdvectionSrcForOpenBC_Normal (const Box& bx,
                              const int& dir,
                              const Array4<      Real>& rhs_arr,
                              const Array4<const Real>& vel_norm_arr,
                              const Array4<const Real>& cell_data_arr,
                              const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                              const bool do_lo)
{
    // NOTE: Klemp, J. B., and R. Wilhelmson, 1978: The simulation of three-dimensional
    //       convective storm dynamics, J. Atmos. Sci., 35, 1070-1096.
    // NOTE: Implementation is for the high bndry side. The low bndry side is obtained
    //       by flipping sgn = -1.
    // NOTE: Indices (i,j,k) correspond to data that is ON the open bdy.
    int sgn = 1; if (do_lo) sgn = -1;
    Real c_o_star = Real(sgn)*30.0;
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        IntVect ivu1(i,j,k); if ( do_lo) ivu1[dir] -= sgn; // Vel indexed into domain for do_lo
        IntVect ivu2(i,j,k); if (!do_lo) ivu2[dir] -= sgn; // Vel indexed into domain for do_hi

        IntVect ivr1(i,j,k); if (!do_lo) ivr1[dir] -= sgn; // Rho indexed into domain for do_hi
        IntVect ivr2(i,j,k); if ( do_lo) ivr2[dir] += sgn; // Rho indexed out  domain for do_lo

        Real rho_face  = 0.5 * ( cell_data_arr(ivr1,Rho_comp) + cell_data_arr(ivr2,Rho_comp) );
        Real mom_star  = rho_face * Real(sgn) * max( Real(sgn)*(vel_norm_arr(ivu1) + c_o_star), 0.0 );
        Real vel_grad  =  ( vel_norm_arr(ivu1) - vel_norm_arr(ivu2) ) * dxInv[dir];
        Real flux      = -( mom_star * vel_grad );
        rhs_arr(i,j,k) = flux;
    });
}

/** Compute advection tendencies for x momentum tangential to BC (2nd order)*/
void
AdvectionSrcForOpenBC_Tangent_Xmom (const Box& bxx,
                                    const int& dir,
                                    const Array4<      Real>& rho_u_rhs,
                                    const Array4<const Real>& u,
                                    const Array4<const Real>& rho_u,
                                    const Array4<const Real>& rho_v,
                                    const Array4<const Real>& Omega,
                                    const Array4<const Real>& ax,
                                    const Array4<const Real>& az,
                                    const Array4<const Real>& detJ,
                                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo)
{
    AMREX_ALWAYS_ASSERT(dir==1);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    ParallelFor(bxx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (ax(i,j,k) > 0.) {
            Real xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (u(i+1,j,k) + u(i,j,k)) * 0.5 * (ax(i,j,k) + ax(i+1,j,k));
            Real xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (u(i-1,j,k) + u(i,j,k)) * 0.5 * (ax(i,j,k) + ax(i-1,j,k));

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k)) * az(i,j,k+1);
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k)) * az(i,j,k  );

            Real x_src = (xflux_hi - xflux_lo) * dxInv;
            Real y_src = AdvectionSrcForOpenBC_Tangent(i, j, k, 0, dir, u, rho_v, dyInv, do_lo);
            Real z_src = (zflux_hi - zflux_lo) * dzInv;
            Real advectionSrc = x_src + y_src + z_src;

            rho_u_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i-1,j,k)));
        } else {
            rho_u_rhs(i, j, k) = 0.0;
        }
    });
}

/** Compute advection tendencies for x momentum tangential to BC (2nd order)*/
void
AdvectionSrcForOpenBC_Tangent_Ymom (const Box& bxy,
                                    const int& dir,
                                    const Array4<      Real>& rho_v_rhs,
                                    const Array4<const Real>& v,
                                    const Array4<const Real>& rho_u,
                                    const Array4<const Real>& rho_v,
                                    const Array4<const Real>& Omega,
                                    const Array4<const Real>& ay,
                                    const Array4<const Real>& az,
                                    const Array4<const Real>& detJ,
                                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo)
{
    AMREX_ALWAYS_ASSERT(dir==0);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    ParallelFor(bxy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (ay(i,j,k) > 0.) {
            Real yflux_hi = 0.25 * (rho_v(i  ,j+1,k) + rho_v(i  ,j  ,k)) * (v(i,j+1,k) + v(i,j,k)) * 0.5 * (ay(i,j,k) + ay(i,j+1,k));
            Real yflux_lo = 0.25 * (rho_v(i  ,j  ,k) + rho_v(i  ,j-1,k)) * (v(i,j-1,k) + v(i,j,k)) * 0.5 * (ay(i,j,k) + ay(i,j-1,k));

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k)) * az(i,j,k+1);
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k)) * az(i,j,k  );

            Real x_src = AdvectionSrcForOpenBC_Tangent(i, j, k, 0, dir, v, rho_u, dxInv, do_lo);
            Real y_src = (yflux_hi - yflux_lo) * dyInv;
            Real z_src = (zflux_hi - zflux_lo) * dzInv;
            Real advectionSrc = x_src + y_src + z_src;

            rho_v_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i,j-1,k)));
        } else {
            rho_v_rhs(i, j, k) = 0.0;
        }
    });
}

/** Compute advection tendencies for x momentum tangential to BC (2nd order)*/
void
AdvectionSrcForOpenBC_Tangent_Zmom (const Box& bxz,
                                    const int& dir,
                                    const Array4<      Real>& rho_w_rhs,
                                    const Array4<const Real>& w,
                                    const Array4<const Real>& rho_u,
                                    const Array4<const Real>& rho_v,
                                    const Array4<const Real>& Omega,
                                    const Array4<const Real>& ax,
                                    const Array4<const Real>& ay,
                                    const Array4<const Real>& az,
                                    const Array4<const Real>& detJ,
                                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const int domhi_z,
                                    const bool do_lo)
{
    AMREX_ALWAYS_ASSERT(dir!=2);

    bool xopen = (dir==0);
    bool yopen = (dir==1);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (az(i,j,k) > 0.) {
            Real xflux_hi = 0.25*(rho_u(i+1,j  ,k) + rho_u(i+1, j, k-1)) * (w(i+1,j,k) + w(i,j,k)) *
                            0.5 *(ax(i+1,j,k) + ax(i+1,j,k-1));

            Real xflux_lo = 0.25*(rho_u(i  ,j  ,k) + rho_u(i  , j, k-1)) * (w(i-1,j,k) + w(i,j,k)) *
                            0.5 *(ax(i,j,k) + ax(i,j,k-1));

            Real yflux_hi = 0.25*(rho_v(i  ,j+1,k) + rho_v(i, j+1, k-1)) * (w(i,j+1,k) + w(i,j,k)) *
                            0.5 *(ay(i,j+1,k) + ay(i,j+1,k-1));

            Real yflux_lo = 0.25*(rho_v(i  ,j  ,k) + rho_v(i, j  , k-1)) * (w(i,j-1,k) + w(i,j,k)) *
                            0.5 *(ay(i,j,k) + ay(i,j,k-1));

            Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1)) *
                            0.5 *(az(i,j,k) + az(i,j,k-1));

            Real zflux_hi = (k == domhi_z+1) ? Omega(i,j,k) * w(i,j,k) * az(i,j,k):
                                               0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1)) *
                                               0.5  * (az(i,j,k) + az(i,j,k-1));

            Real x_src = (xopen) ? AdvectionSrcForOpenBC_Tangent(i, j, k, 0, dir, w, rho_u, dxInv, do_lo) :
                                   (xflux_hi - xflux_lo) * dxInv;
            Real y_src = (yopen) ? AdvectionSrcForOpenBC_Tangent(i, j, k, 0, dir, w, rho_v, dyInv, do_lo) :
                                   (yflux_hi - yflux_lo) * dyInv;
            Real z_src = (zflux_hi - zflux_lo) * dzInv;
            Real advectionSrc = x_src + y_src + z_src;

            rho_w_rhs(i, j, k) = -advectionSrc / (0.5*(detJ(i,j,k) + detJ(i,j,k-1)));
        } else {
            rho_w_rhs(i, j, k) = 0.0;
        }
    });
}

/** Compute advection tendencies for x momentum tangential to BC (2nd order)*/
void
AdvectionSrcForOpenBC_Tangent_Cons (const Box& bx,
                                    const int& dir,
                                    const int& icomp,
                                    const int& ncomp,
                                    const Array4<      Real>& cell_rhs,
                                    const Array4<const Real>& cell_prim,
                                    const Array4<const Real>& avg_xmom,
                                    const Array4<const Real>& avg_ymom,
                                    const Array4<const Real>& avg_zmom,
                                    const Array4<const Real>& detJ,
                                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                                    const bool do_lo)
{
    AMREX_ALWAYS_ASSERT(dir!=2 && icomp>0);

    bool xopen = (dir==0);
    bool yopen = (dir==1);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        if (detJ(i,j,k) > 0.) {
            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            Real prim_xlo = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index));
            Real prim_xhi = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index));
            Real xflux_lo = avg_xmom(i  ,j,k) * prim_xlo;
            Real xflux_hi = avg_xmom(i+1,j,k) * prim_xhi;

            Real prim_ylo = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index));
            Real prim_yhi = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index));
            Real yflux_lo = avg_ymom(i,j  ,k) * prim_ylo;
            Real yflux_hi = avg_ymom(i,j+1,k) * prim_yhi;

            Real prim_zlo = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index));
            Real prim_zhi = 0.5 * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index));
            Real zflux_lo = avg_zmom(i,j,k  ) * prim_zlo;
            Real zflux_hi = avg_zmom(i,j,k+1) * prim_zhi;

            Real x_src = (xopen) ? AdvectionSrcForOpenBC_Tangent(i, j, k, prim_index, dir, cell_prim, avg_xmom, dxInv, do_lo) :
                                              (xflux_hi - xflux_lo) * dxInv;
            Real y_src = (yopen) ? AdvectionSrcForOpenBC_Tangent(i, j, k, prim_index, dir, cell_prim, avg_ymom, dyInv, do_lo) :
                                              (yflux_hi - yflux_lo) * dyInv;
            Real z_src = (zflux_hi - zflux_lo) * dzInv;
            Real advectionSrc = x_src + y_src + z_src;
            cell_rhs(i,j,k,cons_index) = -advectionSrc / detJ(i,j,k);
        }
    });
}


AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
Real
AdvectionSrcForOpenBC_Tangent (const int& i,
                               const int& j,
                               const int& k,
                               const int& nprim,
                               const int& dir,
                               const Array4<const Real>& prim_tang_arr,
                               const Array4<const Real>& mom_norm_arr,
                               const Real& dxInv,
                               const bool do_lo)
{
    // NOTE: Klemp, J. B., and R. Wilhelmson, 1978: The simulation of three-dimensional
    //       convective storm dynamics, J. Atmos. Sci., 35, 1070-1096.
    // NOTE: Implementation is for the high bndry side. The low bndry side is obtained
    //       by flipping sgn = -1.
    // NOTE: Indices (i,j,k) correspond to data that is index 1/2 dx off open bdy.
    int sgn = 1; if (do_lo) sgn = -1;

    IntVect ivm1(i,j,k); if ( do_lo) ivm1[dir] -= sgn; // Mom indexed into domain for do_lo
    IntVect ivm2(i,j,k); if (!do_lo) ivm1[dir] += sgn; // Mom indexed out  domain for do_hi

    IntVect ivs1(i,j,k); if ( do_lo) ivs1[dir] -= sgn; // Scalar indexed into domain for do_hi
    IntVect ivs2(i,j,k); if (!do_lo) ivs2[dir] -= sgn; // Scalar indexed into domain for do_lo

    Real mom_at_cc = 0.5 * (mom_norm_arr(ivm1) + mom_norm_arr(ivm2));
    Real mom_star  =    Real(sgn) * max( Real(sgn)*mom_at_cc, 0.0 );
    Real mom_grad  = ( mom_norm_arr(ivm1) - mom_norm_arr(ivm2) ) * dxInv;
    Real prim_grad = ( prim_tang_arr(ivs1,nprim) - prim_tang_arr(ivs2,nprim) ) * dxInv;

    // NOTE: Negative sign taken care of by wrapper function
    Real src       = ( mom_star*prim_grad + prim_tang_arr(ivs1,nprim)*mom_grad );
    return src;
}
