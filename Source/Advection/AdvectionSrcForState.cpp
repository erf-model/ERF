#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Advection.H>
#include <Interpolation.H>
#include <Interpolation_WENO.H>

using namespace amrex;

void
AdvectionSrcForRhoAndTheta (const Box& bx, const Box& valid_bx,
                            const Array4<Real>& advectionSrc,
                            const Array4<const Real>& rho_u,
                            const Array4<const Real>& rho_v,
                            const Array4<const Real>& Omega, Real fac,
                            const Array4<      Real>& avg_xmom,
                            const Array4<      Real>& avg_ymom,
                            const Array4<      Real>& avg_zmom,
                            const Array4<const Real>& cell_prim,
                            const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                            const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                            const Array4<const Real>& mf_m,
                            const Array4<const Real>& mf_u,
                            const Array4<const Real>& mf_v,
                            const bool all_use_WENO,
                            const int  spatial_order_WENO,
                            const int horiz_spatial_order,
                            const int vert_spatial_order,
                            const int use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForRhoAndTheta", AdvectionSrcForRhoAndTheta);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // We note that valid_bx is the actual grid, while bx may be a tile within that grid
    const auto& vbx_hi = amrex::ubound(valid_bx);

    if ( use_terrain && (std::max(horiz_spatial_order,vert_spatial_order) == 2) && !all_use_WENO) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real invdetJ = 1./ detJ(i,j,k);

            Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
            Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
            Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
            Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k,cellSizeInv,z_nd);
            xflux_lo *= met_h_zeta_xlo;
            Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i+1,j  ,k,cellSizeInv,z_nd);
            xflux_hi *= met_h_zeta_xhi;

            Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k,cellSizeInv,z_nd);
            yflux_lo *= met_h_zeta_ylo;
            Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j+1,k,cellSizeInv,z_nd);
            yflux_hi *= met_h_zeta_yhi;

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
                avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,0) = - invdetJ * (
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - invdetJ * 0.5 * (
                ( xflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index)) -
                  xflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index)) ) * dxInv * mfsq +
                ( yflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index)) -
                  yflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index)) ) * dyInv * mfsq +
                ( zflux_hi * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index)) -
                  zflux_lo * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index)) ) * dzInv);
        });
    } else if ( use_terrain && (std::max(horiz_spatial_order,vert_spatial_order) > 2) && !all_use_WENO) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real invdetJ = 1./ detJ(i,j,k);

            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i+1,j  ,k,cellSizeInv,z_nd);
            Real xflux_hi = rho_u(i+1,j  ,k) * met_h_zeta_xhi * 1./mf_u(i  ,j  ,0);
            Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k,cellSizeInv,z_nd);
            Real xflux_lo = rho_u(i  ,j  ,k) * met_h_zeta_xlo * 1./mf_u(i+1,j  ,0);
            Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j+1,k,cellSizeInv,z_nd);
            Real yflux_hi = rho_v(i  ,j+1,k) * met_h_zeta_yhi * 1./mf_v(i  ,j  ,0);
            Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k,cellSizeInv,z_nd);
            Real yflux_lo = rho_v(i  ,j  ,k) * met_h_zeta_ylo * 1./mf_v(i  ,j+1,0);

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
                avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mf   = mf_m(i,j,0);
            Real mfsq = mf*mf;

            advectionSrc(i,j,k,0) = - invdetJ * (
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - invdetJ * (
                ( xflux_hi * InterpolateInX(i+1,j  , k  , cell_prim, prim_index, rho_u(i+1,j,k), horiz_spatial_order) -
                  xflux_lo * InterpolateInX(i  ,j  , k  , cell_prim, prim_index, rho_u(i  ,j,k), horiz_spatial_order) ) * dxInv * mfsq +
                ( yflux_hi * InterpolateInY(i  ,j+1, k  , cell_prim, prim_index, rho_v(i,j+1,k), horiz_spatial_order) -
                  yflux_lo * InterpolateInY(i  ,j  , k  , cell_prim, prim_index, rho_v(i  ,j,k), horiz_spatial_order) ) * dyInv * mfsq +
                ( zflux_hi * InterpolateInZ(i  ,j  , k+1, cell_prim, prim_index, Omega(i,j,k+1), vert_spatial_order) -
                  zflux_lo * InterpolateInZ(i  ,j  , k  , cell_prim, prim_index, Omega(i  ,j,k), vert_spatial_order) ) * dzInv);
        });
    } else if (use_terrain) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real invdetJ = 1./ detJ(i,j,k);

            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i+1,j  ,k,cellSizeInv,z_nd);
            Real xflux_hi = rho_u(i+1,j  ,k) * met_h_zeta_xhi * 1./mf_u(i  ,j  ,0);
            Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k,cellSizeInv,z_nd);
            Real xflux_lo = rho_u(i  ,j  ,k) * met_h_zeta_xlo * 1./mf_u(i+1,j  ,0);
            Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j+1,k,cellSizeInv,z_nd);
            Real yflux_hi = rho_v(i  ,j+1,k) * met_h_zeta_yhi * 1./mf_v(i  ,j  ,0);
            Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k,cellSizeInv,z_nd);
            Real yflux_lo = rho_v(i  ,j  ,k) * met_h_zeta_ylo * 1./mf_v(i  ,j+1,0);

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
                avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mf   = mf_m(i,j,0);
            Real mfsq = mf*mf;

            advectionSrc(i,j,k,0) = - invdetJ * (
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - invdetJ * (
                ( xflux_hi * InterpolateInX_WENO(i+1,j  , k  , cell_prim, prim_index, spatial_order_WENO) -
                  xflux_lo * InterpolateInX_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dxInv * mfsq +
                ( yflux_hi * InterpolateInY_WENO(i  ,j+1, k  , cell_prim, prim_index, spatial_order_WENO) -
                  yflux_lo * InterpolateInY_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dyInv * mfsq +
                ( zflux_hi * InterpolateInZ_WENO(i  ,j  , k+1, cell_prim, prim_index, spatial_order_WENO) -
                  zflux_lo * InterpolateInZ_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dzInv);
        });
    } else if ( !use_terrain && (std::max(horiz_spatial_order,vert_spatial_order) == 2) && !all_use_WENO) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
            Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
            Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
            Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
               avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,0) = -(
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - 0.5 * (
              ( xflux_hi * (cell_prim(i+1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                xflux_lo * (cell_prim(i-1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dxInv * mfsq +
              ( yflux_hi * (cell_prim(i,j+1,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                yflux_lo * (cell_prim(i,j-1,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dyInv * mfsq +
              ( zflux_hi * (cell_prim(i,j,k+1,prim_index) + cell_prim(i,j,k,prim_index)) -
                zflux_lo * (cell_prim(i,j,k-1,prim_index) + cell_prim(i,j,k,prim_index)) ) * dzInv);
        });
    } else if ( !use_terrain && (std::max(horiz_spatial_order,vert_spatial_order) > 2) && !all_use_WENO) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
            Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
            Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
            Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
               avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mf   = mf_m(i,j,0);
            Real mfsq = mf*mf;

            advectionSrc(i,j,k,0) = -(
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = -(
              ( xflux_hi * InterpolateInX(i+1,j  , k  , cell_prim, prim_index, rho_u(i+1,j,k), horiz_spatial_order) -
                xflux_lo * InterpolateInX(i  ,j  , k  , cell_prim, prim_index, rho_u(i  ,j,k), horiz_spatial_order) ) * dxInv * mfsq +
              ( yflux_hi * InterpolateInY(i  ,j+1, k  , cell_prim, prim_index, rho_v(i,j+1,k), horiz_spatial_order) -
                yflux_lo * InterpolateInY(i  ,j  , k  , cell_prim, prim_index, rho_v(i  ,j,k), horiz_spatial_order) ) * dyInv * mfsq +
              ( zflux_hi * InterpolateInZ(i  ,j  , k+1, cell_prim, prim_index, Omega(i,j,k+1), vert_spatial_order) -
                zflux_lo * InterpolateInZ(i  ,j  , k  , cell_prim, prim_index, Omega(i  ,j,k), vert_spatial_order) ) * dzInv);
        });
    } else {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = rho_u(i  ,j,k) / mf_u(i  ,j  ,0);
            Real xflux_hi = rho_u(i+1,j,k) / mf_u(i+1,j  ,0);
            Real yflux_lo = rho_v(i,j  ,k) / mf_v(i  ,j  ,0);
            Real yflux_hi = rho_v(i,j+1,k) / mf_v(i  ,j+1,0);
            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
               avg_zmom(i,j,k+1) += fac*zflux_hi;

            Real mf   = mf_m(i,j,0);
            Real mfsq = mf*mf;

            advectionSrc(i,j,k,0) = -(
                ( xflux_hi - xflux_lo ) * dxInv * mfsq +
                ( yflux_hi - yflux_lo ) * dyInv * mfsq +
                ( zflux_hi - zflux_lo ) * dzInv);

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = -(
                ( xflux_hi * InterpolateInX_WENO(i+1,j  , k  , cell_prim, prim_index, spatial_order_WENO) -
                  xflux_lo * InterpolateInX_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dxInv * mfsq +
                ( yflux_hi * InterpolateInY_WENO(i  ,j+1, k  , cell_prim, prim_index, spatial_order_WENO) -
                  yflux_lo * InterpolateInY_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dyInv * mfsq +
                ( zflux_hi * InterpolateInZ_WENO(i  ,j  , k+1, cell_prim, prim_index, spatial_order_WENO) -
                  zflux_lo * InterpolateInZ_WENO(i  ,j  , k  , cell_prim, prim_index, spatial_order_WENO) ) * dzInv);
        });
    }
}

void
AdvectionSrcForScalars (const Box& bx, const int &icomp, const int &ncomp,
                        const Array4<const Real>& avg_xmom, const Array4<const Real>& avg_ymom,
                        const Array4<const Real>& avg_zmom,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& advectionSrc,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& mf_m,
                        const bool all_use_WENO,
                        const bool moist_use_WENO,
                        const int  spatial_order_WENO,
                        const int horiz_spatial_order, const int vert_spatial_order,
                        const int use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForScalars", AdvectionSrcForScalars);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    int ncomp_end{ncomp};
    int moist_off = RhoScalar_comp;
#if defined(ERF_USE_MOISTURE)
    moist_off = RhoQt_comp;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    moist_off = RhoQv_comp;
#endif

    // Running with WENO for moisture but not for other vars
    if(moist_use_WENO && ((icomp+ncomp)==NVAR) ) {
        ncomp_end -= 2;
        amrex::ParallelFor(bx, 2, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

            // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
            //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = moist_off + n;
            const int prim_index = cons_index - 1;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - invdetJ * (
            ( avg_xmom(i+1,j,k) *
                InterpolateInX_WENO(i+1,j,k,cell_prim, prim_index, spatial_order_WENO) -
              avg_xmom(i  ,j,k) *
                InterpolateInX_WENO(i  ,j,k,cell_prim, prim_index, spatial_order_WENO) ) * dxInv * mfsq +
            ( avg_ymom(i,j+1,k) *
                InterpolateInY_WENO(i,j+1,k,cell_prim, prim_index, spatial_order_WENO) -
              avg_ymom(i,j  ,k) *
                InterpolateInY_WENO(i,j  ,k,cell_prim, prim_index, spatial_order_WENO) ) * dyInv * mfsq +
            ( avg_zmom(i,j,k+1) *
                InterpolateInZ_WENO(i,j,k+1,cell_prim, prim_index, spatial_order_WENO) -
              avg_zmom(i,j,k  ) *
                InterpolateInZ_WENO(i,j,k  ,cell_prim, prim_index, spatial_order_WENO) ) * dzInv );
        });
    }

    if ((std::max(horiz_spatial_order,vert_spatial_order) == 2) && !all_use_WENO) {

        amrex::ParallelFor(bx, ncomp_end, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

            // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
            //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - 0.5 * invdetJ * (
              ( avg_xmom(i+1,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index)) -
                avg_xmom(i  ,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index)) ) * dxInv * mfsq +
              ( avg_ymom(i,j+1,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index)) -
                avg_ymom(i,j  ,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index)) ) * dyInv * mfsq +
              ( avg_zmom(i,j,k+1) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index)) -
                avg_zmom(i,j,k  ) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index)) ) * dzInv);
        });

    } else if ((std::max(horiz_spatial_order,vert_spatial_order) > 2) && !all_use_WENO) { // order > 2

        amrex::ParallelFor(bx, ncomp_end, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

            // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
            //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - invdetJ * (
            ( avg_xmom(i+1,j,k) *
                InterpolateInX(i+1,j,k,cell_prim, prim_index, avg_xmom(i+1,j,k), horiz_spatial_order) -
              avg_xmom(i  ,j,k) *
                InterpolateInX(i  ,j,k,cell_prim, prim_index, avg_xmom(i  ,j,k), horiz_spatial_order) ) * dxInv * mfsq +
            ( avg_ymom(i,j+1,k) *
                InterpolateInY(i,j+1,k,cell_prim, prim_index, avg_ymom(i,j+1,k), horiz_spatial_order) -
              avg_ymom(i,j  ,k) *
                InterpolateInY(i,j  ,k,cell_prim, prim_index, avg_ymom(i  ,j,k), horiz_spatial_order) ) * dyInv * mfsq +
            ( avg_zmom(i,j,k+1) *
                InterpolateInZ(i,j,k+1,cell_prim, prim_index, avg_zmom(i,j,k+1), vert_spatial_order) -
              avg_zmom(i,j,k  ) *
                InterpolateInZ(i,j,k  ,cell_prim, prim_index, avg_zmom(i  ,j,k), vert_spatial_order) ) * dzInv );
        });

    } else { // all_use_WENO

        amrex::ParallelFor(bx, ncomp_end, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

            // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
            //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            Real mfsq = mf_m(i,j,0) * mf_m(i,j,0);

            advectionSrc(i,j,k,cons_index) = - invdetJ * (
            ( avg_xmom(i+1,j,k) *
                InterpolateInX_WENO(i+1,j,k,cell_prim, prim_index, spatial_order_WENO) -
              avg_xmom(i  ,j,k) *
                InterpolateInX_WENO(i  ,j,k,cell_prim, prim_index, spatial_order_WENO) ) * dxInv * mfsq +
            ( avg_ymom(i,j+1,k) *
                InterpolateInY_WENO(i,j+1,k,cell_prim, prim_index, spatial_order_WENO) -
              avg_ymom(i,j  ,k) *
                InterpolateInY_WENO(i,j  ,k,cell_prim, prim_index, spatial_order_WENO) ) * dyInv * mfsq +
            ( avg_zmom(i,j,k+1) *
                InterpolateInZ_WENO(i,j,k+1,cell_prim, prim_index, spatial_order_WENO) -
              avg_zmom(i,j,k  ) *
                InterpolateInZ_WENO(i,j,k  ,cell_prim, prim_index, spatial_order_WENO) ) * dzInv );
        });

    }
}
