#include <IndexDefines.H>
#include <TerrainMetrics.H>
#include <Advection.H>
#include <Interpolation.H>

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
                            const int &spatial_order, const int& use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForRhoAndTheta", AdvectionSrcForRhoAndTheta);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    // We note that valid_bx is the actual grid, while bx may be a tile within that grid
    const auto& vbx_hi = amrex::ubound(valid_bx);

    if (use_terrain) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real invdetJ = 1./ detJ(i,j,k);

            Real met_h_xi, met_h_eta, met_h_zeta;

            Real zflux_lo = Omega(i,j,k  );
            Real zflux_hi = Omega(i,j,k+1);

            Real met_h_zeta_xhi = Compute_h_zeta_AtIface(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real xflux_hi = rho_u(i+1,k,k) * met_h_zeta_xhi;
            Real met_h_zeta_xlo = Compute_h_zeta_AtIface(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real xflux_lo = rho_u(i  ,k,k) * met_h_zeta_xlo;
            Real met_h_zeta_yhi = Compute_h_zeta_AtJface(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real yflux_hi = rho_v(i,j+1,k) * met_h_zeta_yhi;
            Real met_h_zeta_ylo = Compute_h_zeta_AtJface(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real yflux_lo = rho_v(i,j  ,k) * met_h_zeta_ylo;

            avg_xmom(i  ,j,k) += fac*xflux_lo;
            if (i == vbx_hi.x)
                avg_xmom(i+1,j,k) += fac*xflux_hi;
            avg_ymom(i,j  ,k) += fac*yflux_lo;
            if (j == vbx_hi.y)
                avg_ymom(i,j+1,k) += fac*yflux_hi;
            avg_zmom(i,j,k  ) += fac*zflux_lo;
            if (k == vbx_hi.z)
                avg_zmom(i,j,k+1) += fac*zflux_hi;

            advectionSrc(i,j,k,0) = - invdetJ * (
                ( xflux_hi - xflux_lo ) * dxInv + ( yflux_hi - yflux_lo ) * dyInv + ( zflux_hi - zflux_lo ) * dzInv );

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - invdetJ * (
                ( xflux_hi * InterpolateInX(i+1,j  , k  , cell_prim, prim_index, rho_u(i+1,j,k), spatial_order) -
                  xflux_lo * InterpolateInX(i  ,j  , k  , cell_prim, prim_index, rho_u(i  ,j,k), spatial_order) ) * dxInv +
                ( yflux_hi * InterpolateInY(i  ,j+1, k  , cell_prim, prim_index, rho_v(i,j+1,k), spatial_order) -
                  yflux_lo * InterpolateInY(i  ,j  , k  , cell_prim, prim_index, rho_v(i  ,j,k), spatial_order) ) * dyInv +
                ( zflux_hi * InterpolateInZ(i  ,j  , k+1, cell_prim, prim_index, Omega(i,j,k+1), spatial_order) -
                  zflux_lo * InterpolateInZ(i  ,j  , k  , cell_prim, prim_index, Omega(i  ,j,k), spatial_order) ) * dzInv );
    });
    } else if (spatial_order == 2) {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = rho_u(i  ,j,k);
            Real xflux_hi = rho_u(i+1,j,k);
            Real yflux_lo = rho_v(i,j  ,k);
            Real yflux_hi = rho_v(i,j+1,k);
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

            advectionSrc(i,j,k,0) = -(
                ( xflux_hi - xflux_lo ) * dxInv + ( yflux_hi - yflux_lo ) * dyInv + ( zflux_hi - zflux_lo ) * dzInv );

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = - 0.5 * (
              ( xflux_hi * (cell_prim(i+1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                xflux_lo * (cell_prim(i-1,j,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dxInv +
              ( yflux_hi * (cell_prim(i,j+1,k,prim_index) + cell_prim(i,j,k,prim_index)) -
                yflux_lo * (cell_prim(i,j-1,k,prim_index) + cell_prim(i,j,k,prim_index)) ) * dyInv +
              ( zflux_hi * (cell_prim(i,j,k+1,prim_index) + cell_prim(i,j,k,prim_index)) -
                zflux_lo * (cell_prim(i,j,k-1,prim_index) + cell_prim(i,j,k,prim_index)) ) * dzInv );
        });
    } else { // order > 2
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_lo = rho_u(i  ,j,k);
            Real xflux_hi = rho_u(i+1,j,k);
            Real yflux_lo = rho_v(i,j  ,k);
            Real yflux_hi = rho_v(i,j+1,k);
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

            advectionSrc(i,j,k,0) = -(
                ( xflux_hi - xflux_lo ) * dxInv + ( yflux_hi - yflux_lo ) * dyInv + ( zflux_hi - zflux_lo ) * dzInv );

            const int prim_index = 0;
            advectionSrc(i,j,k,1) = -(
              ( xflux_hi * InterpolateInX(i+1,j  , k  , cell_prim, prim_index, rho_u(i+1,j,k), spatial_order) -
                xflux_lo * InterpolateInX(i  ,j  , k  , cell_prim, prim_index, rho_u(i  ,j,k), spatial_order) ) * dxInv +
              ( yflux_hi * InterpolateInY(i  ,j+1, k  , cell_prim, prim_index, rho_v(i,j+1,k), spatial_order) -
                yflux_lo * InterpolateInY(i  ,j  , k  , cell_prim, prim_index, rho_v(i  ,j,k), spatial_order) ) * dyInv +
              ( zflux_hi * InterpolateInZ(i  ,j  , k+1, cell_prim, prim_index, Omega(i,j,k+1), spatial_order) -
                zflux_lo * InterpolateInZ(i  ,j  , k  , cell_prim, prim_index, Omega(i  ,j,k), spatial_order) ) * dzInv );
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
                        const int &spatial_order, const int& use_terrain)
{
    BL_PROFILE_VAR("AdvectionSrcForScalars", AdvectionSrcForScalars);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    if (spatial_order == 2)
    {
        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

           // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
           //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

              const int cons_index = icomp + n;
              const int prim_index = cons_index - 1;

              advectionSrc(i,j,k,cons_index) = - 0.5 * invdetJ * (
              ( avg_xmom(i+1,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i+1,j,k,prim_index)) -
                avg_xmom(i  ,j,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i-1,j,k,prim_index)) ) * dxInv +
              ( avg_ymom(i,j+1,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j+1,k,prim_index)) -
                avg_ymom(i,j  ,k) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j-1,k,prim_index)) ) * dyInv +
              ( avg_zmom(i,j,k+1) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k+1,prim_index)) -
                avg_zmom(i,j,k  ) * (cell_prim(i,j,k,prim_index) + cell_prim(i,j,k-1,prim_index)) ) * dzInv);
        });

    } else { // order > 2

        amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real invdetJ = (use_terrain) ?  1. / detJ(i,j,k) : 1.;

           // NOTE: we don't need to weight avg_xmom, avg_ymom, avg_zmom with terrain metrics
           //       because that was done when they were constructed in AdvectionSrcForRhoAndTheta

            const int cons_index = icomp + n;
            const int prim_index = cons_index - 1;

            advectionSrc(i,j,k,cons_index) = - invdetJ * (
            ( avg_xmom(i+1,j,k) *
                InterpolateInX(i+1,j,k,cell_prim, prim_index, avg_xmom(i+1,j,k), spatial_order) -
              avg_xmom(i  ,j,k) *
                InterpolateInX(i  ,j,k,cell_prim, prim_index, avg_xmom(i  ,j,k), spatial_order) ) * dxInv +
            ( avg_ymom(i,j+1,k) *
                InterpolateInY(i,j+1,k,cell_prim, prim_index, avg_ymom(i,j+1,k), spatial_order) -
              avg_ymom(i,j  ,k) *
                InterpolateInY(i,j  ,k,cell_prim, prim_index, avg_ymom(i  ,j,k), spatial_order) ) * dyInv +
            ( avg_zmom(i,j,k+1) *
                InterpolateInZ(i,j,k+1,cell_prim, prim_index, avg_zmom(i,j,k+1), spatial_order) -
                avg_zmom(i,j,k  ) *
                InterpolateInZ(i,j,k  ,cell_prim, prim_index, avg_zmom(i  ,j,k), spatial_order) ) * dzInv);
        });
    }
}
