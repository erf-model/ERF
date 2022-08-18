#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

using namespace amrex;

void
AdvectionSrcForState (const Box& bx, const int &icomp, const int &ncomp,
                      const Array4<const Real>& rho_u, const Array4<const Real>& rho_v,
                      const Array4<const Real>& rho_w, const Array4<const Real>& z_t,
                      const Array4<const Real>& cell_data,
                      const Array4<const Real>& cell_prim,
                      const Array4<Real>& advectionSrc,
                      const Array4<Real>& xflux, const Array4<Real>& yflux, const Array4<Real>& zflux,
                      const Array4<const Real>& z_nd, const Array4<const Real>& detJ,
                      const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                      const int &spatial_order, const int& use_terrain,
                      const int &use_deardorff, const int &use_QKE)
{
    BL_PROFILE_VAR("AdvectionSrcForState", AdvectionSrcForState);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Box bxx = surroundingNodes(bx,0);
    const Box bxy = surroundingNodes(bx,1);
    const Box bxz = surroundingNodes(bx,2);

    amrex::ParallelFor(bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real xflux_lo;
        Real met_h_xi, met_h_eta,  met_h_zeta;

        xflux_lo = rho_u(i  ,j,k);

        if (use_terrain) {
            ComputeMetricAtIface(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
            xflux_lo *=  met_h_zeta;
        }

        // ****************************************************************************************
        // Now that we have the correctly weighted vector components, we can multiply by the
        //     scalar (such as theta or C) on the respective faces
        // ****************************************************************************************

        for (int n = icomp; n < icomp+ncomp; n++)
        {
          if ((n != RhoKE_comp && n != RhoQKE_comp) ||
              (  use_deardorff && n == RhoKE_comp) ||
              (  use_QKE       && n == RhoQKE_comp) )
          {
            if (n != Rho_comp)
            {
                const int prim_index = n - RhoTheta_comp;

                xflux(i,j,k,n) = xflux_lo * InterpolateFromCellOrFace(i  , j, k, cell_prim, prim_index, rho_u(i,j,k), Coord::x, spatial_order);

            } else {
                xflux(i,j,k,n) = xflux_lo;
            }
          } else {
                xflux(i,j,k,n) = 0.;
          }
        }
    });

    amrex::ParallelFor(bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real yflux_lo;
        Real met_h_xi, met_h_eta,  met_h_zeta;

        if (use_terrain) {

            ComputeMetricAtJface(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta,cellSizeInv,z_nd,TerrainMet::h_zeta);
            yflux_lo = rho_v(i,j  ,k) * met_h_zeta;

        } else {
            yflux_lo = rho_v(i,j  ,k);
        }

        // ****************************************************************************************
        // Now that we have the correctly weighted vector components, we can multiply by the
        //     scalar (such as theta or C) on the respective faces
        // ****************************************************************************************

        for (int n = icomp; n < icomp+ncomp; n++)
        {
          if ((n != RhoKE_comp && n != RhoQKE_comp) ||
              (  use_deardorff && n == RhoKE_comp) ||
              (  use_QKE       && n == RhoQKE_comp) )
          {
            if (n != Rho_comp)
            {
                const int prim_index = n - RhoTheta_comp;

                yflux(i,j,k,n) = yflux_lo * InterpolateFromCellOrFace(i, j  , k, cell_prim, prim_index, rho_v(i,j,k), Coord::y, spatial_order);

            } else {
                yflux(i,j,k,n) = yflux_lo;
            }
          } else {
                yflux(i,j,k,n) = 0.;
          }
        }
    });

    amrex::ParallelFor(bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real zflux_lo;

        Real z_t_zlo = (z_t) ? z_t(i,j,k) : 0.;

        if (use_terrain) {
            Real rho_zlo = 0.5 * ( cell_data(i,j,k) + cell_data(i,j,k-1) );
            zflux_lo = (k == 0) ? 0.0_rt : (OmegaFromW(i,j,k  ,rho_w(i,j,k  ),rho_u,rho_v,z_nd,cellSizeInv) - z_t_zlo * rho_zlo);
        } else {
            zflux_lo = rho_w(i,j,k  );
        }

        // ****************************************************************************************
        // Now that we have the correctly weighted vector components, we can multiply by the
        //     scalar (such as theta or C) on the respective faces
        // ****************************************************************************************

        for (int n = icomp; n < icomp+ncomp; n++)
        {
          if ((n != RhoKE_comp && n != RhoQKE_comp) ||
              (  use_deardorff && n == RhoKE_comp) ||
              (  use_QKE       && n == RhoQKE_comp) )
          {
            if (n != Rho_comp)
            {
                const int prim_index = n - RhoTheta_comp;

                zflux(i,j,k,n) = zflux_lo * InterpolateFromCellOrFace(i, j, k  , cell_prim, prim_index, rho_w(i,j,k), Coord::z, spatial_order);

            } else {

                zflux(i,j,k,n) = zflux_lo;
            }
          } else {
                zflux(i,j,k,n) = 0.;
          }
        }
    });

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real invdetJ = 1.;
        if (use_terrain) {
            invdetJ /= detJ(i,j,k);
        }

        for (int n = icomp; n < icomp+ncomp; n++)
        {
          if ((n != RhoKE_comp && n != RhoQKE_comp) ||
              (  use_deardorff && n == RhoKE_comp) ||
              (  use_QKE       && n == RhoQKE_comp) )
          {
            advectionSrc(i,j,k,n) = -( (xflux(i+1,j,k,n) - xflux(i,j,k,n)) * dxInv
                                      +(yflux(i,j+1,k,n) - yflux(i,j,k,n)) * dyInv
                                      +(zflux(i,j,k+1,n) - zflux(i,j,k,n)) * dzInv );
            advectionSrc(i,j,k,n) *= invdetJ;

          } else {

              advectionSrc(i,j,k,n) = 0.;
          }
      } // n
    });
}
