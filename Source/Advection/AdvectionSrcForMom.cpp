#include <AdvectionSrcForMom_N.H>
#include <AdvectionSrcForMom_T.H>

using namespace amrex;

void
AdvectionSrcForMom (const Box& bxx, const Box& bxy, const Box& bxz,
                    const Array4<      Real>& rho_u_rhs, const Array4<      Real>& rho_v_rhs,
                    const Array4<      Real>& rho_w_rhs,
                    const Array4<const Real>& u        , const Array4<const Real>& v,
                    const Array4<const Real>& w        ,
                    const Array4<const Real>& rho_u    , const Array4<const Real>& rho_v,
                    const Array4<const Real>& Omega    ,
                    const Array4<const Real>& z_nd     , const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int spatial_order, const int use_terrain, const int domhi_z)
{
    BL_PROFILE_VAR("AdvectionSrcForMom", AdvectionSrcForMom);

    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    if (use_terrain && (spatial_order == 2)) {

        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi, met_h_eta, met_h_zeta_ylo, met_h_zeta_yhi;

            Real met_h_zeta_xhi = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (u(i+1,j,k) + u(i,j,k)) * met_h_zeta_xhi;

            Real met_h_zeta_xlo = Compute_h_zeta_AtCellCenter(i-1,j  ,k  ,cellSizeInv,z_nd);
            Real xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (u(i-1,j,k) + u(i,j,k)) * met_h_zeta_xlo;

            ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta_yhi,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real yflux_hi = 0.25 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k)) * (u(i,j+1,k) + u(i,j,k)) * met_h_zeta_yhi;

            ComputeMetricAtEdgeCenterK(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta_ylo,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real yflux_lo = 0.25 * (rho_v(i, j  , k) + rho_v(i-1, j  , k)) * (u(i,j-1,k) + u(i,j,k)) * met_h_zeta_ylo;

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k));
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;

            rho_u_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i-1,j,k)));

        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi, met_h_eta,  met_h_zeta_lo, met_h_zeta_hi;

            ComputeMetricAtEdgeCenterK(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta_hi,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (v(i+1,j,k) + v(i,j,k)) * met_h_zeta_hi;

            ComputeMetricAtEdgeCenterK(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta_lo,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (v(i-1,j,k) + v(i,j,k)) * met_h_zeta_lo;

            Real met_h_zeta_yhi = Compute_h_zeta_AtCellCenter(i  ,j  ,k  ,cellSizeInv,z_nd);
            Real yflux_hi = 0.25 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k)) * (v(i,j+1,k) + v(i,j,k)) * met_h_zeta_yhi;

            Real met_h_zeta_ylo = Compute_h_zeta_AtCellCenter(i  ,j-1,k  ,cellSizeInv,z_nd);
            Real yflux_lo = 0.25 * (rho_v(i, j  , k) + rho_v(i-1, j  , k)) * (v(i,j-1,k) + v(i,j,k)) * met_h_zeta_ylo;

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k));
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;
            rho_v_rhs(i, j, k) = -advectionSrc / (0.5 * (detJ(i,j,k) + detJ(i,j-1,k)));
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real met_h_xi, met_h_eta,  met_h_zeta_lo, met_h_zeta_hi;

            ComputeMetricAtEdgeCenterJ(i+1,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta_hi,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real xflux_hi = 0.25*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * (w(i+1,j,k) + w(i,j,k)) * met_h_zeta_hi;

            ComputeMetricAtEdgeCenterJ(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta_lo,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real xflux_lo = 0.25*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * (w(i-1,j,k) + w(i,j,k)) * met_h_zeta_lo;

            ComputeMetricAtEdgeCenterI(i  ,j+1,k  ,met_h_xi,met_h_eta,met_h_zeta_hi,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real yflux_hi = 0.25*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * (w(i,j+1,k) + w(i,j,k)) * met_h_zeta_hi;

            ComputeMetricAtEdgeCenterI(i  ,j  ,k  ,met_h_xi,met_h_eta,met_h_zeta_lo,cellSizeInv,z_nd,TerrainMet::h_zeta);
            Real yflux_lo = 0.25*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * (w(i,j-1,k) + w(i,j,k)) * met_h_zeta_lo;

            Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

            Real zflux_hi = (k == domhi_z+1) ? Omega(i,j,k) * w(i,j,k) :
                0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;

            rho_w_rhs(i, j, k) = -advectionSrc / (0.5*(detJ(i,j,k) + detJ(i,j,k-1)));
        });

    } else if (!use_terrain && (spatial_order == 2)) {

        ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_hi = 0.25 * (rho_u(i, j  , k) + rho_u(i+1, j  , k)) * (u(i+1,j,k) + u(i,j,k));
            Real xflux_lo = 0.25 * (rho_u(i, j  , k) + rho_u(i-1, j  , k)) * (u(i-1,j,k) + u(i,j,k));

            Real yflux_hi = 0.25 * (rho_v(i, j+1, k) + rho_v(i-1, j+1, k)) * (u(i,j+1,k) + u(i,j,k));
            Real yflux_lo = 0.25 * (rho_v(i, j  , k) + rho_v(i-1, j  , k)) * (u(i,j-1,k) + u(i,j,k));

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i-1, j, k+1)) * (u(i,j,k+1) + u(i,j,k));
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i-1, j, k  )) * (u(i,j,k-1) + u(i,j,k));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;
            rho_u_rhs(i, j, k) = -advectionSrc;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_hi = 0.25 * (rho_u(i+1, j, k) + rho_u(i+1, j-1, k)) * (v(i+1,j,k) + v(i,j,k));
            Real xflux_lo = 0.25 * (rho_u(i  , j, k) + rho_u(i  , j-1, k)) * (v(i-1,j,k) + v(i,j,k));

            Real yflux_hi = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j+1, k)) * (v(i,j+1,k) + v(i,j,k));
            Real yflux_lo = 0.25 * (rho_v(i, j, k  ) + rho_v(i  , j-1, k)) * (v(i,j-1,k) + v(i,j,k));

            Real zflux_hi = 0.25 * (Omega(i, j, k+1) + Omega(i, j-1, k+1)) * (v(i,j,k+1) + v(i,j,k));
            Real zflux_lo = 0.25 * (Omega(i, j, k  ) + Omega(i, j-1, k  )) * (v(i,j,k-1) + v(i,j,k));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;
            rho_v_rhs(i, j, k) = -advectionSrc;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real xflux_hi = 0.25*(rho_u(i+1, j, k) + rho_u(i+1, j, k-1)) * (w(i+1,j,k) + w(i,j,k));
            Real xflux_lo = 0.25*(rho_u(i  , j, k) + rho_u(i  , j, k-1)) * (w(i-1,j,k) + w(i,j,k));

            Real yflux_hi = 0.25*(rho_v(i, j+1, k) + rho_v(i, j+1, k-1)) * (w(i,j+1,k) + w(i,j,k));
            Real yflux_lo = 0.25*(rho_v(i, j  , k) + rho_v(i, j  , k-1)) * (w(i,j-1,k) + w(i,j,k));

            Real zflux_lo = 0.25 * (Omega(i,j,k) + Omega(i,j,k-1)) * (w(i,j,k) + w(i,j,k-1));

            Real zflux_hi = (k == domhi_z+1) ? Omega(i,j,k) * w(i,j,k) :
                0.25 * (Omega(i,j,k) + Omega(i,j,k+1)) * (w(i,j,k) + w(i,j,k+1));

            Real advectionSrc = (xflux_hi - xflux_lo) * dxInv
                              + (yflux_hi - yflux_lo) * dyInv
                              + (zflux_hi - zflux_lo) * dzInv;
            rho_w_rhs(i, j, k) = -advectionSrc;
        });

    } else if (use_terrain && (spatial_order > 2)) {

        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_u_rhs(i, j, k) = -AdvectionSrcForXMom_T(i, j, k, rho_u, rho_v, Omega, u, z_nd, detJ,
                                                        cellSizeInv, spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_v_rhs(i, j, k) = -AdvectionSrcForYMom_T(i, j, k, rho_u, rho_v, Omega, v, z_nd, detJ,
                                                        cellSizeInv, spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_w_rhs(i, j, k) = -AdvectionSrcForZMom_T(i, j, k, rho_u, rho_v, Omega, w, z_nd, detJ,
                                                        cellSizeInv, spatial_order, domhi_z);
        });

    } else if (!use_terrain && (spatial_order > 2)) {

        amrex::ParallelFor(bxx, bxy, bxz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_u_rhs(i, j, k) = -AdvectionSrcForXMom_N(i, j, k, rho_u, rho_v, Omega, u,
                                                        cellSizeInv, spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_v_rhs(i, j, k) = -AdvectionSrcForYMom_N(i, j, k, rho_u, rho_v, Omega, v,
                                                        cellSizeInv, spatial_order);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            rho_w_rhs(i, j, k) = -AdvectionSrcForZMom_N(i, j, k, rho_u, rho_v, Omega, w,
                                                        cellSizeInv, spatial_order, domhi_z);
        });
    }
}

