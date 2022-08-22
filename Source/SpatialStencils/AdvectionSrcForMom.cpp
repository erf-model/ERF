#include <IndexDefines.H>
#include <SpatialStencils.H>
#include <AdvectionSrcForMom.H>
#include <TerrainMetrics.H>
#include <Interpolation.H>

using namespace amrex;

void
AdvectionSrcForMom (int level, const Box& bx, const Box& valid_bx,
                    const Array4<const int>& mlo_x, const Array4<const int>& mlo_y,
                    const Array4<const int>& mlo_z, const Array4<const int>& mhi_x,
                    const Array4<const int>& mhi_y, const Array4<const int>& mhi_z,
                    const Array4<      Real>& rho_u_rhs, const Array4<      Real>& rho_v_rhs,
                    const Array4<      Real>& rho_w_rhs,
                    const Array4<const Real>& u        , const Array4<const Real>& v,
                    const Array4<const Real>& w        ,
                    const Array4<const Real>& rho_u    , const Array4<const Real>& rho_v,
                    const Array4<const Real>& rho_w    , const Array4<const Real>& z_t,
                    const Array4<const Real>& z_nd     , const Array4<const Real>& detJ,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const int spatial_order, const int use_terrain, const int domhi_z)
{
    BL_PROFILE_VAR("AdvectionSrcForMom", AdvectionSrcForMom);
    auto dxInv = cellSizeInv[0], dyInv = cellSizeInv[1], dzInv = cellSizeInv[2];

    const Box bxx = surroundingNodes(bx,0);
    const Box bxy = surroundingNodes(bx,1);
    const Box bxz = surroundingNodes(bx,2);

    int vlo_x = valid_bx.smallEnd(0);
    int vhi_x = valid_bx.bigEnd(0);
    int vlo_y = valid_bx.smallEnd(1);
    int vhi_y = valid_bx.bigEnd(1);
    int vlo_z = valid_bx.smallEnd(2);
    int vhi_z = valid_bx.bigEnd(2);

    bool left_edge_dirichlet = ( level > 0 && mlo_x(vlo_x,vlo_y,vlo_z) );
    bool rght_edge_dirichlet = ( level > 0 && mhi_x(vhi_x,vhi_y,vhi_z) );
    bool  bot_edge_dirichlet = ( level > 0 && mlo_y(vlo_x,vlo_y,vlo_z) );
    bool  top_edge_dirichlet = ( level > 0 && mhi_y(vhi_x,vhi_y,vhi_z) );

    Box gbx(bx); gbx.growHi(0,1); gbx.growHi(1,1); gbx.growHi(2,1);

    amrex::ParallelFor(
        gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if ( j <= vhi_y && k <= vhi_z &&
                !(i == vlo_x && left_edge_dirichlet) && !(i == vhi_x+1 && rght_edge_dirichlet) )
                rho_u_rhs(i, j, k) = -AdvectionSrcForXMom(i, j, k, rho_u, rho_v, rho_w, z_t, u, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain);
            if (i <= vhi_x && k <= vhi_z &&
                !(j == vlo_y &&  bot_edge_dirichlet) && !(j == vhi_y+1 &&  top_edge_dirichlet) )
                rho_v_rhs(i, j, k) = -AdvectionSrcForYMom(i, j, k, rho_u, rho_v, rho_w, z_t, v, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain);
            if (i <= vhi_x && j <= vhi_y)
                rho_w_rhs(i, j, k) = -AdvectionSrcForZMom(i, j, k, rho_u, rho_v, rho_w, z_t, w, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain, domhi_z);
        }
    );

#if 0
    amrex::ParallelFor(
        bxx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            bool on_coarse_fine_boundary = false;
//          if (level > 0)
//          {
//             on_coarse_fine_boundary =
//               ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
//          }

//          if (!on_coarse_fine_boundary)
//          {
                // Add advective terms
                rho_u_rhs(i, j, k) = -AdvectionSrcForXMom(i, j, k, rho_u, rho_v, rho_w, z_t, u, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain);
//          }
        }
    );
    amrex::ParallelFor(
        bxy, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
//          bool on_coarse_fine_boundary = false;
//          if (level > 0)
//          {
//             on_coarse_fine_boundary =
//               ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
//          }

//          if (!on_coarse_fine_boundary)
//          {
                // Add advective terms
                rho_v_rhs(i, j, k) = -AdvectionSrcForYMom(i, j, k, rho_u, rho_v, rho_w, z_t, v, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain);
 //         }
        }
    );
    amrex::ParallelFor(
        bxz, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
 //         bool on_coarse_fine_boundary = false;
 //         if (level > 0)
 //         {
 //            on_coarse_fine_boundary =
 //              ( (i == vlo_x && mlo_x(i,j,k)) || (i == vhi_x+1 && mhi_x(i,j,k)) );
 //         }

 //         if (!on_coarse_fine_boundary)
 //         {
                // Add advective terms
                rho_w_rhs(i, j, k) = -AdvectionSrcForZMom(i, j, k, rho_u, rho_v, rho_w, z_t, w, z_nd, detJ,
                                                          cellSizeInv, spatial_order, use_terrain, domhi_z);
 //         }
        }
    );
#endif
}

