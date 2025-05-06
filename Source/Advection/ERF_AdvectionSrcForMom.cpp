#include "AMReX_BCRec.H"

#include "ERF_Advection.H"
#include "ERF_AdvectionSrcForMom_N.H"
#include "ERF_AdvectionSrcForMom_T.H"
#include "ERF_EB.H"

using namespace amrex;

/**
 * Function for computing the advective tendency for the momentum equations
 * This routine has explicit expressions for all cases (terrain or not) when
 * the horizontal and vertical spatial orders are <= 2, and calls more specialized
 * functions when either (or both) spatial order(s) is greater than 2.
 *
 * @param[in] mfi MultiFab Iterator
 * @param[in] bxx box over which the x-momentum is updated
 * @param[in] bxy box over which the y-momentum is updated
 * @param[in] bxz box over which the z-momentum is updated
 * @param[in] bxx_grown grown boxes of bxx to loop over the nodal grids of bxx
 * @param[in] bxy_grown grown boxes of bxy to loop over the nodal grids of bxy
 * @param[in] bxz_grown grown boxes of bxz to loop over the nodal grids of bxz
 * @param[out] rho_u_rhs tendency for the x-momentum equation
 * @param[out] rho_v_rhs tendency for the y-momentum equation
 * @param[out] rho_w_rhs tendency for the z-momentum equation
 * @param[in] u x-component of the velocity
 * @param[in] v y-component of the velocity
 * @param[in] w z-component of the velocity
 * @param[in] rho_u x-component of the momentum
 * @param[in] rho_v y-component of the momentum
 * @param[in] Omega component of the momentum normal to the z-coordinate surface
 * @param[in] z_nd height coordinate at nodes
 * @param[in] ax   Area fraction of x-faces
 * @param[in] ay   Area fraction of y-faces
 * @param[in] az   Area fraction of z-faces
 * @param[in] detJ Jacobian of the metric transformation (= 1 if use_terrain_fitted_coords is false)
 * @param[in] cellSizeInv inverse of the mesh spacing
 * @param[in] mf_mx map factor at cell centers
 * @param[in] mf_my map factor at cell centers
 * @param[in] mf_ux map factor at x-faces
 * @param[in] mf_vy map factor at y-faces
 * @param[in] ebfact EB factories for cell- and face-centered variables
 * @param[in] horiz_adv_type sets the spatial order to be used for lateral derivatives
 * @param[in] vert_adv_type  sets the spatial order to be used for vertical derivatives
 * @param[in] ebfact EB factories for cell- and face-centered variables
 * @param[in] flx_u_arr Container of fluxes for x-momentum
 * @param[in] flx_v_arr Container of fluxes for y-momentum
 * @param[in] flx_w_arr Container of fluxes for z-momentum
 * @param[in] physbnd_mask Vector of masks for flux interpolation (=1 otherwise, =0 if physbnd)
 * @param[in] already_on_centroids flag whether flux interpolation is unnecessary
 */
void
AdvectionSrcForMom (const MFIter& mfi,
                    const Box& bx,
                    const Box& bxx, const Box& bxy, const Box& bxz,
                    const Vector<Box>& bxx_grown,
                    const Vector<Box>& bxy_grown,
                    const Vector<Box>& bxz_grown,
                    const Array4<      Real>& rho_u_rhs,
                    const Array4<      Real>& rho_v_rhs,
                    const Array4<      Real>& rho_w_rhs,
                    const Array4<const Real>& cell_data,
                    const Array4<const Real>& u,
                    const Array4<const Real>& v,
                    const Array4<const Real>& w,
                    const Array4<const Real>& rho_u,
                    const Array4<const Real>& rho_v,
                    const Array4<const Real>& omega,
                    const Array4<const Real>& z_nd,
                    const Array4<const Real>& ax,
                    const Array4<const Real>& ay,
                    const Array4<const Real>& az,
                    const Array4<const Real>& detJ,
                          Gpu::DeviceVector<Real>& stretched_dz_d,
                    const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                    const Array4<const Real>& mf_mx,
                    const Array4<const Real>& mf_ux,
                    const Array4<const Real>& mf_vx,
                    const Array4<const Real>& mf_my,
                    const Array4<const Real>& mf_uy,
                    const Array4<const Real>& mf_vy,
                    const AdvType horiz_adv_type,
                    const AdvType vert_adv_type,
                    const Real horiz_upw_frac,
                    const Real vert_upw_frac,
                    MeshType& mesh_type,
                    TerrainType& terrain_type,
                    const eb_& ebfact,
                          GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_u_arr,
                          GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_v_arr,
                          GpuArray<Array4<Real>, AMREX_SPACEDIM>& flx_w_arr,
                    const Vector<iMultiFab>& physbnd_mask,
                    const bool already_on_centroids,
                    const int lo_z_face, const int hi_z_face,
                    const Box& domain,
                    const BCRec* bc_ptr_h)
{
    BL_PROFILE_VAR("AdvectionSrcForMom", AdvectionSrcForMom);

    AMREX_ALWAYS_ASSERT(bxz.smallEnd(2) > 0);

    // compute mapfactor inverses
    Box box2d_u(bxx);   box2d_u.setRange(2,0);   box2d_u.grow({3,3,0});
    Box box2d_v(bxy);   box2d_v.setRange(2,0);   box2d_v.grow({3,3,0});
    FArrayBox mf_u_invFAB(box2d_u,1,The_Async_Arena());
    FArrayBox mf_v_invFAB(box2d_v,1,The_Async_Arena());
    const Array4<Real>& mf_u_inv = mf_u_invFAB.array();
    const Array4<Real>& mf_v_inv = mf_v_invFAB.array();

    const bool use_terrain_fitted_coords = ( terrain_type == TerrainType::StaticFittedMesh ||
                                             terrain_type == TerrainType::MovingFittedMesh);

    ParallelFor(box2d_u, box2d_v,
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_u_inv(i,j,0) = 1. / mf_ux(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        mf_v_inv(i,j,0) = 1. / mf_vy(i,j,0);
    });

    if (mesh_type == MeshType::ConstantDz && terrain_type != TerrainType::EB)
    {
        // amrex::Print() << "ADV:CONSTANT DZ " << std::endl;
        AdvectionSrcForMom_ConstantDz(bxx, bxy, bxz,
                                      rho_u_rhs, rho_v_rhs, rho_w_rhs, u, v, w,
                                      rho_u, rho_v, omega,
                                      cellSizeInv, stretched_dz_d,
                                      mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy,
                                      horiz_adv_type, vert_adv_type,
                                      horiz_upw_frac, vert_upw_frac,
                                      terrain_type, lo_z_face, hi_z_face);
    }
    else if (mesh_type == MeshType::StretchedDz && terrain_type != TerrainType::EB)
    {
        // amrex::Print() << "ADV:STRETCHED DZ " << std::endl;
        AdvectionSrcForMom_StretchedDz(bxx, bxy, bxz,
                                       rho_u_rhs, rho_v_rhs, rho_w_rhs,
                                       u, v, w, rho_u, rho_v, omega,
                                       cellSizeInv, stretched_dz_d,
                                       mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy,
                                       horiz_adv_type, vert_adv_type,
                                       horiz_upw_frac, vert_upw_frac,
                                       lo_z_face, hi_z_face);
    }
    else if ( terrain_type == TerrainType::EB)
    {
        // amrex::Print() << "ADV:EB " << std::endl;
        AdvectionSrcForMom_EB(mfi, bxx, bxy, bxz, bxx_grown, bxy_grown, bxz_grown,
                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                              u, v, w,
                              rho_u, rho_v, omega,
                              cellSizeInv,
                              mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy,
                              horiz_adv_type, vert_adv_type,
                              horiz_upw_frac, vert_upw_frac,
                              ebfact, flx_u_arr, flx_v_arr, flx_w_arr,
                              physbnd_mask, already_on_centroids,
                              lo_z_face, hi_z_face, domain);
    }
    else
    {
        AMREX_ALWAYS_ASSERT(use_terrain_fitted_coords);
        // amrex::Print() << "ADV:TF " << std::endl;
        AdvectionSrcForMom_TF(bxx, bxy, bxz,
                              rho_u_rhs, rho_v_rhs, rho_w_rhs,
                              u, v, w,
                              rho_u, rho_v, omega,
                              z_nd, ax, ay, az, detJ,
                              cellSizeInv,
                              mf_mx, mf_ux, mf_vx, mf_my, mf_uy, mf_vy,
                              horiz_adv_type, vert_adv_type,
                              horiz_upw_frac, vert_upw_frac,
                              lo_z_face, hi_z_face);

    }

    // Open bc will be imposed upon all vars (we only access cons here for simplicity)
    const bool xlo_open = (bc_ptr_h[BCVars::cons_bc].lo(0) == ERFBCType::open);
    const bool xhi_open = (bc_ptr_h[BCVars::cons_bc].hi(0) == ERFBCType::open);
    const bool ylo_open = (bc_ptr_h[BCVars::cons_bc].lo(1) == ERFBCType::open);
    const bool yhi_open = (bc_ptr_h[BCVars::cons_bc].hi(1) == ERFBCType::open);

    // We recreate tbx, tbz, tbz here rather than using bxx, bxy, bxz because those
    //    have already been shrunk by one in the case of open BCs.
    Box tbx(surroundingNodes(bx,0));
    Box tby(surroundingNodes(bx,1));
    Box tbz(surroundingNodes(bx,2)); tbz.growLo(2,-1); tbz.growHi(2,-1);

    const int domhi_z = domain.bigEnd(2);

    // Special advection operator for open BC (bndry normal/tangent operations)
    if (xlo_open)
    {
        Box tbx_xlo, tby_xlo, tbz_xlo;
        if (tbx.smallEnd(0) == domain.smallEnd(0)) { tbx_xlo = makeSlab(tbx,0,domain.smallEnd(0));}
        if (tby.smallEnd(0) == domain.smallEnd(0)) { tby_xlo = makeSlab(tby,0,domain.smallEnd(0));}
        if (tbz.smallEnd(0) == domain.smallEnd(0)) { tbz_xlo = makeSlab(tbz,0,domain.smallEnd(0));}

        bool do_lo = true;

        AdvectionSrcForOpenBC_Normal(tbx_xlo, 0, rho_u_rhs, u, cell_data, cellSizeInv, do_lo);
        AdvectionSrcForOpenBC_Tangent_Ymom(tby_xlo, 0, rho_v_rhs, v,
                                           rho_u, rho_v, omega,
                                           ay, az, detJ, cellSizeInv,
                                           do_lo);
        AdvectionSrcForOpenBC_Tangent_Zmom(tbz_xlo, 0, rho_w_rhs, w,
                                           rho_u, rho_v, omega,
                                           ax, ay, az, detJ, cellSizeInv,
                                           domhi_z, do_lo);
    }
    if (xhi_open)
    {
        Box tbx_xhi, tby_xhi, tbz_xhi;
        if (tbx.bigEnd(0) == domain.bigEnd(0)+1)   { tbx_xhi = makeSlab(tbx,0,domain.bigEnd(0)+1);}
        if (tby.bigEnd(0) == domain.bigEnd(0))     { tby_xhi = makeSlab(tby,0,domain.bigEnd(0)  );}
        if (tbz.bigEnd(0) == domain.bigEnd(0))     { tbz_xhi = makeSlab(tbz,0,domain.bigEnd(0)  );}

        AdvectionSrcForOpenBC_Normal(tbx_xhi, 0, rho_u_rhs, u, cell_data, cellSizeInv);
        AdvectionSrcForOpenBC_Tangent_Ymom(tby_xhi, 0, rho_v_rhs, v,
                                           rho_u, rho_v, omega,
                                           ay, az, detJ, cellSizeInv);
        AdvectionSrcForOpenBC_Tangent_Zmom(tbz_xhi, 0, rho_w_rhs, w,
                                           rho_u, rho_v, omega,
                                           ax, ay, az, detJ, cellSizeInv,
                                           domhi_z);
    }
    if (ylo_open)
    {
        Box tbx_ylo, tby_ylo, tbz_ylo;
        if (tbx.smallEnd(1) == domain.smallEnd(1)) { tbx_ylo = makeSlab(tbx,1,domain.smallEnd(1));}
        if (tby.smallEnd(1) == domain.smallEnd(1)) { tby_ylo = makeSlab(tby,1,domain.smallEnd(1));}
        if (tbz.smallEnd(1) == domain.smallEnd(1)) { tbz_ylo = makeSlab(tbz,1,domain.smallEnd(1));}

        bool do_lo = true;
        AdvectionSrcForOpenBC_Tangent_Xmom(tbx_ylo, 1, rho_u_rhs, u,
                                           rho_u, rho_v, omega,
                                           ax, az, detJ, cellSizeInv,
                                           do_lo);
        AdvectionSrcForOpenBC_Normal(tby_ylo, 1, rho_v_rhs, v, cell_data, cellSizeInv, do_lo);
        AdvectionSrcForOpenBC_Tangent_Zmom(tbz_ylo, 1, rho_w_rhs, w,
                                           rho_u, rho_v, omega,
                                           ax, ay, az, detJ, cellSizeInv,
                                           domhi_z, do_lo);
    }
    if (yhi_open)
    {
        Box tbx_yhi, tby_yhi, tbz_yhi;
        if (tbx.bigEnd(1) == domain.bigEnd(1))     { tbx_yhi = makeSlab(tbx,1,domain.bigEnd(1)  );}
        if (tby.bigEnd(1) == domain.bigEnd(1)+1)   { tby_yhi = makeSlab(tby,1,domain.bigEnd(1)+1);}
        if (tbz.bigEnd(1) == domain.bigEnd(1))     { tbz_yhi = makeSlab(tbz,1,domain.bigEnd(1)  );}

        AdvectionSrcForOpenBC_Tangent_Xmom(tbx_yhi, 1, rho_u_rhs, u,
                                           rho_u, rho_v, omega,
                                           ax, az, detJ, cellSizeInv);
        AdvectionSrcForOpenBC_Normal(tby_yhi, 1, rho_v_rhs, v, cell_data, cellSizeInv);
        AdvectionSrcForOpenBC_Tangent_Zmom(tbz_yhi, 1, rho_w_rhs, w,
                                           rho_u, rho_v, omega,
                                           ax, ay, az, detJ, cellSizeInv,
                                           domhi_z);
    }
}

