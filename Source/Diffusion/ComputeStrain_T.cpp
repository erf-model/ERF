#include <Diffusion.H>
#include <TerrainMetrics.H>

using namespace amrex;

/**
 * Function for computing the strain rates with terrain.
 *
 * @param[in] bxcc cell center box for tau_ii
 * @param[in] tbxxy nodal xy box for tau_12
 * @param[in] tbxxz nodal xz box for tau_13
 * @param[in] tbxyz nodal yz box for tau_23
 * @param[in] u x-direction velocity
 * @param[in] v y-direction velocity
 * @param[in] w z-direction velocity
 * @param[out] tau11 11 strain
 * @param[out] tau22 22 strain
 * @param[out] tau33 33 strain
 * @param[out] tau12 12 strain
 * @param[out] tau13 13 strain
 * @param[out] tau21 21 strain
 * @param[out] tau23 23 strain
 * @param[out] tau31 31 strain
 * @param[out] tau32 32 strain
 * @param[in] z_nd nodal array of physical z heights
 * @param[in] bc_ptr container with boundary condition types
 * @param[in] dxInv inverse cell size array
 * @param[in] mf_u map factor at x-face
 * @param[in] mf_v map factor at y-face
 */
void
ComputeStrain_T (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Box domain,
                 const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                 Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                 Array4<Real>& tau12, Array4<Real>& tau13,
                 Array4<Real>& tau21, Array4<Real>& tau23,
                 Array4<Real>& tau31, Array4<Real>& tau32,
                 const Array4<const Real>& z_nd  ,
                 const BCRec* bc_ptr, const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                 const Array4<const Real>& /*mf_m*/,
                 const Array4<const Real>& mf_u,
                 const Array4<const Real>& mf_v)
{
    // Conver domain to each index type to test if we are on dirichlet boundary
    Box domain_xy = convert(domain, tbxxy.ixType());
    Box domain_xz = convert(domain, tbxxz.ixType());
    Box domain_yz = convert(domain, tbxyz.ixType());

    // Dirichlet on left or right plane
    bool xl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_v_dir = ( xl_v_dir && (tbxxy.smallEnd(0) == domain_xy.smallEnd(0)) );

    bool xh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_v_dir = ( xh_v_dir && (tbxxy.bigEnd(0) == domain_xy.bigEnd(0)) );

    bool xl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
         xl_w_dir = ( xl_w_dir && (tbxxz.smallEnd(0) == domain_xz.smallEnd(0)) );

    bool xh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );
         xh_w_dir = ( xh_w_dir && (tbxxz.bigEnd(0) == domain_xz.bigEnd(0)) );

    // Dirichlet on front or back plane
    bool yl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_u_dir = ( yl_u_dir && (tbxxy.smallEnd(1) == domain_xy.smallEnd(1)) );

    bool yh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );
         yh_u_dir = ( yh_u_dir && (tbxxy.bigEnd(1) == domain_xy.bigEnd(1)) );

    bool yl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
         yl_w_dir = ( yl_w_dir && (tbxyz.smallEnd(1) == domain_yz.smallEnd(1)) );

    bool yh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );
         yh_w_dir = ( yh_w_dir && (tbxyz.bigEnd(1) == domain_yz.bigEnd(1)) );

    // Dirichlet on top or bottom plane
    bool zl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
         zl_u_dir = ( zl_u_dir && (tbxxz.smallEnd(2) == domain_xz.smallEnd(2)) );

    bool zh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );
         zh_u_dir = ( zh_u_dir && (tbxxz.bigEnd(2) == domain_xz.bigEnd(2)) );

    bool zl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
         zl_v_dir = ( zl_v_dir && (tbxyz.smallEnd(2) == domain_yz.smallEnd(2)) );

    bool zh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );
         zh_v_dir = ( zh_v_dir && (tbxyz.bigEnd(2) == domain_yz.bigEnd(2)) );


    // X-Dirichlet
    //***********************************************************************************
    if (xl_v_dir) {
        Box planexy = tbxxy; planexy.setBig(0, planexy.smallEnd(0) );
        tbxxy.growLo(0,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i  ,j-1,k+1)/mf_u(i,j-1,0)
                                             -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i  ,j-1,k-1)/mf_u(i,j-1,0) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i-1,j  ,k+1)/mf_v(i-1,j,0)
                                             -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i-1,j  ,k-1)/mf_v(i-1,j,0) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

            tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i, j-1, k)/mf_u(i,j-1,0))*dxInv[1]
                               + (-(8./3.) * v(i-1,j,k)/mf_v(i-1,j,0) + 3. * v(i,j,k)/mf_v(i,j,0) - (1./3.) * v(i+1,j,k)/mf_v(i+1,j,0))*dxInv[0]
                                 - (met_h_eta/met_h_zeta)*GradUz
                                 - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (xh_v_dir) {
        Box planexy = tbxxy; planexy.setSmall(0, planexy.bigEnd(0) );
        tbxxy.growHi(0,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i  ,j-1,k+1)/mf_u(i,j-1,0)
                                             -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i  ,j-1,k-1)/mf_u(i,j-1,0) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i-1,j  ,k+1)/mf_v(i-1,j,0)
                                             -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i-1,j  ,k-1)/mf_v(i-1,j,0) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

            tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i, j-1, k)/mf_u(i,j-1,0))*dxInv[1]
                               - (-(8./3.) * v(i,j,k)/mf_v(i,j,0) + 3. * v(i-1,j,k)/mf_v(i-1,j,0) - (1./3.) * v(i-2,j,k)/mf_v(i-2,j,0))*dxInv[0]
                               - (met_h_eta/met_h_zeta)*GradUz
                               - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (xl_w_dir) {
        Box planexz = tbxxz; planexz.setBig(0, planexz.smallEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growLo(0,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                 + ( (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dxInv[0]
                                     - (met_h_xi/met_h_zeta)*GradWz ) * mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (xh_w_dir) {
        Box planexz = tbxxz; planexz.setSmall(0, planexz.bigEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growHi(0,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                 - ( (-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dxInv[0]
                                     - (met_h_xi/met_h_zeta)*GradWz ) * mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    // Y-Dirichlet
    //***********************************************************************************
    if (yl_u_dir) {
        Box planexy = tbxxy; planexy.setBig(1, planexy.smallEnd(1) );
        tbxxy.growLo(1,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i  ,j-1,k+1)/mf_u(i,j-1,0)
                                             -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i  ,j-1,k-1)/mf_u(i,j-1,0) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i-1,j  ,k+1)/mf_v(i-1,j,0)
                                             -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i-1,j  ,k-1)/mf_v(i-1,j,0) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

            tau12(i,j,k) = 0.5 * ( (-(8./3.) * u(i,j-1,k)/mf_u(i,j-1,0) + 3. * u(i,j,k)/mf_u(i,j,0) - (1./3.) * u(i,j+1,k)/mf_u(i,j+1,0))*dxInv[1]
                               + (v(i, j, k)/mf_v(i,j,0) - v(i-1, j, k)/mf_v(i,j,0))*dxInv[0]
                               - (met_h_eta/met_h_zeta)*GradUz
                               - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (yh_u_dir) {
        Box planexy = tbxxy; planexy.setSmall(1, planexy.bigEnd(1) );
        tbxxy.growHi(1,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i  ,j-1,k+1)/mf_u(i,j-1,0)
                                             -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i  ,j-1,k-1)/mf_u(i,j-1,0) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i-1,j  ,k+1)/mf_v(i-1,j,0)
                                             -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i-1,j  ,k-1)/mf_v(i-1,j,0) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

            tau12(i,j,k) = 0.5 * ( -(-(8./3.) * u(i,j,k)/mf_u(i,j,0) + 3. * u(i,j-1,k)/mf_u(i,j-1,0) - (1./3.) * u(i,j-2,k)/mf_u(i,j-2,0))*dxInv[1] +
                               + (v(i, j, k)/mf_v(i,j,0) - v(i-1, j, k)/mf_v(i-1,j,0))*dxInv[0]
                               - (met_h_eta/met_h_zeta)*GradUz
                               - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (yl_w_dir) {
        Box planeyz = tbxyz; planeyz.setBig(1, planeyz.smallEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growLo(1,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                 + ( (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dxInv[1]*mf_v(i,j,0)
                                   - (met_h_eta/met_h_zeta)*GradWz ) * mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }
    if (yh_w_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(1, planeyz.bigEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growHi(1,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
             Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                             - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                 - ( (-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dxInv[1]
                                     - (met_h_eta/met_h_zeta)*GradWz ) * mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-Dirichlet
    //***********************************************************************************
    if (zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2]/met_h_zeta
                                 + ( (w(i, j, k) - w(i-1, j, k))*dxInv[0]
                                   - (met_h_xi/met_h_zeta)*GradWz ) * mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dxInv[2]/met_h_zeta
                               + (w(i, j, k) - w(i-1, j, k))*dxInv[0]*mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    if (zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );

            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dxInv[2]/met_h_zeta
                                 + ( (w(i, j, k) - w(i, j-1, k))*dxInv[1]
                                   - (met_h_eta/met_h_zeta)*GradWz ) * mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }
    if (zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dxInv[2]/met_h_zeta
                                 + (w(i, j, k) - w(i, j-1, k))*dxInv[1]*mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // HO derivatives w/ Dirichlet BC (\partival <var> / \partial z from terrain transform)
    if (zl_u_dir && zl_v_dir) {
        Box planecc = bxcc; planecc.setBig(2, planecc.smallEnd(2) );
        bxcc.growLo(2,-1);
        ParallelFor(planecc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.5 * dxInv[2] * ( (-(8./3.) * u(i  ,j,k-1) + 3. * u(i  ,j,k) - (1./3.) * u(i  ,j,k+1))
                                           + (-(8./3.) * u(i-1,j,k-1) + 3. * u(i-1,j,k) - (1./3.) * u(i-1,j,k+1)) );
            Real GradVz = 0.5 * dxInv[2] * ( (-(8./3.) * v(i,j  ,k-1) + 3. * v(i,j  ,k) - (1./3.) * v(i,j  ,k+1))
                                           + (-(8./3.) * v(i,j-1,k-1) + 3. * v(i,j-1,k) - (1./3.) * v(i,j-1,k+1)) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);

            tau11(i,j,k) = ( (u(i+1, j, k)/mf_u(i+1,j,0) - u(i, j, k)/mf_u(i,j,0))*dxInv[0]
                           - (met_h_xi/met_h_zeta)*GradUz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau22(i,j,k) = ( (v(i, j+1, k)/mf_v(i,j+1,0) - v(i, j, k)/mf_v(i,j,0))*dxInv[1]
                           - (met_h_eta/met_h_zeta)*GradVz ) * mf_v(i,j,0)*mf_v(i,j,0);
            tau33(i,j,k) = (w(i, j, k+1) - w(i, j, k))*dxInv[2]/met_h_zeta;
        });

        Box planexy = tbxxy; planexy.setBig(2, planexy.smallEnd(2) );
        tbxxy.growLo(2,-1);
        ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.5 * dxInv[2] * ( (-(8./3.) * u(i,j  ,k-1) + 3. * u(i,j  ,k) - (1./3.) * u(i,j  ,k+1))
                                           + (-(8./3.) * u(i,j-1,k-1) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-1,k+1)) );
            Real GradVz = 0.5 * dxInv[2] * ( (-(8./3.) * v(i  ,j,k-1) + 3. * v(i  ,j,k) - (1./3.) * v(i  ,j,k+1))
                                           + (-(8./3.) * v(i-1,j,k-1) + 3. * v(i-1,j,k) - (1./3.) * v(i-1,j,k+1)) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

            tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i  , j-1, k)/mf_u(i,j-1,0))*dxInv[1]
                                 + (v(i, j, k)/mf_v(i,j,0) - v(i-1, j  , k)/mf_v(i-1,j,0))*dxInv[0]
                                 - (met_h_eta/met_h_zeta)*GradUz
                                 - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    // Z-lo w/out Z-Dirichlet (GradWz extrapolation)
    //***********************************************************************************
    if (!zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i  , j, k-1))*dxInv[2]/met_h_zeta
                                 + ( (w(i, j, k) - w(i-1, j, k  ))*dxInv[0]
                                   - (met_h_xi/met_h_zeta)*GradWz ) * mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );

            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j  , k-1))*dxInv[2]/met_h_zeta
                                 + ( (w(i, j, k) - w(i, j-1, k  ))*dxInv[1]
                                   - (met_h_eta/met_h_zeta)*GradWz ) * mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-hi w/out Z-Dirichlet (h_xi = h_eta = 0)
    //***********************************************************************************
    if (!zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i  , j, k-1))*dxInv[2]/met_h_zeta
                                 + (w(i, j, k) - w(i-1, j, k  ))*dxInv[0]*mf_u(i,j,0) );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j  , k-1))*dxInv[2]/met_h_zeta
                                 + (w(i, j, k) - w(i, j-1, k  ))*dxInv[1]*mf_v(i,j,0) );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Fill the interior cells
    //***********************************************************************************
    // Cell centered strains
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i-1,j  ,k+1)/mf_u(i-1,j,0)
                                         -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i-1,j  ,k-1)/mf_u(i-1,j,0) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i  ,j-1,k+1)/mf_v(i,j-1,0)
                                         -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i  ,j-1,k-1)/mf_v(i,j-1,0) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);

        tau11(i,j,k) = ( (u(i+1, j, k)/mf_u(i+1,j,0) - u(i, j, k)/mf_u(i,j,0))*dxInv[0]
                       - (met_h_xi/met_h_zeta)*GradUz ) * mf_u(i,j,0)*mf_u(i,j,0);
        tau22(i,j,k) = ( (v(i, j+1, k)/mf_v(i,j+1,0) - v(i, j, k)/mf_v(i,j,0))*dxInv[1]
                       - (met_h_eta/met_h_zeta)*GradVz ) * mf_v(i,j,0)*mf_v(i,j,0);
        tau33(i,j,k) = (w(i, j, k+1) - w(i, j, k))*dxInv[2]/met_h_zeta;
    });

    // Off-diagonal strains
    ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1)/mf_u(i,j,0) + u(i  ,j-1,k+1)/mf_u(i,j-1,0)
                                         -u(i  ,j  ,k-1)/mf_u(i,j,0) - u(i  ,j-1,k-1)/mf_u(i,j-1,0) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1)/mf_v(i,j,0) + v(i-1,j  ,k+1)/mf_v(i-1,j,0)
                                         -v(i  ,j  ,k-1)/mf_v(i,j,0) - v(i-1,j  ,k-1)/mf_v(i-1,j,0) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        tau12(i,j,k) = 0.5 * ( (u(i, j, k)/mf_u(i,j,0) - u(i  , j-1, k)/mf_u(i,j-1,0))*dxInv[1]
                             + (v(i, j, k)/mf_v(i,j,0) - v(i-1, j  , k)/mf_v(i-1,j,0))*dxInv[0]
                             - (met_h_eta/met_h_zeta)*GradUz
                             - (met_h_xi /met_h_zeta)*GradVz ) * mf_u(i,j,0)*mf_u(i,j,0);
        tau21(i,j,k) = tau12(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                        - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

        Real met_h_xi,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        tau13(i,j,k) = 0.5 * ( (u(i, j, k) - u(i  , j, k-1))*dxInv[2]/met_h_zeta
                             + ( (w(i, j, k) - w(i-1, j, k  ))*dxInv[0]
                               - (met_h_xi/met_h_zeta)*GradWz ) * mf_u(i,j,0) );
        tau31(i,j,k) = tau13(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                        - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

        Real met_h_eta,met_h_zeta;
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        tau23(i,j,k) = 0.5 * ( (v(i, j, k) - v(i, j  , k-1))*dxInv[2]/met_h_zeta
                             + ( (w(i, j, k) - w(i, j-1, k  ))*dxInv[1]
                               - (met_h_eta/met_h_zeta)*GradWz ) * mf_v(i,j,0) );
        tau32(i,j,k) = tau23(i,j,k);
    });
}
