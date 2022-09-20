#include <Diffusion.H>
#include <TerrainMetrics.H>

using namespace amrex;

void
ComputeStressConsVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                        const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau13,
                        Array4<Real>& tau21, Array4<Real>& tau23,
                        Array4<Real>& tau31, Array4<Real>& tau32,
                        const Array4<const Real>& er_arr,
                        const Array4<const Real>& z_nd  ,
                        const BCRec* bc_ptr, const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{    
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D) = Tau
    //       The boxes are copied here for the second block operations
    //***********************************************************************************   
    Box bxcc2  = bxcc;  Box bxcc3  = bxcc;
    Box tbxxy2 = tbxxy; Box tbxxy3 = tbxxy;
    Box tbxxz2 = tbxxz; Box tbxxz3 = tbxxz;
    Box tbxyz2 = tbxyz; Box tbxyz3 = tbxyz;
    
    Real OneThird   = (1./3.);

    // Dirichlet on left or right plane
    bool xl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
    bool xh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );

    bool xl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
    bool xh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );

    // Dirichlet on front or back plane
    bool yl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
    bool yh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );

    bool yl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
    bool yh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );

    // Dirichlet on top or bottom plane
    bool zl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
    bool zh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );

    bool zl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
    bool zh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );

    //***********************************************************************************
    // First Step: (S-D)
    //***********************************************************************************

    // X-Dirichlet
    //***********************************************************************************
    if (xl_v_dir) {
        Box planexy = tbxxy; planexy.setBig(0, planexy.smallEnd(0) );
        tbxxy.growLo(0,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                          + (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) - (1./3.) * v(i+1,j,k))*dxInv[0]
                                          - (met_h_eta/met_h_zeta)*GradUz
                                          - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (xh_v_dir) {
        Box planexy = tbxxy; planexy.setSmall(0, planexy.bigEnd(0) );
        tbxxy.growHi(0,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                          - (-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) - (1./3.) * v(i-2,j,k))*dxInv[0]
                                          - (met_h_eta/met_h_zeta)*GradUz
                                          - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (xl_w_dir) {
        Box planexz = tbxxz; planexz.setBig(0, planexz.smallEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growLo(0,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                          + (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dxInv[0]
                                          - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (xh_w_dir) {
        Box planexz = tbxxz; planexz.setSmall(0, planexz.bigEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growHi(0,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
            tau13(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                          - (-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dxInv[0]
                                          - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    // Y-Dirichlet
    //***********************************************************************************
    if (yl_u_dir) {
        Box planexy = tbxxy; planexy.setBig(1, planexy.smallEnd(1) );
        tbxxy.growLo(1,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
        
            tau12(i,j,k) = 0.5 * mu_eff * ( (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) - (1./3.) * u(i,j+1,k))*dxInv[1]
                                          + (v(i, j, k) - v(i-1, j, k))*dxInv[0]
                                          - (met_h_eta/met_h_zeta)*GradUz
                                          - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (yh_u_dir) {
        Box planexy = tbxxy; planexy.setSmall(1, planexy.bigEnd(1) );
        tbxxy.growHi(1,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_eff * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-2,k))*dxInv[1] +
                                           + (v(i, j, k) - v(i-1, j, k))*dxInv[0]
                                           - (met_h_eta/met_h_zeta)*GradUz
                                           - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (yl_w_dir) {
        Box planeyz = tbxyz; planeyz.setBig(1, planeyz.smallEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growLo(1,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );
            
             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
             
            tau23(i,j,k) = 0.5 * mu_eff * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                          + (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dxInv[1]
                                          - (met_h_eta/met_h_zeta)*GradWz);
            tau32(i,j,k) = tau23(i,j,k);
        });
    }
    if (yh_w_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(1, planeyz.bigEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growHi(1,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
             Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                             - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
             
            tau23(i,j,k) = 0.5 * mu_eff * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                          - (-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dxInv[1]
                                          - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-Dirichlet
    //***********************************************************************************
    if (zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_eff * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i-1, j, k))*dxInv[0]
                                          - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_eff * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dxInv[2]/met_h_zeta
                                           + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    if (zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );

            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_eff * ( (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i, j-1, k))*dxInv[1]
                                          - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }     
    if (zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_eff * ( -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dxInv[2]/met_h_zeta
                                           + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-lo w/out Z-Dirichlet (GradWz extrapolation)
    //***********************************************************************************
    if (!zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );
            
            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_eff * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i-1, j, k))*dxInv[0]
                                          - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );
            
            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
            tau23(i,j,k) = 0.5 * mu_eff * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i, j-1, k))*dxInv[1]
                                          - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-hi w/out Z-Dirichlet (h_xi = h_eta = 0)
    //***********************************************************************************
    if (!zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_eff * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Fill the remaining cells
    //***********************************************************************************
    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i-1,j  ,k+1)
                                         -u(i  ,j  ,k-1) - u(i-1,j  ,k-1) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i  ,j-1,k+1)
                                         -v(i  ,j  ,k-1) - v(i  ,j-1,k-1) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);
        
        tau11(i,j,k) = mu_eff * ( (u(i+1, j, k) - u(i, j, k))*dxInv[0]
                               - (met_h_xi/met_h_zeta) * GradUz
                               - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = mu_eff * ( (v(i, j+1, k) - v(i, j, k))*dxInv[1]
                               - (met_h_eta/met_h_zeta) * GradVz
                               - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = mu_eff * ( (w(i, j, k+1) - w(i, j, k))*dxInv[2]/met_h_zeta
                               - OneThird*er_arr(i,j,k) );
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                         -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                         -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
  
        tau12(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i  , j-1, k))*dxInv[1]
                                      + (v(i, j, k) - v(i-1, j  , k))*dxInv[0]
                                      - (met_h_eta/met_h_zeta)*GradUz
                                      - (met_h_xi /met_h_zeta)*GradVz );
        tau21(i,j,k) = tau12(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                        - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

        Real met_h_xi,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
        tau13(i,j,k) = 0.5 * mu_eff * ( (u(i, j, k) - u(i  , j, k-1))*dxInv[2]/met_h_zeta
                                      + (w(i, j, k) - w(i-1, j, k  ))*dxInv[0]
                                      - (met_h_xi/met_h_zeta)*GradWz );
        tau31(i,j,k) = tau13(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                        - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

        Real met_h_eta,met_h_zeta;
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        tau23(i,j,k) = 0.5 * mu_eff * ( (v(i, j, k) - v(i, j  , k-1))*dxInv[2]/met_h_zeta
                                      + (w(i, j, k) - w(i, j-1, k  ))*dxInv[1]
                                      - (met_h_eta/met_h_zeta)*GradWz );
        tau32(i,j,k) = tau23(i,j,k);
    });
    

    //***********************************************************************************
    // Second Step: JT*2mu*(S-D) = Tau
    //***********************************************************************************
    
    // Must fill  tau33, tau13, tau23 first (linear combinations)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 at bottom
    {
        Box planexz2 = tbxxz2; planexz2.setBig(2, planexz2.smallEnd(2) );
        tbxxz2.growLo(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = 1.5*tau11lo - 0.5*tau11hi;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = 1.5*tau12lo - 0.5*tau12hi;
            
            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        });

        Box planeyz2 = tbxyz2; planeyz2.setBig(2, planeyz2.smallEnd(2) );
        tbxyz2.growLo(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = 1.5*tau21lo - 0.5*tau21hi;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = 1.5*tau22lo - 0.5*tau22hi;
            
            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        });
    }
    // Extrapolate tau13 & tau23 at top
    {
        Box planexz2 = tbxxz2; planexz2.setSmall(2, planexz2.bigEnd(2) );
        tbxxz2.growHi(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = 1.5*tau11hi - 0.5*tau11lo;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = 1.5*tau12hi - 0.5*tau12lo;
            
            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        });
        
        Box planeyz2 = tbxyz2; planeyz2.setSmall(2, planeyz2.bigEnd(2) );
        tbxyz2.growHi(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = 1.5*tau21hi - 0.5*tau21lo;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = 1.5*tau22hi - 0.5*tau22lo;
            
            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        });
    }

    // We don't need x/y ghost cells for tau33 (avoids linear comb issues)
    bxcc2.growLo(0,-1); bxcc2.growLo(1,-1);
    bxcc2.growHi(0,-1); bxcc2.growHi(1,-1);

    // Standard operations
    amrex::ParallelFor(bxcc2,tbxxz2,tbxyz2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {                   
        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = 0.25 * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                               + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = 0.25 * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                               + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );

        tau33(i,j,k) -= met_h_xi*tau31bar + met_h_eta*tau32bar;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi,met_h_eta;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = 0.25 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                               + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                               + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );

        tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi,met_h_eta;
        met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = 0.25 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                               + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = 0.25 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                               + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
        
        tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
    });
    
    // Fill the remaining components (h_zeta multiplication)
    //-----------------------------------------------------------------------------------
    // Cell centered strains
    amrex::ParallelFor(bxcc3, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);
        
        tau11(i,j,k) *= met_h_zeta; 
        tau22(i,j,k) *= met_h_zeta;
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy3,tbxxz3,tbxyz3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        tau12(i,j,k) *= met_h_zeta;
        tau21(i,j,k) *= met_h_zeta; 
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
        tau31(i,j,k) *= met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        tau32(i,j,k) *= met_h_zeta;
    });
    
    /*
    // AML DEBUG
    Box tmpcc = bxcc3;
    tmpcc.growLo(1,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC1 test: " << IntVect(i,j,k) << ' '
                       << tau11(i,j,k) << "\n";
    });
    amrex::Print() << "CC1 CLEARED" << "\n";
    amrex::Print() << "\n";
    
    tmpcc = bxcc3;
    tmpcc.growLo(0,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC2 test: " << IntVect(i,j,k) << ' '
                       << tau22(i,j,k) << "\n";
    });
    amrex::Print() << "CC2 CLEARED" << "\n";
    amrex::Print() << "\n";

    tmpcc = bxcc3;
    tmpcc.growLo(0,-1); tmpcc.growLo(1,-1);
    tmpcc.growHi(0,-1); tmpcc.growHi(1,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC3 test: " << IntVect(i,j,k) << ' '
                       << tau33(i,j,k) << "\n";
    });
    amrex::Print() << "CC3 CLEARED" << "\n";
    amrex::Print() << "\n";

    amrex::ParallelFor(tbxxy3,tbxxz3,tbxyz3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "12 test: " << IntVect(i,j,k) << ' '
                       << tau12(i,j,k) << ' '
                       << tau21(i,j,k) <<"\n";
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "13 test: " << IntVect(i,j,k) << ' '
                       << tau13(i,j,k) << ' '
                       << tau31(i,j,k) << "\n";
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "23 test: " << IntVect(i,j,k) << ' '
                       << tau23(i,j,k) << ' '
                       << tau32(i,j,k) << "\n";
    });
    amrex::Print() << "ALL CLEARED" << "\n";
    amrex::Print() << "\n";
    exit(0);
    */
    
}


void
ComputeStressVarVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                       const Array4<const Real>& K_turb,
                       const Array4<const Real>& u, const Array4<const Real>& v, const Array4<const Real>& w,
                       Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                       Array4<Real>& tau12, Array4<Real>& tau13,
                       Array4<Real>& tau21, Array4<Real>& tau23,
                       Array4<Real>& tau31, Array4<Real>& tau32,
                       const Array4<const Real>& er_arr,
                       const Array4<const Real>& z_nd  ,
                       const BCRec* bc_ptr, const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D) = Tau
    //       The boxes are copied here for the second block operations
    //***********************************************************************************   
    Box bxcc2  = bxcc;  Box bxcc3  = bxcc;
    Box tbxxy2 = tbxxy; Box tbxxy3 = tbxxy;
    Box tbxxz2 = tbxxz; Box tbxxz3 = tbxxz;
    Box tbxyz2 = tbxyz; Box tbxyz3 = tbxyz;
    
    Real OneThird   = (1./3.);

    // Dirichlet on left or right plane
    bool xl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
    bool xh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );

    bool xl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(0) == ERFBCType::ext_dir_ingested) );
    bool xh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(0) == ERFBCType::ext_dir_ingested) );

    // Dirichlet on front or back plane
    bool yl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
    bool yh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );

    bool yl_w_dir = ( (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].lo(1) == ERFBCType::ext_dir_ingested) );
    bool yh_w_dir = ( (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::zvel_bc].hi(1) == ERFBCType::ext_dir_ingested) );

    // Dirichlet on top or bottom plane
    bool zl_u_dir = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
    bool zh_u_dir = ( (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::xvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );

    bool zl_v_dir = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir_ingested) );
    bool zh_v_dir = ( (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir)          ||
                      (bc_ptr[BCVars::yvel_bc].hi(2) == ERFBCType::ext_dir_ingested) );

    //***********************************************************************************
    // First Step: (S-D)
    //***********************************************************************************

    // X-Dirichlet
    //***********************************************************************************
    if (xl_v_dir) {
        Box planexy = tbxxy; planexy.setBig(0, planexy.smallEnd(0) );
        tbxxy.growLo(0,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                       + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );

            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_12 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                         + (-(8./3.) * v(i-1,j,k) + 3. * v(i,j,k) - (1./3.) * v(i+1,j,k))*dxInv[0]
                                         - (met_h_eta/met_h_zeta)*GradUz
                                         - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (xh_v_dir) {
        Box planexy = tbxxy; planexy.setSmall(0, planexy.bigEnd(0) );
        tbxxy.growHi(0,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                       + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );

            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_12 * ( (u(i, j, k) - u(i, j-1, k))*dxInv[1]
                                         - (-(8./3.) * v(i,j,k) + 3. * v(i-1,j,k) - (1./3.) * v(i-2,j,k))*dxInv[0]
                                         - (met_h_eta/met_h_zeta)*GradUz
                                         - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (xl_w_dir) {
        Box planexz = tbxxz; planexz.setBig(0, planexz.smallEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growLo(0,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_13 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                         + (-(8./3.) * w(i-1,j,k) + 3. * w(i,j,k) - (1./3.) * w(i+1,j,k))*dxInv[0]
                                         - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (xh_w_dir) {
        Box planexz = tbxxz; planexz.setSmall(0, planexz.bigEnd(0) );
        planexz.setSmall(2, planexz.smallEnd(2)+1 ); planexz.setBig(2, planexz.bigEnd(2)-1 );
        tbxxz.growHi(0,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
            tau13(i,j,k) = 0.5 * mu_13 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                         - (-(8./3.) * w(i,j,k) + 3. * w(i-1,j,k) - (1./3.) * w(i-2,j,k))*dxInv[0]
                                         - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    // Y-Dirichlet
    //***********************************************************************************
    if (yl_u_dir) {
        Box planexy = tbxxy; planexy.setBig(1, planexy.smallEnd(1) );
        tbxxy.growLo(1,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                       + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );

            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
        
            tau12(i,j,k) = 0.5 * mu_12 * ( (-(8./3.) * u(i,j-1,k) + 3. * u(i,j,k) - (1./3.) * u(i,j+1,k))*dxInv[1]
                                         + (v(i, j, k) - v(i-1, j, k))*dxInv[0]
                                         - (met_h_eta/met_h_zeta)*GradUz
                                         - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }
    if (yh_u_dir) {
        Box planexy = tbxxy; planexy.setSmall(1, planexy.bigEnd(1) );
        tbxxy.growHi(1,-1);
        amrex::ParallelFor(planexy,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                       + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );

            Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                             -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
            Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                             -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
            
            tau12(i,j,k) = 0.5 * mu_12 * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j-1,k) - (1./3.) * u(i,j-2,k))*dxInv[1] +
                                          + (v(i, j, k) - v(i-1, j, k))*dxInv[0]
                                          - (met_h_eta/met_h_zeta)*GradUz
                                          - (met_h_xi /met_h_zeta)*GradVz );
            tau21(i,j,k) = tau12(i,j,k);
        });
    }

    if (yl_w_dir) {
        Box planeyz = tbxyz; planeyz.setBig(1, planeyz.smallEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growLo(1,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );
            
             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
             
            tau23(i,j,k) = 0.5 * mu_23 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                         + (-(8./3.) * w(i,j-1,k) + 3. * w(i,j  ,k) - (1./3.) * w(i,j+1,k))*dxInv[1]
                                         - (met_h_eta/met_h_zeta)*GradWz);
            tau32(i,j,k) = tau23(i,j,k);
        });
    }
    if (yh_w_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(1, planeyz.bigEnd(1) );
        planeyz.setSmall(2, planeyz.smallEnd(2)+1 ); planeyz.setBig(2, planeyz.bigEnd(2)-1 );
        tbxyz.growHi(1,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

             Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                             - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

             Real met_h_eta,met_h_zeta;
             met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
             met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
             
            tau23(i,j,k) = 0.5 * mu_23 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                         - (-(8./3.) * w(i,j  ,k) + 3. * w(i,j-1,k) - (1./3.) * w(i,j-2,k))*dxInv[1]
                                         - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-Dirichlet
    //***********************************************************************************
    if (zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_13 * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i-1, j, k))*dxInv[0]
                                         - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_13 * ( -(-(8./3.) * u(i,j,k) + 3. * u(i,j,k-1) - (1./3.) * u(i,j,k-2))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }

    if (zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );

            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_23 * ( (-(8./3.) * v(i,j,k-1) + 3. * v(i,j,k  ) - (1./3.) * v(i,j,k+1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i, j-1, k))*dxInv[1]
                                         - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }     
    if (zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_23 * ( -(-(8./3.) * v(i,j,k  ) + 3. * v(i,j,k-1) - (1./3.) * v(i,j,k-2))*dxInv[2]/met_h_zeta
                                          + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-lo w/out Z-Dirichlet (GradWz extrapolation)
    //***********************************************************************************
    if (!zl_u_dir) {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                            - w(i  ,j  ,k  ) - w(i-1,j  ,k  ) );
            
            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_13 * ( (-(8./3.) * u(i,j,k-1) + 3. * u(i,j,k) - (1./3.) * u(i,j,k+1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i-1, j, k))*dxInv[0]
                                         - (met_h_xi/met_h_zeta)*GradWz );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zl_v_dir) {
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real GradWz = 0.5  * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                            - w(i  ,j  ,k  ) - w(i  ,j-1,k  ) );
            
            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
            tau23(i,j,k) = 0.5 * mu_23 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i, j-1, k))*dxInv[1]
                                         - (met_h_eta/met_h_zeta)*GradWz );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Z-hi w/out Z-Dirichlet (h_xi = h_eta = 0)
    //***********************************************************************************
    if (!zh_u_dir) {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
            
            tau13(i,j,k) = 0.5 * mu_13 * ( (u(i, j, k) - u(i, j, k-1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i-1, j, k))*dxInv[0] );
            tau31(i,j,k) = tau13(i,j,k);
        });
    }
    if (!zh_v_dir) {
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                       + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
            
            tau23(i,j,k) = 0.5 * mu_23 * ( (v(i, j, k) - v(i, j, k-1))*dxInv[2]/met_h_zeta
                                         + (w(i, j, k) - w(i, j-1, k))*dxInv[1] );
            tau32(i,j,k) = tau23(i,j,k);
        });
    }

    // Fill the remaining cells
    //***********************************************************************************
    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_11 = mu_eff + K_turb(i, j, k, EddyDiff::Mom_h);
        Real mu_22 = mu_11;
        Real mu_33 = mu_eff + K_turb(i, j, k, EddyDiff::Mom_v);

        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i-1,j  ,k+1)
                                         -u(i  ,j  ,k-1) - u(i-1,j  ,k-1) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i  ,j-1,k+1)
                                         -v(i  ,j  ,k-1) - v(i  ,j-1,k-1) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);
        
        tau11(i,j,k) = mu_11 * ( (u(i+1, j, k) - u(i, j, k))*dxInv[0]
                               - (met_h_xi/met_h_zeta) * GradUz
                               - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = mu_22 * ( (v(i, j+1, k) - v(i, j, k))*dxInv[1]
                               - (met_h_eta/met_h_zeta) * GradVz
                               - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = mu_33 * ( (w(i, j, k+1) - w(i, j, k))*dxInv[2]/met_h_zeta
                               - OneThird*er_arr(i,j,k) );
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                   + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );

        Real GradUz = 0.25 * dxInv[2] * ( u(i  ,j  ,k+1) + u(i  ,j-1,k+1)
                                         -u(i  ,j  ,k-1) - u(i  ,j-1,k-1) );
        Real GradVz = 0.25 * dxInv[2] * ( v(i  ,j  ,k+1) + v(i-1,j  ,k+1)
                                         -v(i  ,j  ,k-1) - v(i-1,j  ,k-1) );

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterK  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterK (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);
  
        tau12(i,j,k) = 0.5 * mu_12 * ( (u(i, j, k) - u(i  , j-1, k))*dxInv[1]
                                     + (v(i, j, k) - v(i-1, j  , k))*dxInv[0]
                                     - (met_h_eta/met_h_zeta)*GradUz
                                     - (met_h_xi /met_h_zeta)*GradVz );
        tau21(i,j,k) = tau12(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                   + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i-1,j  ,k+1)
                                        - w(i  ,j  ,k-1) - w(i-1,j  ,k-1) );

        Real met_h_xi,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
        tau13(i,j,k) = 0.5 * mu_13 * ( (u(i, j, k) - u(i  , j, k-1))*dxInv[2]/met_h_zeta
                                     + (w(i, j, k) - w(i-1, j, k  ))*dxInv[0]
                                     - (met_h_xi/met_h_zeta)*GradWz );
        tau31(i,j,k) = tau13(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                   + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );

        Real GradWz = 0.25 * dxInv[2] * ( w(i  ,j  ,k+1) + w(i  ,j-1,k+1)
                                        - w(i  ,j  ,k-1) - w(i  ,j-1,k-1) );

        Real met_h_eta,met_h_zeta;
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        tau23(i,j,k) = 0.5 * mu_23 * ( (v(i, j, k) - v(i, j  , k-1))*dxInv[2]/met_h_zeta
                                     + (w(i, j, k) - w(i, j-1, k  ))*dxInv[1]
                                     - (met_h_eta/met_h_zeta)*GradWz );
        tau32(i,j,k) = tau23(i,j,k);
    });

    //***********************************************************************************
    // Second Step: JT*2mu*(S-D) = Tau
    //***********************************************************************************
    
    // Must fill  tau33, tau13, tau23 first (linear combinations)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 at bottom
    {
        Box planexz2 = tbxxz2; planexz2.setBig(2, planexz2.smallEnd(2) );
        tbxxz2.growLo(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = 1.5*tau11lo - 0.5*tau11hi;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = 1.5*tau12lo - 0.5*tau12hi;
            
            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        });

        Box planeyz2 = tbxyz2; planeyz2.setBig(2, planeyz2.smallEnd(2) );
        tbxyz2.growLo(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = 1.5*tau21lo - 0.5*tau21hi;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = 1.5*tau22lo - 0.5*tau22hi;
            
            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        });
    }
    // Extrapolate tau13 & tau23 at top
    {
        Box planexz2 = tbxxz2; planexz2.setSmall(2, planexz2.bigEnd(2) );
        tbxxz2.growHi(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = 0.5 * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = 0.5 * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = 1.5*tau11hi - 0.5*tau11lo;

            Real tau12lo  = 0.5 * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = 0.5 * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = 1.5*tau12hi - 0.5*tau12lo;
            
            tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
        });
        
        Box planeyz2 = tbxyz2; planeyz2.setSmall(2, planeyz2.bigEnd(2) );
        tbxyz2.growHi(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi,met_h_eta;
            met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = 0.5 * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = 0.5 * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = 1.5*tau21hi - 0.5*tau21lo;

            Real tau22lo  = 0.5 * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = 0.5 * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = 1.5*tau22hi - 0.5*tau22lo;
            
            tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
        });
    }

    // We don't need x/y ghost cells for tau33 (avoids linear comb issues)
    bxcc2.growLo(0,-1); bxcc2.growLo(1,-1);
    bxcc2.growHi(0,-1); bxcc2.growHi(1,-1);

    // Standard operations
    amrex::ParallelFor(bxcc2,tbxxz2,tbxyz2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {                   
        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = 0.25 * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                               + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = 0.25 * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                               + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );

        tau33(i,j,k) -= met_h_xi*tau31bar + met_h_eta*tau32bar;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi,met_h_eta;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_eta = Compute_h_eta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = 0.25 * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                               + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                               + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );

        tau13(i,j,k) -= met_h_xi*tau11bar + met_h_eta*tau12bar;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi,met_h_eta;
        met_h_xi  = Compute_h_xi_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = 0.25 * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                               + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = 0.25 * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                               + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
        
        tau23(i,j,k) -= met_h_xi*tau21bar + met_h_eta*tau22bar;
    });
    
    // Fill the remaining components (h_zeta multiplication)
    //-----------------------------------------------------------------------------------
    // Cell centered strains
    amrex::ParallelFor(bxcc3, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);
        
        tau11(i,j,k) *= met_h_zeta; 
        tau22(i,j,k) *= met_h_zeta;
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy3,tbxxz3,tbxyz3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        tau12(i,j,k) *= met_h_zeta;
        tau21(i,j,k) *= met_h_zeta; 
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);
        
        tau31(i,j,k) *= met_h_zeta;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        tau32(i,j,k) *= met_h_zeta;
    });

    /*
    // AML DEBUG
    Box tmpcc = bxcc3;
    tmpcc.growLo(1,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC1 test: " << IntVect(i,j,k) << ' '
                       << tau11(i,j,k) << "\n";
    });
    amrex::Print() << "CC1 CLEARED" << "\n";
    amrex::Print() << "\n";
    
    tmpcc = bxcc3;
    tmpcc.growLo(0,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC2 test: " << IntVect(i,j,k) << ' '
                       << tau22(i,j,k) << "\n";
    });
    amrex::Print() << "CC2 CLEARED" << "\n";
    amrex::Print() << "\n";

    tmpcc = bxcc3;
    tmpcc.growLo(0,-1); tmpcc.growLo(1,-1);
    tmpcc.growHi(0,-1); tmpcc.growHi(1,-1);
    amrex::ParallelFor(tmpcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "CC3 test: " << IntVect(i,j,k) << ' '
                       << tau33(i,j,k) << "\n";
    });
    amrex::Print() << "CC3 CLEARED" << "\n";
    amrex::Print() << "\n";

    amrex::ParallelFor(tbxxy3,tbxxz3,tbxyz3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "12 test: " << IntVect(i,j,k) << ' '
                       << tau12(i,j,k) << ' '
                       << tau21(i,j,k) <<"\n";
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "13 test: " << IntVect(i,j,k) << ' '
                       << tau13(i,j,k) << ' '
                       << tau31(i,j,k) << "\n";
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        amrex::Print() << "23 test: " << IntVect(i,j,k) << ' '
                       << tau23(i,j,k) << ' '
                       << tau32(i,j,k) << "\n";
    });
    amrex::Print() << "ALL CLEARED" << "\n";
    amrex::Print() << "\n";
    exit(0);
    */
}
