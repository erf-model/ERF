#include <Diffusion.H>
#include <TerrainMetrics.H>

using namespace amrex;

void
ComputeStressConsVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau13,
                        Array4<Real>& tau21, Array4<Real>& tau23,
                        Array4<Real>& tau31, Array4<Real>& tau32,
                        const Array4<const Real>& er_arr,
                        const Array4<const Real>& z_nd  ,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D) = Tau
    //       The boxes are copied here for the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;  Box bxcc3  = bxcc;
                        Box tbxxy3 = tbxxy;
    Box tbxxz2 = tbxxz; Box tbxxz3 = tbxxz;
    Box tbxyz2 = tbxyz; Box tbxyz3 = tbxyz;

    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);

        tau11(i,j,k) -= met_h_zeta*OneThird*er_arr(i,j,k);
        tau11(i,j,k) *= mu_eff;
        
        tau22(i,j,k) -= met_h_zeta*OneThird*er_arr(i,j,k);
        tau22(i,j,k) *= mu_eff;
        
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) *= mu_eff;
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau12(i,j,k) *= mu_eff;
        tau21(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau31(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau32(i,j,k) *= mu_eff;
    });

    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz2 = tbxxz2; planexz2.setBig(2, planexz2.smallEnd(2) );
        tbxxz2.growLo(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i-1, j  , k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau13(i,j,k) -= met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff;
        });

        Box planeyz2 = tbxyz2; planeyz2.setBig(2, planeyz2.smallEnd(2) );
        tbxyz2.growLo(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i  , j-1, k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau23(i,j,k) -= met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz2 = tbxxz2; planexz2.setSmall(2, planexz2.bigEnd(2) );
        tbxxz2.growHi(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i-1, j  , k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau13(i,j,k) -= met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff;
        });

        Box planeyz2 = tbxyz2; planeyz2.setSmall(2, planeyz2.bigEnd(2) );
        tbxyz2.growHi(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i  , j-1, k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau23(i,j,k) -= met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff;
        });
    }
    
    // Standard operations
    amrex::ParallelFor(tbxxz2,tbxyz2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  )
                             + er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );

        tau13(i,j,k) -= met_h_xi*errbar;
        tau13(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_eta;
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  )
                             + er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );

        tau23(i,j,k) -= met_h_eta*OneThird*errbar;
        tau23(i,j,k) *= mu_eff;
    });
}


void
ComputeStressVarVisc_T(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                       const Array4<const Real>& K_turb,
                       Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                       Array4<Real>& tau12, Array4<Real>& tau13,
                       Array4<Real>& tau21, Array4<Real>& tau23,
                       Array4<Real>& tau31, Array4<Real>& tau32,
                       const Array4<const Real>& er_arr,
                       const Array4<const Real>& z_nd  ,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D) = Tau
    //       The boxes are copied here for the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;  Box bxcc3  = bxcc;
                        Box tbxxy3 = tbxxy;
    Box tbxxz2 = tbxxz; Box tbxxz3 = tbxxz;
    Box tbxyz2 = tbxyz; Box tbxyz3 = tbxyz;

    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_zeta;
        met_h_zeta = Compute_h_zeta_AtCellCenter(i,j,k,dxInv,z_nd);

        tau11(i,j,k) -= met_h_zeta*OneThird*er_arr(i,j,k);
        tau11(i,j,k) *= mu_eff + K_turb(i, j, k, EddyDiff::Mom_h);;
        
        tau22(i,j,k) -= met_h_zeta*OneThird*er_arr(i,j,k);
        tau22(i,j,k) *= mu_eff + K_turb(i, j, k, EddyDiff::Mom_h);;
        
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) *= mu_eff + K_turb(i, j, k, EddyDiff::Mom_v);;
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_turb = 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                            + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );
        tau12(i,j,k) *= mu_eff + mu_turb;
        tau21(i,j,k) *= mu_eff + mu_turb;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau31(i,j,k) *= mu_eff+ 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                     + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau32(i,j,k) *= mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                      + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
    });

    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz2 = tbxxz2; planexz2.setBig(2, planexz2.smallEnd(2) );
        tbxxz2.growLo(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i-1, j  , k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau13(i,j,k) -= met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                          + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        });

        Box planeyz2 = tbxyz2; planeyz2.setBig(2, planeyz2.smallEnd(2) );
        tbxyz2.growLo(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i  , j-1, k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau23(i,j,k) -= met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                          + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz2 = tbxxz2; planexz2.setSmall(2, planexz2.bigEnd(2) );
        tbxxz2.growHi(2,-1);
        amrex::ParallelFor(planexz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i-1, j  , k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau13(i,j,k) -= met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                          + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        });

        Box planeyz2 = tbxyz2; planeyz2.setSmall(2, planeyz2.bigEnd(2) );
        tbxyz2.growHi(2,-1);
        amrex::ParallelFor(planeyz2,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i  , j-1, k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau23(i,j,k) -= met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                          + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        });
    }
    
    // Standard operations
    amrex::ParallelFor(tbxxz2,tbxyz2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  )
                             + er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );

        tau13(i,j,k) -= met_h_xi*errbar;
        tau13(i,j,k) *= mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                      + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_eta;
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);
        
        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  )
                             + er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );

        tau23(i,j,k) -= met_h_eta*OneThird*errbar;
        tau23(i,j,k) *= mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                      + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
    });
}
