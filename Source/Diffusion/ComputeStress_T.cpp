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

    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i-1, j  , k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau13(i,j,k) += met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff;
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i  , j-1, k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            tau23(i,j,k) += met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i-1, j  , k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau13(i,j,k) += met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff;
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i  , j-1, k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            tau23(i,j,k) += met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff;
        });
    }

    // Standard operations
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau12(i,j,k) *= mu_eff;
        tau21(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  )
                             + er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );

        tau13(i,j,k) += met_h_xi*OneThird*errbar;
        tau13(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_eta;
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  )
                             + er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );

        tau23(i,j,k) += met_h_eta*OneThird*errbar;
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

    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i-1, j  , k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            Real varlo  = 0.5 * ( K_turb(i  , j  , k  , EddyDiff::Mom_v) + K_turb(i-1, j  , k  , EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k+1, EddyDiff::Mom_v) + K_turb(i-1, j  , k+1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varlo - 0.5*varhi;

            tau13(i,j,k) += met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff + varbar;

            tau31(i,j,k) *= mu_eff + varbar;
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  ) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k+1) + er_arr(i  , j-1, k+1) );
            Real errbar = 1.5*errlo - 0.5*errhi;

            Real varlo  = 0.5 * ( K_turb(i  , j  , k  , EddyDiff::Mom_v) + K_turb(i  , j-1, k  , EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k+1, EddyDiff::Mom_v) + K_turb(i  , j-1, k+1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varlo - 0.5*varhi;

            tau23(i,j,k) += met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff + varbar;

            tau32(i,j,k) *= mu_eff + varbar;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_xi;
            met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i-1, j  , k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            Real varlo  = 0.5 * ( K_turb(i  , j  , k-2, EddyDiff::Mom_v) + K_turb(i-1, j  , k-2, EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k-1, EddyDiff::Mom_v) + K_turb(i-1, j  , k-1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varhi - 0.5*varlo;

            tau13(i,j,k) += met_h_xi*OneThird*errbar;
            tau13(i,j,k) *= mu_eff + varbar;

            tau31(i,j,k) *= mu_eff + varbar;
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real met_h_eta;
            met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real errlo  = 0.5 * ( er_arr(i  , j  , k-2) + er_arr(i  , j-1, k-2) );
            Real errhi  = 0.5 * ( er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );
            Real errbar = 1.5*errhi - 0.5*errlo;

            Real varlo  = 0.5 * ( K_turb(i  , j  , k-2, EddyDiff::Mom_v) + K_turb(i  , j-1, k-2, EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k-1, EddyDiff::Mom_v) + K_turb(i  , j-1, k-1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varhi - 0.5*varlo;

            tau23(i,j,k) += met_h_eta*OneThird*errbar;
            tau23(i,j,k) *= mu_eff + varbar;

            tau32(i,j,k) *= mu_eff + varbar;
        });
    }

    // Standard operations
    amrex::ParallelFor(tbxxy, tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_turb = 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                            + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );
        tau12(i,j,k) *= mu_eff + mu_turb;
        tau21(i,j,k) *= mu_eff + mu_turb;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_xi;
        met_h_xi  = Compute_h_xi_AtEdgeCenterJ (i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i-1, j  , k  )
                             + er_arr(i  , j  , k-1) + er_arr(i-1, j  , k-1) );
        Real varbar = 0.25 * ( K_turb(i-1, j  , k  , EddyDiff::Mom_v) + K_turb(i  , j  , k  , EddyDiff::Mom_v)
                             + K_turb(i-1, j  , k-1, EddyDiff::Mom_v) + K_turb(i  , j  , k-1, EddyDiff::Mom_v) )

        tau13(i,j,k) += met_h_xi*OneThird*errbar;
        tau13(i,j,k) *= mu_eff + varbar;

        tau31(i,j,k) *= mu_eff + varbar;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real met_h_eta;
        met_h_eta = Compute_h_eta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real errbar = 0.25 * ( er_arr(i  , j  , k  ) + er_arr(i  , j-1, k  )
                             + er_arr(i  , j  , k-1) + er_arr(i  , j-1, k-1) );
        Real varbar = 0.25 * ( K_turb(i  , j-1, k  , EddyDiff::Mom_v) + K_turb(i  , j  , k  , EddyDiff::Mom_v)
                             + K_turb(i  , j-1, k-1, EddyDiff::Mom_v) + K_turb(i  , j  , k-1, EddyDiff::Mom_v) )

        tau23(i,j,k) += met_h_eta*OneThird*errbar;
        tau23(i,j,k) *= mu_eff + varbar;

        tau32(i,j,k) *= mu_eff + varbar;
    });
}
