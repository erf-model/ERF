#include <Diffusion.H>

using namespace amrex;

void
ComputeStressConsVisc_N(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                        const Array4<const Real>& er_arr,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        tau11(i,j,k) = mu_eff * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = mu_eff * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = mu_eff * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau12(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau13(i,j,k) *= mu_eff;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau23(i,j,k) *= mu_eff;
    });

}


void
ComputeStressVarVisc_N(Box& bxcc, Box& tbxxy, Box& tbxxz, Box& tbxyz, Real mu_eff,
                       const Array4<const Real>& K_turb,
                       Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                       Array4<Real>& tau12, Array4<Real>& tau13, Array4<Real>& tau23,
                       const Array4<const Real>& er_arr,
                       const GpuArray<Real, AMREX_SPACEDIM>& dxInv)
{
    Real OneThird   = (1./3.);

    // Cell centered strains
    amrex::ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_11 = mu_eff + K_turb(i, j, k, EddyDiff::Mom_h);
        Real mu_22 = mu_11;
        Real mu_33 = mu_eff + K_turb(i, j, k, EddyDiff::Mom_v);
        tau11(i,j,k) = mu_11 * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = mu_22 * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = mu_33 * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Extrapolate K_turb to top and bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real varlo  = 0.5 * ( K_turb(i  , j  , k  , EddyDiff::Mom_v) + K_turb(i-1, j  , k  , EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k+1, EddyDiff::Mom_v) + K_turb(i-1, j  , k+1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varlo - 0.5*varhi;

            tau13(i,j,k) *= varbar;
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real varlo  = 0.5 * ( K_turb(i  , j  , k  , EddyDiff::Mom_v) + K_turb(i  , j-1, k  , EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k+1, EddyDiff::Mom_v) + K_turb(i  , j-1, k+1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varlo - 0.5*varhi;

            tau23(i,j,k) *= varbar;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        amrex::ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real varlo  = 0.5 * ( K_turb(i  , j  , k-2, EddyDiff::Mom_v) + K_turb(i-1, j  , k-2, EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k-1, EddyDiff::Mom_v) + K_turb(i-1, j  , k-1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varhi - 0.5*varlo;

            tau13(i,j,k) *= varbar;
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        amrex::ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            Real varlo  = 0.5 * ( K_turb(i  , j  , k-2, EddyDiff::Mom_v) + K_turb(i  , j-1, k-2, EddyDiff::Mom_v) );
            Real varhi  = 0.5 * ( K_turb(i  , j  , k-1, EddyDiff::Mom_v) + K_turb(i  , j-1, k-1, EddyDiff::Mom_v) );
            Real varbar = 1.5*varhi - 0.5*varlo;

            tau23(i,j,k) *= varbar;
        });
    }

    // Off-diagonal strains
    amrex::ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_12 = mu_eff + 0.25*( K_turb(i-1, j  , k, EddyDiff::Mom_h) + K_turb(i, j  , k, EddyDiff::Mom_h)
                                   + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );
        tau12(i,j,k) *= mu_12;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_13 = mu_eff + 0.25*( K_turb(i-1, j, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                   + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        tau13(i,j,k) *= mu_13;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_23 = mu_eff + 0.25*( K_turb(i, j-1, k  , EddyDiff::Mom_v) + K_turb(i, j, k  , EddyDiff::Mom_v)
                                   + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
        tau23(i,j,k) *= mu_23;
    });
}
