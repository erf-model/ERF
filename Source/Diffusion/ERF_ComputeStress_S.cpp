#include <ERF_Diffusion.H>
#include <ERF_TerrainMetrics.H>

using namespace amrex;

/**
 * Function for computing the stress with constant viscosity and with terrain.
 *
 * @param[in]  bxcc cell center box for tau_ii
 * @param[in]  tbxxy nodal xy box for tau_12
 * @param[in]  tbxxz nodal xz box for tau_13
 * @param[in]  tbxyz nodal yz box for tau_23
 * @param[in]  mu_eff constant molecular viscosity
 * @param[in]  cell_data to access rho if ConstantAlpha
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau21 21 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in,out] tau31 31 strain -> stress
 * @param[in,out] tau32 32 strain -> stress
 * @param[in]  er_arr expansion rate
 * @param[in]  z_nd nodal array of physical z heights
 * @param[in]  dxInv inverse cell size array
 */
void
ComputeStressConsVisc_S (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                         const Array4<const Real>& cell_data,
                         Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                         Array4<Real>& tau12, Array4<Real>& tau21,
                         Array4<Real>& tau13, Array4<Real>& tau31,
                         Array4<Real>& tau23, Array4<Real>& tau32,
                         const Array4<const Real>& er_arr,
                         const Gpu::DeviceVector<Real>& stretched_dz_d,
                         const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                         const Array4<const Real>& mf_mx,
                         const Array4<const Real>& mf_ux,
                         const Array4<const Real>& mf_vx,
                         const Array4<const Real>& mf_my,
                         const Array4<const Real>& mf_uy,
                         const Array4<const Real>& mf_vy)
{
    // Handle constant alpha case, in which the provided mu_eff is actually
    // "alpha" and the viscosity needs to be scaled by rho. This can be further
    // optimized with if statements below instead of creating a new FAB,
    // but this is implementation is cleaner.
    FArrayBox temp;
    Box gbx = bxcc; // Note: bxcc have been grown in x/y only.
    gbx.grow(IntVect(0,0,1));
    temp.resize(gbx,1, The_Async_Arena());
    Array4<Real> rhoAlpha = temp.array();
    if (cell_data) {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = cell_data(i, j, k, Rho_comp) * mu_eff;
        });
    } else {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = mu_eff;
        });
    }

    auto dz_ptr = stretched_dz_d.data();

    // First block: cell centered stresses
    //***********************************************************************************
    Real OneThird   = (1./3.);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real mu_tot     = rhoAlpha(i,j,k);

        tau11(i,j,k) = -mu_tot / mfy * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = -mu_tot / mfx * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = -mu_tot       * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Second block: off diagonal stresses
    //***********************************************************************************
    ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = 0.5 * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = 0.5 * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real mu_tot = 0.25*( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                           + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );

        tau12(i,j,k) *= -mu_tot / mfx;
        tau21(i,j,k) *= -mu_tot / mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfy = mf_uy(i,j,0);

        Real mu_tot = 0.25 * ( rhoAlpha(i-1, j  , k  ) + rhoAlpha(i  , j  , k  )
                             + rhoAlpha(i-1, j  , k-1) + rhoAlpha(i  , j  , k-1) );

        tau13(i,j,k) *= -mu_tot;
        tau31(i,j,k) *= -mu_tot / mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_vx(i,j,0);

        Real mu_tot = 0.25 * ( rhoAlpha(i  , j-1, k  ) + rhoAlpha(i  , j  , k  )
                             + rhoAlpha(i  , j-1, k-1) + rhoAlpha(i  , j  , k-1) );

        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot / mfx;
    });
}

/**
 * Function for computing the stress with constant viscosity and with terrain.
 *
 * @param[in]  bxcc cell center box for tau_ii
 * @param[in]  tbxxy nodal xy box for tau_12
 * @param[in]  tbxxz nodal xz box for tau_13
 * @param[in]  tbxyz nodal yz box for tau_23
 * @param[in]  mu_eff constant molecular viscosity
 * @param[in]  mu_turb variable turbulent viscosity
 * @param[in]  cell_data to access rho if ConstantAlpha
 * @param[in,out] tau11 11 strain -> stress
 * @param[in,out] tau22 22 strain -> stress
 * @param[in,out] tau33 33 strain -> stress
 * @param[in,out] tau12 12 strain -> stress
 * @param[in,out] tau13 13 strain -> stress
 * @param[in,out] tau21 21 strain -> stress
 * @param[in,out] tau23 23 strain -> stress
 * @param[in,out] tau31 31 strain -> stress
 * @param[in,out] tau32 32 strain -> stress
 * @param[in]  er_arr expansion rate
 * @param[in]  z_nd nodal array of physical z heights
 * @param[in]  dxInv inverse cell size array
 */
void
ComputeStressVarVisc_S (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                        const Array4<const Real>& mu_turb,
                        const Array4<const Real>& cell_data,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau21,
                        Array4<Real>& tau13, Array4<Real>& tau31,
                        Array4<Real>& tau23, Array4<Real>& tau32,
                        const Array4<const Real>& er_arr,
                        const Gpu::DeviceVector<Real>& stretched_dz_d,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                        const Array4<const Real>& mf_mx,
                        const Array4<const Real>& mf_ux,
                        const Array4<const Real>& mf_vx,
                        const Array4<const Real>& mf_my,
                        const Array4<const Real>& mf_uy,
                        const Array4<const Real>& mf_vy)
{
    // Handle constant alpha case, in which the provided mu_eff is actually
    // "alpha" and the viscosity needs to be scaled by rho. This can be further
    // optimized with if statements below instead of creating a new FAB,
    // but this is implementation is cleaner.
    FArrayBox temp;
    Box gbx = bxcc; // Note: bxcc have been grown in x/y only.
    gbx.grow(IntVect(0,0,1));
    temp.resize(gbx,1, The_Async_Arena());
    Array4<Real> rhoAlpha = temp.array();
    if (cell_data) {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = cell_data(i, j, k, Rho_comp) * mu_eff;
        });
    } else {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = mu_eff;
        });
    }

    auto dz_ptr = stretched_dz_d.data();

    // First block: cell centered stresses
    //***********************************************************************************
    Real OneThird   = (1./3.);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real mu_tot     = rhoAlpha(i,j,k) + 2.0*mu_turb(i, j, k, EddyDiff::Mom_h);

        tau11(i,j,k) = -mu_tot / mfy * ( tau11(i,j,k) - OneThird*er_arr(i,j,k) );
        tau22(i,j,k) = -mu_tot / mfx * ( tau22(i,j,k) - OneThird*er_arr(i,j,k) );
        tau33(i,j,k) = -mu_tot       * ( tau33(i,j,k) - OneThird*er_arr(i,j,k) );
    });

    // Second block: off diagonal stresses
    //***********************************************************************************
    ParallelFor(tbxxy,tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = 0.5 * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = 0.5 * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real mu_bar = 0.25*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                           + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real rhoAlpha_bar = 0.25*( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                                 + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau12(i,j,k) *= -mu_tot / mfx;
        tau21(i,j,k) *= -mu_tot / mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfy = mf_uy(i,j,0);

        Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                           + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = 0.25*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau13(i,j,k) *= -mu_tot;
        tau31(i,j,k) *= -mu_tot / mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_vx(i,j,0);

        Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                           + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = 0.25*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot / mfx;
    });
}
