#include <ERF_Diffusion.H>
#include <ERF_TerrainMetrics.H>

using namespace amrex;

/**
 * Function for computing the stress with constant viscosity and with stretched dz
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
 * @param[in]  stretched_dz_d  vertical mesh spacing
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_mx mapfactor in x-direction at cell centers
 * @param[in]  mf_ux mapfactor in x-direction at x-faces
 * @param[in]  mf_vx mapfactor in x-direction at y-faces
 * @param[in]  mf_my mapfactor in y-direction at cell centers
 * @param[in]  mf_uy mapfactor in y-direction at x-faces
 * @param[in]  mf_vy mapfactor in y-direction at y-faces
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
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = cell_data(i, j, k, Rho_comp) * mu_eff;
        });
    } else {
        ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = mu_eff;
        });
    }

    auto dz_ptr = stretched_dz_d.data();

    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;
    bxcc2.grow(IntVect(-1,-1,0));

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (1./3.);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_tot   = rhoAlpha(i,j,k);

        tau33(i,j,k) *= -mu_tot;
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_tot = 0.25*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );

            tau13(i,j,k) *= -mu_tot;
            tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vy(i,j,0);
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

            Real mu_tot = 0.25*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );

            tau23(i,j,k) *= -mu_tot;
            tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_tot = 0.25*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );

            tau13(i,j,k) *= -mu_tot;
            tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vy(i,j,0);
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_tot = 0.25*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );

            tau23(i,j,k) *= -mu_tot;
            tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
        });
    }

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
    //-----------------------------------------------------------------------------------
    ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_tot = 0.25 * ( rhoAlpha(i-1, j  , k  ) + rhoAlpha(i  , j  , k  )
                             + rhoAlpha(i-1, j  , k-1) + rhoAlpha(i  , j  , k-1) );

        tau13(i,j,k) *= -mu_tot;
        tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_uy(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_tot = 0.25 * ( rhoAlpha(i  , j-1, k  ) + rhoAlpha(i  , j  , k  )
                             + rhoAlpha(i  , j-1, k-1) + rhoAlpha(i  , j  , k-1) );

        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_tot = rhoAlpha(i,j,k);

        tau11(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_my(i,j,0);
        tau22(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_mx(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = 0.5 * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = 0.5 * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real dz0 = 0.5 * (dz_ptr[k] + dz_ptr[k-1]) * dxInv[2];

        Real mu_tot = 0.25*( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                           + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );

        tau12(i,j,k) *= -mu_tot*dz0/mfx;
        tau21(i,j,k) *= -mu_tot*dz0/mfy;
    });
}

/**
 * Function for computing the stress with variable viscosity and with stretched dz
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
 * @param[in]  stretched_dz_d  vertical mesh spacing
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_mx mapfactor in x-direction at cell centers
 * @param[in]  mf_ux mapfactor in x-direction at x-faces
 * @param[in]  mf_vx mapfactor in x-direction at y-faces
 * @param[in]  mf_my mapfactor in y-direction at cell centers
 * @param[in]  mf_uy mapfactor in y-direction at x-faces
 * @param[in]  mf_vy mapfactor in y-direction at y-faces
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

    //***********************************************************************************
    // NOTE: The first  block computes (S-D).
    //       The second block computes 2mu*JT*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;
    bxcc2.grow(IntVect(-1,-1,0));

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (1./3.);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_tot   = rhoAlpha(i,j,k) + 2.0*mu_turb(i, j, k, EddyDiff::Mom_v);

        tau33(i,j,k) *= -mu_tot;
    });

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = 0.25*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

            tau13(i,j,k) *= -mu_tot;

            tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_uy(i,j,0);
        });

        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = 0.25*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

            tau23(i,j,k) *= -mu_tot;

            tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        ParallelFor(planexz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_bar = 0.25*( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = 0.25*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

            tau13(i,j,k) *= -mu_tot;

            tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_uy(i,j,0);
        });

        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planeyz,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mu_bar = 0.25*( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = 0.25*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

            tau23(i,j,k) *= -mu_tot;

            tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
        });
    }

    // Second block: compute 2mu*JT*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
    //-----------------------------------------------------------------------------------
    ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_bar = 0.25 * ( mu_turb(i-1, j  , k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                             + mu_turb(i-1, j  , k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = 0.25 * ( rhoAlpha(i-1, j  , k  ) + rhoAlpha(i  , j  , k  )
                                   + rhoAlpha(i-1, j  , k-1) + rhoAlpha(i  , j  , k-1) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau13(i,j,k) *= -mu_tot;

        tau31(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_uy(i,j,0);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_bar = 0.25 * ( mu_turb(i  , j-1, k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                             + mu_turb(i  , j-1, k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = 0.25 * ( rhoAlpha(i  , j-1, k  ) + rhoAlpha(i  , j  , k  )
                                   + rhoAlpha(i  , j-1, k-1) + rhoAlpha(i  , j  , k-1) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mf_vx(i,j,0);
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real mu_tot = rhoAlpha(i,j,k) + 2.0*mu_turb(i, j, k, EddyDiff::Mom_h);

        tau11(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mfy;
        tau22(i,j,k) *= -mu_tot*dz_ptr[k]*dxInv[2]/mfx;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = 0.5 * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = 0.5 * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real dz0 = 0.5 * (dz_ptr[k] + dz_ptr[k-1]) * dxInv[2];

        Real mu_bar = 0.25*( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                           + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real rhoAlpha_bar = 0.25*( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                                 + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        Real mu_tot = rhoAlpha_bar + 2.0*mu_bar;

        tau12(i,j,k) *= -mu_tot*dz0/mfx;
        tau21(i,j,k) *= -mu_tot*dz0/mfy;
    });
}
