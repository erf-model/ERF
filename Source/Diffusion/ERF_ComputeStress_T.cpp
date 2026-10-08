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
 * @param[in]  detJ Jacobian determinant
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[in,out] tau13i contribution to stress from du/dz
 * @param[in,out] tau23i contribution to stress from dv/dz
 * @param[in,out] tau33i contribution to stress from dw/dz
 */
void
ComputeStressConsVisc_T (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                         const Array4<const Real>& cell_data,
                         Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                         Array4<Real>& tau12, Array4<Real>& tau21,
                         Array4<Real>& tau13, Array4<Real>& tau31,
                         Array4<Real>& tau23, Array4<Real>& tau32,
                         const Array4<const Real>& er_arr,
                         const Array4<const Real>& z_nd,
                         const Array4<const Real>& detJ,
                         const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                         const Array4<const Real>& mf_mx,
                         const Array4<const Real>& mf_ux,
                         const Array4<const Real>& mf_vx,
                         const Array4<const Real>& mf_my,
                         const Array4<const Real>& mf_uy,
                         const Array4<const Real>& mf_vy,
                         Array4<Real>& tau13i,
                         Array4<Real>& tau23i,
                         Array4<Real>& tau33i)
{
    // NOTE: mu_eff includes factor of 2

    // Handle constant alpha case, in which the provided mu_eff is actually
    // "alpha" and the viscosity needs to be scaled by rho. This can be further
    // optimized with if statements below instead of creating a new FAB,
    // but this is implementation is cleaner.
    FArrayBox temp;
    Box gbx = bxcc; // Note: bxcc have been grown in x/y only.
    gbx.grow(IntVect(0,0,1));
    temp.resize(gbx,1, The_Async_Arena());
    Array4<Real> rhoAlpha = temp.array();

    if (cell_data)
    // constant alpha (stored in mu_eff)
    {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = cell_data(i, j, k, Rho_comp) * mu_eff;
        });
    }
    else
    // constant mu_eff
    {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = mu_eff;
        });
    }

    //***********************************************************************************
    // NOTE: The first  block computes Tau = K*(S-D).
    //       The second block computes the projection JT*Tau
    //       The implicit parts are not projected: tau13i and tau23i carry the whole
    //       second vertical derivative of u and v, -(K_v + K_h M) du_i/dz, where the
    //       slope factor M (Compute_TerrainVertDiffFac) is what the projection of
    //       tau11, tau12 and tau22 adds; tau33i is -K_v dw/dz, without the expansion rate.
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;            // Grown by 1 in x and y directions
    bxcc2.grow(IntVect(-1,-1,0)); // CC box without halo cells

    // First block: compute Tau = K*(S-D)
    //***********************************************************************************
    Real OneThird   = (one/three);
    ParallelFor(bxcc, tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_tot = rhoAlpha(i,j,k);
        if (tau33i) tau33i(i,j,k) = -mu_tot * tau33(i,j,k);

        tau11(i,j,k) = -mu_tot * (tau11(i,j,k) - OneThird*er_arr(i,j,k));
        tau22(i,j,k) = -mu_tot * (tau22(i,j,k) - OneThird*er_arr(i,j,k));
        tau33(i,j,k) = -mu_tot * (tau33(i,j,k) - OneThird*er_arr(i,j,k));
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_tot = fourth * ( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                               + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        tau12(i,j,k) *= -mu_tot;
        tau21(i,j,k) *= -mu_tot;
    });
    ParallelFor(tbxxz, tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_tot = fourth * ( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
        tau13(i,j,k) *= -mu_tot;
        tau31(i,j,k) *= -mu_tot;
        if (tau13i) {
            Real met_fac = Compute_TerrainVertDiffFac<0>(i,j,k,dxInv,z_nd,mf_ux,mf_uy);
            tau13i(i,j,k) *= -mu_tot * (one + met_fac);
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_tot = fourth * ( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                               + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot;
        if (tau23i) {
            Real met_fac = Compute_TerrainVertDiffFac<1>(i,j,k,dxInv,z_nd,mf_vx,mf_vy);
            tau23i(i,j,k) *= -mu_tot * (one + met_fac);
        }
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = fourth * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                                 + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = fourth * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                                 + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );

        tau33(i,j,k) -= met_h_xi*mfx*tau31bar + met_h_eta*mfy*tau32bar;
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom (tau31 & tau32 normal operations)
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);

        ParallelFor(planexz, planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = myhalf * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = Real(1.5)*tau11lo - myhalf*tau11hi;

            Real tau12lo  = myhalf * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = myhalf * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = Real(1.5)*tau12lo - myhalf*tau12hi;

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

            tau31(i,j,k) *= met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = myhalf * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = Real(1.5)*tau21lo - myhalf*tau21hi;

            Real tau22lo  = myhalf * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = myhalf * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = Real(1.5)*tau22lo - myhalf*tau22hi;

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

            tau32(i,j,k) *= met_h_zeta/mfx;
        });
    }
    // Extrapolate tau13 & tau23 to top (tau31 & tau32 normal operations)
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);

        ParallelFor(planexz, planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = myhalf * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = Real(1.5)*tau11hi - myhalf*tau11lo;

            Real tau12lo  = myhalf * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = myhalf * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = Real(1.5)*tau12hi - myhalf*tau12lo;

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

            tau31(i,j,k) *= met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = myhalf * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = Real(1.5)*tau21hi - myhalf*tau21lo;

            Real tau22lo  = myhalf * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = myhalf * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = Real(1.5)*tau22hi - myhalf*tau22lo;

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

            tau32(i,j,k) *= met_h_zeta/mfx;
        });
    }

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau13 and tau23 in valid averaging region (tau31 & tau32 normal operations)
    //-----------------------------------------------------------------------------------
    ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_ux(i,j,0);
        Real mfy = mf_uy(i,j,0);

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = fourth * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                                 + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = fourth * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                                 + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );

        tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

        tau31(i,j,k) *= met_h_zeta/mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_vx(i,j,0);
        Real mfy = mf_vy(i,j,0);

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = fourth * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                                 + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = fourth * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                                 + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );

        tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

        tau32(i,j,k) *= met_h_zeta/mfx;
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Finally project tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_zeta = detJ(i,j,k);

        tau11(i,j,k) *= met_h_zeta/mfy;
        tau22(i,j,k) *= met_h_zeta/mfx;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = myhalf * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = myhalf * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        tau12(i,j,k) *= met_h_zeta/mfx;
        tau21(i,j,k) *= met_h_zeta/mfy;
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
 * @param[in]  detJ Jacobian determinant
 * @param[in]  dxInv inverse cell size array
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[in,out] tau13i contribution to stress from du/dz
 * @param[in,out] tau23i contribution to stress from dv/dz
 * @param[in,out] tau33i contribution to stress from dw/dz
 */
void
ComputeStressVarVisc_T (Box bxcc, Box tbxxy, Box tbxxz, Box tbxyz, Real mu_eff,
                        const Array4<const Real>& mu_turb,
                        const Array4<const Real>& cell_data,
                        Array4<Real>& tau11, Array4<Real>& tau22, Array4<Real>& tau33,
                        Array4<Real>& tau12, Array4<Real>& tau21,
                        Array4<Real>& tau13, Array4<Real>& tau31,
                        Array4<Real>& tau23, Array4<Real>& tau32,
                        const Array4<const Real>& er_arr,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& dxInv,
                        const Array4<const Real>& mf_mx,
                        const Array4<const Real>& mf_ux,
                        const Array4<const Real>& mf_vx,
                        const Array4<const Real>& mf_my,
                        const Array4<const Real>& mf_uy,
                        const Array4<const Real>& mf_vy,
                        Array4<Real>& tau13i,
                        Array4<Real>& tau23i,
                        Array4<Real>& tau33i)
{
    // NOTE: mu_eff includes factor of 2

    // Handle constant alpha case, in which the provided mu_eff is actually
    // "alpha" and the viscosity needs to be scaled by rho. This can be further
    // optimized with if statements below instead of creating a new FAB,
    // but this is implementation is cleaner.
    FArrayBox temp;
    Box gbx = bxcc; // Note: bxcc have been grown in x/y only.
    gbx.grow(IntVect(0,0,1));
    temp.resize(gbx,1, The_Async_Arena());
    Array4<Real> rhoAlpha = temp.array();

    if (cell_data)
    // constant alpha (stored in mu_eff)
    {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = cell_data(i, j, k, Rho_comp) * mu_eff;
        });
    }
    else
    // constant mu_eff
    {
        ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            rhoAlpha(i,j,k) = mu_eff;
        });
    }

    //***********************************************************************************
    // NOTE: The first  block computes Tau = K*(S-D).
    //       The second block computes the projection JT*Tau
    //       The implicit parts are not projected: tau13i and tau23i carry the whole
    //       second vertical derivative of u and v, -(K_v + K_h M) du_i/dz, where the
    //       slope factor M (Compute_TerrainVertDiffFac) is what the projection of
    //       tau11, tau12 and tau22 adds; tau33i is -K_v dw/dz, without the expansion rate.
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;            // Grown by 1 in x and y directions
    bxcc2.grow(IntVect(-1,-1,0)); // CC box without halo cells

    // First block: compute Tau = K*(S-D)
    //***********************************************************************************
    Real OneThird   = (one/three);
    ParallelFor(bxcc, tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_tot_h = rhoAlpha(i,j,k) + two*mu_turb(i, j, k, EddyDiff::Mom_h);
        Real mu_tot_v = rhoAlpha(i,j,k) + two*mu_turb(i, j, k, EddyDiff::Mom_v);
        if (tau33i) tau33i(i,j,k) = -mu_tot_v * tau33(i,j,k);

        tau11(i,j,k) = -mu_tot_h * (tau11(i,j,k) - OneThird*er_arr(i,j,k));
        tau22(i,j,k) = -mu_tot_h * (tau22(i,j,k) - OneThird*er_arr(i,j,k));
        tau33(i,j,k) = -mu_tot_v * (tau33(i,j,k) - OneThird*er_arr(i,j,k));
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = fourth * ( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                               + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                                     + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;
        tau12(i,j,k) *= -mu_tot;
        tau21(i,j,k) *= -mu_tot;
    });
    ParallelFor(tbxxz, tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = fourth * ( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;
        tau13(i,j,k) *= -mu_tot;
        tau31(i,j,k) *= -mu_tot;
        if (tau13i) {
            Real mu_bar_h = fourth * ( mu_turb(i-1, j, k  , EddyDiff::Mom_h) + mu_turb(i, j, k  , EddyDiff::Mom_h)
                                     + mu_turb(i-1, j, k-1, EddyDiff::Mom_h) + mu_turb(i, j, k-1, EddyDiff::Mom_h) );
            Real mu_tot_h = rhoAlpha_bar + two*mu_bar_h;
            Real met_fac  = Compute_TerrainVertDiffFac<0>(i,j,k,dxInv,z_nd,mf_ux,mf_uy);
            tau13i(i,j,k) *= -(mu_tot + mu_tot_h * met_fac);
        }
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        Real mu_bar = fourth * ( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                               + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                     + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;
        tau23(i,j,k) *= -mu_tot;
        tau32(i,j,k) *= -mu_tot;
        if (tau23i) {
            Real mu_bar_h = fourth * ( mu_turb(i, j-1, k  , EddyDiff::Mom_h) + mu_turb(i, j, k  , EddyDiff::Mom_h)
                                     + mu_turb(i, j-1, k-1, EddyDiff::Mom_h) + mu_turb(i, j, k-1, EddyDiff::Mom_h) );
            Real mu_tot_h = rhoAlpha_bar + two*mu_bar_h;
            Real met_fac  = Compute_TerrainVertDiffFac<1>(i,j,k,dxInv,z_nd,mf_vx,mf_vy);
            tau23i(i,j,k) *= -(mu_tot + mu_tot_h * met_fac);
        }
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau33 first (no linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc2,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_xi,met_h_eta;
        met_h_xi   = Compute_h_xi_AtCellCenter  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtCellCenter (i,j,k,dxInv,z_nd);

        Real tau31bar = fourth * ( tau31(i  , j  , k  ) + tau31(i+1, j  , k  )
                                 + tau31(i  , j  , k+1) + tau31(i+1, j  , k+1) );
        Real tau32bar = fourth * ( tau32(i  , j  , k  ) + tau32(i  , j+1, k  )
                                 + tau32(i  , j  , k+1) + tau32(i  , j+1, k+1) );

        tau33(i,j,k) -= met_h_xi*mfx*tau31bar + met_h_eta*mfy*tau32bar;
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom (tau31 & tau32 normal operations)
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);

        ParallelFor(planexz, planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  ) );
            Real tau11hi  = myhalf * ( tau11(i  , j  , k+1) + tau11(i-1, j  , k+1) );
            Real tau11bar = Real(1.5)*tau11lo - myhalf*tau11hi;

            Real tau12lo  = myhalf * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  ) );
            Real tau12hi  = myhalf * ( tau12(i  , j  , k+1) + tau12(i  , j+1, k+1) );
            Real tau12bar = Real(1.5)*tau12lo - myhalf*tau12hi;

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

            tau31(i,j,k) *= met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  ) );
            Real tau21hi  = myhalf * ( tau21(i  , j  , k+1) + tau21(i+1, j  , k+1) );
            Real tau21bar = Real(1.5)*tau21lo - myhalf*tau21hi;

            Real tau22lo  = myhalf * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  ) );
            Real tau22hi  = myhalf * ( tau22(i  , j  , k+1) + tau22(i  , j-1, k+1) );
            Real tau22bar = Real(1.5)*tau22lo - myhalf*tau22hi;

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

            tau32(i,j,k) *= met_h_zeta/mfx;
        });
    }
    // Extrapolate tau13 & tau23 to top (tau31 & tau32 normal operations)
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);

        ParallelFor(planexz, planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( tau11(i  , j  , k-2) + tau11(i-1, j  , k-2) );
            Real tau11hi  = myhalf * ( tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
            Real tau11bar = Real(1.5)*tau11hi - myhalf*tau11lo;

            Real tau12lo  = myhalf * ( tau12(i  , j  , k-2) + tau12(i  , j+1, k-2) );
            Real tau12hi  = myhalf * ( tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );
            Real tau12bar = Real(1.5)*tau12hi - myhalf*tau12lo;

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

            tau31(i,j,k) *= met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( tau21(i  , j  , k-2) + tau21(i+1, j  , k-2) );
            Real tau21hi  = myhalf * ( tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
            Real tau21bar = Real(1.5)*tau21hi - myhalf*tau21lo;

            Real tau22lo  = myhalf * ( tau22(i  , j  , k-2) + tau22(i  , j-1, k-2) );
            Real tau22hi  = myhalf * ( tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );
            Real tau22bar = Real(1.5)*tau22hi - myhalf*tau22lo;

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

            tau32(i,j,k) *= met_h_zeta/mfx;
        });
    }

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Fill tau13 and tau23 in valid averaging region (tau31 & tau32 normal operations)
    //-----------------------------------------------------------------------------------
    ParallelFor(tbxxz,tbxyz,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_ux(i,j,0);
        Real mfy = mf_uy(i,j,0);

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

        Real tau11bar = fourth * ( tau11(i  , j  , k  ) + tau11(i-1, j  , k  )
                                 + tau11(i  , j  , k-1) + tau11(i-1, j  , k-1) );
        Real tau12bar = fourth * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                                 + tau12(i  , j  , k-1) + tau12(i  , j+1, k-1) );

        tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;

        tau31(i,j,k) *= met_h_zeta/mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_vx(i,j,0);
        Real mfy = mf_vy(i,j,0);

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = fourth * ( tau21(i  , j  , k  ) + tau21(i+1, j  , k  )
                                 + tau21(i  , j  , k-1) + tau21(i+1, j  , k-1) );
        Real tau22bar = fourth * ( tau22(i  , j  , k  ) + tau22(i  , j-1, k  )
                                 + tau22(i  , j  , k-1) + tau22(i  , j-1, k-1) );

        tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;

        tau32(i,j,k) *= met_h_zeta/mfx;
    });

    // Second block: compute JT*Tau
    //***********************************************************************************
    // Finally project tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_zeta = detJ(i,j,k);

        tau11(i,j,k) *= met_h_zeta/mfy;
        tau22(i,j,k) *= met_h_zeta/mfx;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = myhalf * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = myhalf * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        tau12(i,j,k) *= met_h_zeta/mfx;
        tau21(i,j,k) *= met_h_zeta/mfy;
    });
}
