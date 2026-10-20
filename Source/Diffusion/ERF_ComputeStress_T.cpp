#include <ERF_Diffusion.H>
#include <ERF_TerrainMetrics.H>
#include <ERF_TerrainImplicitMetric.H>

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
 * @param[in]  implicit_metric erf.implicit_terrain_metric: tau13i/tau23i also hold the compact
 *             terrain-metric term on the interior faces (see ERF_TerrainImplicitMetric.H)
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
                         Array4<Real>& tau33i,
                         const bool implicit_metric)
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
    // NOTE: The first  block computes (S-D).
    //       The second block computes JT*2mu*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;
    bxcc2.grow(IntVect(-1,-1,0));

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (one/three);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (tau33i) tau33i(i,j,k) = tau33(i,j,k);

        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute JT*2mu*(S-D)
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
        Real mu_tot   = rhoAlpha(i,j,k);

        tau33(i,j,k) -= met_h_xi*mfx*tau31bar + met_h_eta*mfy*tau32bar;
        tau33(i,j,k) *= -mu_tot;

        if (tau33i) { tau33i(i,j,k) *= -mu_tot; }
    });

    // Second block: compute JT*2mu*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planexz,planeyz,
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

            Real mu_tot = fourth*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
            tau13(i,j,k) *= -mu_tot;
            if (tau13i) { tau13i(i,j,k) *= -mu_tot; }

            tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
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

            Real mu_tot = fourth*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
            tau23(i,j,k) *= -mu_tot;
            if (tau23i) { tau23i(i,j,k) *= -mu_tot; }

            tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planexz,planeyz,
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

            Real mu_tot = fourth*( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );

            tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
            tau13(i,j,k) *= -mu_tot;
            if (tau13i) { tau13i(i,j,k) *= -mu_tot; }

            tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
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

            Real mu_tot = fourth*( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                 + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );

            tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
            tau23(i,j,k) *= -mu_tot;
            if (tau23i) { tau23i(i,j,k) *= -mu_tot; }

            tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        });
    }

    // Second block: compute JT*2mu*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
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
        Real mu_tot = fourth * ( rhoAlpha(i-1, j  , k  ) + rhoAlpha(i  , j  , k  )
                               + rhoAlpha(i-1, j  , k-1) + rhoAlpha(i  , j  , k-1) );

        tau13(i,j,k) -= met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
        tau13(i,j,k) *= -mu_tot;
        if (tau13i) {
            // erf.implicit_terrain_metric: the implicit solve also takes the compact
            // K_h*M*d/dz metric term on these interior faces, so tau13i holds it too
            const Real metric = (implicit_metric) ?
                TerrainMetricMomFactor(0, met_h_xi, met_h_eta, mfx, mfy) : zero;
            tau13i(i,j,k) *= -(mu_tot + mu_tot*metric);
        }

        tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
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
        Real mu_tot = fourth * ( rhoAlpha(i  , j-1, k  ) + rhoAlpha(i  , j  , k  )
                               + rhoAlpha(i  , j-1, k-1) + rhoAlpha(i  , j  , k-1) );

        tau23(i,j,k) -= met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
        tau23(i,j,k) *= -mu_tot;
        if (tau23i) {
            // erf.implicit_terrain_metric: the implicit solve also takes the compact
            // K_h*M*d/dz metric term on these interior faces, so tau23i holds it too
            const Real metric = (implicit_metric) ?
                TerrainMetricMomFactor(1, met_h_xi, met_h_eta, mfx, mfy) : zero;
            tau23i(i,j,k) *= -(mu_tot + mu_tot*metric);
        }

        tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_zeta = detJ(i,j,k);
        Real mu_tot = rhoAlpha(i,j,k);

        tau11(i,j,k) *= -mu_tot*met_h_zeta/mfy;
        tau22(i,j,k) *= -mu_tot*met_h_zeta/mfx;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = myhalf * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = myhalf * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        Real mu_tot = fourth*( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                             + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );

        tau12(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        tau21(i,j,k) *= -mu_tot*met_h_zeta/mfy;
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
 * @param[in]  implicit_metric erf.implicit_terrain_metric: tau13i/tau23i also hold the compact
 *             terrain-metric term on the interior faces (see ERF_TerrainImplicitMetric.H)
 *
 * NOTE: The zeta-face stresses tau13/tau23 are the terrain-normal combination
 *       tau_i3 - h_xi*tau_i1 - h_eta*tau_i2. tau_i3 is a vertical stress and takes K_v
 *       (EddyDiff::Mom_v); tau_i1 and tau_i2 are horizontal stresses and take K_h
 *       (EddyDiff::Mom_h), as the h*Fx terms of the scalar fluxes do. Applying K_v to the
 *       projected terms makes the operator non-symmetric and, when K_h*h^2 > 2*K_v,
 *       anti-diffusive. The projected stresses K_h*S are formed at the cells and xy edges
 *       and averaged to the zeta edges. On a uniform slope with uniform dz this is the
 *       transpose of the metric term in S11/S22, so the strain part of the operator
 *       dissipates energy for any K_h and K_v (an edge-averaged K_h times averaged strains
 *       does not, where K_h varies from cell to cell). tau13i/tau23i (the part the
 *       implicit solve takes) are pure K_v terms.
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
                        Array4<Real>& tau33i,
                        const bool implicit_metric)
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
    // NOTE: The first  block computes (S-D).
    //       The second block computes JT*2K*(S-D)
    //       Boxes are copied here for extrapolations in the second block operations
    //***********************************************************************************
    Box bxcc2  = bxcc;
    bxcc2.grow(IntVect(-1,-1,0));

    // First block: compute S-D
    //***********************************************************************************
    Real OneThird   = (one/three);
    ParallelFor(bxcc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        if (tau33i) tau33i(i,j,k) = tau33(i,j,k);

        tau11(i,j,k) -= OneThird*er_arr(i,j,k);
        tau22(i,j,k) -= OneThird*er_arr(i,j,k);
        tau33(i,j,k) -= OneThird*er_arr(i,j,k);
    });

    // Second block: compute JT*2K*(S-D)
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

        Real mu_tot   = rhoAlpha(i,j,k) + two*mu_turb(i, j, k, EddyDiff::Mom_v);

        tau33(i,j,k) -= met_h_xi*mfx*tau31bar + met_h_eta*mfy*tau32bar;
        tau33(i,j,k) *= -mu_tot;

        if (tau33i) { tau33i(i,j,k) *= -mu_tot; }
    });

    // The zeta-face stresses project the horizontal stresses K_h*S11, K_h*S12, K_h*S21 and
    // K_h*S22 (see the NOTE above).  Form them where they live -- at the cells and the xy
    // edges, with the coefficients the tau11/tau22/tau12 kernels below apply -- and average
    // the stresses to the zeta edges, as the scalar fluxes average K_h*grad.  On a uniform
    // slope with uniform dz that is the transpose of the metric term in S11/S22, which keeps
    // the strain part of the operator dissipative for any K_h and K_v (see the NOTE above).
    //-----------------------------------------------------------------------------------
    // kh12 is read at (i,j+1) over tbxxz and kh21 at (i+1,j) over tbxyz.  On a tile that is not
    // the last in x or y the nodal tile box tbxxy stops one node short of those reads, so form
    // them on the nodes the reads need; the strains are valid there (they are computed on the
    // halo-grown boxes).
    Box khxy = tbxxy;
    khxy.setSmall(0, amrex::min(tbxxz.smallEnd(0), tbxyz.smallEnd(0)));
    khxy.setSmall(1, amrex::min(tbxxz.smallEnd(1), tbxyz.smallEnd(1)));
    khxy.setBig  (0, amrex::max(tbxxz.bigEnd(0)  , tbxyz.bigEnd(0)+1));
    khxy.setBig  (1, amrex::max(tbxxz.bigEnd(1)+1, tbxyz.bigEnd(1)  ));
    FArrayBox kh11_fab(bxcc,1,The_Async_Arena()), kh22_fab(bxcc,1,The_Async_Arena());
    FArrayBox kh12_fab(khxy,1,The_Async_Arena()), kh21_fab(khxy,1,The_Async_Arena());
    const Array4<Real> kh11 = kh11_fab.array();
    const Array4<Real> kh22 = kh22_fab.array();
    const Array4<Real> kh12 = kh12_fab.array();
    const Array4<Real> kh21 = kh21_fab.array();
    ParallelFor(bxcc, khxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_h_tot = rhoAlpha(i,j,k) + two*mu_turb(i, j, k, EddyDiff::Mom_h);
        kh11(i,j,k) = mu_h_tot*tau11(i,j,k);
        kh22(i,j,k) = mu_h_tot*tau22(i,j,k);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mu_bar = fourth * ( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                               + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                                     + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        Real mu_h_tot = rhoAlpha_bar + two*mu_bar;
        kh12(i,j,k) = mu_h_tot*tau12(i,j,k);
        kh21(i,j,k) = mu_h_tot*tau21(i,j,k);
    });

    // Second block: compute JT*2K*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (linear combination extrapolation)
    //-----------------------------------------------------------------------------------
    // Extrapolate tau13 & tau23 to bottom
    {
        Box planexz = tbxxz; planexz.setBig(2, planexz.smallEnd(2) );
        tbxxz.growLo(2,-1);
        Box planeyz = tbxyz; planeyz.setBig(2, planeyz.smallEnd(2) );
        tbxyz.growLo(2,-1);
        ParallelFor(planexz,planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( kh11(i  , j  , k  ) + kh11(i-1, j  , k  ) );
            Real tau11hi  = myhalf * ( kh11(i  , j  , k+1) + kh11(i-1, j  , k+1) );
            Real tau11bar = Real(1.5)*tau11lo - myhalf*tau11hi;

            Real tau12lo  = myhalf * ( kh12(i  , j  , k  ) + kh12(i  , j+1, k  ) );
            Real tau12hi  = myhalf * ( kh12(i  , j  , k+1) + kh12(i  , j+1, k+1) );
            Real tau12bar = Real(1.5)*tau12lo - myhalf*tau12hi;

            Real mu_bar = fourth * ( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                                   + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                         + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + two*mu_bar;

            // K_v on S13; the projected horizontal stresses carry K_h (see the NOTE above)
            tau13(i,j,k) *= -mu_tot;
            tau13(i,j,k) += met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
            if (tau13i) { tau13i(i,j,k) *= -mu_tot; }

            tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( kh21(i  , j  , k  ) + kh21(i+1, j  , k  ) );
            Real tau21hi  = myhalf * ( kh21(i  , j  , k+1) + kh21(i+1, j  , k+1) );
            Real tau21bar = Real(1.5)*tau21lo - myhalf*tau21hi;

            Real tau22lo  = myhalf * ( kh22(i  , j  , k  ) + kh22(i  , j-1, k  ) );
            Real tau22hi  = myhalf * ( kh22(i  , j  , k+1) + kh22(i  , j-1, k+1) );
            Real tau22bar = Real(1.5)*tau22lo - myhalf*tau22hi;

            Real mu_bar = fourth * ( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                                   + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = fourth * ( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                         + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + two*mu_bar;

            // K_v on S23; the projected horizontal stresses carry K_h (see the NOTE above)
            tau23(i,j,k) *= -mu_tot;
            tau23(i,j,k) += met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
            if (tau23i) { tau23i(i,j,k) *= -mu_tot; }

            tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        });
    }
    // Extrapolate tau13 & tau23 to top
    {
        Box planexz = tbxxz; planexz.setSmall(2, planexz.bigEnd(2) );
        tbxxz.growHi(2,-1);
        Box planeyz = tbxyz; planeyz.setSmall(2, planeyz.bigEnd(2) );
        tbxyz.growHi(2,-1);
        ParallelFor(planexz,planeyz,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_ux(i,j,0);
            Real mfy = mf_uy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterJ  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterJ (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterJ(i,j,k,dxInv,z_nd);

            Real tau11lo  = myhalf * ( kh11(i  , j  , k-2) + kh11(i-1, j  , k-2) );
            Real tau11hi  = myhalf * ( kh11(i  , j  , k-1) + kh11(i-1, j  , k-1) );
            Real tau11bar = Real(1.5)*tau11hi - myhalf*tau11lo;

            Real tau12lo  = myhalf * ( kh12(i  , j  , k-2) + kh12(i  , j+1, k-2) );
            Real tau12hi  = myhalf * ( kh12(i  , j  , k-1) + kh12(i  , j+1, k-1) );
            Real tau12bar = Real(1.5)*tau12hi - myhalf*tau12lo;

            Real mu_bar = fourth * ( mu_turb(i-1, j, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                                   + mu_turb(i-1, j, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j, k  ) + rhoAlpha(i, j, k  )
                                         + rhoAlpha(i-1, j, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + two*mu_bar;

            // K_v on S13; the projected horizontal stresses carry K_h (see the NOTE above)
            tau13(i,j,k) *= -mu_tot;
            tau13(i,j,k) += met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
            if (tau13i) { tau13i(i,j,k) *= -mu_tot; }

            tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real mfx = mf_vx(i,j,0);
            Real mfy = mf_vy(i,j,0);

            Real met_h_xi,met_h_eta,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
            met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

            Real tau21lo  = myhalf * ( kh21(i  , j  , k-2) + kh21(i+1, j  , k-2) );
            Real tau21hi  = myhalf * ( kh21(i  , j  , k-1) + kh21(i+1, j  , k-1) );
            Real tau21bar = Real(1.5)*tau21hi - myhalf*tau21lo;

            Real tau22lo  = myhalf * ( kh22(i  , j  , k-2) + kh22(i  , j-1, k-2) );
            Real tau22hi  = myhalf * ( kh22(i  , j  , k-1) + kh22(i  , j-1, k-1) );
            Real tau22bar = Real(1.5)*tau22hi - myhalf*tau22lo;

            Real mu_bar = fourth * ( mu_turb(i, j-1, k  , EddyDiff::Mom_v) + mu_turb(i, j, k  , EddyDiff::Mom_v)
                                   + mu_turb(i, j-1, k-1, EddyDiff::Mom_v) + mu_turb(i, j, k-1, EddyDiff::Mom_v) );
            Real rhoAlpha_bar = fourth * ( rhoAlpha(i, j-1, k  ) + rhoAlpha(i, j, k  )
                                         + rhoAlpha(i, j-1, k-1) + rhoAlpha(i, j, k-1) );
            Real mu_tot = rhoAlpha_bar + two*mu_bar;

            // K_v on S23; the projected horizontal stresses carry K_h (see the NOTE above)
            tau23(i,j,k) *= -mu_tot;
            tau23(i,j,k) += met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
            if (tau23i) { tau23i(i,j,k) *= -mu_tot; }

            tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        });
    }

    // Second block: compute JT*2K*(S-D)
    //***********************************************************************************
    // Fill tau13, tau23 next (valid averaging region)
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

        Real tau11bar = fourth * ( kh11(i  , j  , k  ) + kh11(i-1, j  , k  )
                               + kh11(i  , j  , k-1) + kh11(i-1, j  , k-1) );
        Real tau12bar = fourth * ( kh12(i  , j  , k  ) + kh12(i  , j+1, k  )
                               + kh12(i  , j  , k-1) + kh12(i  , j+1, k-1) );

        Real mu_bar = fourth * ( mu_turb(i-1, j  , k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                               + mu_turb(i-1, j  , k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j  , k  ) + rhoAlpha(i  , j  , k  )
                                     + rhoAlpha(i-1, j  , k-1) + rhoAlpha(i  , j  , k-1) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;

        // K_v on S13; the projected horizontal stresses carry K_h (see the NOTE above)
        tau13(i,j,k) *= -mu_tot;
        tau13(i,j,k) += met_h_xi*mfx*tau11bar + met_h_eta*mfy*tau12bar;
        if (tau13i) {
            // erf.implicit_terrain_metric: the implicit solve also takes the compact
            // K_h*M*d/dz metric term on these interior faces, with K_h averaged to the edge
            // as getRhoAlphaForFaces averages it there, so tau13i holds it too
            Real metric_coef = zero;
            if (implicit_metric) {
                Real mu_h_bar = fourth * ( mu_turb(i-1, j  , k  , EddyDiff::Mom_h) + mu_turb(i  , j  , k  , EddyDiff::Mom_h)
                                         + mu_turb(i-1, j  , k-1, EddyDiff::Mom_h) + mu_turb(i  , j  , k-1, EddyDiff::Mom_h) );
                metric_coef = (rhoAlpha_bar + two*mu_h_bar)
                            * TerrainMetricMomFactor(0, met_h_xi, met_h_eta, mfx, mfy);
            }
            tau13i(i,j,k) *= -(mu_tot + metric_coef);
        }

        tau31(i,j,k) *= -mu_tot*met_h_zeta/mfy;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_vx(i,j,0);
        Real mfy = mf_vy(i,j,0);

        Real met_h_xi,met_h_eta,met_h_zeta;
        met_h_xi   = Compute_h_xi_AtEdgeCenterI  (i,j,k,dxInv,z_nd);
        met_h_eta  = Compute_h_eta_AtEdgeCenterI (i,j,k,dxInv,z_nd);
        met_h_zeta = Compute_h_zeta_AtEdgeCenterI(i,j,k,dxInv,z_nd);

        Real tau21bar = fourth * ( kh21(i  , j  , k  ) + kh21(i+1, j  , k  )
                                 + kh21(i  , j  , k-1) + kh21(i+1, j  , k-1) );
        Real tau22bar = fourth * ( kh22(i  , j  , k  ) + kh22(i  , j-1, k  )
                                 + kh22(i  , j  , k-1) + kh22(i  , j-1, k-1) );

        Real mu_bar = fourth * ( mu_turb(i  , j-1, k  , EddyDiff::Mom_v) + mu_turb(i  , j  , k  , EddyDiff::Mom_v)
                               + mu_turb(i  , j-1, k-1, EddyDiff::Mom_v) + mu_turb(i  , j  , k-1, EddyDiff::Mom_v) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i  , j-1, k  ) + rhoAlpha(i  , j  , k  )
                                     + rhoAlpha(i  , j-1, k-1) + rhoAlpha(i  , j  , k-1) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;

        // K_v on S23; the projected horizontal stresses carry K_h (see the NOTE above)
        tau23(i,j,k) *= -mu_tot;
        tau23(i,j,k) += met_h_xi*mfx*tau21bar + met_h_eta*mfy*tau22bar;
        if (tau23i) {
            // erf.implicit_terrain_metric: the implicit solve also takes the compact
            // K_h*M*d/dz metric term on these interior faces, with K_h averaged to the edge
            // as getRhoAlphaForFaces averages it there, so tau23i holds it too
            Real metric_coef = zero;
            if (implicit_metric) {
                Real mu_h_bar = fourth * ( mu_turb(i  , j-1, k  , EddyDiff::Mom_h) + mu_turb(i  , j  , k  , EddyDiff::Mom_h)
                                         + mu_turb(i  , j-1, k-1, EddyDiff::Mom_h) + mu_turb(i  , j  , k-1, EddyDiff::Mom_h) );
                metric_coef = (rhoAlpha_bar + two*mu_h_bar)
                            * TerrainMetricMomFactor(1, met_h_xi, met_h_eta, mfx, mfy);
            }
            tau23i(i,j,k) *= -(mu_tot + metric_coef);
        }

        tau32(i,j,k) *= -mu_tot*met_h_zeta/mfx;
    });

    // Fill the remaining components: tau11, tau22, tau12/21
    //-----------------------------------------------------------------------------------
    ParallelFor(bxcc,tbxxy,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = mf_mx(i,j,0);
        Real mfy = mf_my(i,j,0);

        Real met_h_zeta = detJ(i,j,k);

        Real mu_tot = rhoAlpha(i,j,k) + two*mu_turb(i, j, k, EddyDiff::Mom_h);

        tau11(i,j,k) *= -mu_tot*met_h_zeta/mfy;
        tau22(i,j,k) *= -mu_tot*met_h_zeta/mfx;
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfx = myhalf * (mf_ux(i,j,0) + mf_ux(i,j-1,0));
        Real mfy = myhalf * (mf_vy(i,j,0) + mf_vy(i-1,j,0));

        Real met_h_zeta = Compute_h_zeta_AtEdgeCenterK(i,j,k,dxInv,z_nd);

        Real mu_bar = fourth * ( mu_turb(i-1, j  , k, EddyDiff::Mom_h) + mu_turb(i, j  , k, EddyDiff::Mom_h)
                               + mu_turb(i-1, j-1, k, EddyDiff::Mom_h) + mu_turb(i, j-1, k, EddyDiff::Mom_h) );
        Real rhoAlpha_bar = fourth * ( rhoAlpha(i-1, j  , k) + rhoAlpha(i, j  , k)
                                     + rhoAlpha(i-1, j-1, k) + rhoAlpha(i, j-1, k) );
        Real mu_tot = rhoAlpha_bar + two*mu_bar;

        tau12(i,j,k) *= -mu_tot*met_h_zeta/mfx;
        tau21(i,j,k) *= -mu_tot*met_h_zeta/mfy;
    });
}
