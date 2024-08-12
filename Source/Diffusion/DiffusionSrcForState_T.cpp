#include <Diffusion.H>
#include <EddyViscosity.H>
#include <TerrainMetrics.H>
#include <PBLModels.H>

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator without terrain.
 *
 * @param[in]  bx cell center box to loop over
 * @param[in]  domain box of the whole domain
 * @param[in]  start_comp starting component index
 * @param[in]  num_comp number of components
 * @param[in]  u velocity in x-dir
 * @param[in]  v velocity in y-dir
 * @param[in]  cell_data conserved cell center vars
 * @param[in]  cell_prim primitive cell center vars
 * @param[out] cell_rhs RHS for cell center vars
 * @param[in]  xflux flux in x-dir
 * @param[in]  yflux flux in y-dir
 * @param[in]  zflux flux in z-dir
 * @param[in]  z_nd physical z height
 * @param[in]  detJ Jacobian determinant
 * @param[in]  cellSizeInv inverse cell size array
 * @param[in]  SmnSmn_a strain rate magnitude
 * @param[in]  mf_m map factor at cell center
 * @param[in]  mf_u map factor at x-face
 * @param[in]  mf_v map factor at y-face
 * @param[inout]  hfx_z heat flux in z-dir
 * @param[inout]  qfx1_z heat flux in z-dir
 * @param[out]    qfx2_z heat flux in z-dir
 * @param[in]  diss dissipation of TKE
 * @param[in]  mu_turb turbulent viscosity
 * @param[in]  diffChoice container of diffusion parameters
 * @param[in]  turbChoice container of turbulence parameters
 * @param[in]  tm_arr theta mean array
 * @param[in]  grav_gpu gravity vector
 * @param[in]  bc_ptr container with boundary conditions
 * @param[in]  use_most whether we have turned on MOST BCs
 */
void
DiffusionSrcForState_T (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const bool& exp_most,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& ax,
                        const Array4<const Real>& ay,
                        const Array4<const Real>& az,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& SmnSmn_a,
                        const Array4<const Real>& mf_m,
                        const Array4<const Real>& mf_u,
                        const Array4<const Real>& mf_v,
                              Array4<      Real>& hfx_z,
                              Array4<      Real>& /*qfx1_z*/,
                              Array4<      Real>& /*qfx2_z*/,
                              Array4<      Real>& diss,
                        const Array4<const Real>& mu_turb,
                        const DiffChoice &diffChoice,
                        const TurbChoice &turbChoice,
                        const Array4<const Real>& tm_arr,
                        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu,
                        const BCRec* bc_ptr,
                        const bool use_most)
{
    BL_PROFILE_VAR("DiffusionSrcForState_T()",DiffusionSrcForState_T);

    amrex::ignore_unused(use_most);

    const Real dx_inv = cellSizeInv[0];
    const Real dy_inv = cellSizeInv[1];
    const Real dz_inv = cellSizeInv[2];

    const auto& dom_hi = ubound(domain);

    bool l_use_QKE       = turbChoice.use_QKE && turbChoice.advect_QKE;
    bool l_use_deardorff = (turbChoice.les_type == LESType::Deardorff);
    Real l_inv_theta0    = 1.0 / turbChoice.theta_ref;
    Real l_abs_g         = std::abs(grav_gpu[2]);

    bool l_consA  = (diffChoice.molec_diff_type == MolecDiffType::ConstantAlpha);
    bool l_turb   = ( (turbChoice.les_type == LESType::Smagorinsky) ||
                      (turbChoice.les_type == LESType::Deardorff  ) ||
                      (turbChoice.pbl_type == PBLType::MYNN25     ) ||
                      (turbChoice.pbl_type == PBLType::YSU        ) );

    const Box xbx = surroundingNodes(bx,0);
    const Box ybx = surroundingNodes(bx,1);
    const Box zbx = surroundingNodes(bx,2);
    Box zbx3 = zbx;

    const int end_comp   = start_comp + num_comp - 1;
    const int qty_offset = RhoTheta_comp;

    // Theta, KE, QKE, Scalar
    Vector<Real> alpha_eff(NVAR_max, 0.0);
    if (l_consA) {
        for (int i = 0; i < NVAR_max; ++i) {
           switch (i) {
               case PrimTheta_comp:
                    alpha_eff[PrimTheta_comp] = diffChoice.alpha_T;
                    break;
               case PrimScalar_comp:
                    alpha_eff[PrimScalar_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ1_comp:
                    alpha_eff[PrimQ1_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ2_comp:
                    alpha_eff[PrimQ2_comp] = diffChoice.alpha_C;
                    break;
               case PrimQ3_comp:
                   alpha_eff[PrimQ3_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ4_comp:
                   alpha_eff[PrimQ4_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ5_comp:
                   alpha_eff[PrimQ5_comp] = diffChoice.alpha_C;
                   break;
               case PrimQ6_comp:
                   alpha_eff[PrimQ6_comp] = diffChoice.alpha_C;
                   break;
               default:
                    alpha_eff[i] = 0.0;
                    break;
          }
       }
    } else {
        for (int i = 0; i < NVAR_max; ++i) {
           switch (i) {
               case PrimTheta_comp:
                    alpha_eff[PrimTheta_comp] = diffChoice.rhoAlpha_T;
                    break;
               case PrimScalar_comp:
                    alpha_eff[PrimScalar_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ1_comp:
                    alpha_eff[PrimQ1_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ2_comp:
                    alpha_eff[PrimQ2_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ3_comp:
                    alpha_eff[PrimQ3_comp] = diffChoice.rhoAlpha_C;
                    break;
               case PrimQ4_comp:
                   alpha_eff[PrimQ4_comp] = diffChoice.rhoAlpha_C;
                   break;
               case PrimQ5_comp:
                   alpha_eff[PrimQ5_comp] = diffChoice.rhoAlpha_C;
                   break;
               case PrimQ6_comp:
                   alpha_eff[PrimQ6_comp] = diffChoice.rhoAlpha_C;
                   break;
               default:
                    alpha_eff[i] = 0.0;
                    break;
          }
       }
    }

    Vector<int> eddy_diff_idx{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::QKE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h, EddyDiff::Q_h ,
                              EddyDiff::Q_h    , EddyDiff::Q_h, EddyDiff::Q_h };
    Vector<int> eddy_diff_idy{EddyDiff::Theta_h, EddyDiff::KE_h, EddyDiff::QKE_h, EddyDiff::Scalar_h,
                              EddyDiff::Q_h    , EddyDiff::Q_h, EddyDiff::Q_h ,
                              EddyDiff::Q_h    , EddyDiff::Q_h, EddyDiff::Q_h };
    Vector<int> eddy_diff_idz{EddyDiff::Theta_v, EddyDiff::KE_v, EddyDiff::QKE_v, EddyDiff::Scalar_v,
                              EddyDiff::Q_v    , EddyDiff::Q_v, EddyDiff::Q_v ,
                              EddyDiff::Q_v    , EddyDiff::Q_v, EddyDiff::Q_v };

    // Device vectors
    Gpu::AsyncVector<Real> alpha_eff_d;
    Gpu::AsyncVector<int>  eddy_diff_idx_d,eddy_diff_idy_d,eddy_diff_idz_d;
    alpha_eff_d.resize(alpha_eff.size());
    eddy_diff_idx_d.resize(eddy_diff_idx.size());
    eddy_diff_idy_d.resize(eddy_diff_idy.size());
    eddy_diff_idz_d.resize(eddy_diff_idz.size());

    Gpu::copy(Gpu::hostToDevice, alpha_eff.begin()    , alpha_eff.end()    , alpha_eff_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idx.begin(), eddy_diff_idx.end(), eddy_diff_idx_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idy.begin(), eddy_diff_idy.end(), eddy_diff_idy_d.begin());
    Gpu::copy(Gpu::hostToDevice, eddy_diff_idz.begin(), eddy_diff_idz.end(), eddy_diff_idz_d.begin());

    // Capture pointers for device code
    Real* d_alpha_eff     = alpha_eff_d.data();
    int*  d_eddy_diff_idx = eddy_diff_idx_d.data();
    int*  d_eddy_diff_idy = eddy_diff_idy_d.data();
    int*  d_eddy_diff_idz = eddy_diff_idz_d.data();

    // Constant alpha & Turb model
    if (l_consA && l_turb) {
        ParallelFor(xbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                              + mu_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ax(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i-1, j, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i-1, j, k-1, prim_index) );
            Real GradCx =        dx_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k,qty_index) = rhoAlpha * mf_u(i,j,0) * ( GradCx - (met_h_xi/met_h_zeta)*GradCz );
        });
        ParallelFor(ybx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                              + mu_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ay(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i, j-1, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i, j-1, k-1, prim_index) );
            Real GradCy =        dy_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k,qty_index) = rhoAlpha * mf_v(i,j,0) * ( GradCy - (met_h_eta/met_h_zeta)*GradCz );
        });
        ParallelFor(zbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                              + mu_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            Real met_h_zeta = az(i,j,k);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == 0);
            bool ext_dir_on_zhi = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir_prim) )
                                    && k == dom_hi.z+1);
            bool most_on_zlo    = ( use_most && exp_most &&
                                  (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::foextrap) && k == 0);

            if (ext_dir_on_zlo) {
                GradCz = dz_inv * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                        + 3. * cell_prim(i, j, k  , prim_index)
                                   - (1./3.) * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                GradCz = dz_inv * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                        - 3. * cell_prim(i, j, k-1, prim_index)
                                   + (1./3.) * cell_prim(i, j, k-2, prim_index) );
            } else {
                GradCz = dz_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (most_on_zlo && (qty_index == RhoTheta_comp)) {
                // set the exact value from MOST, don't need finite diff
                zflux(i,j,k,qty_index) = -rhoFace * hfx_z(i,j,0);
            } else {
                zflux(i,j,k,qty_index) = rhoAlpha * GradCz / met_h_zeta;
            }
        });
    // Constant rho*alpha & Turb model
    } else if (l_turb) {
        ParallelFor(xbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                              + mu_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ax(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i-1, j, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i-1, j, k-1, prim_index) );
            Real GradCx =        dx_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k,qty_index) = rhoAlpha * mf_u(i,j,0) * ( GradCx - (met_h_xi/met_h_zeta)*GradCz );
        });
        ParallelFor(ybx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                              + mu_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ay(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i, j-1, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i, j-1, k-1, prim_index) );
            Real GradCy =        dy_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k,qty_index) = rhoAlpha * mf_v(i,j,0) * ( GradCy - (met_h_eta/met_h_zeta)*GradCz );
        });
        ParallelFor(zbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];

            rhoAlpha += 0.5 * ( mu_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                              + mu_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            Real met_h_zeta = az(i,j,k);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == 0);
            bool ext_dir_on_zhi = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool most_on_zlo    = ( use_most && exp_most &&
                                  (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::foextrap) && k == 0);

            if (ext_dir_on_zlo) {
                GradCz = dz_inv * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                        + 3. * cell_prim(i, j, k  , prim_index)
                                   - (1./3.) * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                GradCz = dz_inv * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                        - 3. * cell_prim(i, j, k-1, prim_index)
                                   + (1./3.) * cell_prim(i, j, k-2, prim_index) );
            } else {
                GradCz = dz_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (most_on_zlo && (qty_index == RhoTheta_comp)) {
                // set the exact value from MOST, don't need finite diff
                Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
                zflux(i,j,k,qty_index) = -rhoFace * hfx_z(i,j,0);
            } else {
                zflux(i,j,k,qty_index) = rhoAlpha * GradCz / met_h_zeta;
            }
        });
    // Constant alpha & no LES/PBL model
    } else if(l_consA) {
        ParallelFor(xbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ax(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i-1, j, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i-1, j, k-1, prim_index) );
            Real GradCx =        dx_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k,qty_index) = rhoAlpha * mf_u(i,j,0) * ( GradCx - (met_h_xi/met_h_zeta)*GradCz );
        });
        ParallelFor(ybx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
            Real met_h_zeta = ay(i,j,k);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i, j-1, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i, j-1, k-1, prim_index) );
            Real GradCy =        dy_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k,qty_index) = rhoAlpha * mf_v(i,j,0) * ( GradCy - (met_h_eta/met_h_zeta)*GradCz );
        });
        ParallelFor(zbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            Real met_h_zeta = az(i,j,k);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == 0);
            bool ext_dir_on_zhi = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool most_on_zlo    = ( use_most && exp_most &&
                                  (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::foextrap) && k == 0);

            if (ext_dir_on_zlo) {
                GradCz = dz_inv * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                        + 3. * cell_prim(i, j, k  , prim_index)
                                   - (1./3.) * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                GradCz = dz_inv * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                        - 3. * cell_prim(i, j, k-1, prim_index)
                                   + (1./3.) * cell_prim(i, j, k-2, prim_index) );
            } else {
                GradCz = dz_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (most_on_zlo && (qty_index == RhoTheta_comp)) {
                // set the exact value from MOST, don't need finite diff
                zflux(i,j,k,qty_index) = -rhoFace * hfx_z(i,j,0);
            } else {
                zflux(i,j,k,qty_index) = rhoAlpha * GradCz / met_h_zeta;
            }
        });
    // Constant rho*alpha & no LES/PBL model
    } else {
        ParallelFor(xbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_xi,met_h_zeta;
            met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtIface(i,j,k,cellSizeInv,z_nd);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i-1, j, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i-1, j, k-1, prim_index) );
            Real GradCx =        dx_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k,qty_index) = rhoAlpha * mf_u(i,j,0) * ( GradCx - (met_h_xi/met_h_zeta)*GradCz );
        });
        ParallelFor(ybx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_eta,met_h_zeta;
            met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
            met_h_zeta = Compute_h_zeta_AtJface(i,j,k,cellSizeInv,z_nd);

            Real GradCz = 0.25 * dz_inv * ( cell_prim(i, j, k+1, prim_index) + cell_prim(i, j-1, k+1, prim_index)
                                          - cell_prim(i, j, k-1, prim_index) - cell_prim(i, j-1, k-1, prim_index) );
            Real GradCy =        dy_inv * ( cell_prim(i, j, k  , prim_index) - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k,qty_index) = rhoAlpha * mf_v(i,j,0) * ( GradCy - (met_h_eta/met_h_zeta)*GradCz );
        });
        ParallelFor(zbx, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            const int  qty_index = start_comp + n;
            const int prim_index = qty_index - qty_offset;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == 0);
            bool ext_dir_on_zhi = ( ((bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[BCVars::cons_bc+qty_index].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool most_on_zlo    = ( use_most && exp_most &&
                                  (bc_ptr[BCVars::cons_bc+qty_index].lo(2) == ERFBCType::foextrap) && k == 0);

            if (ext_dir_on_zlo) {
                GradCz = dz_inv * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                        + 3. * cell_prim(i, j, k  , prim_index)
                                   - (1./3.) * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                GradCz = dz_inv * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                        - 3. * cell_prim(i, j, k-1, prim_index)
                                   + (1./3.) * cell_prim(i, j, k-2, prim_index) );
            } else {
                GradCz = dz_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (most_on_zlo && (qty_index == RhoTheta_comp)) {
                // set the exact value from MOST, don't need finite diff
                Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
                zflux(i,j,k,qty_index) = -rhoFace * hfx_z(i,j,0);
            } else {
                zflux(i,j,k,qty_index) = rhoAlpha * GradCz / met_h_zeta;
            }
        });
    }

    // Linear combinations for z-flux with terrain
    //-----------------------------------------------------------------------------------
    // Extrapolate top and bottom cells
    {
      Box planexy = zbx; planexy.setBig(2, planexy.smallEnd(2) );
      int k_lo = zbx.smallEnd(2); int k_hi = zbx.bigEnd(2);
      zbx3.growLo(2,-1); zbx3.growHi(2,-1);
      ParallelFor(planexy, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int , int n) noexcept
      {
          const int  qty_index = start_comp + n;
          Real met_h_xi,met_h_eta;

          { // Bottom face
            met_h_xi  = Compute_h_xi_AtKface (i,j,k_lo,cellSizeInv,z_nd);
            met_h_eta = Compute_h_eta_AtKface(i,j,k_lo,cellSizeInv,z_nd);

            Real xfluxlo  = 0.5 * ( xflux(i  , j  , k_lo  , qty_index) + xflux(i+1, j  , k_lo  , qty_index) );
            Real xfluxhi  = 0.5 * ( xflux(i  , j  , k_lo+1, qty_index) + xflux(i+1, j  , k_lo+1, qty_index) );
            Real xfluxbar = 1.5*xfluxlo - 0.5*xfluxhi;

            Real yfluxlo  = 0.5 * ( yflux(i  , j  , k_lo  , qty_index) + yflux(i  , j+1, k_lo  , qty_index) );
            Real yfluxhi  = 0.5 * ( yflux(i  , j  , k_lo+1, qty_index) + yflux(i  , j+1, k_lo+1, qty_index) );
            Real yfluxbar = 1.5*yfluxlo - 0.5*yfluxhi;

            zflux(i,j,k_lo,qty_index) -= met_h_xi*xfluxbar + met_h_eta*yfluxbar;
          }

          { // Top face
            met_h_xi  = Compute_h_xi_AtKface (i,j,k_hi,cellSizeInv,z_nd);
            met_h_eta = Compute_h_eta_AtKface(i,j,k_hi,cellSizeInv,z_nd);

            Real xfluxlo  = 0.5 * ( xflux(i  , j  , k_hi-2, qty_index) + xflux(i+1, j  , k_hi-2, qty_index) );
            Real xfluxhi  = 0.5 * ( xflux(i  , j  , k_hi-1, qty_index) + xflux(i+1, j  , k_hi-1, qty_index) );
            Real xfluxbar = 1.5*xfluxhi - 0.5*xfluxlo;

            Real yfluxlo  = 0.5 * ( yflux(i  , j  , k_hi-2, qty_index) + yflux(i  , j+1, k_hi-2, qty_index) );
            Real yfluxhi  = 0.5 * ( yflux(i  , j  , k_hi-1, qty_index) + yflux(i  , j+1, k_hi-1, qty_index) );
            Real yfluxbar = 1.5*yfluxhi - 0.5*yfluxlo;

            zflux(i,j,k_hi,qty_index) -= met_h_xi*xfluxbar + met_h_eta*yfluxbar;
          }
      });
    }
    // Average interior cells
    ParallelFor(zbx3, num_comp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int  qty_index = start_comp + n;

      Real met_h_xi,met_h_eta;
      met_h_xi  = Compute_h_xi_AtKface (i,j,k,cellSizeInv,z_nd);
      met_h_eta = Compute_h_eta_AtKface(i,j,k,cellSizeInv,z_nd);

      Real xfluxbar = 0.25 * ( xflux(i  , j  , k  , qty_index) + xflux(i+1, j  , k  , qty_index)
                             + xflux(i  , j  , k-1, qty_index) + xflux(i+1, j  , k-1, qty_index) );
      Real yfluxbar = 0.25 * ( yflux(i  , j  , k  , qty_index) + yflux(i  , j+1, k  , qty_index)
                             + yflux(i  , j  , k-1, qty_index) + yflux(i  , j+1, k-1, qty_index) );

      zflux(i,j,k,qty_index) -= met_h_xi*xfluxbar + met_h_eta*yfluxbar;
    });
    // Multiply h_zeta by x/y-fluxes
    ParallelFor(xbx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int  qty_index = start_comp + n;
      Real met_h_zeta      = Compute_h_zeta_AtIface(i,j,k,cellSizeInv,z_nd);
      xflux(i,j,k,qty_index) *= met_h_zeta;
    });
    ParallelFor(ybx, num_comp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      const int  qty_index = start_comp + n;
      Real met_h_zeta      = Compute_h_zeta_AtJface(i,j,k,cellSizeInv,z_nd);
      yflux(i,j,k,qty_index) *= met_h_zeta;
    });


    // Use fluxes to compute RHS
    //-----------------------------------------------------------------------------------
    for (int n(0); n < num_comp; ++n)
    {
        int qty_index = start_comp + n;
        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real stateContrib = (xflux(i+1,j  ,k  ,qty_index) - xflux(i, j, k, qty_index)) * dx_inv * mf_m(i,j,0)  // Diffusive flux in x-dir
                               +(yflux(i  ,j+1,k  ,qty_index) - yflux(i, j, k, qty_index)) * dy_inv * mf_m(i,j,0)  // Diffusive flux in y-dir
                               +(zflux(i  ,j  ,k+1,qty_index) - zflux(i, j, k, qty_index)) * dz_inv;  // Diffusive flux in z-dir

            stateContrib /= detJ(i,j,k);

            cell_rhs(i,j,k,qty_index) += stateContrib;
        });
    }

    // Using Deardorff (see Sullivan et al 1994)
    //
    // Note: At this point, the thermal diffusivity ("Khv" field in ERF), the
    //       subgrid heat flux ("hfx_z" here), and the subgrid dissipation
    //       ("diss" here) have been updated by ComputeTurbulentViscosityLES --
    //       at the beginning of each timestep.
    //       The strain rate magnitude is updated at the beginning of the first
    //       RK stage only, therefore the shear production term also does not
    //       change between RK stages.
    //       The surface heat flux hfx_z(i,j,-1) is updated in MOSTStress at
    //       each RK stage if using the ERF_EXPLICIT_MOST_STRESS path, but that
    //       does not change the buoyancy production term here.
    if (l_use_deardorff && start_comp <= RhoKE_comp && end_comp >=RhoKE_comp) {
        int qty_index = RhoKE_comp;
        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Add Buoyancy Source
            // where the SGS buoyancy flux tau_{theta,i} = -KH * dtheta/dx_i,
            // such that for dtheta/dz < 0, there is a positive (upward) heat
            // flux; the TKE buoyancy production is then
            //   B = g/theta_0 * tau_{theta,w}
            // for a dry atmosphere.
            // TODO: To account for moisture, the Brunt-Vaisala frequency,
            //   N^2 = g[1/theta * dtheta/dz + ...]
            // **should** be a function of the water vapor and total water
            // mixing ratios, depending on whether conditions are saturated or
            // not (see the WRF model description, Skamarock et al 2019).
            cell_rhs(i,j,k,qty_index) += l_abs_g * l_inv_theta0 * hfx_z(i,j,k);

            // TKE shear production
            //   P = -tau_ij * S_ij = 2 * mu_turb * S_ij * S_ij
            // Note: This assumes that the horizontal and vertical diffusivities
            // of momentum are equal
            cell_rhs(i,j,k,qty_index) += 2.0*mu_turb(i,j,k,EddyDiff::Mom_v) * SmnSmn_a(i,j,k);

            // TKE dissipation
            cell_rhs(i,j,k,qty_index) -= diss(i,j,k);
        });
    }

    // Using PBL
    if (l_use_QKE && start_comp <= RhoQKE_comp && end_comp >=RhoQKE_comp) {
        int qty_index = RhoQKE_comp;
        auto pbl_mynn_B1_l = turbChoice.pbl_mynn_B1;
        bool c_ext_dir_on_zlo = ( (bc_ptr[BCVars::cons_bc].lo(2) == ERFBCType::ext_dir) );
        bool c_ext_dir_on_zhi = ( (bc_ptr[BCVars::cons_bc].lo(5) == ERFBCType::ext_dir) );
        bool u_ext_dir_on_zlo = ( (bc_ptr[BCVars::xvel_bc].lo(2) == ERFBCType::ext_dir) );
        bool u_ext_dir_on_zhi = ( (bc_ptr[BCVars::xvel_bc].lo(5) == ERFBCType::ext_dir) );
        bool v_ext_dir_on_zlo = ( (bc_ptr[BCVars::yvel_bc].lo(2) == ERFBCType::ext_dir) );
        bool v_ext_dir_on_zhi = ( (bc_ptr[BCVars::yvel_bc].lo(5) == ERFBCType::ext_dir) );
        ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const Real met_h_zeta = detJ(i,j,k);
            cell_rhs(i, j, k, qty_index) += ComputeQKESourceTerms(i,j,k,u,v,cell_data,cell_prim,
                                                                  mu_turb,cellSizeInv,domain,pbl_mynn_B1_l,tm_arr(i,j,0),
                                                                  c_ext_dir_on_zlo, c_ext_dir_on_zhi,
                                                                  u_ext_dir_on_zlo, u_ext_dir_on_zhi,
                                                                  v_ext_dir_on_zlo, v_ext_dir_on_zhi,
                                                                  met_h_zeta);
        });
    }
}
