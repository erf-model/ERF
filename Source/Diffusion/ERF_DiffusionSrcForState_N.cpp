#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_PBLModels.H"

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
 * @param[in]  use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in]  implicit_fac -- factor of implicitness for vertical differences only
 */
void
DiffusionSrcForState_N (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& SmnSmn_a,
                        const Array4<const Real>& mf_mx,
                        const Array4<const Real>& mf_ux,
                        const Array4<const Real>& mf_vx,
                        const Array4<const Real>& mf_my,
                        const Array4<const Real>& mf_uy,
                        const Array4<const Real>& mf_vy,
                              Array4<      Real>& hfx_z,
                              Array4<      Real>& qfx1_z,
                              Array4<      Real>& qfx2_z,
                              Array4<      Real>& diss,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const int level,
                        const Array4<const Real>& tm_arr,
                        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("DiffusionSrcForState_N()",DiffusionSrcForState_N);

    const Real explicit_fac = 1.0 - implicit_fac;

#include "ERF_SetupDiff.H"
    Real l_abs_g      = std::abs(grav_gpu[2]);

    const Real dz_inv = cellSizeInv[2];

    for (int n(0); n<num_comp; ++n) {
        const int qty_index = start_comp + n;

    // Constant alpha & Turb model
    if (l_consA && l_turb) {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;
            const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_scal_index];
            rhoAlpha += 0.5 * ( mu_turb(i  , j, k, d_eddy_diff_idx[prim_scal_index])
                              + mu_turb(i-1, j, k, d_eddy_diff_idx[prim_scal_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_xlo = ( (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_upwind && u(dom_lo.x,j,k) >= 0.) );
            ext_dir_on_xlo &= (i == dom_lo.x);

            bool ext_dir_on_xhi = ( (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_upwind && u(dom_hi.x+1,j,k) <= 0.) );
            ext_dir_on_xlo &= (i == dom_hi.x+1);

            if (ext_dir_on_xlo) {
                xflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i-1, j, k, prim_index)
                                                 + 3. * cell_prim(i  , j, k, prim_index)
                                            - (1./3.) * cell_prim(i+1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else if (ext_dir_on_xhi) {
                xflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i  , j, k, prim_index)
                                                 - 3. * cell_prim(i-1, j, k, prim_index)
                                            + (1./3.) * cell_prim(i-2, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * (  cell_prim(i  , j, k, prim_index)
                                            - cell_prim(i-1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;
            const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_scal_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j  , k, d_eddy_diff_idy[prim_scal_index])
                              + mu_turb(i, j-1, k, d_eddy_diff_idy[prim_scal_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_ylo = ( (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_upwind && v(i,dom_lo.y,k) >= 0.) );
            ext_dir_on_ylo &= (j == dom_lo.y);
            bool ext_dir_on_yhi = ( (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_upwind && v(i,dom_hi.y+1,k) <= 0.) );
            ext_dir_on_yhi &= (j == dom_hi.y+1);

            if (ext_dir_on_ylo) {
                yflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j-1, k, prim_index)
                                                 + 3. * cell_prim(i, j  , k, prim_index)
                                            - (1./3.) * cell_prim(i, j+1, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else if (ext_dir_on_yhi) {
                yflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j  , k, prim_index)
                                                 - 3. * cell_prim(i, j-1, k, prim_index)
                                            + (1./3.) * cell_prim(i, j-2, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else {
                yflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;
            const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_scal_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j, k  , d_eddy_diff_idz[prim_scal_index])
                              + mu_turb(i, j, k-1, d_eddy_diff_idz[prim_scal_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k==dom_lo.z);

            if (ext_dir_on_zlo) {
                zflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                                 + 3. * cell_prim(i, j, k  , prim_index)
                                            - (1./3.) * cell_prim(i, j, k+1, prim_index) ) * dz_inv;
            } else if (ext_dir_on_zhi) {
                zflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                                 - 3. * cell_prim(i, j, k-1, prim_index)
                                            + (1./3.) * cell_prim(i, j, k-2, prim_index) ) * dz_inv;
            } else if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
                }
            } else  if (qty_index == RhoQ1_comp) {
                if (!SurfLayer_on_zlo) {
                    qfx1_z(i,j,k) = zflux(i,j,k);
                }
            } else  if (qty_index == RhoQ2_comp) {
                qfx2_z(i,j,k) = zflux(i,j,k);
            }
        });
    // Constant rho*alpha & Turb model
    } else if (l_turb) {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                              + mu_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_xlo = ( (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_upwind && u(dom_lo.x,j,k) >= 0.) );
            ext_dir_on_xlo &= (i == dom_lo.x);

            bool ext_dir_on_xhi = ( (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_upwind && u(dom_hi.x+1,j,k) <= 0.) );
            ext_dir_on_xhi &= (i == dom_hi.x+1);

            if (ext_dir_on_xlo) {
                xflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i-1, j, k, prim_index)
                                                 + 3. * cell_prim(i  , j, k, prim_index)
                                            - (1./3.) * cell_prim(i+1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else if (ext_dir_on_xhi) {
                xflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i  , j, k, prim_index)
                                                 - 3. * cell_prim(i-1, j, k, prim_index)
                                            + (1./3.) * cell_prim(i-2, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * ( cell_prim(i  , j, k, prim_index)
                                           - cell_prim(i-1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                              + mu_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_ylo = ( (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_upwind && v(i,dom_lo.y,k) >= 0.) );
            ext_dir_on_ylo &= (j == dom_lo.y);

            bool ext_dir_on_yhi = ( (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_upwind && v(i,dom_hi.y+1,k) <= 0.) );
            ext_dir_on_yhi &= (j == dom_hi.y+1);

            if (ext_dir_on_ylo) {
                yflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j-1, k, prim_index)
                                                 + 3. * cell_prim(i, j  , k, prim_index)
                                            - (1./3.) * cell_prim(i, j+1, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else if (ext_dir_on_yhi) {
                yflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j  , k, prim_index)
                                                 - 3. * cell_prim(i, j-1, k, prim_index)
                                            + (1./3.) * cell_prim(i, j-2, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else {
              yflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                              + mu_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                zflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                                 + 3. * cell_prim(i, j, k  , prim_index)
                                            - (1./3.) * cell_prim(i, j, k+1, prim_index) ) * dz_inv;
            } else if (ext_dir_on_zhi) {
                zflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                                 - 3. * cell_prim(i, j, k-1, prim_index)
                                            + (1./3.) * cell_prim(i, j, k-2, prim_index) ) * dz_inv;
            } else if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
                }
            } else  if (qty_index == RhoQ1_comp) {
                if (!SurfLayer_on_zlo) {
                    qfx1_z(i,j,k) = zflux(i,j,k);
                }
            } else  if (qty_index == RhoQ2_comp) {
                qfx2_z(i,j,k) = zflux(i,j,k);
            }
        });
    // Constant alpha & no LES/PBL model
    } else if(l_consA) {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_xlo = ( (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_upwind && u(dom_lo.x,j,k) >= 0.) );
            ext_dir_on_xlo &= (i == dom_lo.x);

            bool ext_dir_on_xhi = ( (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_upwind && u(dom_hi.x+1,j,k) <= 0.) );
            ext_dir_on_xhi &= (i == dom_hi.x+1);

            if (ext_dir_on_xlo) {
                xflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i-1, j, k, prim_index)
                                                 + 3. * cell_prim(i  , j, k, prim_index)
                                            - (1./3.) * cell_prim(i+1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else if (ext_dir_on_xhi) {
                xflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i  , j, k, prim_index)
                                                 - 3. * cell_prim(i-1, j, k, prim_index)
                                            + (1./3.) * cell_prim(i-2, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * ( cell_prim(i  , j, k, prim_index)
                                           - cell_prim(i-1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_ylo = ( (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_upwind && v(i,dom_lo.y,k) >= 0.) );
            ext_dir_on_ylo &= (j == dom_lo.y);

            bool ext_dir_on_yhi = ( (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_upwind && v(i,dom_hi.y+1,k) <= 0.) );
            ext_dir_on_yhi &= (j == dom_hi.y+1);

            if (ext_dir_on_ylo) {
                yflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j-1, k, prim_index)
                                                 + 3. * cell_prim(i, j  , k, prim_index)
                                            - (1./3.) * cell_prim(i, j+1, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else if (ext_dir_on_yhi) {
                yflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j  , k, prim_index)
                                                 - 3. * cell_prim(i, j-1, k, prim_index)
                                            + (1./3.) * cell_prim(i, j-2, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else {
              yflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                zflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                                 + 3. * cell_prim(i, j, k  , prim_index)
                                            - (1./3.) * cell_prim(i, j, k+1, prim_index) ) * dz_inv;
            } else if (ext_dir_on_zhi) {
                zflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                                 - 3. * cell_prim(i, j, k-1, prim_index)
                                            + (1./3.) * cell_prim(i, j, k-2, prim_index) ) * dz_inv;
            } else if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
                }
            } else  if (qty_index == RhoQ1_comp) {
                if (!SurfLayer_on_zlo) {
                    qfx1_z(i,j,k) = zflux(i,j,k);
                }
            } else  if (qty_index == RhoQ2_comp) {
                qfx2_z(i,j,k) = zflux(i,j,k);
            }
        });
    // Constant rho*alpha & no LES/PBL model
    } else {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_xlo = ( (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(0) == ERFBCType::ext_dir_upwind && u(dom_lo.x,j,k) >= 0.) );
            ext_dir_on_xlo &= (i == dom_lo.x);

            bool ext_dir_on_xhi = ( (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(0) == ERFBCType::ext_dir_upwind && u(dom_hi.x+1,j,k) <= 0.) );
            ext_dir_on_xhi &= (i == dom_hi.x+1);

            if (ext_dir_on_xlo) {
                xflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i-1, j, k, prim_index)
                                                 + 3. * cell_prim(i  , j, k, prim_index)
                                            - (1./3.) * cell_prim(i+1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else if (ext_dir_on_xhi) {
                xflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i  , j, k, prim_index)
                                                 - 3. * cell_prim(i-1, j, k, prim_index)
                                            + (1./3.) * cell_prim(i-2, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * ( cell_prim(i  , j, k, prim_index)
                                           - cell_prim(i-1, j, k, prim_index) ) * dx_inv * mf_ux(i,j,0)/mf_uy(i,j,0);
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool ext_dir_on_ylo = ( (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].lo(1) == ERFBCType::ext_dir_upwind && v(i,dom_lo.y,k) >= 0.) );
            ext_dir_on_ylo &= (j == dom_lo.y);

            bool ext_dir_on_yhi = ( (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir)      ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_prim) ||
                                    (bc_ptr[bc_comp].hi(1) == ERFBCType::ext_dir_upwind && v(i,dom_hi.y+1,k) <= 0.) );
            ext_dir_on_yhi &= (j == dom_hi.y+1);

            if (ext_dir_on_ylo) {
                yflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j-1, k, prim_index)
                                                 + 3. * cell_prim(i, j  , k, prim_index)
                                            - (1./3.) * cell_prim(i, j+1, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else if (ext_dir_on_yhi) {
                yflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j  , k, prim_index)
                                                 - 3. * cell_prim(i, j-1, k, prim_index)
                                            + (1./3.) * cell_prim(i, j-2, k, prim_index) ) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            } else {
              yflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j-1, k, prim_index)) * dy_inv * mf_vy(i,j,0)/mf_vx(i,j,0);
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                zflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                                 + 3. * cell_prim(i, j, k  , prim_index)
                                            - (1./3.) * cell_prim(i, j, k+1, prim_index) ) * dz_inv;
            } else if (ext_dir_on_zhi) {
                zflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                                 - 3. * cell_prim(i, j, k-1, prim_index)
                                            + (1./3.) * cell_prim(i, j, k-2, prim_index) ) * dz_inv;
            } else if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
                }
            } else  if (qty_index == RhoQ1_comp) {
                if (!SurfLayer_on_zlo) {
                    qfx1_z(i,j,k) = zflux(i,j,k);
                }
            } else  if (qty_index == RhoQ2_comp) {
                qfx2_z(i,j,k) = zflux(i,j,k);
            }
        });
    }

    // This allows us to do semi-implicit discretization of the vertical diffusive terms
    if (qty_index == RhoTheta_comp) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            zflux(i,j,k) *= explicit_fac;
        });
    }

    // Use fluxes to compute RHS
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
        cell_rhs(i,j,k,qty_index) -= (xflux(i+1,j  ,k  ) - xflux(i, j, k)) * dx_inv * mfsq  // Diffusive flux in x-dir
                                    +(yflux(i  ,j+1,k  ) - yflux(i, j, k)) * dy_inv * mfsq  // Diffusive flux in y-dir
                                    +(zflux(i  ,j  ,k+1) - zflux(i, j, k)) * dz_inv;        // Diffusive flux in z-dir
    });
    } // n

#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
