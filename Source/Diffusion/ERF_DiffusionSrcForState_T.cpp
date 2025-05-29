#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_PBLModels.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator with terrain-fitted coordinates.
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
 * @param[in]  use_SurfLayer whether we have turned on subgrid diffusion
 */
void
DiffusionSrcForState_T (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const bool& rotate,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& z_cc,
                        const Array4<const Real>& ax,
                        const Array4<const Real>& ay,
                        const Array4<const Real>& /*az*/,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& SmnSmn_a,
                        const Array4<const Real>& mf_mx,
                        const Array4<const Real>& mf_ux,
                        const Array4<const Real>& mf_vx,
                        const Array4<const Real>& mf_my,
                        const Array4<const Real>& mf_uy,
                        const Array4<const Real>& mf_vy,
                              Array4<      Real>& hfx_x,
                              Array4<      Real>& hfx_y,
                              Array4<      Real>& hfx_z,
                              Array4<      Real>& qfx1_x,
                              Array4<      Real>& qfx1_y,
                              Array4<      Real>& qfx1_z,
                              Array4<      Real>& qfx2_z,
                              Array4<      Real>& diss,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const int level,
                        const Array4<const Real>& tm_arr,
                        const GpuArray<Real,AMREX_SPACEDIM> grav_gpu,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer)
{
    BL_PROFILE_VAR("DiffusionSrcForState_T()",DiffusionSrcForState_T);

#include "ERF_DiffSetup.H"

    const Real dz_inv = cellSizeInv[2];

    Box zbx3 = zbx;

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

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = 1.0 / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i-1, j, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i-1, j, k-1, prim_index)*idz_lo );
            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                xflux(i,j,k) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                xflux(i,j,k) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * ( GradCx - met_h_xi*GradCz );
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

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = 1.0 / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i, j-1, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i, j-1, k-1, prim_index)*idz_lo );
            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                yflux(i,j,k) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                yflux(i,j,k) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * ( GradCy - met_h_eta*GradCz );
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

            Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);

            Real GradCz;
            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir_prim) )
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * ( c1 * cell_prim(i, j, k-1, prim_index)
                                + c2 * cell_prim(i, j, k  , prim_index)
                                + c3 * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * (  -( c1 * cell_prim(i, j, k  , prim_index)
                                    + c2 * cell_prim(i, j, k-1, prim_index)
                                    + c3 * cell_prim(i, j, k-2, prim_index) ) );
            } else {
                GradCz = (dz_inv/met_h_zeta) * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * GradCz;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k);
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

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = 1.0 / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i-1, j, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i-1, j, k-1, prim_index)*idz_lo );
            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                xflux(i,j,k) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                xflux(i,j,k) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * ( GradCx - met_h_xi*GradCz );
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i, j  , k, d_eddy_diff_idy[prim_index])
                              + mu_turb(i, j-1, k, d_eddy_diff_idy[prim_index]) );

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = 1.0 / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i, j-1, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i, j-1, k-1, prim_index)*idz_lo );
            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                yflux(i,j,k) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                yflux(i,j,k) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * ( GradCy - met_h_eta*GradCz );
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            rhoAlpha += 0.5 * ( mu_turb(i, j, k  , d_eddy_diff_idz[prim_index])
                              + mu_turb(i, j, k-1, d_eddy_diff_idz[prim_index]) );

            Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);

            Real GradCz;
            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * ( c1 * cell_prim(i, j, k-1, prim_index)
                                + c2 * cell_prim(i, j, k  , prim_index)
                                + c3 * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * (  -( c1 * cell_prim(i, j, k  , prim_index)
                                    + c2 * cell_prim(i, j, k-1, prim_index)
                                    + c3 * cell_prim(i, j, k-2, prim_index) ) );
            } else {
                GradCz = (dz_inv/met_h_zeta) * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * GradCz;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k);
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

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = 1.0 / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i-1, j, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i-1, j, k-1, prim_index)*idz_lo );
            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                xflux(i,j,k) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                xflux(i,j,k) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * ( GradCx - met_h_xi*GradCz );
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = 1.0 / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i, j-1, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i, j-1, k-1, prim_index)*idz_lo );
            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                yflux(i,j,k) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                yflux(i,j,k) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * ( GradCy - met_h_eta*GradCz );
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);

            Real GradCz;
            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * ( c1 * cell_prim(i, j, k-1, prim_index)
                                + c2 * cell_prim(i, j, k  , prim_index)
                                + c3 * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * (  -( c1 * cell_prim(i, j, k  , prim_index)
                                    + c2 * cell_prim(i, j, k-1, prim_index)
                                    + c3 * cell_prim(i, j, k-2, prim_index) ) );
            } else {
                GradCz = (dz_inv/met_h_zeta) * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * GradCz;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k);
                }
            } else if (qty_index == RhoQ1_comp) {
                if (!SurfLayer_on_zlo) {
                    qfx1_z(i,j,k) = zflux(i,j,k);
                }
            } else if (qty_index == RhoQ2_comp) {
                qfx2_z(i,j,k) = zflux(i,j,k);
            }
        });
    // Constant rho*alpha & no LES/PBL model
    } else {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_xi = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = 1.0 / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i-1, j, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i-1, j, k-1, prim_index)*idz_lo );
            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                xflux(i,j,k) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                xflux(i,j,k) = qfx1_x(i,j,0);
            } else {
              xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * ( GradCx - met_h_xi*GradCz );
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && rotate && k == dom_lo.z);

            Real idz_hi = 1.0 / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = 1.0 / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    0.5 * ( cell_prim(i, j, k+1, prim_index)*idz_hi + cell_prim(i, j-1, k+1, prim_index)*idz_lo
                                   - cell_prim(i, j, k-1, prim_index)*idz_hi - cell_prim(i, j-1, k-1, prim_index)*idz_lo );
            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                yflux(i,j,k) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                yflux(i,j,k) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * ( GradCy - met_h_eta*GradCz );
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            Real met_h_zeta;
            met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);

            Real GradCz;
            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
            bool ext_dir_on_zlo = ( ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir) ||
                                     (bc_ptr[bc_comp].lo(5) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * ( c1 * cell_prim(i, j, k-1, prim_index)
                                + c2 * cell_prim(i, j, k  , prim_index)
                                + c3 * cell_prim(i, j, k+1, prim_index) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = 1.0 / dz0;
                Real f    = (dz1 / dz0) + 2.0;
                Real f2   = f*f;
                Real c3   = 2.0 / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(1.0-f2)*c3;
                GradCz = idz0 * (  -( c1 * cell_prim(i, j, k  , prim_index)
                                    + c2 * cell_prim(i, j, k-1, prim_index)
                                    + c3 * cell_prim(i, j, k-2, prim_index) ) );
            } else {
                GradCz = (dz_inv/met_h_zeta) * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            if (SurfLayer_on_zlo && (qty_index == RhoTheta_comp)) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else if (SurfLayer_on_zlo && (qty_index == RhoQ1_comp)) {
                zflux(i,j,k) = qfx1_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * GradCz;
            }

            if (qty_index == RhoTheta_comp) {
                if (!SurfLayer_on_zlo) {
                    hfx_z(i,j,k) = zflux(i,j,k);
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

    // Linear combinations for z-flux with terrain
    //-----------------------------------------------------------------------------------
    // Extrapolate top and bottom cells
    {
      Box planexy = zbx; planexy.setBig(2, planexy.smallEnd(2) );
      int k_lo = zbx.smallEnd(2); int k_hi = zbx.bigEnd(2);
      zbx3.growLo(2,-1); zbx3.growHi(2,-1);
      ParallelFor(planexy, [=] AMREX_GPU_DEVICE (int i, int j, int ) noexcept
      {
          Real met_h_xi,met_h_eta;

          { // Bottom face
              met_h_xi  = Compute_h_xi_AtKface (i,j,k_lo,cellSizeInv,z_nd);
              met_h_eta = Compute_h_eta_AtKface(i,j,k_lo,cellSizeInv,z_nd);

              Real xfluxlo  = 0.5 * ( xflux(i  , j  , k_lo  ) + xflux(i+1, j  , k_lo  ) );
              Real xfluxhi  = 0.5 * ( xflux(i  , j  , k_lo+1) + xflux(i+1, j  , k_lo+1) );
              Real xfluxbar = 1.5*xfluxlo - 0.5*xfluxhi;

              Real yfluxlo  = 0.5 * ( yflux(i  , j  , k_lo  ) + yflux(i  , j+1, k_lo  ) );
              Real yfluxhi  = 0.5 * ( yflux(i  , j  , k_lo+1) + yflux(i  , j+1, k_lo+1) );
              Real yfluxbar = 1.5*yfluxlo - 0.5*yfluxhi;

              zflux(i,j,k_lo) -= met_h_xi*mf_mx(i,j,0)*xfluxbar + met_h_eta*mf_my(i,j,0)*yfluxbar;
          }

          { // Top face
              met_h_xi  = Compute_h_xi_AtKface (i,j,k_hi,cellSizeInv,z_nd);
              met_h_eta = Compute_h_eta_AtKface(i,j,k_hi,cellSizeInv,z_nd);

              Real xfluxlo  = 0.5 * ( xflux(i  , j  , k_hi-2) + xflux(i+1, j  , k_hi-2) );
              Real xfluxhi  = 0.5 * ( xflux(i  , j  , k_hi-1) + xflux(i+1, j  , k_hi-1) );
              Real xfluxbar = 1.5*xfluxhi - 0.5*xfluxlo;

              Real yfluxlo  = 0.5 * ( yflux(i  , j  , k_hi-2) + yflux(i  , j+1, k_hi-2) );
              Real yfluxhi  = 0.5 * ( yflux(i  , j  , k_hi-1) + yflux(i  , j+1, k_hi-1) );
              Real yfluxbar = 1.5*yfluxhi - 0.5*yfluxlo;

              zflux(i,j,k_hi) -= met_h_xi*mf_mx(i,j,0)*xfluxbar + met_h_eta*mf_my(i,j,0)*yfluxbar;
          }
      });
    }
    // Average interior cells
    ParallelFor(zbx3, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi,met_h_eta;
        met_h_xi  = Compute_h_xi_AtKface (i,j,k,cellSizeInv,z_nd);
        met_h_eta = Compute_h_eta_AtKface(i,j,k,cellSizeInv,z_nd);

        Real xfluxbar = 0.25 * ( xflux(i  , j  , k  ) + xflux(i+1, j  , k  )
                               + xflux(i  , j  , k-1) + xflux(i+1, j  , k-1) );
        Real yfluxbar = 0.25 * ( yflux(i  , j  , k  ) + yflux(i  , j+1, k  )
                               + yflux(i  , j  , k-1) + yflux(i  , j+1, k-1) );

        zflux(i,j,k) -= met_h_xi*mf_mx(i,j,0)*xfluxbar + met_h_eta*mf_my(i,j,0)*yfluxbar;
    });
    // Multiply h_zeta by x/y-fluxes
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        xflux(i,j,k) *= ax(i,j,k)/mf_uy(i,j,0);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        yflux(i,j,k) *= ay(i,j,k)/mf_vx(i,j,0);
    });


    // Use fluxes to compute RHS
    //-----------------------------------------------------------------------------------
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
        Real stateContrib = (xflux(i+1,j  ,k  ) - xflux(i, j, k)) * dx_inv * mfsq  // Diffusive flux in x-dir
                           +(yflux(i  ,j+1,k  ) - yflux(i, j, k)) * dy_inv * mfsq  // Diffusive flux in y-dir
                           +(zflux(i  ,j  ,k+1) - zflux(i, j, k)) * dz_inv;        // Diffusive flux in z-dir

        stateContrib /= detJ(i,j,k);

        cell_rhs(i,j,k,qty_index) -= stateContrib;
    });
    } // n

#include "ERF_DiffTKEAdjustment.H"
#include "ERF_DiffQKEAdjustment.H"
}
