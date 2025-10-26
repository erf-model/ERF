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
 * @param[in]  implicit_fac -- factor of implicitness for vertical differences only
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
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("DiffusionSrcForState_T()",DiffusionSrcForState_T);

    const Real explicit_fac = 1.0 - implicit_fac;

#include "ERF_SetupDiff.H"
    Real l_abs_g      = std::abs(grav_gpu[2]);

    const Real dz_inv = cellSizeInv[2];

    // We need to grow these boxes in the vertical direction when tiling so that we can access xflux and yflux
    //    to modify zflux
    Box xbx_g1(xbx); Box ybx_g1(ybx);
    if (xbx_g1.smallEnd(2) != dom_lo.z) xbx_g1.growLo(2,1);
    if (ybx_g1.smallEnd(2) != dom_lo.z) ybx_g1.growLo(2,1);
    if (xbx_g1.bigEnd(2)   != dom_hi.z) xbx_g1.growHi(2,1);
    if (ybx_g1.bigEnd(2)   != dom_hi.z) ybx_g1.growHi(2,1);

    for (int n(0); n<num_comp; ++n) {
        const int qty_index = start_comp + n;

    // Constant alpha & Turb model
    if (l_consA && l_turb) {
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
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
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
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
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
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
                    hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
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
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
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
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
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

    //-----------------------------------------------------------------------------------
    //
    // Modify fluxes by terrain and use fluxes to compute RHS
    //
    // Note that we combine all of these operations in order to keep this section
    //      of the loop tiling-safe.
    //-----------------------------------------------------------------------------------
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real met_h_xi_lo  = Compute_h_xi_AtKface (i,j,k  ,cellSizeInv,z_nd);
        Real met_h_xi_hi  = Compute_h_xi_AtKface (i,j,k+1,cellSizeInv,z_nd);

        Real met_h_eta_lo = Compute_h_eta_AtKface(i,j,k  ,cellSizeInv,z_nd);
        Real met_h_eta_hi = Compute_h_eta_AtKface(i,j,k+1,cellSizeInv,z_nd);

        Real xfluxbar_lo, yfluxbar_lo;
        if (k == dom_lo.z) {
            Real xfluxlo  = 0.5 * ( xflux(i,j,k  ) + xflux(i+1,j,k  ) );
            Real xfluxhi  = 0.5 * ( xflux(i,j,k+1) + xflux(i+1,j,k+1) );
            xfluxbar_lo = 1.5*xfluxlo - 0.5*xfluxhi;

            Real yfluxlo  = 0.5 * ( yflux(i,j,k  ) + yflux(i,j+1,k  ) );
            Real yfluxhi  = 0.5 * ( yflux(i,j,k+1) + yflux(i,j+1,k+1) );
            yfluxbar_lo = 1.5*yfluxlo - 0.5*yfluxhi;
        } else {
            xfluxbar_lo = 0.25 * ( xflux(i,j,k  ) + xflux(i+1,j  ,k  )
                                 + xflux(i,j,k-1) + xflux(i+1,j  ,k-1) );
            yfluxbar_lo = 0.25 * ( yflux(i,j,k  ) + yflux(i  ,j+1,k  )
                                 + yflux(i,j,k-1) + yflux(i  ,j+1,k-1) );
        }

        Real xfluxbar_hi, yfluxbar_hi;
        if (k == dom_hi.z) {
            Real xfluxlo  = 0.5 * ( xflux(i,j,k-1) + xflux(i+1,j,k-1) );
            Real xfluxhi  = 0.5 * ( xflux(i,j,k  ) + xflux(i+1,j,k  ) );
            xfluxbar_hi = 1.5*xfluxhi - 0.5*xfluxlo;

            Real yfluxlo  = 0.5 * ( yflux(i,j,k-1) + yflux(i,j+1,k-1) );
            Real yfluxhi  = 0.5 * ( yflux(i,j,k  ) + yflux(i,j+1,k  ) );
            yfluxbar_hi = 1.5*yfluxhi - 0.5*yfluxlo;
        } else {
            xfluxbar_hi = 0.25 * ( xflux(i,j,k+1) + xflux(i+1,j  ,k+1)
                                 + xflux(i,j,k  ) + xflux(i+1,j  ,k  ) );
            yfluxbar_hi = 0.25 * ( yflux(i,j,k+1) + yflux(i  ,j+1,k+1)
                                 + yflux(i,j,k  ) + yflux(i  ,j+1,k  ) );
        }

        // Allow semi-implicit discretization of the vertical diffusive terms
        Real zflux_lo = explicit_fac * zflux(i,j,k  )
                      - met_h_xi_lo  * mf_mx(i,j,0) * xfluxbar_lo
                      - met_h_eta_lo * mf_my(i,j,0) * yfluxbar_lo;
        Real zflux_hi = explicit_fac * zflux(i,j,k+1)
                      - met_h_xi_hi  * mf_mx(i,j,0) * xfluxbar_hi
                      - met_h_eta_hi * mf_my(i,j,0) * yfluxbar_hi;

        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
        Real stateContrib = ( xflux(i+1,j  ,k  ) * ax(i+1,j,k) / mf_uy(i+1,j,0)
                             -xflux(i  ,j  ,k  ) * ax(i  ,j,k) / mf_uy(i  ,j,0) ) * dx_inv * mfsq  // Diffusive flux in x-dir
                           +( yflux(i  ,j+1,k  ) * ay(i,j+1,k) / mf_vx(i,j+1,0)
                             -yflux(i  ,j  ,k  ) * ay(i,j  ,k) / mf_vx(i,j  ,0) ) * dy_inv * mfsq  // Diffusive flux in y-dir
                           +( zflux_hi - zflux_lo)                                * dz_inv;        // Diffusive flux in z-dir

        stateContrib /= detJ(i,j,k);

        cell_rhs(i,j,k,qty_index) -= stateContrib;
    });
    } // n

#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
