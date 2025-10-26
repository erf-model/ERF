#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_PBLModels.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator with stretched dz
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
DiffusionSrcForState_S (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Gpu::DeviceVector<Real>& stretched_dz_d,
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
    BL_PROFILE_VAR("DiffusionSrcForState_S()",DiffusionSrcForState_S);

    const Real explicit_fac = 1.0 - implicit_fac;

#include "ERF_SetupDiff.H"
    Real l_abs_g      = std::abs(grav_gpu[2]);

    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    auto dz_ptr = stretched_dz_d.data();

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

            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * GradCx;
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

            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * GradCy;
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
                Real dz0  = dz_ptr[k];
                Real dz1  = dz_ptr[k+1];
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
                Real dz0  = dz_ptr[k-1];
                Real dz1  = dz_ptr[k-2];
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
                Real dzk_inv;
                if (k==klo) {
                    dzk_inv =  1.0 / dz_ptr[k];
                } else if (k==(khi+1)) {
                    dzk_inv =  1.0 / dz_ptr[k-1];
                } else {
                    dzk_inv =  2.0 / (dz_ptr[k] + dz_ptr[k-1]);
                }
                GradCz = dzk_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];
            rhoAlpha += 0.5 * ( mu_turb(i  , j, k, d_eddy_diff_idx[prim_index])
                              + mu_turb(i-1, j, k, d_eddy_diff_idx[prim_index]) );

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * GradCx;
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

            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * GradCy;
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

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real dz0  = dz_ptr[k];
                Real dz1  = dz_ptr[k+1];
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
                Real dz0  = dz_ptr[k-1];
                Real dz1  = dz_ptr[k-2];
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
                Real dzk_inv;
                if (k==klo) {
                    dzk_inv =  1.0 / dz_ptr[k];
                } else if (k==(khi+1)) {
                    dzk_inv =  1.0 / dz_ptr[k-1];
                } else {
                    dzk_inv =  2.0 / (dz_ptr[k] + dz_ptr[k-1]);
                }
                GradCz = dzk_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * GradCx;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoFace  = 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) );
            Real rhoAlpha = rhoFace * d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * GradCy;
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

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real dz0  = dz_ptr[k];
                Real dz1  = dz_ptr[k+1];
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
                Real dz0  = dz_ptr[k-1];
                Real dz1  = dz_ptr[k-2];
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
                Real dzk_inv;
                if (k==klo) {
                    dzk_inv =  1.0 / dz_ptr[k];
                } else if (k==(khi+1)) {
                    dzk_inv =  1.0 / dz_ptr[k-1];
                } else {
                    dzk_inv =  2.0 / (dz_ptr[k] + dz_ptr[k-1]);
                }
                GradCz = dzk_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            Real GradCx = dx_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i-1, j, k  , prim_index) );

            xflux(i,j,k) = -rhoAlpha * mf_ux(i,j,0) * GradCx;
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int prim_index = qty_index - 1;

            Real rhoAlpha = d_alpha_eff[prim_index];

            int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                           BCVars::RhoScalar_bc_comp : qty_index;
            if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);

            Real GradCy = dy_inv * ( cell_prim(i, j, k  , prim_index)        - cell_prim(i, j-1, k  , prim_index) );

            yflux(i,j,k) = -rhoAlpha * mf_vy(i,j,0) * GradCy;
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

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real dz0  = dz_ptr[k];
                Real dz1  = dz_ptr[k+1];
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
                Real dz0  = dz_ptr[k-1];
                Real dz1  = dz_ptr[k-2];
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
                Real dzk_inv;
                if (k==klo) {
                    dzk_inv =  1.0 / dz_ptr[k];
                } else if (k==(khi+1)) {
                    dzk_inv =  1.0 / dz_ptr[k-1];
                } else {
                    dzk_inv =  2.0 / (dz_ptr[k] + dz_ptr[k-1]);
                }
                GradCz = dzk_inv * ( cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index) );
            }

            bool SurfLayer_on_zlo = ( use_SurfLayer && k == dom_lo.z);

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

    // Adjust with map factors
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        xflux(i,j,k) /= mf_uy(i,j,0);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        yflux(i,j,k) /= mf_vx(i,j,0);
    });

    // This allows us to do semi-implicit discretization of the vertical diffusive terms
    if (qty_index == RhoTheta_comp) {
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            zflux(i,j,k) *= explicit_fac;
        });
    }

    // Use fluxes to compute RHS
    //-----------------------------------------------------------------------------------
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
        Real dzk_inv = 1.0 / dz_ptr[k];
        Real stateContrib = (xflux(i+1,j  ,k  ) - xflux(i, j, k)) * dx_inv * mfsq  // Diffusive flux in x-dir
                           +(yflux(i  ,j+1,k  ) - yflux(i, j, k)) * dy_inv * mfsq  // Diffusive flux in y-dir
                           +(zflux(i  ,j  ,k+1) - zflux(i, j, k)) * dzk_inv;       // Diffusive flux in z-dir

        cell_rhs(i,j,k,qty_index) -= stateContrib;
    });

    } // n

#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
