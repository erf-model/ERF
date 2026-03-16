#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"

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
 * @param[inout]  hfx_z heat flux in z-dir
 * @param[inout]  qfx1_z heat flux in z-dir
 * @param[out]    qfx2_z heat flux in z-dir
 * @param[in]  mu_turb turbulent viscosity
 * @param[in]  diffChoice container of diffusion parameters
 * @param[in]  turbChoice container of turbulence parameters
 * @param[in]  grav_gpu gravity vector
 * @param[in]  bc_ptr container with boundary conditions
 * @param[in]  use_SurfLayer whether we have turned on subgrid diffusion
 */
void
DiffusionSrcForState_EB (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Array4<const EBCellFlag>& cfg_arr,
                        const Array4<const Real>& ax_arr,
                        const Array4<const Real>& ay_arr,
                        const Array4<const Real>& az_arr,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        [[maybe_unused]] Array4<Real>& hfx_z,
                        [[maybe_unused]] Array4<Real>& qfx1_z,
                        [[maybe_unused]] Array4<Real>& qfx2_z,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const int level,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer)
{
    BL_PROFILE_VAR("DiffusionSrcForState_EB()",DiffusionSrcForState_EB);

#include "ERF_SetupDiff.H"

    const Real dz_inv = cellSizeInv[2];

    for (int n(0); n<num_comp; ++n) {
        const int qty_index = start_comp + n;
        const int prim_index = qty_index - 1;
        const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;
        const int eff_index = (l_consA && l_turb) ? prim_scal_index : prim_index;
        int bc_comp = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ?
                       BCVars::RhoScalar_bc_comp : qty_index;
        if (bc_comp > BCVars::RhoScalar_bc_comp) bc_comp -= (NSCALARS-1);
        const Real alpha_mol = d_alpha_eff[eff_index];
        const int  eddy_x    = d_eddy_diff_idx[eff_index];
        const int  eddy_y    = d_eddy_diff_idy[eff_index];
        const int  eddy_z    = d_eddy_diff_idz[eff_index];

        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoFace  = l_consA ? 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i-1, j, k, Rho_comp) ) : 1.0;
            Real rhoAlpha = rhoFace * alpha_mol;
            if (l_turb) {
                rhoAlpha += 0.5 * ( mu_turb(i  , j, k, eddy_x)
                                  + mu_turb(i-1, j, k, eddy_x) );
            }

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
                                            - (1./3.) * cell_prim(i+1, j, k, prim_index) ) * dx_inv;
            } else if (ext_dir_on_xhi) {
                xflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i  , j, k, prim_index)
                                                 - 3. * cell_prim(i-1, j, k, prim_index)
                                            + (1./3.) * cell_prim(i-2, j, k, prim_index) ) * dx_inv;
            } else {
                if (cfg_arr(i,j,k).isCovered()) {
                    xflux(i,j,k) = -rhoAlpha * ( cell_prim(i-3, j, k, prim_index)
                                            - 3.*cell_prim(i-2, j, k, prim_index)
                                            + 2.*cell_prim(i-1, j, k, prim_index) ) * dx_inv;
                } else if (cfg_arr(i-1,j,k).isCovered()) {
                    xflux(i,j,k) = -rhoAlpha * ( 3.*cell_prim(i+1, j, k, prim_index)
                                            -    cell_prim(i+2, j, k, prim_index)
                                            - 2.*cell_prim(i, j, k, prim_index) ) * dx_inv;
                } else {
                    xflux(i,j,k) = -rhoAlpha * ( cell_prim(i  , j, k, prim_index)
                                               - cell_prim(i-1, j, k, prim_index) ) * dx_inv;
                }
            }
        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoFace  = l_consA ? 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j-1, k, Rho_comp) ) : 1.0;
            Real rhoAlpha = rhoFace * alpha_mol;
            if (l_turb) {
                rhoAlpha += 0.5 * ( mu_turb(i, j  , k, eddy_y)
                                  + mu_turb(i, j-1, k, eddy_y) );
            }

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
                                            - (1./3.) * cell_prim(i, j+1, k, prim_index) ) * dy_inv;
            } else if (ext_dir_on_yhi) {
                yflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j  , k, prim_index)
                                                 - 3. * cell_prim(i, j-1, k, prim_index)
                                            + (1./3.) * cell_prim(i, j-2, k, prim_index) ) * dy_inv;
            } else {
                if (cfg_arr(i,j,k).isCovered()) {
                    yflux(i,j,k) = -rhoAlpha * ( cell_prim(i, j-3, k, prim_index)
                                            - 3.*cell_prim(i, j-2, k, prim_index)
                                            + 2.*cell_prim(i, j-1, k, prim_index) ) * dy_inv;
                } else if (cfg_arr(i,j-1,k).isCovered()) {
                    yflux(i,j,k) = -rhoAlpha * ( 3.*cell_prim(i, j+1, k, prim_index)
                                            -    cell_prim(i, j+2, k, prim_index)
                                            - 2.*cell_prim(i, j, k, prim_index) ) * dy_inv;
                } else {
                    yflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index)
                                               - cell_prim(i, j-1, k, prim_index)) * dy_inv;
                }
            }
        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoFace  = l_consA ? 0.5 * ( cell_data(i, j, k, Rho_comp) + cell_data(i, j, k-1, Rho_comp) ) : 1.0;
            Real rhoAlpha = rhoFace * alpha_mol;
            if (l_turb) {
                rhoAlpha += 0.5 * ( mu_turb(i, j, k  , eddy_z)
                                  + mu_turb(i, j, k-1, eddy_z) );
            }

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
                if (cfg_arr(i,j,k).isCovered()) {
                    zflux(i,j,k) = -rhoAlpha * ( cell_prim(i, j, k-3, prim_index)
                                            - 3.*cell_prim(i, j, k-2, prim_index)
                                            + 2.*cell_prim(i, j, k-1, prim_index) ) * dz_inv;
                } else if (cfg_arr(i,j,k-1).isCovered()) {
                    zflux(i,j,k) = -rhoAlpha * ( 3.*cell_prim(i, j, k+1, prim_index)
                                            -    cell_prim(i, j, k+2, prim_index)
                                            - 2.*cell_prim(i, j, k, prim_index) ) * dz_inv;
                } else {
                    zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index)
                                            - cell_prim(i, j, k-1, prim_index)) * dz_inv;
                }
            }

            // Store z-boundary fluxes.
            // if (qty_index == RhoTheta_comp) {
            //     if (!SurfLayer_on_zlo) {
            //         hfx_z(i,j,k) = zflux(i,j,k) * explicit_fac;
            //     }
            // } else  if (qty_index == RhoQ1_comp) {
            //     if (!SurfLayer_on_zlo) {
            //         qfx1_z(i,j,k) = zflux(i,j,k);
            //     }
            // } else  if (qty_index == RhoQ2_comp) {
            //     qfx2_z(i,j,k) = zflux(i,j,k);
            // }
        });

    // Use fluxes to compute RHS
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!cfg_arr(i,j,k).isCovered()) {
            cell_rhs(i,j,k,qty_index) -= ((ax_arr(i+1,j,k) * xflux(i+1,j  ,k  ) - ax_arr(i,j,k) * xflux(i, j, k)) * dx_inv
                                        +(ay_arr(i,j+1,k) * yflux(i  ,j+1,k  ) - ay_arr(i,j,k) * yflux(i, j, k)) * dy_inv
                                        +(az_arr(i,j,k+1) * zflux(i  ,j  ,k+1) - az_arr(i,j,k) * zflux(i, j, k)) * dz_inv)
                                        / detJ(i,j,k);
        }
    });
    } // n
}
