#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_PBLModels.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator without terrain.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] start_comp starting component index
 * @param[in   ] num_comp number of components
 * @param[in   ] u velocity in x-dir
 * @param[in   ] v velocity in y-dir
 * @param[inout] cell_data conserved cell center vars
 * @param[in   ] zflux flux in z-dir
 * @param[in   ] detJ Jacobian determinant
 * @param[in   ] cellSizeInv inverse cell size array
 * @param[in   ] SmnSmn_a strain rate magnitude in general; here just empty pointer
 * @param[inout] hfx_z heat flux in z-dir
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] diffChoice container of diffusion parameters
 * @param[in   ] turbChoice container of turbulence parameters
 * @param[in   ] tm_arr theta mean array
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 */
void
ImplicitDiffForState_N (const Box& bx, const Box& domain, const Real dt,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<      Real>& cell_data,
                        const Array4<      Real>& zflux,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<const Real>& SmnSmn_a,
                        const Array4<      Real>& hfx_z,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const int level,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForState_N()",ImplicitDiffForState_N);

#include "ERF_DiffSetup.H"

    const int         n = RhoTheta_comp;
    const int qty_index = RhoTheta_comp;
    const int prim_index = qty_index - 1;
    const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;

    // Box bounds
    int ilo = bx.smallEnd(0);
    int ihi = bx.bigEnd(0);
    int jlo = bx.smallEnd(1);
    int jhi = bx.bigEnd(1);
    int klo = bx.smallEnd(2);
    int khi = bx.bigEnd(2);

 // Temporary FABs for tridiagonal solve (allocated on column)
 // A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]
    amrex::FArrayBox RHS_fab, soln_fab, coeffA_fab, coeffB_fab, inv_coeffB_fab, coeffC_fab;
           RHS_fab.resize(bx,1, amrex::The_Async_Arena());
          soln_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffA_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffB_fab.resize(bx,1, amrex::The_Async_Arena());
    inv_coeffB_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffC_fab.resize(bx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& coeffA_a     =     coeffA_fab.array(); // lower diagonal
    auto const& coeffB_a     =     coeffB_fab.array(); // diagonal
    auto const& inv_coeffB_a = inv_coeffB_fab.array();
    auto const& coeffC_a     =     coeffC_fab.array(); // upper diagonal

    int bc_comp = qty_index;

    Real rhoAlpha_lo;
    Real rhoAlpha_hi;

    Real dz_inv        = cellSizeInv[2];

    bool ext_dir_on_zlo = ((bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir) ||
                           (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim));
    bool ext_dir_on_zhi = ((bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir) ||
                           (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim));

    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {

        // Dirichlet top/bottom
//             RHS_a(i,j,klo) = 0.0;
//             RHS_a(i,j,khi) = 0.0;
//          coeffA_a(i,j,klo) = 0.0;
//          coeffA_a(i,j,khi) = 0.0;
//          coeffC_a(i,j,klo) = 0.0;
//          coeffC_a(i,j,khi) = 0.0;
//      inv_coeffB_a(i,j,klo) = 1.0;
//      inv_coeffB_a(i,j,khi) = 1.0;

        // Build the coefficients and RHS
        for (int k(klo); k <= khi; k++)
        {
            if (l_consA && l_turb) {
                rhoAlpha_lo = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                            + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                                    + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_scal_index]) );
                rhoAlpha_hi = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_scal_index]
                            + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_scal_index])
                                    + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_scal_index]) );
            }
            else if (l_turb) // with MolecDiffType::Constant or None
            {
                rhoAlpha_lo = d_alpha_eff[prim_index]
                            + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                                    + mu_turb(i,j,k-1, d_eddy_diff_idz[prim_index]) );
                rhoAlpha_hi =  d_alpha_eff[prim_index]
                            + 0.5 * ( mu_turb(i,j,k  , d_eddy_diff_idz[prim_index])
                                    + mu_turb(i,j,k+1, d_eddy_diff_idz[prim_index]) );
            }
            else if (l_consA) // without an LES/PBL model
            {
                rhoAlpha_lo = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) ) * d_alpha_eff[prim_index];
                rhoAlpha_hi = 0.5 * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k+1,Rho_comp) ) * d_alpha_eff[prim_index];
            }
            else // with MolecDiffType::Constant or None - without an LES/PBL model
            {
                rhoAlpha_lo = d_alpha_eff[prim_index];
                rhoAlpha_hi = d_alpha_eff[prim_index];
            }

#if 0
            if (ext_dir_on_zlo) {
                zflux(i,j,k) = -rhoAlpha * ( -(8./3.) * cell_prim(i, j, k-1, prim_index)
                                                 + 3. * cell_prim(i, j, k  , prim_index)
                                            - (1./3.) * cell_prim(i, j, k+1, prim_index) ) * dz_inv;
            } else if (ext_dir_on_zhi) {
                zflux(i,j,k) = -rhoAlpha * (  (8./3.) * cell_prim(i, j, k  , prim_index)
                                                 - 3. * cell_prim(i, j, k-1, prim_index)
                                        + (1./3.) * cell_prim(i, j, k-2, prim_index) ) * dz_inv;
            } else if (SurfLayer_on_zlo) {
                zflux(i,j,k) = hfx_z(i,j,0);
            } else {
                zflux(i,j,k) = -rhoAlpha * (cell_prim(i, j, k, prim_index) - cell_prim(i, j, k-1, prim_index)) * dz_inv;
            }
#endif

            coeffA_a(i,j,k) = -implicit_fac * rhoAlpha_lo * dt * dz_inv * dz_inv;
            coeffC_a(i,j,k) = -implicit_fac * rhoAlpha_hi * dt * dz_inv * dz_inv;

            // TODO: inhomogeneous BCs
            if (k == klo) {
                coeffA_a(i,j,klo) = 0.; // Zero gradient at bottom boundary
            }
            if (k == khi) {
                coeffC_a(i,j,khi) = 0.; // Zero gradient at top boundary
            }

            coeffB_a(i,j,k) = cell_data(i,j,k,Rho_comp) - coeffA_a(i,j,k) - coeffC_a(i,j,k);

            //amrex::Print() <<" A B C " << k << " " <<
            //            coeffA_a(i,j,k) << " " << coeffB_a(i,j,k) << " " << coeffC_a(i,j,k) << std::endl;

            RHS_a(i,j,k)  = cell_data(i,j,k,n); // Note this is rho*theta, whereas solution will be theta
        } // k

        // Forward sweep

        Real bet = coeffB_a(i,j,klo);

        for (int k(klo+1); k<=khi; ++k) {
            Real gam = coeffC_a(i,j,k-1) / bet;
            bet = coeffB_a(i,j,k) - coeffA_a(i,j,k)*gam;
            AMREX_ASSERT(bet != 0.0);
            coeffB_a(i,j,k) = bet;
        }

        for (int k(klo); k<=khi; ++k) {
            inv_coeffB_a(i,j,k) = 1.0 / coeffB_a(i,j,k);
        }

#if 0
        if (use_SurfLayer) {
            RHS_a(i,j,klo) += implicit_fac * dt * hfx_z(i,j,0);
        }
#endif

        //
        // Tridiagonal solve
        //
        soln_a(i,j,klo) = RHS_a(i,j,klo) * inv_coeffB_a(i,j,klo);

        for (int k(klo+1); k<=khi; ++k) {
            soln_a(i,j,k) = (RHS_a(i,j,k)-coeffA_a(i,j,k)*soln_a(i,j,k-1)) * inv_coeffB_a(i,j,k);
        }

        for (int k(khi-1); k>=klo; --k) {
            soln_a(i,j,k) -= ( coeffC_a(i,j,k) * inv_coeffB_a(i,j,k) ) * soln_a(i,j,k+1);
        }

#if 0
        zflux(i,j,klo  ) = 0.0;  // We assume zero gradient at bottom (except if prescribed gradient)
        if (!use_SurfLayer) {
            hfx_z(i,j,0) -= zflux(i,j,klo);
        }
        zflux(i,j,khi+1) = 0.0;  // We assume zero gradient at top

        for (int k(klo+1); k<=khi; ++k) {
            zflux(i,j,k) = (coeffA_a(i,j,k)/dt) * (soln_a(i,j,k,0) - soln_a(i,j,k-1,0));
        }
#endif

        //
        // Transfer back to original array
        //
        for (int k(klo); k<=khi; ++k) {
            cell_data(i,j,k,n) = soln_a(i,j,k) * cell_data(i,j,k,Rho_comp);
        }

      } // i
    } // j
}
