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
 * @param[inout] cell_data conserved cell center vars
 * @param[in   ] cellSizeInv inverse cell size array
 * @param[inout] hfx_z heat flux in z-dir
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] diffChoice container of diffusion parameters
 * @param[in   ] turbChoice container of turbulence parameters
 * @param[in   ] tm_arr theta mean array
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 */
void
ImplicitDiffForState_N (const Box& bx, const Box& domain,
                        const int level,
                        const Real dt,
                        /*int start_comp, int num_comp,*/
                        const GpuArray<Real, AMREX_SPACEDIM*2>& bc_neumann_vals,
                        const Array4<      Real>& cell_data,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                        const Array4<      Real>& hfx_z,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForState_N()",ImplicitDiffForState_N);

    // this uses domain, level, start_comp, num_comp
#include "ERF_SetupVertDiff.H"

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
    //   A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]
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

//  bool ext_dir_on_zlo = (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir ||
//                         bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim)
//  bool ext_dir_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir ||
//                         bc_ptr[bc_comp].hi(2) == ERFBCType::ext_dir_prim)
    bool neumann_on_zlo = (bc_ptr[bc_comp].lo(2) == ERFBCType::neumann);
    bool neumann_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::neumann);

    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {

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

            RHS_a(i,j,k)  = cell_data(i,j,k,n); // Note this is rho*theta, whereas solution will be theta

            coeffA_a(i,j,k) = -implicit_fac * rhoAlpha_lo * dt * dz_inv * dz_inv;
            coeffC_a(i,j,k) = -implicit_fac * rhoAlpha_hi * dt * dz_inv * dz_inv;

            if (k == dom_lo.z) {
                if (use_SurfLayer) {
                    RHS_a(i,j,klo) += implicit_fac * dt * dz_inv * hfx_z(i,j,0);
                } else if (neumann_on_zlo) {
                    RHS_a(i,j,klo) += coeffA_a(i,j,klo) * bc_neumann_vals[2] / dz_inv;
                }

                coeffA_a(i,j,klo) = 0.;
            }
            if (k == dom_hi.z) {
                if (neumann_on_zhi) {
                    RHS_a(i,j,khi) -= coeffC_a(i,j,khi) * bc_neumann_vals[5] / dz_inv;
                }

                coeffC_a(i,j,khi) = 0.;
            }

            coeffB_a(i,j,k) = cell_data(i,j,k,Rho_comp) - coeffA_a(i,j,k) - coeffC_a(i,j,k);
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

        //
        // Transfer back to original array
        //
        for (int k(klo); k<=khi; ++k) {
            cell_data(i,j,k,n) = soln_a(i,j,k) * cell_data(i,j,k,Rho_comp);
        }

      } // i
    } // j
}
