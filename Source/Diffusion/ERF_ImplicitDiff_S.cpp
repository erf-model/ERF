#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_GetRhoAlpha.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator without terrain.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] bc_neumann_vals values of derivatives if bc_type == Neumann
 * @param[inout] cell_data conserved cell center vars
 * @param[in   ] stretched_dz_d array over z of dz[k]
 * @param[inout] hfx_z heat flux in z-dir
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] solverChoice container of parameters
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in   ] implicit_fac if 1 then fully implicit; if 0 then fully explicit
 */
void
ImplicitDiffForState_S (const Box& bx, const Box& domain,
                        const int level,
                        const Real dt,
                        const GpuArray<Real, AMREX_SPACEDIM*2>& bc_neumann_vals,
                        const Array4<      Real>& cell_data,
                        const Gpu::DeviceVector<Real>& stretched_dz_d,
                        const Array4<      Real>& hfx_z,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForState_S()",ImplicitDiffForState_T);

    // setup quantities for getRhoAlpha()
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
    amrex::ignore_unused(ilo, ihi, jlo, jhi);

    AMREX_ALWAYS_ASSERT( (klo == dom_lo.z) && (khi == dom_hi.z) );

    // Temporary FABs for tridiagonal solve (allocated on column)
    //   A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]
    amrex::FArrayBox RHS_fab, soln_fab, inv_coeffB_fab, coeffC_fab;
           RHS_fab.resize(bx,1, amrex::The_Async_Arena());
          soln_fab.resize(bx,1, amrex::The_Async_Arena());
    inv_coeffB_fab.resize(bx,1, amrex::The_Async_Arena());
        coeffC_fab.resize(bx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& inv_coeffB_a = inv_coeffB_fab.array();
    auto const& coeffC_a     =     coeffC_fab.array(); // upper diagonal

    int bc_comp = qty_index;

    auto dz_ptr = stretched_dz_d.data();

    bool neumann_on_zlo = (bc_ptr[bc_comp].lo(2) == ERFBCType::neumann);
    bool neumann_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::neumann);

#ifdef AMREX_USE_GPU
    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
#else
    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {
#endif
          Real rhoAlpha_lo, rhoAlpha_hi;
          Real dz_inv, dz_inv_lo, dz_inv_hi;

          // Bottom boundary coefficients and RHS for L decomp
          //===================================================
          Real a_tmp, b_tmp;
          {
              getRhoAlpha(i, j, klo, rhoAlpha_lo, rhoAlpha_hi,
                           cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                           prim_index, prim_scal_index, l_consA, l_turb);

              dz_inv    = 1.0 / dz_ptr[klo];
              dz_inv_lo = dz_inv;
              dz_inv_hi = 2.0 / (dz_ptr[klo] + dz_ptr[klo+1]);

              a_tmp             = 0.;
              coeffC_a(i,j,klo) = -implicit_fac * rhoAlpha_hi * dt * dz_inv * dz_inv_hi;

              b_tmp                 = cell_data(i,j,klo,Rho_comp) - a_tmp - coeffC_a(i,j,klo);
              inv_coeffB_a(i,j,klo) = 1.;

              RHS_a(i,j,klo) = cell_data(i,j,klo,n); // Note this is rho*theta, whereas solution will be theta
              if (use_SurfLayer) {
                  RHS_a(i,j,klo) +=  implicit_fac * dt * dz_inv * hfx_z(i,j,klo);
              } else if (neumann_on_zlo) {
                  RHS_a(i,j,klo) += -implicit_fac * dt * dz_inv * rhoAlpha_lo * bc_neumann_vals[2];
              }

              RHS_a(i,j,klo)    /= b_tmp; // NOTE: this is now "rho"
              coeffC_a(i,j,klo) /= b_tmp; // NOTE: this is now "gamma"
          }


          // Build the coefficients and RHS for L decomp
          //===================================================
          for (int k(klo+1); k < khi; k++)
          {
              getRhoAlpha(i, j, k, rhoAlpha_lo, rhoAlpha_hi,
                           cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                           prim_index, prim_scal_index, l_consA, l_turb);

              dz_inv    = 1.0 / dz_ptr[k];
              dz_inv_lo = 2.0 / (dz_ptr[k] + dz_ptr[k-1]);
              dz_inv_hi = 2.0 / (dz_ptr[k] + dz_ptr[k+1]);

              a_tmp           = -implicit_fac * rhoAlpha_lo * dt * dz_inv * dz_inv_lo;
              coeffC_a(i,j,k) = -implicit_fac * rhoAlpha_hi * dt * dz_inv * dz_inv_hi;

              b_tmp               = cell_data(i,j,k,Rho_comp) - a_tmp - coeffC_a(i,j,k);
              inv_coeffB_a(i,j,k) = 1. / (b_tmp - a_tmp * coeffC_a(i,j,k-1));

              RHS_a(i,j,k)     = cell_data(i,j,k,n); // NOTE: this is rho*theta, whereas solution will be theta

              RHS_a(i,j,k)     = (RHS_a(i,j,k) - a_tmp * RHS_a(i,j,k-1)) * inv_coeffB_a(i,j,k); // NOTE: This is now "rho"
              coeffC_a(i,j,k) *= inv_coeffB_a(i,j,k); // NOTE: this is now "gamma"
          } // k


          // Top boundary coefficients and RHS for L decomp
          //===================================================
          {
              getRhoAlpha(i, j, khi, rhoAlpha_lo, rhoAlpha_hi,
                           cell_data, mu_turb, d_alpha_eff, d_eddy_diff_idz,
                           prim_index, prim_scal_index, l_consA, l_turb);

              dz_inv    = 1.0 / dz_ptr[khi];
              dz_inv_lo = 2.0 / (dz_ptr[khi] + dz_ptr[khi-1]);
              dz_inv_hi = dz_inv;

              a_tmp             = -implicit_fac * rhoAlpha_lo * dt * dz_inv * dz_inv_lo;
              coeffC_a(i,j,khi) = 0.;

              b_tmp                 = cell_data(i,j,khi,Rho_comp) - a_tmp - coeffC_a(i,j,khi);
              inv_coeffB_a(i,j,khi) = 1. / (b_tmp - a_tmp * coeffC_a(i,j,khi-1));

              RHS_a(i,j,khi)    = cell_data(i,j,khi,n); // Note this is rho*theta, whereas solution will be theta
              if (neumann_on_zhi) {
                  RHS_a(i,j,khi) -= -implicit_fac * dt * dz_inv * rhoAlpha_hi * bc_neumann_vals[5];
              }

              soln_a(i,j,khi) = (RHS_a(i,j,khi) - a_tmp * RHS_a(i,j,khi-1)) * inv_coeffB_a(i,j,khi);
          }


          // Back sweep the U decomp solution
          //===================================================
          for (int k(khi-1); k>=klo; --k) {
              soln_a(i,j,k) = RHS_a(i,j,k) - coeffC_a(i,j,k) * soln_a(i,j,k+1);
          }


          // Transfer back to original array
          //===================================================
          for (int k(klo); k<=khi; ++k) {
              cell_data(i,j,k,n) = soln_a(i,j,k) * cell_data(i,j,k,Rho_comp);
          }

#ifdef AMREX_USE_GPU
    });
#else
      } // i
    } // j
#endif
}
