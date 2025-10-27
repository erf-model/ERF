#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_PBLModels.H"

using namespace amrex;

/**
 * Function for computing the implicit contribution to the vertical diffusion
 * of theta, with a vertically stretched grid over flat terrain.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] bc_neumann_vals values of derivatives if bc_type == Neumann
 * @param[inout] cell_data conserved cell-centered rho, rho theta
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
                        const Array4<const Real>& hfx_z,
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const BCRec* bc_ptr,
                        const bool use_SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForState_S()",ImplicitDiffForState_S);

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

    auto dz_ptr = stretched_dz_d.data();

    bool foextrap_on_zlo = (bc_ptr[bc_comp].lo(2) == ERFBCType::foextrap);
    bool foextrap_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::foextrap);
    bool neumann_on_zlo  = (bc_ptr[bc_comp].lo(2) == ERFBCType::neumann);
    bool neumann_on_zhi  = (bc_ptr[bc_comp].hi(2) == ERFBCType::neumann);

    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zlo || neumann_on_zlo || use_SurfLayer,
                              "Unexpected lower BC used with implicit vertical diffusion");
    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zhi || neumann_on_zhi,
                              "Unexpected upper BC used with implicit vertical diffusion");

#ifdef AMREX_USE_GPU
    ParallelFor(makeSlab(bx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
#else
    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {
#endif
        // Build the coefficients and RHS
        for (int k(klo); k <= khi; k++)
        {
#include "ERF_GetRhoAlpha.H"

            Real dz_inv = 1.0 / dz_ptr[k];
            Real dz_inv_lo = (k == dom_lo.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k-1]);
            Real dz_inv_hi = (k == dom_hi.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k+1]);

            RHS_a(i,j,k)  = cell_data(i,j,k,n); // Note this is rho*theta, whereas solution will be theta

            coeffA_a(i,j,k) = -implicit_fac * rhoAlpha_lo * dt * dz_inv * dz_inv_lo;
            coeffC_a(i,j,k) = -implicit_fac * rhoAlpha_hi * dt * dz_inv * dz_inv_hi;

            if (k == dom_lo.z) {
                if (use_SurfLayer) {
                    RHS_a(i,j,klo) += implicit_fac * dt * dz_inv * hfx_z(i,j,0);
                } else if (neumann_on_zlo) {
                    RHS_a(i,j,klo) += coeffA_a(i,j,klo) * bc_neumann_vals[2] / dz_inv_lo;
                }

                coeffA_a(i,j,klo) = 0.;
            }
            if (k == dom_hi.z) {
                if (neumann_on_zhi) {
                    RHS_a(i,j,khi) -= coeffC_a(i,j,khi) * bc_neumann_vals[5] / dz_inv_hi;
                }

                coeffC_a(i,j,khi) = 0.;
            }

            coeffB_a(i,j,k) = cell_data(i,j,k,Rho_comp) - coeffA_a(i,j,k) - coeffC_a(i,j,k);
        } // k

#include "ERF_SolveTridiag.H"
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
