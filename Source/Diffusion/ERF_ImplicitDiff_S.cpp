#include "ERF_Diffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_SolveTridiag.H"

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

#include "ERF_SetupVertDiff.H"

    // define get_rhoAlpha()
    const int         n = RhoTheta_comp;
    const int qty_index = RhoTheta_comp;
    const int prim_index = qty_index - 1;
    const int prim_scal_index = (qty_index >= RhoScalar_comp && qty_index < RhoScalar_comp+NSCALARS) ? PrimScalar_comp : prim_index;
#include "ERF_GetRhoAlpha.H"

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
            Real rhoAlpha_lo, rhoAlpha_hi;
            get_rhoAlpha(i, j, k, rhoAlpha_lo, rhoAlpha_hi);

            Real dz_inv = 1.0 / dz_ptr[k];
            Real dz_inv_lo = (k == dom_lo.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k-1]);
            Real dz_inv_hi = (k == dom_hi.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k+1]);

            RHS_a(i,j,k)  = cell_data(i,j,k,n); // Note this is rho*theta, whereas solution will be theta

            // This represents the cell-centered finite difference of two
            // face-centered finite differences (hi and lo)
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

        SolveTridiag(i,j,klo,khi,soln_a,coeffA_a,coeffB_a,inv_coeffB_a,coeffC_a,RHS_a);
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

/**
 * Function for computing the implicit contribution to the vertical diffusion
 * of momentum, with a vertically stretched grid over flat terrain.
 *
 * This function (explicitly instantiated below) handles staggering in x or y
 * through the template parameter, stagdir.
 *
 * @param[in   ] bx cell-centered box to loop over
 * @param[in   ] domain box of the whole domain
 * @param[in   ] dt time step
 * @param[in   ] cell_data conserved cell-centered rho, rho theta
 * @param[inout] face_data conserved momentum
 * @param[in   ] stretched_dz_d array over z of dz[k]
 * @param[inout] hfx_z heat flux in z-dir
 * @param[in   ] mu_turb turbulent viscosity
 * @param[in   ] solverChoice container of parameters
 * @param[in   ] bc_ptr container with boundary conditions
 * @param[in   ] use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in   ] implicit_fac if 1 then fully implicit; if 0 then fully explicit
 */
template <int stagdir>
void
ImplicitDiffForMom_S (const Box& bx,
                      const Box& domain,
                      const int level,
                      const Real dt,
                      const Array4<const Real>& cell_data,
                      const Array4<      Real>& face_data,
                      const Array4<const Real>& tau_corr,
                      const Gpu::DeviceVector<Real>& stretched_dz_d,
                      const Array4<const Real>& mu_turb,
                      const SolverChoice &solverChoice,
                      const BCRec* bc_ptr,
                      const bool use_SurfLayer,
                      const Real implicit_fac)
{
    BL_PROFILE_VAR("ImplicitDiffForMom_S()",ImplicitDiffForMom_S);

#include "ERF_SetupVertDiff.H"

    // g(S*) coefficient
    constexpr Real gfac = (stagdir == 2) ? 4.0/3.0 : 0.5;

    // offsets used to average to faces
    constexpr int ioff = (stagdir == 0) ? 1 : 0;
    constexpr int joff = (stagdir == 1) ? 1 : 0;

    // define get_rhoAlpha()
#include "ERF_GetRhoAlphaAtFaces.H"

    // Box bounds
    const Box bxx = convert(bx, IntVect(ioff, joff, 0));
    int ilo = bxx.smallEnd(0);
    int ihi = bxx.bigEnd(0);
    int jlo = bxx.smallEnd(1);
    int jhi = bxx.bigEnd(1);
    int klo = bxx.smallEnd(2);
    int khi = bxx.bigEnd(2);

    // Temporary FABs for tridiagonal solve (allocated on column)
    //   A[k] * x[k-1] + B[k] * x[k] + C[k+1] = RHS[k]
    amrex::FArrayBox RHS_fab, soln_fab, coeffA_fab, coeffB_fab, inv_coeffB_fab, coeffC_fab;
           RHS_fab.resize(bxx,1, amrex::The_Async_Arena());
          soln_fab.resize(bxx,1, amrex::The_Async_Arena());
        coeffA_fab.resize(bxx,1, amrex::The_Async_Arena());
        coeffB_fab.resize(bxx,1, amrex::The_Async_Arena());
    inv_coeffB_fab.resize(bxx,1, amrex::The_Async_Arena());
        coeffC_fab.resize(bxx,1, amrex::The_Async_Arena());
    auto const& RHS_a        =        RHS_fab.array();
    auto const& soln_a       =       soln_fab.array();
    auto const& coeffA_a     =     coeffA_fab.array(); // lower diagonal
    auto const& coeffB_a     =     coeffB_fab.array(); // diagonal
    auto const& inv_coeffB_a = inv_coeffB_fab.array();
    auto const& coeffC_a     =     coeffC_fab.array(); // upper diagonal

    auto dz_ptr = stretched_dz_d.data();

    int bc_comp = BCVars::xvel_bc + stagdir;


  //bool ext_dir_on_zlo  = (bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir ||
  //                        bc_ptr[bc_comp].lo(2) == ERFBCType::ext_dir_prim);
    bool foextrap_on_zhi = (bc_ptr[bc_comp].hi(2) == ERFBCType::foextrap);

  //AMREX_ASSERT_WITH_MESSAGE(ext_dir_on_zlo || use_SurfLayer,
    AMREX_ASSERT_WITH_MESSAGE(use_SurfLayer,
                              "Unexpected lower BC used with implicit vertical diffusion");
    AMREX_ASSERT_WITH_MESSAGE(foextrap_on_zhi,
                              "Unexpected upper BC used with implicit vertical diffusion");

#ifdef AMREX_USE_GPU
    ParallelFor(makeSlab(bxx,2,0), [=] AMREX_GPU_DEVICE (int i, int j, int)
    {
#else
    for (int j(jlo); j<=jhi; ++j) {
      for (int i(ilo); i<=ihi; ++i) {
#endif
        // Build the coefficients and RHS
        for (int k(klo); k <= khi; k++)
        {
            // Note: either ioff or joff are 1
            Real rhoface = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i-ioff,j-joff,k,Rho_comp));

            Real rhoAlpha_lo, rhoAlpha_hi;
            get_rhoAlpha(i, j, k, rhoAlpha_lo, rhoAlpha_hi);

            Real dz_inv = 1.0 / dz_ptr[k];
            Real dz_inv_lo = (k == dom_lo.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k-1]);
            Real dz_inv_hi = (k == dom_hi.z) ? dz_inv
                                             : 2.0 / (dz_ptr[k] + dz_ptr[k+1]);

            // Face data currently holds the _fully_ explicit solution, which
            // will be used to determine the velocity gradient for the bottom
            // BC
            RHS_a(i,j,k) = face_data(i,j,k); // Note this is momentum but solution will be velocity

            // Notes:
            //
            // - In DiffusionSrcForMom (e.g., for x-mom)
            //
            //     Real diffContrib = ...
            //                      + (tau13(i,j,k+1) - tau13(i,j,k)) / dz_ptr[k]
            //     rho_u_rhs(i,j,k) -= diffContrib;  // note the negative sign
            //
            // - We need to scale the explicit _part_ of `tau13` (for x-mom) by (1 - implicit_fac)
            //   The part that needs to be scaled is stored in `tau_corr`.
            //   E.g., tau13 = 0.5 * (du/dz + dw/dx)
            //         tau13_corr = 0.5 * du/dz
            //
            // - The momentum (`face_data`) was set to `S_old + S_rhs * dt`
            //   prior to including "ERF_Implicit.H". Recall that S_rhs includes
            //   sources from advection and other forcings, not just diffusion.
            //
            // - To correct momentum, we need to subtract `implicit_fac * diffContrib_corr`
            //   from S_rhs to recover `(1 - implicit_fac) * diffContrib_corr`,
            //   where `diffContrib_corr = -d(tau_corr)/dz`. The negative sign
            //   comes from our convention for the RHS diffusion source.
            //
            //   Subtracting a negative gives the += below; multiply by dt to
            //   get the intermediate momentum on the RHS of the tridiagonal
            //   system.
            RHS_a(i,j,k) += implicit_fac * (tau_corr(i,j,k+1) - tau_corr(i,j,k))*dz_inv * dt;

            // This represents the face-centered finite difference of two
            // edge-centered finite differences (hi and lo)
            coeffA_a(i,j,k) = -implicit_fac * gfac * rhoAlpha_lo * dt * dz_inv * dz_inv_lo;
            coeffC_a(i,j,k) = -implicit_fac * gfac * rhoAlpha_hi * dt * dz_inv * dz_inv_hi;

            // Setup BCs
            if (k == dom_lo.z) {
                if (use_SurfLayer) {
                    //RHS_a(i,j,klo) += implicit_fac * dt * dz_inv * hfx_z(i,j,0);
                    //RHS_a(i,j,klo) += coeffA_a(i,j,klo) * bc_neumann_vals[2] / dz_inv_lo;
                    Real uhi = 2.0 * face_data(i,j,klo  ) / (cell_data(i,j,klo  ,Rho_comp) + cell_data(i,j,klo-1,Rho_comp));
                    Real ulo = 2.0 * face_data(i,j,klo-1) / (cell_data(i,j,klo-1,Rho_comp) + cell_data(i,j,klo-2,Rho_comp));
                    RHS_a(i,j,klo) += coeffA_a(i,j,klo) * (uhi - ulo);
                }
                coeffA_a(i,j,klo) = 0.;
            }
            if (k == dom_hi.z) {
                coeffC_a(i,j,khi) = 0.; // foextrap
            }

            coeffB_a(i,j,k) = rhoface - coeffA_a(i,j,k) - coeffC_a(i,j,k);
        } // k

        SolveTridiag(i,j,klo,khi,soln_a,coeffA_a,coeffB_a,inv_coeffB_a,coeffC_a,RHS_a);
        for (int k(klo); k<=khi; ++k) {
            Real rhoface = 0.5 * (cell_data(i,j,k,Rho_comp) + cell_data(i-ioff,j-joff,k,Rho_comp));
            face_data(i,j,k) = rhoface * soln_a(i,j,k);
        }

#ifdef AMREX_USE_GPU
    });
#else
      } // i
    } // j
#endif
}


#define INSTANTIATE_IMPLICIT_DIFF_FOR_MOM(STAGDIR) \
    template void ImplicitDiffForMom_S<STAGDIR> ( \
        const Box&, \
        const Box&, \
        const int, \
        const Real, \
        const Array4<const Real>&, \
        const Array4<      Real>&, \
        const Array4<const Real>&, \
        const Gpu::DeviceVector<Real>&, \
        const Array4<const Real>&, \
        const SolverChoice&, \
        const BCRec*, \
        const bool, \
        const Real);
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM(0)
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM(1)
INSTANTIATE_IMPLICIT_DIFF_FOR_MOM(2)
#undef INSTANTIATE_IMPLICIT_DIFF_FOR_MOM
