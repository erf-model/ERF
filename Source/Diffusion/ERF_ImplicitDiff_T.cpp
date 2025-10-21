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
ImplicitDiffForState_T (const Box& bx, const Box& domain, const Real dt,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<      Real>& cell_data,
                        const Array4<      Real>& zflux,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& detJ,
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
    BL_PROFILE_VAR("ImplicitDiffForState_T()",ImplicitDiffForState_T);

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


      } // i
    } // j
}
