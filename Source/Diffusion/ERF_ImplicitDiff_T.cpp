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
}
