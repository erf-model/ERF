#include "ERF_Diffusion.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_NativeScalarDiffusion.H"
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
 * @param[in]  stretched_dz_d array of vertical grid spacings
 * @param[in]  cellSizeInv inverse cell size array
 * @param[in]  SmnSmn_a strain rate magnitude
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[inout]  hfx_z heat flux in z-dir
 * @param[inout]  qfx1_z heat flux in z-dir
 * @param[out]    qfx2_z heat flux in z-dir
 * @param[in]  diss dissipation of TKE
 * @param[in]  mu_turb turbulent viscosity
 * @param[in]  solverChoice container of solver and diffusion parameters
 * @param[in]  level AMR level
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
                        const Vector<std::unique_ptr<SurfaceLayer>>& SurfLayer,
                        const Real implicit_fac)
{
    BL_PROFILE_VAR("DiffusionSrcForState_S()",DiffusionSrcForState_S);

    const Real explicit_fac = one - implicit_fac;

#include "ERF_SetupScalarDiffusion.H"
    Real l_abs_g      = std::abs(grav_gpu[2]);

    int klo = domain.smallEnd(2);
    int khi = domain.bigEnd(2);
    auto dz_ptr = stretched_dz_d.data();

    for (int n(0); n<num_comp; ++n) {
        const int qty_index = start_comp + n;
        const NativeScalarDiffusionPolicy native_policy =
            ResolveNativeScalarDiffusionPolicy(qty_index, diffChoice);

        ScalarDiffusionFieldViews field;
        field.scalar = cell_prim;
        field.scalar_comp = native_policy.scalar_comp;
        field.density = cell_data;
        field.rho_comp = Rho_comp;
        field.mu_turb = mu_turb;
        field.xflux = xflux;
        field.yflux = yflux;
        field.zflux = zflux;
        field.flux_comp = native_policy.flux_comp;
        field.rhs = cell_rhs;
        field.rhs_comp = native_policy.rhs_comp;

        ScalarDiffusionFluxPolicy flux_policy;
        flux_policy.coefficients = native_policy.coefficients;
        flux_policy.coefficient_mode = {l_consA, l_turb};
        flux_policy.surface = {
            SurfLayer_xlo, SurfLayer_xhi, SurfLayer_ylo, SurfLayer_yhi,
            SurfLayer_zlo, SurfLayer_zhi
        };
        if (native_policy.is_theta) {
            flux_policy.prescribed.use_x_on_side = true;
            flux_policy.prescribed.use_y_on_side = true;
            flux_policy.prescribed.use_z_on_side = true;
            flux_policy.prescribed.x = Array4<const Real>(hfx_x);
            flux_policy.prescribed.y = Array4<const Real>(hfx_y);
            flux_policy.prescribed.z = Array4<const Real>(hfx_z);
            flux_policy.diagnostic.enabled = true;
            flux_policy.diagnostic.z = hfx_z;
        } else if (native_policy.is_q1) {
            flux_policy.prescribed.use_x_on_side = true;
            flux_policy.prescribed.use_y_on_side = true;
            flux_policy.prescribed.use_z_on_side = true;
            flux_policy.prescribed.x = Array4<const Real>(qfx1_x);
            flux_policy.prescribed.y = Array4<const Real>(qfx1_y);
            flux_policy.prescribed.z = Array4<const Real>(qfx1_z);
            flux_policy.diagnostic.enabled = true;
            flux_policy.diagnostic.z = qfx1_z;
        } else if (native_policy.is_q2) {
            flux_policy.diagnostic.enabled = true;
            flux_policy.diagnostic.write_on_surface = true;
            flux_policy.diagnostic.z = qfx2_z;
        }

        const Real dx_inv = cellSizeInv[0];
        const Real dy_inv = cellSizeInv[1];
        BuildScalarDiffusionFluxes_S(
            bx, domain, field, flux_policy, dx_inv, dy_inv, dz_ptr, klo, khi,
            mf_ux, mf_uy, mf_vx, mf_vy, bc_ptr, native_policy.bc_comp);
        if (native_policy.scale_raw_vertical_flux) {
            ScaleScalarDiffusionVerticalFlux(
                zbx, zflux, field.flux_comp, explicit_fac);
        }
        ApplyScalarDiffusionFluxDivergence_S(
            bx, xflux, yflux, zflux, field.flux_comp, cell_rhs, field.rhs_comp,
            mf_mx, mf_my, dx_inv, dy_inv, dz_ptr);
    } // n

    const PBLDerivativeDzInv_S pbl_derivative_dz_inv{dz_ptr, klo, khi};
#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
