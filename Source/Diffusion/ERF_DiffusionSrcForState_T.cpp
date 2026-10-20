#include "ERF_Diffusion.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_NativeScalarDiffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_TerrainMetrics.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator with terrain-fitted coordinates.
 *
 * @param[in]  bx cell center box to loop over
 * @param[in]  domain box of the whole domain
 * @param[in]  start_comp starting component index
 * @param[in]  num_comp number of components
 * @param[in]  rotate flag to rotate terrain-aligned fluxes
 * @param[in]  cell_data conserved cell center vars
 * @param[in]  cell_prim primitive cell center vars
 * @param[out] cell_rhs RHS for cell center vars
 * @param[in]  xflux flux in x-dir
 * @param[in]  yflux flux in y-dir
 * @param[in]  zflux flux in z-dir
 * @param[in]  z_nd physical z height
 * @param[in]  z_cc cell-centered physical z height
 * @param[in]  ax area fraction of x-faces
 * @param[in]  ay area fraction of y-faces
 * @param[in]  detJ Jacobian determinant
 * @param[in]  cellSizeInv inverse cell size array
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[inout]  hfx_x heat flux in x-dir
 * @param[inout]  hfx_y heat flux in y-dir
 * @param[inout]  hfx_z heat flux in z-dir
 * @param[inout]  qfx1_x heat flux in x-dir
 * @param[inout]  qfx1_y heat flux in y-dir
 * @param[inout]  qfx1_z heat flux in z-dir
 * @param[out]    qfx2_z heat flux in z-dir
 * @param[in]  mu_turb turbulent viscosity
 * @param[in]  solverChoice container of solver and diffusion parameters
 * @param[in]  level AMR level
 * @param[in]  bc_ptr container with boundary conditions (device, read in kernels)
 * @param[in]  bc_ptr_h the same boundary conditions on the host, read here to set up the kernels
 * @param[in]  use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in]  implicit_fac -- factor of implicitness for vertical differences only
 * @param[in]  reflux_flx   if not null, face views (indexed by state component) to which the
 *                           fluxes of each component are added for refluxing
 */
void
DiffusionSrcForState_T (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const bool& rotate,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
                        const Array4<const Real>& z_nd,
                        const Array4<const Real>& z_cc,
                        const Array4<const Real>& ax,
                        const Array4<const Real>& ay,
                        const Array4<const Real>& /*az*/,
                        const Array4<const Real>& detJ,
                        const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
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
                        const Array4<const Real>& mu_turb,
                        const SolverChoice &solverChoice,
                        const int level,
                        const BCRec* bc_ptr,
                        const BCRec* bc_ptr_h,
                        const bool use_SurfLayer,
                        const Vector<std::unique_ptr<SurfaceLayer>>& SurfLayer,
                        const Real implicit_fac,
                        const ScalarDiffusionRefluxViews* reflux_flx)
{
    BL_PROFILE_VAR("DiffusionSrcForState_T()",DiffusionSrcForState_T);

    const Real explicit_fac = one - implicit_fac;

#include "ERF_SetupScalarDiffusion.H"

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
            flux_policy.prescribed.use_x_on_bottom_z = SurfLayer_zlo && rotate;
            flux_policy.prescribed.use_y_on_bottom_z = SurfLayer_zlo && rotate;
            flux_policy.prescribed.x = Array4<const Real>(hfx_x);
            flux_policy.prescribed.y = Array4<const Real>(hfx_y);
            flux_policy.prescribed.z = Array4<const Real>(hfx_z);
            flux_policy.diagnostic.enabled = true;
            flux_policy.diagnostic.z = hfx_z;
        } else if (native_policy.is_q1) {
            flux_policy.prescribed.use_x_on_side = true;
            flux_policy.prescribed.use_y_on_side = true;
            flux_policy.prescribed.use_z_on_side = true;
            flux_policy.prescribed.use_x_on_bottom_z = SurfLayer_zlo && rotate;
            flux_policy.prescribed.use_y_on_bottom_z = SurfLayer_zlo && rotate;
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

        BuildScalarDiffusionFluxes_T(
            bx, domain, field, flux_policy, rotate, z_nd, z_cc, cellSizeInv,
            mf_ux, mf_vy, bc_ptr, native_policy.bc_comp);
        if (native_policy.scale_raw_vertical_flux) {
            ScaleScalarDiffusionVerticalFlux(
                zbx, zflux, field.flux_comp, explicit_fac);
        }
        // erf.implicit_terrain_metric: the implicit solve of this quantity also takes the
        // compact terrain-metric term, so remove that part of it from the explicit flux
        if (native_policy.scale_raw_vertical_flux && solverChoice.implicit_terrain_metric &&
            implicit_fac > zero) {
            SubtractTerrainMetricImplicitPart_T(
                bx, domain, field, flux_policy, z_nd, cellSizeInv, mf_mx, mf_my, implicit_fac);
        }

        // No diffusive transfer through a bottom or top domain face where the surface layer
        // supplies the flux of this quantity separately, or where its boundary condition is
        // zero gradient (foextrap, as at a slip or no-slip wall without a prescribed value).
        // There the condition is no flux through the face along its normal, which is the whole
        // of Gz; the terrain-adjusted ghost cells only approximate that, so on sloping terrain
        // building Gz from them moved the scalar through the wall (erf-model/ERF#4330).
        const bool suppress_lower_z_flux = SurfLayer_zlo
            ? !(native_policy.is_theta || native_policy.is_q1)
            : (bc_ptr_h[native_policy.bc_comp].lo(2) == ERFBCType::foextrap);
        const bool suppress_upper_z_flux =
            !SurfLayer_zhi && (bc_ptr_h[native_policy.bc_comp].hi(2) == ERFBCType::foextrap);
        ApplyScalarDiffusionMappedDivergence_T(
            bx, domain, field, z_nd, ax, ay, detJ, cellSizeInv, mf_mx, mf_uy,
            mf_my, mf_vx, suppress_lower_z_flux, suppress_upper_z_flux);
        if (reflux_flx) {
            AddScalarDiffusionRefluxFluxes_T(
                bx, domain, xflux, yflux, zflux, field.flux_comp, *reflux_flx, field.rhs_comp,
                z_nd, ax, ay, cellSizeInv, mf_mx, mf_my, mf_uy, mf_vx,
                suppress_lower_z_flux, suppress_upper_z_flux);
        }
    } // n
}
