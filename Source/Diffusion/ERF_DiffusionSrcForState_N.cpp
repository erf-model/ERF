#include "ERF_Diffusion.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_NativeScalarDiffusion.H"
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
 * @param[in]  mf_mx x map factor at cell centers
 * @param[in]  mf_ux x map factor at x-faces
 * @param[in]  mf_vx x map factor at y-faces
 * @param[in]  mf_my y map factor at cell centers
 * @param[in]  mf_uy y map factor at x-faces
 * @param[in]  mf_vy y map factor at y-faces
 * @param[inout]  hfx_z heat flux in z-dir
 * @param[inout]  qfx1_z heat flux in z-dir
 * @param[out]    qfx2_z heat flux in z-dir
 * @param[in]  mu_turb turbulent viscosity
 * @param[in]  solverChoice container of solver and diffusion parameters
 * @param[in]  level AMR level
 * @param[in]  bc_ptr container with boundary conditions
 * @param[in]  use_SurfLayer whether we have turned on subgrid diffusion
 * @param[in]  implicit_fac -- factor of implicitness for vertical differences only
 * @param[in]  reflux_flx   if not null, face views (indexed by state component) to which the
 *                           fluxes of each component are added for refluxing
 */
void
DiffusionSrcForState_N (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
                        const Array4<const Real>& cell_data,
                        const Array4<const Real>& cell_prim,
                        const Array4<Real>& cell_rhs,
                        const Array4<Real>& xflux,
                        const Array4<Real>& yflux,
                        const Array4<Real>& zflux,
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
                        const bool use_SurfLayer,
                        const Vector<std::unique_ptr<SurfaceLayer>>& SurfLayer,
                        const Real implicit_fac,
                        const ScalarDiffusionRefluxViews* reflux_flx,
                        const Array4<const Real>& wall_face33,
                        const Array4<const Real>& wall_hfx,
                        const Array4<const Real>& wall_qfx)
{
    BL_PROFILE_VAR("DiffusionSrcForState_N()",DiffusionSrcForState_N);

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

        BuildScalarDiffusionFluxes_N(
            bx, domain, field, flux_policy, cellSizeInv, u, v, mf_ux, mf_uy,
            mf_vx, mf_vy, bc_ptr, native_policy.bc_comp);
        if (native_policy.scale_raw_vertical_flux) {
            ScaleScalarDiffusionVerticalFlux(
                zbx, zflux, field.flux_comp, explicit_fac);
        }
        // Fraction-stress immersed wall law: no diffusive flux of any scalar (theta, TKE, the
        // passive scalars, the moisture species) through the bottom face of the wall cell, as
        // through zlo on flat ground (ERF_ImmersedWallCell.H). The surface carries the wall heat
        // flux and the wall moisture flux of the wall law (below), both applied as sources in the
        // wall cell by the immersed scalar forcing.
        if (wall_face33) {
            const int fc = field.flux_comp;
            const Box mbx = zbx & amrex::convert(Box(wall_face33), zbx.ixType());
            ParallelFor(mbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                zflux(i,j,k,fc) *= wall_face33(i,j,k);
            });
            // the qv and qc face-flux diagnostics were filled from the open face: close them the
            // same way (the TKE buoyancy sources read qfx1_z); on the wall face the water vapour
            // flux is the wall moisture flux of the wall law (rho times the kinematic flux, as the
            // surface layer writes it at zlo), zero without a surface moisture condition; the
            // cloud-water face is zero. wall_qfx has no ghost cells (see wall_hfx below).
            if ((native_policy.is_q1 && qfx1_z) || (native_policy.is_q2 && qfx2_z)) {
                const Array4<Real>& qfx_z = native_policy.is_q1 ? qfx1_z : qfx2_z;
                const Array4<const Real> q_wall_arr = native_policy.is_q1 ? wall_qfx : Array4<const Real>{};
                const Box qbx = mbx & amrex::convert(Box(qfx_z), zbx.ixType());
                ParallelFor(qbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const Real f_wall = (q_wall_arr && q_wall_arr.contains(i,j,k))
                                      ? cell_data(i,j,k,Rho_comp) * q_wall_arr(i,j,k) : zero;
                    qfx_z(i,j,k) = (wall_face33(i,j,k) == zero) ? f_wall : wall_face33(i,j,k) * qfx_z(i,j,k);
                });
            }
            if (qty_index == RhoTheta_comp && hfx_z) {
                // the TKE buoyancy reads hfx_z: on the wall face it is the wall heat flux of the
                // wall law (rho times the kinematic flux, as the surface layer writes it at zlo).
                // wall_hfx has no ghost cells: on the top face of a box split in z (never with
                // amr.no_box_split_dir = 2) the face reads zero here and its value from the box above
                const Box hbx = mbx & amrex::convert(Box(hfx_z), zbx.ixType());
                ParallelFor(hbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const Real q_wall = (wall_hfx && wall_hfx.contains(i,j,k))
                                      ? cell_data(i,j,k,Rho_comp) * wall_hfx(i,j,k) : zero;
                    hfx_z(i,j,k) = (wall_face33(i,j,k) == zero) ? q_wall : wall_face33(i,j,k) * hfx_z(i,j,k);
                });
            }
        }
        ApplyScalarDiffusionFluxDivergence_N(
            bx, xflux, yflux, zflux, field.flux_comp, cell_rhs, field.rhs_comp,
            cellSizeInv, mf_mx, mf_my);
        if (reflux_flx) {
            AddScalarDiffusionRefluxFluxes_N(
                bx, xflux, yflux, zflux, field.flux_comp, *reflux_flx, field.rhs_comp, mf_mx, mf_my);
        }
    } // n
}
