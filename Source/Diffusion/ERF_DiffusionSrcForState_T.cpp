#include "ERF_Diffusion.H"
#include "ERF_ScalarDiffusion.H"
#include "ERF_EddyViscosity.H"
#include "ERF_TerrainMetrics.H"
#include "ERF_PBLModels.H"

using namespace amrex;

/**
 * Function for computing the scalar RHS for diffusion operator with terrain-fitted coordinates.
 *
 * @param[in]  bx cell center box to loop over
 * @param[in]  domain box of the whole domain
 * @param[in]  start_comp starting component index
 * @param[in]  num_comp number of components
 * @param[in]  rotate flag to rotate terrain-aligned fluxes
 * @param[in]  u velocity in x-dir
 * @param[in]  v velocity in y-dir
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
 * @param[in]  SmnSmn_a strain rate magnitude
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
DiffusionSrcForState_T (const Box& bx, const Box& domain,
                        int start_comp, int num_comp,
                        const bool& rotate,
                        const Array4<const Real>& u,
                        const Array4<const Real>& v,
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
    BL_PROFILE_VAR("DiffusionSrcForState_T()",DiffusionSrcForState_T);

    const Real explicit_fac = one - implicit_fac;

#include "ERF_SetupScalarDiffusion.H"
    Real l_abs_g      = std::abs(grav_gpu[2]);

    const Real dz_inv = cellSizeInv[2];

    // We need to grow these boxes in the vertical direction when tiling so that we can access xflux and yflux
    //    to modify zflux
    Box xbx_g1(xbx); Box ybx_g1(ybx);
    if (xbx_g1.smallEnd(2) != dom_lo.z) xbx_g1.growLo(2,1);
    if (ybx_g1.smallEnd(2) != dom_lo.z) ybx_g1.growLo(2,1);
    if (xbx_g1.bigEnd(2)   != dom_hi.z) xbx_g1.growHi(2,1);
    if (ybx_g1.bigEnd(2)   != dom_hi.z) ybx_g1.growHi(2,1);

    for (int n(0); n<num_comp; ++n) {
        const int qty_index = start_comp + n;
        const NativeScalarDiffusionPolicy native_policy =
            ResolveNativeScalarDiffusionPolicy(qty_index, diffChoice);
        const int scalar_comp = native_policy.scalar_comp;
        const int rhs_comp = native_policy.rhs_comp;
        constexpr int flux_comp = 0;
        const int rho_comp = Rho_comp;
        const Array4<const Real>& scalar = cell_prim;

    // Constant alpha & Turb model
    if (l_consA && l_turb) {
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,native_policy.coefficients.eddy_h_comp);

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = one / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i-1, j, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i-1, j, k-1, scalar_comp)*idz_lo );
            Real GradCx = dx_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i-1, j, k  , scalar_comp) );

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Tx(rhoAlpha,mf_ux(i,j,0),
                                                               GradCx,met_h_xi,GradCz);
            }

        });
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,native_policy.coefficients.eddy_h_comp);

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = one / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i, j-1, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i, j-1, k-1, scalar_comp)*idz_lo );
            Real GradCy = dy_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i, j-1, k  , scalar_comp) );

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Ty(rhoAlpha,mf_vy(i,j,0),
                                                               GradCy,met_h_eta,GradCz);
            }

        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,0,1,native_policy.coefficients.eddy_v_comp);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir_prim) )
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && k == dom_lo.z);
            bool SurfLayer_on_zhi = ( SurfLayer_zhi && k == dom_hi.z + 1);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * ( c1 * scalar(i, j, k-1, scalar_comp)
                                + c2 * scalar(i, j, k  , scalar_comp)
                                + c3 * scalar(i, j, k+1, scalar_comp) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * (  -( c1 * scalar(i, j, k  , scalar_comp)
                                    + c2 * scalar(i, j, k-1, scalar_comp)
                                    + c3 * scalar(i, j, k-2, scalar_comp) ) );
            } else {
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
                GradCz = (dz_inv/met_h_zeta) * ( scalar(i, j, k, scalar_comp) - scalar(i, j, k-1, scalar_comp) );
            }

            if (SurfLayer_on_zlo || SurfLayer_on_zhi) {
                if (native_policy.is_theta) {
                    zflux(i,j,k,flux_comp) = hfx_z(i,j,k);
                } else if (native_policy.is_q1) {
                    zflux(i,j,k,flux_comp) = qfx1_z(i,j,k);
                } else {
                    zflux(i,j,k,flux_comp) = zero;
                }
            } else {
                zflux(i,j,k,flux_comp) = -rhoAlpha * GradCz;
            }

            if (native_policy.is_theta) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    hfx_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q1) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    qfx1_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q2) {
                qfx2_z(i,j,k) = zflux(i,j,k,flux_comp);
            }
        });
    // Constant rho*alpha & Turb model
    } else if (l_turb) {
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,native_policy.coefficients.eddy_h_comp);

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = one / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i-1, j, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i-1, j, k-1, scalar_comp)*idz_lo );
            Real GradCx = dx_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i-1, j, k  , scalar_comp) );

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Tx(rhoAlpha,mf_ux(i,j,0),
                                                               GradCx,met_h_xi,GradCz);
            }

        });
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,native_policy.coefficients.eddy_h_comp);

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = one / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i, j-1, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i, j-1, k-1, scalar_comp)*idz_lo );
            Real GradCy = dy_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i, j-1, k  , scalar_comp) );

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Ty(rhoAlpha,mf_vy(i,j,0),
                                                               GradCy,met_h_eta,GradCz);
            }

        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,0,1,native_policy.coefficients.eddy_v_comp);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && k == dom_lo.z);
            bool SurfLayer_on_zhi = ( SurfLayer_zhi && k == dom_hi.z + 1);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * ( c1 * scalar(i, j, k-1, scalar_comp)
                                + c2 * scalar(i, j, k  , scalar_comp)
                                + c3 * scalar(i, j, k+1, scalar_comp) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * (  -( c1 * scalar(i, j, k  , scalar_comp)
                                    + c2 * scalar(i, j, k-1, scalar_comp)
                                    + c3 * scalar(i, j, k-2, scalar_comp) ) );
            } else {
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
                GradCz = (dz_inv/met_h_zeta) * ( scalar(i, j, k, scalar_comp) - scalar(i, j, k-1, scalar_comp) );
            }

            if (SurfLayer_on_zlo || SurfLayer_on_zhi) {
                if (native_policy.is_theta) {
                    zflux(i,j,k,flux_comp) = hfx_z(i,j,k);
                } else if (native_policy.is_q1) {
                    zflux(i,j,k,flux_comp) = qfx1_z(i,j,k);
                } else {
                    zflux(i,j,k,flux_comp) = zero;
                }
            } else {
                zflux(i,j,k,flux_comp) = -rhoAlpha * GradCz;
            }

            if (native_policy.is_theta) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    hfx_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q1) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    qfx1_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q2) {
                qfx2_z(i,j,k) = zflux(i,j,k,flux_comp);
            }
        });
    // Constant alpha & no LES/PBL model
    } else if(l_consA) {
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,0);

            Real met_h_xi   = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = one / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i-1, j, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i-1, j, k-1, scalar_comp)*idz_lo );
            Real GradCx = dx_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i-1, j, k  , scalar_comp) );

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Tx(rhoAlpha,mf_ux(i,j,0),
                                                               GradCx,met_h_xi,GradCz);
            }

        });
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,0);

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = one / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i, j-1, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i, j-1, k-1, scalar_comp)*idz_lo );
            Real GradCy = dy_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i, j-1, k  , scalar_comp) );

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Ty(rhoAlpha,mf_vy(i,j,0),
                                                               GradCy,met_h_eta,GradCz);
            }

        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,0,1,0);

            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && k == dom_lo.z);
            bool SurfLayer_on_zhi = ( SurfLayer_zhi && k == dom_hi.z + 1);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * ( c1 * scalar(i, j, k-1, scalar_comp)
                                + c2 * scalar(i, j, k  , scalar_comp)
                                + c3 * scalar(i, j, k+1, scalar_comp) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * (  -( c1 * scalar(i, j, k  , scalar_comp)
                                    + c2 * scalar(i, j, k-1, scalar_comp)
                                    + c3 * scalar(i, j, k-2, scalar_comp) ) );
            } else {
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
                GradCz = (dz_inv/met_h_zeta) * ( scalar(i, j, k, scalar_comp) - scalar(i, j, k-1, scalar_comp) );
            }

            if (SurfLayer_on_zlo || SurfLayer_on_zhi) {
                if (native_policy.is_theta) {
                    zflux(i,j,k,flux_comp) = hfx_z(i,j,k);
                } else if (native_policy.is_q1) {
                    zflux(i,j,k,flux_comp) = qfx1_z(i,j,k);
                } else {
                    zflux(i,j,k,flux_comp) = zero;
                }
            } else {
                zflux(i,j,k,flux_comp) = -rhoAlpha * GradCz;
            }

            if (native_policy.is_theta) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    hfx_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q1) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    qfx1_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else if (native_policy.is_q2) {
                qfx2_z(i,j,k) = zflux(i,j,k,flux_comp);
            }
        });
    // Constant rho*alpha & no LES/PBL model
    } else {
        ParallelFor(xbx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients, i,j,k,0,0,0,0);

            Real met_h_xi = Compute_h_xi_AtIface  (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i  ,j,k+1) - z_cc(i  ,j,k-1));
            Real idz_lo = one / (z_cc(i-1,j,k+1) - z_cc(i-1,j,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i-1, j, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i-1, j, k-1, scalar_comp)*idz_lo );
            Real GradCx = dx_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i-1, j, k  , scalar_comp) );

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,0);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Tx(rhoAlpha,mf_ux(i,j,0),
                                                               GradCx,met_h_xi,GradCz);
            }

        });
        ParallelFor(ybx_g1, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients, i,j,k,0,0,0,0);

            Real met_h_eta  = Compute_h_eta_AtJface (i,j,k,cellSizeInv,z_nd);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && rotate && k == dom_lo.z);

            Real idz_hi = one / (z_cc(i,j  ,k+1) - z_cc(i,j  ,k-1));
            Real idz_lo = one / (z_cc(i,j-1,k+1) - z_cc(i,j-1,k-1));
            Real GradCz =    myhalf * ( scalar(i, j, k+1, scalar_comp)*idz_hi + scalar(i, j-1, k+1, scalar_comp)*idz_lo
                                   - scalar(i, j, k-1, scalar_comp)*idz_hi - scalar(i, j-1, k-1, scalar_comp)*idz_lo );
            Real GradCy = dy_inv * ( scalar(i, j, k  , scalar_comp)        - scalar(i, j-1, k  , scalar_comp) );

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else if (SurfLayer_on_zlo && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,0);
            } else if (SurfLayer_on_zlo && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,0);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_Ty(rhoAlpha,mf_vy(i,j,0),
                                                               GradCy,met_h_eta,GradCz);
            }

        });
        ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients, i,j,k,0,0,0,0);


            Real GradCz;
            bool ext_dir_on_zlo = ( ((bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].lo(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_lo.z);
            bool ext_dir_on_zhi = ( ((bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir) ||
                                     (bc_ptr[native_policy.bc_comp].hi(2) == ERFBCType::ext_dir_prim))
                                    && k == dom_hi.z+1);
            bool SurfLayer_on_zlo = ( SurfLayer_zlo && k == dom_lo.z);
            bool SurfLayer_on_zhi = ( SurfLayer_zhi && k == dom_hi.z + 1);

            if (ext_dir_on_zlo) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k+1,z_nd);
                Real dz0  = zm - Compute_Z_AtWFace(i,j,k,z_nd);
                Real dz1  = Compute_Z_AtWFace(i,j,k+2,z_nd) - zm;
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * ( c1 * scalar(i, j, k-1, scalar_comp)
                                + c2 * scalar(i, j, k  , scalar_comp)
                                + c3 * scalar(i, j, k+1, scalar_comp) );
            } else if (ext_dir_on_zhi) {
                // Third order stencil with variable dz
                Real zm   = Compute_Z_AtWFace(i,j,k-1,z_nd);
                Real dz0  = Compute_Z_AtWFace(i,j,k,z_nd) - zm;
                Real dz1  = zm - Compute_Z_AtWFace(i,j,k-2,z_nd);
                Real idz0 = one / dz0;
                Real f    = (dz1 / dz0) + two;
                Real f2   = f*f;
                Real c3   = two / (f - f2);
                Real c2   = -f2*c3;
                Real c1   = -(one-f2)*c3;
                GradCz = idz0 * (  -( c1 * scalar(i, j, k  , scalar_comp)
                                    + c2 * scalar(i, j, k-1, scalar_comp)
                                    + c3 * scalar(i, j, k-2, scalar_comp) ) );
            } else {
                Real met_h_zeta = Compute_h_zeta_AtKface(i,j,k,cellSizeInv,z_nd);
                GradCz = (dz_inv/met_h_zeta) * ( scalar(i, j, k, scalar_comp) - scalar(i, j, k-1, scalar_comp) );
            }

            if (SurfLayer_on_zlo || SurfLayer_on_zhi) {
                if (native_policy.is_theta) {
                    zflux(i,j,k,flux_comp) = hfx_z(i,j,k);
                } else if (native_policy.is_q1) {
                    zflux(i,j,k,flux_comp) = qfx1_z(i,j,k);
                } else {
                    zflux(i,j,k,flux_comp) = zero;
                }
            } else {
                zflux(i,j,k,flux_comp) = -rhoAlpha * GradCz;
            }

            if (native_policy.is_theta) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    hfx_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q1) {
                if (!(SurfLayer_on_zlo || SurfLayer_on_zhi)) {
                    qfx1_z(i,j,k) = zflux(i,j,k,flux_comp);
                }
            } else  if (native_policy.is_q2) {
                qfx2_z(i,j,k) = zflux(i,j,k,flux_comp);
            }
        });
    }

    // NOTE: With terrain, we implicitly treat the leading order vertical gradient (no metric terms)
    // This allows us to do semi-implicit discretization of the vertical diffusive terms
    if (native_policy.scale_raw_vertical_flux) {
        ScaleScalarDiffusionVerticalFlux(zbx, zflux, flux_comp, explicit_fac);
    }

    //-----------------------------------------------------------------------------------
    //
    // Modify fluxes by terrain and use fluxes to compute RHS
    //
    // Note that we combine all of these operations in order to keep this section
    //      of the loop tiling-safe.
    //-----------------------------------------------------------------------------------
    ParallelFor(bx,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real xfluxbar_lo, yfluxbar_lo;
        Real met_h_xi_lo  = Compute_h_xi_AtKface (i,j,k  ,cellSizeInv,z_nd);
        Real met_h_eta_lo = Compute_h_eta_AtKface(i,j,k  ,cellSizeInv,z_nd);
        if (k == dom_lo.z) {
            Real xfluxlo  = myhalf * ( xflux(i,j,k  ,flux_comp) + xflux(i+1,j,k  ,flux_comp) );
            Real xfluxhi  = myhalf * ( xflux(i,j,k+1,flux_comp) + xflux(i+1,j,k+1,flux_comp) );
            xfluxbar_lo = Real(1.5)*xfluxlo - myhalf*xfluxhi;

            Real yfluxlo  = myhalf * ( yflux(i,j,k  ,flux_comp) + yflux(i,j+1,k  ,flux_comp) );
            Real yfluxhi  = myhalf * ( yflux(i,j,k+1,flux_comp) + yflux(i,j+1,k+1,flux_comp) );
            yfluxbar_lo = Real(1.5)*yfluxlo - myhalf*yfluxhi;
        } else {
            xfluxbar_lo = fourth * ( xflux(i,j,k  ,flux_comp) + xflux(i+1,j  ,k  ,flux_comp)
                                   + xflux(i,j,k-1,flux_comp) + xflux(i+1,j  ,k-1,flux_comp) );
            yfluxbar_lo = fourth * ( yflux(i,j,k  ,flux_comp) + yflux(i  ,j+1,k  ,flux_comp)
                                   + yflux(i,j,k-1,flux_comp) + yflux(i  ,j+1,k-1,flux_comp) );
        }

        Real xfluxbar_hi, yfluxbar_hi;
        Real met_h_xi_hi  = Compute_h_xi_AtKface (i,j,k+1,cellSizeInv,z_nd);
        Real met_h_eta_hi = Compute_h_eta_AtKface(i,j,k+1,cellSizeInv,z_nd);
        if (k == dom_hi.z) {
            Real xfluxlo  = myhalf * ( xflux(i,j,k-1,flux_comp) + xflux(i+1,j,k-1,flux_comp) );
            Real xfluxhi  = myhalf * ( xflux(i,j,k  ,flux_comp) + xflux(i+1,j,k  ,flux_comp) );
            xfluxbar_hi = Real(1.5)*xfluxhi - myhalf*xfluxlo;

            Real yfluxlo  = myhalf * ( yflux(i,j,k-1,flux_comp) + yflux(i,j+1,k-1,flux_comp) );
            Real yfluxhi  = myhalf * ( yflux(i,j,k  ,flux_comp) + yflux(i,j+1,k  ,flux_comp) );
            yfluxbar_hi = Real(1.5)*yfluxhi - myhalf*yfluxlo;
        } else {
            xfluxbar_hi = fourth * ( xflux(i,j,k+1,flux_comp) + xflux(i+1,j  ,k+1,flux_comp)
                                   + xflux(i,j,k  ,flux_comp) + xflux(i+1,j  ,k  ,flux_comp) );
            yfluxbar_hi = fourth * ( yflux(i,j,k+1,flux_comp) + yflux(i  ,j+1,k+1,flux_comp)
                                   + yflux(i,j,k  ,flux_comp) + yflux(i  ,j+1,k  ,flux_comp) );
        }

        // Allow semi-implicit discretization of the vertical diffusive terms
        Real zflux_lo;
        if ( SurfLayer_zlo &&
             k == dom_lo.z &&
             !native_policy.is_theta &&
             !native_policy.is_q1 ) {
            zflux_lo = zero;
        } else {
            zflux_lo = TerrainDiffusionGz(zflux(i,j,k,flux_comp), mf_mx(i,j,0),
                                           met_h_xi_lo, xfluxbar_lo, mf_my(i,j,0),
                                           met_h_eta_lo, yfluxbar_lo);
        }
        Real zflux_hi = TerrainDiffusionGz(zflux(i,j,k+1,flux_comp), mf_mx(i,j,0),
                                           met_h_xi_hi, xfluxbar_hi, mf_my(i,j,0),
                                           met_h_eta_hi, yfluxbar_hi);

        Real stateContrib = TerrainDiffusionDivergence_T(
            xflux(i+1,j,k,flux_comp), ax(i+1,j,k), mf_uy(i+1,j,0),
            xflux(i,j,k,flux_comp), ax(i,j,k), mf_uy(i,j,0),
            yflux(i,j+1,k,flux_comp), ay(i,j+1,k), mf_vx(i,j+1,0),
            yflux(i,j,k,flux_comp), ay(i,j,k), mf_vx(i,j,0),
            zflux_hi, zflux_lo, dx_inv, dy_inv, dz_inv,
            mf_mx(i,j,0), mf_my(i,j,0), detJ(i,j,k));

        cell_rhs(i,j,k,rhs_comp) -= stateContrib;
    });
    } // n

    const PBLDerivativeDzInv_T pbl_derivative_dz_inv{z_cc};
#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
