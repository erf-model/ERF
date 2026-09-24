#include "ERF_Diffusion.H"
#include "ERF_ScalarDiffusion.H"
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
        const int scalar_comp = native_policy.scalar_comp;
        const int rhs_comp = native_policy.rhs_comp;
        constexpr int flux_comp = 0;
        const int rho_comp = Rho_comp;
        const Array4<const Real>& scalar = cell_prim;

    // Constant alpha & Turb model
    if (l_consA && l_turb) {
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,native_policy.coefficients.eddy_h_comp);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<0>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dx_inv,mf_ux(i,j,0));
            }

        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,native_policy.coefficients.eddy_h_comp);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<1>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dy_inv,mf_vy(i,j,0));
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

            GradCz = StretchedScalarGradient(scalar, scalar_comp, i,j,k,dz_ptr,klo,khi,
                                              ext_dir_on_zlo,ext_dir_on_zhi);

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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,native_policy.coefficients.eddy_h_comp);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<0>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dx_inv,mf_ux(i,j,0));
            }

        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,true>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,native_policy.coefficients.eddy_h_comp);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<1>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dy_inv,mf_vy(i,j,0));
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

            GradCz = StretchedScalarGradient(scalar, scalar_comp, i,j,k,dz_ptr,klo,khi,
                                              ext_dir_on_zlo,ext_dir_on_zhi);

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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,1,0,0,0);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else {
                xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<0>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dx_inv,mf_ux(i,j,0));
            }

        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<true,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients,
                i,j,k,0,1,0,0);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<1>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dy_inv,mf_vy(i,j,0));
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

            GradCz = StretchedScalarGradient(scalar, scalar_comp, i,j,k,dz_ptr,klo,khi,
                                              ext_dir_on_zlo,ext_dir_on_zhi);

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
        ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients, i,j,k,0,0,0,0);
bool SurfLayer_on_xlo = ( SurfLayer_xlo && i == dom_lo.x);
            bool SurfLayer_on_xhi = ( SurfLayer_xhi && i == dom_hi.x + 1);

            if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_theta)) {
                xflux(i,j,k,flux_comp) = hfx_x(i,j,k);
            } else if ((SurfLayer_on_xlo || SurfLayer_on_xhi) && (native_policy.is_q1)) {
                xflux(i,j,k,flux_comp) = qfx1_x(i,j,k);
            } else {
              xflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<0>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dx_inv,mf_ux(i,j,0));
            }

        });
        ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real rhoAlpha = ScalarDiffusionFaceCoefficient<false,false>(
                cell_data, rho_comp, mu_turb, native_policy.coefficients, i,j,k,0,0,0,0);
bool SurfLayer_on_ylo = ( SurfLayer_ylo && j == dom_lo.y);
            bool SurfLayer_on_yhi = ( SurfLayer_yhi && j == dom_hi.y + 1);

            if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_theta)) {
                yflux(i,j,k,flux_comp) = hfx_y(i,j,k);
            } else if ((SurfLayer_on_ylo || SurfLayer_on_yhi) && (native_policy.is_q1)) {
                yflux(i,j,k,flux_comp) = qfx1_y(i,j,k);
            } else {
                yflux(i,j,k,flux_comp) = ScalarDiffusionFlux_S_Horizontal<1>(
                    scalar, scalar_comp, i,j,k,rhoAlpha,dy_inv,mf_vy(i,j,0));
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

            GradCz = StretchedScalarGradient(scalar, scalar_comp, i,j,k,dz_ptr,klo,khi,
                                              ext_dir_on_zlo,ext_dir_on_zhi);

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

    // Adjust with map factors
    ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        xflux(i,j,k,flux_comp) /= mf_uy(i,j,0);
    });
    ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        yflux(i,j,k,flux_comp) /= mf_vx(i,j,0);
    });

    // This allows us to do semi-implicit discretization of the vertical diffusive terms
    if (native_policy.scale_raw_vertical_flux) {
        ScaleScalarDiffusionVerticalFlux(zbx, zflux, flux_comp, explicit_fac);
    }

    // Apply the component-explicit conservative divergence.
    ApplyScalarDiffusionFluxDivergence_S(bx, xflux, yflux, zflux, flux_comp,
                                         cell_rhs, rhs_comp, mf_mx, mf_my,
                                         dx_inv, dy_inv, dz_ptr);

    } // n

    const PBLDerivativeDzInv_S pbl_derivative_dz_inv{dz_ptr, klo, khi};
#include "ERF_AddTKESources.H"
#include "ERF_AddQKESources.H"
}
