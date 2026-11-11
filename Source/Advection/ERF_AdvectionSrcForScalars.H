#ifndef ERF_ADVECTION_SRC_FOR_SCALARS_H_
#define ERF_ADVECTION_SRC_FOR_SCALARS_H_

#include <AMReX.H>

#include <ERF_IndexDefines.H>
#include <ERF_Interpolation.H>

/**
 * Build scalar face fluxes with a selected pair of interpolation policies.
 *
 * @tparam InterpType_H horizontal interpolation type
 * @tparam InterpType_V vertical interpolation type
 * @tparam FaceValue write face values instead of fluxes
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] upw_x x-face field whose sign picks the upwind side (the mass flux, or a direction)
 * @param[in] upw_y y-face field whose sign picks the upwind side
 * @param[in] upw_z z-face field whose sign picks the upwind side
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
template<typename InterpType_H, typename InterpType_V, bool FaceValue = false>
inline void
BuildScalarAdvectionFluxesWithPolicy (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom,
    const amrex::Array4<const amrex::Real>& upw_x,
    const amrex::Array4<const amrex::Real>& upw_y,
    const amrex::Array4<const amrex::Real>& upw_z,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac)
{
    InterpType_H interp_h(scalar, horiz_upw_frac);
    InterpType_V interp_v(scalar, vert_upw_frac);

    const amrex::Box xbx = amrex::surroundingNodes(bx,0);
    const amrex::Box ybx = amrex::surroundingNodes(bx,1);
    const amrex::Box zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real interpx(zero);
        interp_h.InterpolateInX(i,j,k,scalar_comp,interpx,upw_x(i,j,k));
        flx_arr[0](i,j,k,flux_comp) = FaceValue ? interpx : avg_xmom(i,j,k) * interpx;
    });
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real interpy(zero);
        interp_h.InterpolateInY(i,j,k,scalar_comp,interpy,upw_y(i,j,k));
        flx_arr[1](i,j,k,flux_comp) = FaceValue ? interpy : avg_ymom(i,j,k) * interpy;
    });
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real interpz(zero);
        interp_v.InterpolateInZ(i,j,k,scalar_comp,interpz,upw_z(i,j,k));
        flx_arr[2](i,j,k,flux_comp) = FaceValue ? interpz : avg_zmom(i,j,k) * interpz;
    });
}

/**
 * Build centered-second-order scalar face fluxes without a reconstruction temporary.
 *
 * @tparam FaceValue write face values instead of fluxes
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 */
template<bool FaceValue = false>
inline void
BuildScalarAdvectionFluxesCentered2 (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom)
{
    const amrex::Box xbx = amrex::surroundingNodes(bx,0);
    const amrex::Box ybx = amrex::surroundingNodes(bx,1);
    const amrex::Box zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real scalar_on_face = myhalf *
            (scalar(i,j,k,scalar_comp) + scalar(i-1,j,k,scalar_comp));
        flx_arr[0](i,j,k,flux_comp) = FaceValue ? scalar_on_face : avg_xmom(i,j,k) * scalar_on_face;
    });
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real scalar_on_face = myhalf *
            (scalar(i,j,k,scalar_comp) + scalar(i,j-1,k,scalar_comp));
        flx_arr[1](i,j,k,flux_comp) = FaceValue ? scalar_on_face : avg_ymom(i,j,k) * scalar_on_face;
    });
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real scalar_on_face = myhalf *
            (scalar(i,j,k,scalar_comp) + scalar(i,j,k-1,scalar_comp));
        flx_arr[2](i,j,k,flux_comp) = FaceValue ? scalar_on_face : avg_zmom(i,j,k) * scalar_on_face;
    });
}

/**
 * Dispatch scalar face-flux construction for a selected vertical scheme.
 *
 * @tparam InterpType_H horizontal interpolation type
 * @tparam FaceValue write face values instead of fluxes
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] upw_x x-face field whose sign picks the upwind side (the mass flux, or a direction)
 * @param[in] upw_y y-face field whose sign picks the upwind side
 * @param[in] upw_z z-face field whose sign picks the upwind side
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 * @param[in] vert_adv_type vertical advection scheme
 */
template<typename InterpType_H, bool FaceValue = false>
inline void
BuildScalarAdvectionFluxesVert (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom,
    const amrex::Array4<const amrex::Real>& upw_x,
    const amrex::Array4<const amrex::Real>& upw_y,
    const amrex::Array4<const amrex::Real>& upw_z,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac,
    const AdvType vert_adv_type)
{
    switch (vert_adv_type) {
    case AdvType::Centered_2nd:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED2,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND3,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd_SL:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND3SL,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_4th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED4,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_5th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND5,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_6th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED6,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    default:
        // Not an assert: under NDEBUG that would compile away and leave the
        // face fluxes unwritten, so the divergence below would silently
        // consume whatever stale data is in this flux component.  Note that
        // the WENO schemes have no vertical-only variant -- they must be
        // selected for both the horizontal and the vertical direction.
        amrex::Abort("Unknown or unsupported vertical advection scheme!");
    }
}

/**
 * Dispatch shared by BuildScalarAdvectionFluxes (FaceValue = false: write the mass flux times
 * the reconstructed face value) and BuildScalarAdvectionFaceValues (FaceValue = true: write the
 * reconstructed face value itself).
 *
 * @tparam FaceValue write face values instead of fluxes
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux (or face-value) views
 * @param[in] flux_comp component of the supplied face views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] upw_x x-face field whose sign picks the upwind side (the mass flux, or a direction)
 * @param[in] upw_y y-face field whose sign picks the upwind side
 * @param[in] upw_z z-face field whose sign picks the upwind side
 * @param[in] horiz_adv_type horizontal advection scheme
 * @param[in] vert_adv_type vertical advection scheme
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
template<bool FaceValue>
inline void
BuildScalarAdvectionFluxesImpl (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom,
    const amrex::Array4<const amrex::Real>& upw_x,
    const amrex::Array4<const amrex::Real>& upw_y,
    const amrex::Array4<const amrex::Real>& upw_z,
    const AdvType horiz_adv_type,
    const AdvType vert_adv_type,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac)
{
    if (horiz_adv_type == AdvType::Centered_2nd &&
        vert_adv_type == AdvType::Centered_2nd) {
        BuildScalarAdvectionFluxesCentered2<FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom);
        return;
    }

    switch (horiz_adv_type) {
    case AdvType::Centered_2nd:
        BuildScalarAdvectionFluxesVert<CENTERED2,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_3rd:
        BuildScalarAdvectionFluxesVert<UPWIND3,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_3rd_SL:
        BuildScalarAdvectionFluxesVert<UPWIND3SL,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Centered_4th:
        BuildScalarAdvectionFluxesVert<CENTERED4,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_5th:
        BuildScalarAdvectionFluxesVert<UPWIND5,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Centered_6th:
        BuildScalarAdvectionFluxesVert<CENTERED6,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Weno_3:
        BuildScalarAdvectionFluxesWithPolicy<WENO3,WENO3,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_5:
        BuildScalarAdvectionFluxesWithPolicy<WENO5,WENO5,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_7:
        BuildScalarAdvectionFluxesWithPolicy<WENO7,WENO7,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_3Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z3,WENO_Z3,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_3MZQ:
        BuildScalarAdvectionFluxesWithPolicy<WENO_MZQ3,WENO_MZQ3,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_5Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z5,WENO_Z5,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_7Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z7,WENO_Z7,FaceValue>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, upw_x, upw_y, upw_z, horiz_upw_frac, vert_upw_frac);
        break;
    default:
        // See the note in BuildScalarAdvectionFluxesVert: an assert here would
        // vanish under NDEBUG and leave stale fluxes in place of an error.
        amrex::Abort("Unknown horizontal advection scheme!");
    }
}

/**
 * Build scalar face fluxes using ERF's compile-time interpolation policies.
 * WENO schemes intentionally use the same policy horizontally and vertically,
 * matching the existing native scalar-advection dispatch.
 *
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] horiz_adv_type horizontal advection scheme
 * @param[in] vert_adv_type vertical advection scheme
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
inline void
BuildScalarAdvectionFluxes (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom,
    const AdvType horiz_adv_type,
    const AdvType vert_adv_type,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac)
{
    BuildScalarAdvectionFluxesImpl<false>(bx, scalar, scalar_comp, flx_arr, flux_comp,
                                          avg_xmom, avg_ymom, avg_zmom,
                                          avg_xmom, avg_ymom, avg_zmom,
                                          horiz_adv_type, vert_adv_type,
                                          horiz_upw_frac, vert_upw_frac);
}

/**
 * As above, but with the upwind side on each face picked from the sign of upw_x, upw_y and
 * upw_z instead of the sign of the mass flux.  A pseudo-incompressible level passes the upwind
 * direction its projection set (see ERF::pi_upw_dir), so that theta is reconstructed on each face
 * exactly as the projection reconstructed rho0 theta0 / rho.
 *
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] upw_x x-face field whose sign picks the upwind side
 * @param[in] upw_y y-face field whose sign picks the upwind side
 * @param[in] upw_z z-face field whose sign picks the upwind side
 * @param[in] horiz_adv_type horizontal advection scheme
 * @param[in] vert_adv_type vertical advection scheme
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
inline void
BuildScalarAdvectionFluxes (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<const amrex::Real>& avg_xmom,
    const amrex::Array4<const amrex::Real>& avg_ymom,
    const amrex::Array4<const amrex::Real>& avg_zmom,
    const amrex::Array4<const amrex::Real>& upw_x,
    const amrex::Array4<const amrex::Real>& upw_y,
    const amrex::Array4<const amrex::Real>& upw_z,
    const AdvType horiz_adv_type,
    const AdvType vert_adv_type,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac)
{
    BuildScalarAdvectionFluxesImpl<false>(bx, scalar, scalar_comp, flx_arr, flux_comp,
                                          avg_xmom, avg_ymom, avg_zmom,
                                          upw_x, upw_y, upw_z,
                                          horiz_adv_type, vert_adv_type,
                                          horiz_upw_frac, vert_upw_frac);
}

/**
 * The face values of a cell-centered field that BuildScalarAdvectionFluxes would multiply by
 * the mass flux: the same reconstruction, chosen by the same scheme and upwind fraction, with
 * the upwind side picked from the sign of the given mass flux (every scheme uses only its sign,
 * and a zero flux selects the centered blend, as in advection).  Writes the face values instead
 * of the fluxes.
 *
 * @param[in] bx cell-centered box over whose faces the values are computed
 * @param[in] scalar supplied cell-centered field view
 * @param[in] scalar_comp component of the supplied field view to reconstruct
 * @param[out] face_arr face-value views
 * @param[in] face_comp component of the supplied face-value views to write
 * @param[in] xflux x-face mass flux (its sign picks the upwind side)
 * @param[in] yflux y-face mass flux (its sign picks the upwind side)
 * @param[in] zflux z-face mass flux (its sign picks the upwind side)
 * @param[in] horiz_adv_type horizontal advection scheme
 * @param[in] vert_adv_type vertical advection scheme
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
inline void
BuildScalarAdvectionFaceValues (
    const amrex::Box& bx,
    const amrex::Array4<const amrex::Real>& scalar,
    const int scalar_comp,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& face_arr,
    const int face_comp,
    const amrex::Array4<const amrex::Real>& xflux,
    const amrex::Array4<const amrex::Real>& yflux,
    const amrex::Array4<const amrex::Real>& zflux,
    const AdvType horiz_adv_type,
    const AdvType vert_adv_type,
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac)
{
    BuildScalarAdvectionFluxesImpl<true>(bx, scalar, scalar_comp, face_arr, face_comp,
                                         xflux, yflux, zflux,
                                         xflux, yflux, zflux,
                                         horiz_adv_type, vert_adv_type,
                                         horiz_upw_frac, vert_upw_frac);

}

/**
 * Apply the mapped divergence of component-selected scalar face fluxes.
 *
 * @param[in] bx cell-centered box over which the tendency is computed
 * @param[in] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to read
 * @param[out] rhs supplied tendency view
 * @param[in] rhs_comp component of the supplied tendency view to write
 * @param[in] detJ metric Jacobian
 * @param[in] cellSizeInv inverse grid spacing
 * @param[in] mf_mx x map factor at cell centers
 * @param[in] mf_my y map factor at cell centers
 */
inline void
ApplyScalarAdvectionFluxDivergence (
    const amrex::Box& bx,
    const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
    const int flux_comp,
    const amrex::Array4<amrex::Real>& rhs,
    const int rhs_comp,
    const amrex::Array4<const amrex::Real>& detJ,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
    const amrex::Array4<const amrex::Real>& mf_mx,
    const amrex::Array4<const amrex::Real>& mf_my)
{
    const auto dxInv = cellSizeInv[0];
    const auto dyInv = cellSizeInv[1];
    const auto dzInv = cellSizeInv[2];

    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (detJ(i,j,k) > zero) {
            const amrex::Real invdetJ = one / detJ(i,j,k);
            const amrex::Real mfsq = mf_mx(i,j,0) * mf_my(i,j,0);
            rhs(i,j,k,rhs_comp) = - invdetJ * mfsq * (
                (flx_arr[0](i+1,j,k,flux_comp) - flx_arr[0](i,j,k,flux_comp)) * dxInv +
                (flx_arr[1](i,j+1,k,flux_comp) - flx_arr[1](i,j,k,flux_comp)) * dyInv +
                (flx_arr[2](i,j,k+1,flux_comp) - flx_arr[2](i,j,k,flux_comp)) * dzInv);
        } else {
            rhs(i,j,k,rhs_comp) = zero;
        }
    });
}

#endif // ERF_ADVECTION_SRC_FOR_SCALARS_H_
