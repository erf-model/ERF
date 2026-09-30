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
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 */
template<typename InterpType_H, typename InterpType_V>
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
        interp_h.InterpolateInX(i,j,k,scalar_comp,interpx,avg_xmom(i,j,k));
        flx_arr[0](i,j,k,flux_comp) = avg_xmom(i,j,k) * interpx;
    });
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real interpy(zero);
        interp_h.InterpolateInY(i,j,k,scalar_comp,interpy,avg_ymom(i,j,k));
        flx_arr[1](i,j,k,flux_comp) = avg_ymom(i,j,k) * interpy;
    });
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real interpz(zero);
        interp_v.InterpolateInZ(i,j,k,scalar_comp,interpz,avg_zmom(i,j,k));
        flx_arr[2](i,j,k,flux_comp) = avg_zmom(i,j,k) * interpz;
    });
}

/**
 * Build centered-second-order scalar face fluxes without a reconstruction temporary.
 *
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied cell-centered scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 */
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
        flx_arr[0](i,j,k,flux_comp) = avg_xmom(i,j,k) * scalar_on_face;
    });
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real scalar_on_face = myhalf *
            (scalar(i,j,k,scalar_comp) + scalar(i,j-1,k,scalar_comp));
        flx_arr[1](i,j,k,flux_comp) = avg_ymom(i,j,k) * scalar_on_face;
    });
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const amrex::Real scalar_on_face = myhalf *
            (scalar(i,j,k,scalar_comp) + scalar(i,j,k-1,scalar_comp));
        flx_arr[2](i,j,k,flux_comp) = avg_zmom(i,j,k) * scalar_on_face;
    });
}

/**
 * Dispatch scalar face-flux construction for a selected vertical scheme.
 *
 * @tparam InterpType_H horizontal interpolation type
 * @param[in] bx cell-centered box over which scalar fluxes are computed
 * @param[in] scalar supplied scalar view
 * @param[in] scalar_comp component of the supplied scalar view to reconstruct
 * @param[out] flx_arr face-flux views
 * @param[in] flux_comp component of the supplied face-flux views to write
 * @param[in] avg_xmom x-component of time-averaged mass flux
 * @param[in] avg_ymom y-component of time-averaged mass flux
 * @param[in] avg_zmom z-component of time-averaged mass flux
 * @param[in] horiz_upw_frac horizontal upwind blending fraction
 * @param[in] vert_upw_frac vertical upwind blending fraction
 * @param[in] vert_adv_type vertical advection scheme
 */
template<typename InterpType_H>
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
    const amrex::Real horiz_upw_frac,
    const amrex::Real vert_upw_frac,
    const AdvType vert_adv_type)
{
    switch (vert_adv_type) {
    case AdvType::Centered_2nd:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED2>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND3>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd_SL:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND3SL>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_4th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED4>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_5th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,UPWIND5>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_6th:
        BuildScalarAdvectionFluxesWithPolicy<InterpType_H,CENTERED6>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
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
    if (horiz_adv_type == AdvType::Centered_2nd &&
        vert_adv_type == AdvType::Centered_2nd) {
        BuildScalarAdvectionFluxesCentered2(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom);
        return;
    }

    switch (horiz_adv_type) {
    case AdvType::Centered_2nd:
        BuildScalarAdvectionFluxesVert<CENTERED2>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_3rd:
        BuildScalarAdvectionFluxesVert<UPWIND3>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_3rd_SL:
        BuildScalarAdvectionFluxesVert<UPWIND3SL>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Centered_4th:
        BuildScalarAdvectionFluxesVert<CENTERED4>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Upwind_5th:
        BuildScalarAdvectionFluxesVert<UPWIND5>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Centered_6th:
        BuildScalarAdvectionFluxesVert<CENTERED6>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac, vert_adv_type);
        break;
    case AdvType::Weno_3:
        BuildScalarAdvectionFluxesWithPolicy<WENO3,WENO3>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_5:
        BuildScalarAdvectionFluxesWithPolicy<WENO5,WENO5>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_7:
        BuildScalarAdvectionFluxesWithPolicy<WENO7,WENO7>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_3Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z3,WENO_Z3>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_3MZQ:
        BuildScalarAdvectionFluxesWithPolicy<WENO_MZQ3,WENO_MZQ3>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_5Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z5,WENO_Z5>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Weno_7Z:
        BuildScalarAdvectionFluxesWithPolicy<WENO_Z7,WENO_Z7>(
            bx, scalar, scalar_comp, flx_arr, flux_comp,
            avg_xmom, avg_ymom, avg_zmom, horiz_upw_frac, vert_upw_frac);
        break;
    default:
        // See the note in BuildScalarAdvectionFluxesVert: an assert here would
        // vanish under NDEBUG and leave stale fluxes in place of an error.
        amrex::Abort("Unknown horizontal advection scheme!");
    }
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
