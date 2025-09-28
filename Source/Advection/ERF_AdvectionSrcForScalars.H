#include <ERF_IndexDefines.H>
#include <ERF_Interpolation.H>

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H, typename InterpType_V>
void
AdvectionSrcForScalarsWrapper (const amrex::Box& bx,
                               const int& cons_index,
                               const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                               const amrex::Array4<const amrex::Real>& cell_prim,
                               const amrex::Array4<const amrex::Real>& avg_xmom,
                               const amrex::Array4<const amrex::Real>& avg_ymom,
                               const amrex::Array4<const amrex::Real>& avg_zmom,
                               const amrex::Real horiz_upw_frac,
                               const amrex::Real vert_upw_frac)
{
    // Instantiate structs for vert/horiz interp
    InterpType_H interp_prim_h(cell_prim, horiz_upw_frac);
    InterpType_V interp_prim_v(cell_prim, vert_upw_frac);

    const amrex::Box xbx = amrex::surroundingNodes(bx,0);
    const amrex::Box ybx = amrex::surroundingNodes(bx,1);
    const amrex::Box zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int prim_index = cons_index - 1;
        amrex::Real interpx(0.);
        interp_prim_h.InterpolateInX(i,j,k,prim_index,interpx,avg_xmom(i  ,j  ,k  ));
        (flx_arr[0])(i,j,k) = avg_xmom(i,j,k) * interpx;
    });
    amrex::ParallelFor(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int prim_index = cons_index - 1;
        amrex::Real interpy(0.);
        interp_prim_h.InterpolateInY(i,j,k,prim_index,interpy,avg_ymom(i  ,j  ,k  ));
        (flx_arr[1])(i,j,k) = avg_ymom(i,j,k) * interpy;
    });
    amrex::ParallelFor(zbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int prim_index = cons_index - 1;
        amrex::Real interpz(0.);
        interp_prim_v.InterpolateInZ(i,j,k,prim_index,interpz,avg_zmom(i  ,j  ,k  ));
        (flx_arr[2])(i,j,k) = avg_zmom(i,j,k) * interpz;
    });
}

/**
 * Wrapper function for templating the vertical advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
AdvectionSrcForScalarsVert (const amrex::Box& bx,
                            const int& cons_index,
                            const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                            const amrex::Array4<const amrex::Real>& cell_prim,
                            const amrex::Array4<const amrex::Real>& avg_xmom,
                            const amrex::Array4<const amrex::Real>& avg_ymom,
                            const amrex::Array4<const amrex::Real>& avg_zmom,
                            const amrex::Real horiz_upw_frac,
                            const amrex::Real vert_upw_frac,
                            const AdvType vert_adv_type)
{
    switch(vert_adv_type) {
    case AdvType::Centered_2nd:
        AdvectionSrcForScalarsWrapper<InterpType_H,CENTERED2>(bx, cons_index,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom,
                                                              horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd:
        AdvectionSrcForScalarsWrapper<InterpType_H,UPWIND3>(bx, cons_index,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd_SL:
        AdvectionSrcForScalarsWrapper<InterpType_H,UPWIND3SL>(bx, cons_index,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom,
                                                              horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_4th:
        AdvectionSrcForScalarsWrapper<InterpType_H,CENTERED4>(bx, cons_index,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom,
                                                              horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_5th:
        AdvectionSrcForScalarsWrapper<InterpType_H,UPWIND5>(bx, cons_index,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_6th:
        AdvectionSrcForScalarsWrapper<InterpType_H,CENTERED6>(bx, cons_index,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom,
                                                              horiz_upw_frac, vert_upw_frac);
        break;
    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown vertical advection scheme!");
    }
}
