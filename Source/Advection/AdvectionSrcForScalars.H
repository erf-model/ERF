#include <IndexDefines.H>
#include <Interpolation.H>

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H, typename InterpType_V>
void
AdvectionSrcForScalarsWrapper_N(const amrex::Box& bx,
                                const int& ncomp, const int& icomp,
                                const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                                const amrex::Array4<const amrex::Real>& cell_prim,
                                const amrex::Array4<const amrex::Real>& avg_xmom,
                                const amrex::Array4<const amrex::Real>& avg_ymom,
                                const amrex::Array4<const amrex::Real>& avg_zmom)
{
    // Instantiate structs for vert/horiz interp
    InterpType_H interp_prim_h(cell_prim);
    InterpType_V interp_prim_v(cell_prim);

    const amrex::Box xbx = amrex::surroundingNodes(bx,0);
    const amrex::Box ybx = amrex::surroundingNodes(bx,1);
    const amrex::Box zbx = amrex::surroundingNodes(bx,2);

    amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;
        const int prim_index = cons_index - 1;

        amrex::Real interpx_hi(0.), interpx_lo(0.);

        interp_prim_h.InterpolateInX(i,j,k,prim_index,interpx_lo,avg_xmom(i  ,j  ,k  ));
        (flx_arr[0])(i,j,k,n) = avg_xmom(i,j,k) * interpx_lo;
    });
    amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;
        const int prim_index = cons_index - 1;

        amrex::Real interpy_hi(0.), interpy_lo(0.);
        interp_prim_h.InterpolateInY(i,j,k,prim_index,interpy_lo,avg_ymom(i  ,j  ,k  ));

        (flx_arr[1])(i,j,k,n) = avg_ymom(i,j,k) * interpy_lo;
    });
    amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;
        const int prim_index = cons_index - 1;

        amrex::Real interpz_hi(0.), interpz_lo(0.);

        interp_prim_v.InterpolateInZ(i,j,k,prim_index,interpz_lo,avg_zmom(i  ,j  ,k  ));

        (flx_arr[2])(i,j,k,n) = avg_zmom(i,j,k) * interpz_lo;
    });
}

/**
 * Wrapper function for templating the vertical advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
AdvectionSrcForScalarsVert_N(const amrex::Box& bx,
                             const int& ncomp, const int& icomp,
                             const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                             const amrex::Array4<const amrex::Real>& cell_prim,
                             const amrex::Array4<const amrex::Real>& avg_xmom,
                             const amrex::Array4<const amrex::Real>& avg_ymom,
                             const amrex::Array4<const amrex::Real>& avg_zmom,
                             const AdvType vert_adv_type)
{
    switch(vert_adv_type) {
    case AdvType::Centered_2nd:
        AdvectionSrcForScalarsWrapper_N<InterpType_H,CENTERED2>(bx, ncomp, icomp,
                                                                flx_arr, cell_prim,
                                                                avg_xmom, avg_ymom, avg_zmom);
    case AdvType::Upwind_3rd:
        AdvectionSrcForScalarsWrapper_N<InterpType_H,UPWIND3>(bx, ncomp, icomp,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom);
    case AdvType::Centered_4th:
        AdvectionSrcForScalarsWrapper_N<InterpType_H,CENTERED4>(bx, ncomp, icomp,
                                                                flx_arr, cell_prim,
                                                                avg_xmom, avg_ymom, avg_zmom);
    case AdvType::Upwind_5th:
        AdvectionSrcForScalarsWrapper_N<InterpType_H,UPWIND5>(bx, ncomp, icomp,
                                                              flx_arr, cell_prim,
                                                              avg_xmom, avg_ymom, avg_zmom);
    case AdvType::Centered_6th:
        AdvectionSrcForScalarsWrapper_N<InterpType_H,CENTERED6>(bx, ncomp, icomp,
                                                                flx_arr, cell_prim,
                                                                avg_xmom, avg_ymom, avg_zmom);
    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown vertical advection scheme!");
    }
}
