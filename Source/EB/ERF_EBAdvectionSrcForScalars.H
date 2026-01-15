#include <ERF_IndexDefines.H>
#include <ERF_Interpolation.H>

#include <iostream>
#include <fstream>
#include <AMReX_Vector.H>

/**
 * When terrain_type == TerrainType::EB,
 * Wrapper function for computing the advective tendency w/ spatial order > 2
 */
template<typename InterpType_H, typename InterpType_V>
void
EBAdvectionSrcForScalarsWrapper (const amrex::Box& bx,
                                const int& ncomp, const int& icomp,
                                const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                                const amrex::Array4<const amrex::Real>& cell_prim,
                                const amrex::Array4<const amrex::Real>& avg_xmom,
                                const amrex::Array4<const amrex::Real>& avg_ymom,
                                const amrex::Array4<const amrex::Real>& avg_zmom,
                                const amrex::Array4<const amrex::EBCellFlag>& cellflag,
                                const amrex::Array4<const amrex::Real>& ax_arr,
                                const amrex::Array4<const amrex::Real>& ay_arr,
                                const amrex::Array4<const amrex::Real>& az_arr,
                                const amrex::Real horiz_upw_frac,
                                const amrex::Real vert_upw_frac)
{
    // Instantiate structs for vert/horiz interp
    InterpType_H interp_prim_h(cell_prim, horiz_upw_frac);
    InterpType_V interp_prim_v(cell_prim, vert_upw_frac);

    const int ncells_h = interp_prim_h.GetUpwindCellNumber();
    const int ncells_v = interp_prim_v.GetUpwindCellNumber();

    // Instantiate structs for lower-order vert/horiz interp
    CENTERED2 interp_prim_CEN2(cell_prim, 0);
    UPWIND3   interp_prim_UPW3(cell_prim, horiz_upw_frac);

    const amrex::Box xbx = amrex::surroundingNodes(bx,0).grow(amrex::IntVect(0, 1, 1));
    const amrex::Box ybx = amrex::surroundingNodes(bx,1).grow(amrex::IntVect(1, 0, 1));
    const amrex::Box zbx = amrex::surroundingNodes(bx,2).grow(amrex::IntVect(1, 1, 0));

    amrex::ParallelFor(xbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;

        if ( ax_arr(i,j,k) > 0. ) // Do we need a threshold value for the area fraction to avoid infinitesimal fractions?
        {
            // Find the highest upwind order based on the cell flag

            int icell = 0;

            if (avg_xmom(i,j,k)>0.)
            {
                for (int ii=0; ii<ncells_h; ++ii) {
                    if (cellflag(i-ii-1,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }
            else if (avg_xmom(i,j,k)<0.)
            {
                for (int ii=0; ii<ncells_h; ++ii) {
                    if (cellflag(i+ii,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }

            // Interpolate the scalar variable (cell_prim) using the highest-order scheme

            const int prim_index = cons_index - 1;
            amrex::Real interpx(0.);

            if (icell==ncells_h) {
                interp_prim_h.InterpolateInX(i,j,k,prim_index,interpx,avg_xmom(i,j,k));
            } else {
                if (icell==1) {
                    interp_prim_CEN2.InterpolateInX(i,j,k,prim_index,interpx,avg_xmom(i,j,k));
                } else if (icell==2) {
                    interp_prim_UPW3.InterpolateInX(i,j,k,prim_index,interpx,avg_xmom(i,j,k));
                }
            }
            (flx_arr[0])(i,j,k,cons_index) = avg_xmom(i,j,k) * interpx;
        }
        else
        {
            (flx_arr[0])(i,j,k,cons_index) = 0.;
        }
    });

    amrex::ParallelFor(ybx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;

        if ( ay_arr(i,j,k) > 0. )
        {
            // Find the highest upwind order based on the cell flag

            int jcell = 0;

            if (avg_ymom(i,j,k)>0.)
            {
                for (int jj=0; jj<ncells_h; ++jj) {
                    if (cellflag(i,j-jj-1,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }
            else if (avg_ymom(i,j,k)<0.)
            {
                for (int jj=0; jj<ncells_h; ++jj) {
                    if (cellflag(i,j+jj,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }

            // Interpolate the scalar variable (cell_prim) using the highest-order scheme

            const int prim_index = cons_index - 1;
            amrex::Real interpy(0.);

            if (jcell==ncells_h) {
                interp_prim_h.InterpolateInY(i,j,k,prim_index,interpy,avg_ymom(i,j,k));
            } else {
                if (jcell==1) {
                    interp_prim_CEN2.InterpolateInY(i,j,k,prim_index,interpy,avg_ymom(i,j,k));
                } else if (jcell==2) {
                    interp_prim_UPW3.InterpolateInY(i,j,k,prim_index,interpy,avg_ymom(i,j,k));
                }
            }
            (flx_arr[1])(i,j,k,cons_index) = avg_ymom(i,j,k) * interpy;
        }
        else
        {
            (flx_arr[1])(i,j,k,cons_index) = 0.;
        }
    });

    interp_prim_UPW3.SetUpwinding(vert_upw_frac);

    amrex::ParallelFor(zbx, ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        const int cons_index = icomp + n;

        if ( az_arr(i,j,k) > 0. )
        {
            // Find the highest upwind order based on the cell flag

            int kcell = 0;

            if (avg_zmom(i,j,k)>0.)
            {
                for (int kk=0; kk<ncells_v; ++kk) {
                    if (cellflag(i,j,k-kk-1).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }
            else if (avg_zmom(i,j,k)<0.)
            {
                for (int kk=0; kk<ncells_v; ++kk) {
                    if (cellflag(i,j,k+kk).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }

            // Interpolate the scalar variable (cell_prim) using the highest-order scheme

            const int prim_index = cons_index - 1;
            amrex::Real interpz(0.);

            if (kcell==ncells_v) {
                interp_prim_v.InterpolateInZ(i,j,k,prim_index,interpz,avg_zmom(i,j,k));
            } else {
                if (kcell==1) {
                    interp_prim_CEN2.InterpolateInZ(i,j,k,prim_index,interpz,avg_zmom(i,j,k));
                } else if (kcell==2) {
                    interp_prim_UPW3.InterpolateInZ(i,j,k,prim_index,interpz,avg_zmom(i,j,k));
                }
            }
            (flx_arr[2])(i,j,k,cons_index) = avg_zmom(i,j,k) * interpz;
        }
        else
        {
            (flx_arr[2])(i,j,k,cons_index) = 0.;
        }
    });
}


/**
 * When terrain_type == TerrainType::EB,
 * Wrapper function for templating the vertical advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
EBAdvectionSrcForScalarsVert (const amrex::Box& bx,
                            const int& ncomp, const int& icomp,
                            const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx_arr,
                            const amrex::Array4<const amrex::Real>& cell_prim,
                            const amrex::Array4<const amrex::Real>& avg_xmom,
                            const amrex::Array4<const amrex::Real>& avg_ymom,
                            const amrex::Array4<const amrex::Real>& avg_zmom,
                            const amrex::Array4<const amrex::EBCellFlag>& cellflag,
                            const amrex::Array4<const amrex::Real>& ax_arr,
                            const amrex::Array4<const amrex::Real>& ay_arr,
                            const amrex::Array4<const amrex::Real>& az_arr,
                            const amrex::Real horiz_upw_frac,
                            const amrex::Real vert_upw_frac,
                            const AdvType vert_adv_type)
{
    switch(vert_adv_type) {
    case AdvType::Centered_2nd:
        EBAdvectionSrcForScalarsWrapper<InterpType_H,CENTERED2>(bx, ncomp, icomp,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            cellflag, ax_arr, ay_arr, az_arr,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_3rd:
        EBAdvectionSrcForScalarsWrapper<InterpType_H,UPWIND3>(bx, ncomp, icomp,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            cellflag, ax_arr, ay_arr, az_arr,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_4th:
        EBAdvectionSrcForScalarsWrapper<InterpType_H,CENTERED4>(bx, ncomp, icomp,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            cellflag, ax_arr, ay_arr, az_arr,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Upwind_5th:
        EBAdvectionSrcForScalarsWrapper<InterpType_H,UPWIND5>(bx, ncomp, icomp,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            cellflag, ax_arr, ay_arr, az_arr,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    case AdvType::Centered_6th:
        EBAdvectionSrcForScalarsWrapper<InterpType_H,CENTERED6>(bx, ncomp, icomp,
                                                            flx_arr, cell_prim,
                                                            avg_xmom, avg_ymom, avg_zmom,
                                                            cellflag, ax_arr, ay_arr, az_arr,
                                                            horiz_upw_frac, vert_upw_frac);
        break;
    default:
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown vertical advection scheme!");
    }
}

