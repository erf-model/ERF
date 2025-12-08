#include <ERF_IndexDefines.H>
#include <ERF_Interpolation.H>

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H, typename InterpType_V, typename WallInterpType>
void
EBAdvectionSrcForMomWrapper (const amrex::Vector<amrex::Box>& bxx_grown,
                            const amrex::Vector<amrex::Box>& bxy_grown,
                            const amrex::Vector<amrex::Box>& bxz_grown,
                            const amrex::Array4<const amrex::Real>& rho_u,
                            const amrex::Array4<const amrex::Real>& rho_v,
                            const amrex::Array4<const amrex::Real>& rho_w,
                            const amrex::Array4<const amrex::Real>& u,
                            const amrex::Array4<const amrex::Real>& v,
                            const amrex::Array4<const amrex::Real>& w,
                            const amrex::Array4<const amrex::EBCellFlag>& u_cflag,
                            const amrex::Array4<const amrex::Real>& u_afrac_x,
                            const amrex::Array4<const amrex::Real>& u_afrac_y,
                            const amrex::Array4<const amrex::Real>& u_afrac_z,
                            const amrex::Array4<const amrex::EBCellFlag>& v_cflag,
                            const amrex::Array4<const amrex::Real>& v_afrac_x,
                            const amrex::Array4<const amrex::Real>& v_afrac_y,
                            const amrex::Array4<const amrex::Real>& v_afrac_z,
                            const amrex::Array4<const amrex::EBCellFlag>& w_cflag,
                            const amrex::Array4<const amrex::Real>& w_afrac_x,
                            const amrex::Array4<const amrex::Real>& w_afrac_y,
                            const amrex::Array4<const amrex::Real>& w_afrac_z,
                            const amrex::Array4<const amrex::Real>& mf_ux_inv,
                            const amrex::Array4<const amrex::Real>& mf_vx_inv,
                            const amrex::Array4<const amrex::Real>& mf_uy_inv,
                            const amrex::Array4<const amrex::Real>& mf_vy_inv,
                            const amrex::Real upw_frac_h,
                            const amrex::Real upw_frac_v,
                            const AdvType /*vert_adv_type*/,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                                  amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                            const int lo_z_face, const int hi_z_face)
{
    // Instantiate the appropriate structs
    InterpType_H interp_u_horz(u, upw_frac_h); InterpType_V interp_u_vert(u, upw_frac_v); // X-MOM
    InterpType_H interp_v_horz(v, upw_frac_h); InterpType_V interp_v_vert(v, upw_frac_v); // Y-MOM
    InterpType_H interp_w_horz(w, upw_frac_h); InterpType_V interp_w_vert(w, upw_frac_v); // Z-MOM
    // WallInterpType interp_w_wall(w, upw_frac_v); // Z-MOM @ wall

    const int ncells_u_horz = interp_u_horz.GetUpwindCellNumber();
    const int ncells_u_vert = interp_u_vert.GetUpwindCellNumber();
    const int ncells_v_horz = interp_v_horz.GetUpwindCellNumber();
    const int ncells_v_vert = interp_v_vert.GetUpwindCellNumber();
    const int ncells_w_horz = interp_w_horz.GetUpwindCellNumber();
    const int ncells_w_vert = interp_w_vert.GetUpwindCellNumber();

    // Instantiate structs for lower-order vert/horiz interp (commonly using upw_frac_h)
    CENTERED2 interp_u_CEN2(u, 0);
    CENTERED2 interp_v_CEN2(v, 0);
    CENTERED2 interp_w_CEN2(w, 0);

    UPWIND3   interp_u_UPW3(u, upw_frac_h);
    UPWIND3   interp_v_UPW3(v, upw_frac_h);
    UPWIND3   interp_w_UPW3(w, upw_frac_h);

    // ---------------------------------------------------- //
    // x-momentum
    // ---------------------------------------------------- //

    amrex::ParallelFor(bxx_grown[0], bxx_grown[1], bxx_grown[2],
    // x-momentum: x-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( u_afrac_x(i,j,k) > 0. )
        {
            amrex::Real rho_u_avg = 0.5 * (rho_u(i-1,j,k) * mf_ux_inv(i-1,j,0) + rho_u(i,j,k) * mf_ux_inv(i,j,0));

            int icell = 0;

            if (rho_u_avg>0.)
            {
                for (int ii=0; ii<ncells_u_horz; ++ii) {
                    if (u_cflag(i-ii-1,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }
            else if (rho_u_avg<0.)
            {
                for (int ii=0; ii<ncells_u_horz; ++ii) {
                    if (u_cflag(i+ii,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }

            // Interpolate u using the highest-order scheme

            amrex::Real interpx = 0.;
            if (icell==ncells_u_horz) {
                interp_u_horz.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
            } else {
                if (icell<=1) {
                    interp_u_CEN2.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                } else if (icell==2) {
                    interp_u_UPW3.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                }
            }
            flx_u_arr[0](i,j,k) = rho_u_avg * interpx;
        }
        else
        {
            flx_u_arr[0](i,j,k) = 0.;
        }
    },
    // x-momentum: y-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( u_afrac_y(i,j,k) > 0. )
        {
            amrex::Real rho_v_avg = 0.5 * (rho_v(i,j,k) * mf_vx_inv(i,j,0) + rho_v(i-1,j,k) * mf_vx_inv(i-1,j,0));

            int jcell = 0;

            if (rho_v_avg>0.)
            {
                for (int jj=0; jj<ncells_u_horz; ++jj) {
                    if (u_cflag(i,j-jj-1,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }
            else if (rho_v_avg<0.)
            {
                for (int jj=0; jj<ncells_u_horz; ++jj) {
                    if (u_cflag(i,j+jj,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }

            // Interpolate u using the highest-order scheme

            amrex::Real interpy = 0.;
            if (jcell==ncells_u_horz) {
                interp_u_horz.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
            } else {
                if (jcell<=1) {
                    interp_u_CEN2.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                } else if (jcell==2) {
                    interp_u_UPW3.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                }
            }
            flx_u_arr[1](i,j,k) = rho_v_avg * interpy;
        }
        else
        {
            flx_u_arr[1](i,j,k) = 0.;
        }
    },
    // x-momentum: z-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( u_afrac_z(i,j,k) > 0. )
        {
            amrex::Real rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i-1,j,k));

            int kcell = 0;

            if (rho_w_avg>0.)
            {
                for (int kk=0; kk<ncells_u_vert; ++kk) {
                    if (u_cflag(i,j,k-kk-1).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }
            else if (rho_w_avg<0.)
            {
                for (int kk=0; kk<ncells_u_vert; ++kk) {
                    if (u_cflag(i,j,k+kk).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }

            // Interpolate u using the highest-order scheme

            amrex::Real interpz = 0.;
            if (kcell==ncells_u_vert) {
                interp_u_vert.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
            } else {
                if (kcell<=1) {
                    interp_u_CEN2.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                } else if (kcell==2) {
                    interp_u_UPW3.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                }
            }
            flx_u_arr[2](i,j,k) = rho_w_avg * interpz;
        }
        else
        {
            flx_u_arr[2](i,j,k) = 0.;
        }
    });

    // ---------------------------------------------------- //
    // y-momentum
    // ---------------------------------------------------- //

    amrex::ParallelFor(bxy_grown[0], bxy_grown[1], bxy_grown[2],
    // y-momentum: x-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( v_afrac_x(i,j,k) > 0. )
        {
            amrex::Real rho_u_avg = 0.5 * (rho_u(i,j,k) * mf_uy_inv(i,j,0) + rho_u(i,j-1,k) * mf_uy_inv(i,j-1,0));

            int icell = 0;

            if (rho_u_avg>0.)
            {
                for (int ii=0; ii<ncells_v_horz; ++ii) {
                    if (v_cflag(i-ii-1,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }
            else if (rho_u_avg<0.)
            {
                for (int ii=0; ii<ncells_v_horz; ++ii) {
                    if (v_cflag(i+ii,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }

            // Interpolate v using the highest-order scheme

            amrex::Real interpx = 0.;
            if (icell==ncells_v_horz) {
                interp_v_horz.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
            } else {
                if (icell<=1) {
                    interp_v_CEN2.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                } else if (icell==2) {
                    interp_v_UPW3.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                }
            }
            flx_v_arr[0](i,j,k) = rho_u_avg * interpx;
        }
        else
        {
            flx_v_arr[0](i,j,k) = 0.;
        }
    },
    // y-momentum: y-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( v_afrac_y(i,j,k) > 0. )
        {
            amrex::Real rho_v_avg = 0.5 * (rho_v(i,j,k) * mf_vy_inv(i,j,0) + rho_v(i,j-1,k) * mf_vy_inv(i,j-1,0));

            int jcell = 0;

            if (rho_v_avg>0.)
            {
                for (int jj=0; jj<ncells_v_horz; ++jj) {
                    if (v_cflag(i,j-jj-1,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }
            else if (rho_v_avg<0.)
            {
                for (int jj=0; jj<ncells_v_horz; ++jj) {
                    if (v_cflag(i,j+jj,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }

            // Interpolate v using the highest-order scheme

            amrex::Real interpy = 0.;
            if (jcell==ncells_v_horz) {
                interp_v_horz.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
            } else {
                if (jcell<=1) {
                    interp_v_CEN2.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                } else if (jcell==2) {
                    interp_v_UPW3.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                }
            }
            flx_v_arr[1](i,j,k) = rho_v_avg * interpy;
        }
        else
        {
            flx_v_arr[1](i,j,k) = 0.;
        }
    },
    // y-momentum: z-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( v_afrac_z(i,j,k) > 0. )
        {
            amrex::Real rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j-1,k));

            int kcell = 0;

            if (rho_w_avg>0.)
            {
                for (int kk=0; kk<ncells_v_vert; ++kk) {
                    if (v_cflag(i,j,k-kk-1).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }
            else if (rho_w_avg<0.)
            {
                for (int kk=0; kk<ncells_v_vert; ++kk) {
                    if (v_cflag(i,j,k+kk).isCovered()) {
                        break;
                    }
                    kcell++;
                }
            }

            // Interpolate v using the highest-order scheme

            amrex::Real interpz = 0.;
            if (kcell==ncells_v_vert) {
                interp_v_vert.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
            } else {
                if (kcell<=1) {
                    interp_v_CEN2.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                } else if (kcell==2) {
                    interp_v_UPW3.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                }
            }
            flx_v_arr[2](i,j,k) = rho_w_avg * interpz;
        }
        else
        {
            flx_v_arr[2](i,j,k) = 0.;
        }
    });

    // ---------------------------------------------------- //
    // z-momentum
    // ---------------------------------------------------- //

    amrex::ParallelFor(bxz_grown[0], bxz_grown[1], bxz_grown[2],
    // z-momentum: x-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( w_afrac_x(i,j,k) > 0. )
        {
            amrex::Real rho_u_avg = 0.5 * (rho_u(i,j,k) + rho_u(i,j,k-1)) * mf_ux_inv(i,j,0);

            int icell = 0;

            if (rho_u_avg>0.)
            {
                for (int ii=0; ii<ncells_w_horz; ++ii) {
                    if (w_cflag(i-ii-1,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }
            else if (rho_u_avg<0.)
            {
                for (int ii=0; ii<ncells_w_horz; ++ii) {
                    if (w_cflag(i+ii,j,k).isCovered()) {
                        break;
                    }
                    icell++;
                }
            }

            // Interpolate w using the highest-order scheme

            amrex::Real interpx = 0.;
            if (icell==ncells_w_horz) {
                interp_w_horz.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
            } else {
                if (icell<=1) {
                    interp_w_CEN2.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                } else if (icell==2) {
                    interp_w_UPW3.InterpolateInX(i,j,k,0,interpx,rho_u_avg);
                }
            }
            flx_w_arr[0](i,j,k) = rho_u_avg * interpx;
        }
        else
        {
            flx_w_arr[0](i,j,k) = 0.;
        }
    },
    // z-momentum: y-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( w_afrac_y(i,j,k) > 0. )
        {
            amrex::Real rho_v_avg = 0.5 * (rho_v(i,j,k) + rho_v(i,j,k-1)) * mf_vy_inv(i,j,0);

            int jcell = 0;

            if (rho_v_avg>0.)
            {
                for (int jj=0; jj<ncells_w_horz; ++jj) {
                    if (w_cflag(i,j-jj-1,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }
            else if (rho_v_avg<0.)
            {
                for (int jj=0; jj<ncells_w_horz; ++jj) {
                    if (w_cflag(i,j+jj,k).isCovered()) {
                        break;
                    }
                    jcell++;
                }
            }

            // Interpolate w using the highest-order scheme

            amrex::Real interpy = 0.;
            if (jcell==ncells_w_horz) {
                interp_w_horz.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
            } else {
                if (jcell<=1) {
                    interp_w_CEN2.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                } else if (jcell==2) {
                    interp_w_UPW3.InterpolateInY(i,j,k,0,interpy,rho_v_avg);
                }
            }
            flx_w_arr[1](i,j,k) = rho_v_avg * interpy;
        }
        else
        {
            flx_w_arr[1](i,j,k) = 0.;
        }
    },
    // z-momentum: z-face
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if ( w_afrac_z(i,j,k) > 0. )
        {
            if (k==lo_z_face || k==hi_z_face+1) {

                flx_w_arr[2](i,j,k) = rho_w(i,j,k) * w(i,j,k);

            } else {
                amrex::Real rho_w_avg = 0.5 * (rho_w(i,j,k) + rho_w(i,j,k-1));

                int kcell = 0;

                if (rho_w_avg>0.)
                {
                    for (int kk=0; kk<ncells_w_vert; ++kk) {
                        if (w_cflag(i,j,k-kk-1).isCovered()) {
                            break;
                        }
                        kcell++;
                    }
                }
                else if (rho_w_avg<0.)
                {
                    for (int kk=0; kk<ncells_w_vert; ++kk) {
                        if (w_cflag(i,j,k+kk).isCovered()) {
                            break;
                        }
                        kcell++;
                    }
                }

                // Adjust available cell numbers near the top and bottom boundaries

                if (k == hi_z_face || k == lo_z_face+1) {
                    kcell=std::min(kcell,1);
                } else if (k == hi_z_face-1 || k == lo_z_face+2 ) {
                    kcell=std::min(kcell,2);
                }

                // Interpolate w using the highest-order scheme

                amrex::Real interpz = 0.;
                if (kcell==ncells_w_vert) {
                    interp_w_vert.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                } else {
                    if (kcell<=1) {
                        interp_w_CEN2.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                    } else if (kcell==2) {
                        interp_w_UPW3.InterpolateInZ(i,j,k,0,interpz,rho_w_avg);
                    }
                }
                flx_w_arr[2](i,j,k) = rho_w_avg * interpz;
            } // lo_z_face, hi_z_face
        }
        else
        {
            flx_w_arr[2](i,j,k) = 0.;
        }
    });

}

/**
 * Wrapper function for computing the advective tendency w/ spatial order > 2.
 */
template<typename InterpType_H>
void
EBAdvectionSrcForMomVert (const amrex::Vector<amrex::Box>& bxx_grown,
                        const amrex::Vector<amrex::Box>& bxy_grown,
                        const amrex::Vector<amrex::Box>& bxz_grown,
                        const amrex::Array4<const amrex::Real>& rho_u,
                        const amrex::Array4<const amrex::Real>& rho_v,
                        const amrex::Array4<const amrex::Real>& rho_w,
                        const amrex::Array4<const amrex::Real>& u,
                        const amrex::Array4<const amrex::Real>& v,
                        const amrex::Array4<const amrex::Real>& w,
                        const amrex::Array4<const amrex::EBCellFlag>& u_cflag,
                        const amrex::Array4<const amrex::Real>& u_afrac_x,
                        const amrex::Array4<const amrex::Real>& u_afrac_y,
                        const amrex::Array4<const amrex::Real>& u_afrac_z,
                        const amrex::Array4<const amrex::EBCellFlag>& v_cflag,
                        const amrex::Array4<const amrex::Real>& v_afrac_x,
                        const amrex::Array4<const amrex::Real>& v_afrac_y,
                        const amrex::Array4<const amrex::Real>& v_afrac_z,
                        const amrex::Array4<const amrex::EBCellFlag>& w_cflag,
                        const amrex::Array4<const amrex::Real>& w_afrac_x,
                        const amrex::Array4<const amrex::Real>& w_afrac_y,
                        const amrex::Array4<const amrex::Real>& w_afrac_z,
                        const amrex::Array4<const amrex::Real>& mf_ux_inv,
                        const amrex::Array4<const amrex::Real>& mf_vx_inv,
                        const amrex::Array4<const amrex::Real>& mf_uy_inv,
                        const amrex::Array4<const amrex::Real>& mf_vy_inv,
                        const amrex::Real upw_frac_h,
                        const amrex::Real upw_frac_v,
                        const AdvType vert_adv_type,
                              amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_u_arr,
                              amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_v_arr,
                              amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_w_arr,
                        const int lo_z_face, const int hi_z_face)
{
    if (vert_adv_type == AdvType::Centered_2nd) {
        EBAdvectionSrcForMomWrapper<InterpType_H,CENTERED2,UPWINDALL>(
                                bxx_grown, bxy_grown, bxz_grown,
                                rho_u, rho_v, rho_w, u, v, w,
                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                mf_ux_inv, mf_vx_inv,
                                mf_uy_inv, mf_vy_inv,
                                upw_frac_h, upw_frac_v, vert_adv_type,
                                flx_u_arr, flx_v_arr, flx_w_arr,
                                lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_3rd) {
        EBAdvectionSrcForMomWrapper<InterpType_H,UPWIND3,UPWINDALL>(
                                bxx_grown, bxy_grown, bxz_grown,
                                rho_u, rho_v, rho_w, u, v, w,
                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                mf_ux_inv, mf_vx_inv,
                                mf_uy_inv, mf_vy_inv,
                                upw_frac_h, upw_frac_v, vert_adv_type,
                                flx_u_arr, flx_v_arr, flx_w_arr,
                                lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_4th) {
        EBAdvectionSrcForMomWrapper<InterpType_H,CENTERED4,UPWINDALL>(
                                bxx_grown, bxy_grown, bxz_grown,
                                rho_u, rho_v, rho_w, u, v, w,
                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                mf_ux_inv, mf_vx_inv,
                                mf_uy_inv, mf_vy_inv,
                                upw_frac_h, upw_frac_v, vert_adv_type,
                                flx_u_arr, flx_v_arr, flx_w_arr,
                                lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Upwind_5th) {
        EBAdvectionSrcForMomWrapper<InterpType_H,UPWIND5,UPWINDALL>(
                                bxx_grown, bxy_grown, bxz_grown,
                                rho_u, rho_v, rho_w, u, v, w,
                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                mf_ux_inv, mf_vx_inv,
                                mf_uy_inv, mf_vy_inv,
                                upw_frac_h, upw_frac_v, vert_adv_type,
                                flx_u_arr, flx_v_arr, flx_w_arr,
                                lo_z_face, hi_z_face);
    } else if (vert_adv_type == AdvType::Centered_6th) {
        EBAdvectionSrcForMomWrapper<InterpType_H,CENTERED6,UPWINDALL>(
                                bxx_grown, bxy_grown, bxz_grown,
                                rho_u, rho_v, rho_w, u, v, w,
                                u_cflag, u_afrac_x, u_afrac_y, u_afrac_z,
                                v_cflag, v_afrac_x, v_afrac_y, v_afrac_z,
                                w_cflag, w_afrac_x, w_afrac_y, w_afrac_z,
                                mf_ux_inv, mf_vx_inv,
                                mf_uy_inv, mf_vy_inv,
                                upw_frac_h, upw_frac_v, vert_adv_type,
                                flx_u_arr, flx_v_arr, flx_w_arr,
                                lo_z_face, hi_z_face);
    } else {
        AMREX_ASSERT_WITH_MESSAGE(false, "Unknown advection scheme!");
    }
}
