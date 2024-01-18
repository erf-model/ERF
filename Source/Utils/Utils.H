#ifndef _UTILS_H_
#define _UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>
#include <ERF_FillPatcher.H>

/*
 * Create the Jacobian for the metric transformation when use_terrain is true
 */
void make_J (const amrex::Geometry& geom,
             amrex::MultiFab& z_phys_nd,
             amrex::MultiFab& detJ_cc);

/*
 * Average z_phys_nd on nodes to cell centers
 */
void make_zcc (const amrex::Geometry& geom,
               amrex::MultiFab& z_phys_nd,
               amrex::MultiFab& z_phys_cc);

/*
 * Convert momentum to velocity by dividing by density averaged onto faces
 */
void MomentumToVelocity (amrex::MultiFab& xvel_out,
                         amrex::MultiFab& yvel_out,
                         amrex::MultiFab& zvel_out,
                         const amrex::MultiFab& cons_in,
                         const amrex::MultiFab& xmom_in,
                         const amrex::MultiFab& ymom_in,
                         const amrex::MultiFab& zmom_in);

/*
 * Convert velocity to momentum by multiplying by density averaged onto faces
 */
void VelocityToMomentum (const amrex::MultiFab& xvel_in,
                         const amrex::IntVect& xvel_ngrow,
                         const amrex::MultiFab& yvel_in,
                         const amrex::IntVect& yvel_ngrow,
                         const amrex::MultiFab& zvel_in,
                         const amrex::IntVect& zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         amrex::MultiFab& xmom_out,
                         amrex::MultiFab& ymom_out,
                         amrex::MultiFab& zmom_out,
                         bool l_use_ndiff);

/*
 * Compute boxes for looping over interior/exterior ghost cells
 * for use by fillpatch, erf_slow_rhs_pre, and erf_slow_rhs_post
 */
void compute_interior_ghost_bxs_xy (const amrex::Box& bx,
                                    const amrex::Box& domain,
                                    const int& width,
                                    const int& set_width,
                                    amrex::Box& bx_xlo,
                                    amrex::Box& bx_xhi,
                                    amrex::Box& bx_ylo,
                                    amrex::Box& bx_yhi,
                                    const amrex::IntVect& ng_vect=amrex::IntVect(0,0,0),
                                    const bool get_int_ng=false);

/*
 * Compute relaxation region RHS with wrfbdy
 */
void wrfbdy_compute_interior_ghost_rhs (const std::string& init_type,
                                        const amrex::Real& bdy_time_interval,
                                        const amrex::Real& start_bdy_time,
                                        const amrex::Real& time,
                                        const amrex::Real& delta_t,
                                        const int& width,
                                        const int& set_width,
                                        const amrex::Geometry& geom,
                                        amrex::Vector<amrex::MultiFab>& S_rhs,
                                        amrex::Vector<amrex::MultiFab>& S_data,
                                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi);

/*
 * Compute relaxation region RHS at fine-crse interface
 */
void
fine_compute_interior_ghost_rhs (const amrex::Real& time,
                                 const amrex::Real& delta_t,
                                 const int& width,
                                 const int& set_width,
                                 const amrex::Geometry& geom,
                                 ERFFillPatcher* FPr_c,
                                 ERFFillPatcher* FPr_u,
                                 ERFFillPatcher* FPr_v,
                                 ERFFillPatcher* FPr_w,
                                 amrex::Vector<amrex::BCRec>& domain_bcs_type,
                                 amrex::Vector<amrex::MultiFab>& S_rhs_f,
                                 amrex::Vector<amrex::MultiFab>& S_data_f);

/**
 * Zero RHS in the set region
 *
 * @param[in] icomp component offset
 * @param[in] num_var number of variables to loop
 * @param[in] bx_xlo box for low x relaxation
 * @param[in] bx_xhi box for high x relaxation
 * @param[in] bx_ylo box for low y relaxation
 * @param[in] bx_yhi box for high y relaxation
 * @param[out] rhs_arr RHS array
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
wrfbdy_zero_rhs_in_set_region (const int& icomp,
                               const int& num_var,
                               const amrex::Box& bx_xlo,
                               const amrex::Box& bx_xhi,
                               const amrex::Box& bx_ylo,
                               const amrex::Box& bx_yhi,
                               const amrex::Array4<amrex::Real>& rhs_arr)
{
    amrex::ParallelFor(bx_xlo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
      rhs_arr(i,j,k,n+icomp) = 0.0;
    },
    bx_xhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        rhs_arr(i,j,k,n+icomp) = 0.0;
    });

    amrex::ParallelFor(bx_ylo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        rhs_arr(i,j,k,n+icomp) = 0.0;
    },
    bx_yhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        rhs_arr(i,j,k,n+icomp) = 0.0;
    });
}

/**
 * Compute the Laplacian RHS in the relaxation zone
 *
 * @param[in] delta_t time step
 * @param[in] icomp component offset
 * @param[in] num_var number of variables to loop
 * @param[in] width width of wrf bdy file
 * @param[in] set_width width the set region
 * @param[in] dom_lo low bound of domain
 * @param[in] dom_hi high bound of domain
 * @param[in] F1 drift relaxation parameter
 * @param[in] F2 Laplacian relaxation parameter
 * @param[in] bx_xlo box for low x relaxation
 * @param[in] bx_xhi box for high x relaxation
 * @param[in] bx_ylo box for low y relaxation
 * @param[in] bx_yhi box for high y relaxation
 * @param[in] arr_xlo array for low x relaxation
 * @param[in] arr_xhi array for high x relaxation
 * @param[in] arr_ylo array for low y relaxation
 * @param[in] arr_yhi array for high y relaxation
 * @param[in] data_arr data array
 * @param[out] rhs_arr RHS array
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
wrfbdy_compute_laplacian_relaxation (const amrex::Real& delta_t,
                                     const int& icomp,
                                     const int& num_var,
                                     const int& width,
                                     const int& set_width,
                                     const amrex::Dim3& dom_lo,
                                     const amrex::Dim3& dom_hi,
                                     const amrex::Real& F1,
                                     const amrex::Real& F2,
                                     const amrex::Box& bx_xlo,
                                     const amrex::Box& bx_xhi,
                                     const amrex::Box& bx_ylo,
                                     const amrex::Box& bx_yhi,
                                     const amrex::Array4<const amrex::Real>& arr_xlo,
                                     const amrex::Array4<const amrex::Real>& arr_xhi,
                                     const amrex::Array4<const amrex::Real>& arr_ylo,
                                     const amrex::Array4<const amrex::Real>& arr_yhi,
                                     const amrex::Array4<const amrex::Real>& data_arr,
                                     const amrex::Array4<amrex::Real>& rhs_arr)
{
    // RHS computation
    int Spec_z  = set_width;
    int Relax_z = width - Spec_z + 1;
    amrex::Real num    = amrex::Real(Spec_z + Relax_z);
    amrex::Real denom  = amrex::Real(Relax_z - 1);
    amrex::ParallelFor(bx_xlo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // Corners with x boxes
        int j_lo  = std::min(j-dom_lo.y,width-1);
        int j_hi  = std::min(dom_hi.y-j,width-1);
        int jj    = std::min(j_lo,j_hi);
        int n_ind = std::min(i-dom_lo.x,jj) + 1;
        if (n_ind <= Spec_z) {
            rhs_arr(i,j,k,n) = 0.0;
        } else {
            amrex::Real Factor   = (num - amrex::Real(n_ind))/denom;
            amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
            amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
            amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
            amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
            amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
            amrex::Real delta    = arr_xlo(i  ,j  ,k,n) - d;
            amrex::Real delta_xp = arr_xlo(i+1,j  ,k,n) - d_ip1;
            amrex::Real delta_xm = arr_xlo(i-1,j  ,k,n) - d_im1;
            amrex::Real delta_yp = arr_xlo(i  ,j+1,k,n) - d_jp1;
            amrex::Real delta_ym = arr_xlo(i  ,j-1,k,n) - d_jm1;
            amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
            rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;

            /*
            if ((std::abs((arr_xlo(i  ,j  ,k,n)-data_arr(i  ,j  ,k  ,n+icomp))/data_arr(i  ,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xlo(i+1,j  ,k,n)-data_arr(i+1,j  ,k  ,n+icomp))/data_arr(i+1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xlo(i-1,j  ,k,n)-data_arr(i-1,j  ,k  ,n+icomp))/data_arr(i-1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xlo(i  ,j+1,k,n)-data_arr(i  ,j+1,k  ,n+icomp))/data_arr(i  ,j+1,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xlo(i  ,j-1,k,n)-data_arr(i  ,j-1,k  ,n+icomp))/data_arr(i  ,j-1,k  ,n+icomp)) > 0.5) ) {
                amrex::Print() << "XLO error:\n";
                amrex::Print() << "RELAX x: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j,k) << ' '
                               << data_arr(i  ,j  ,k  ,n+icomp) << ' ' << arr_xlo(i  ,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i+1,j,k) << ' '
                               << data_arr(i+1,j  ,k  ,n+icomp) << ' ' << arr_xlo(i+1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xm: " << icomp << ' ' << n << ' ' << amrex::IntVect(i-1,j,k) << ' '
                               << data_arr(i-1,j  ,k  ,n+icomp) << ' ' << arr_xlo(i-1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX yp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j+1,k) << ' '
                               << data_arr(i  ,j+1,k  ,n+icomp) << ' ' << arr_xlo(i  ,j+1,k,n) << "\n";
                amrex::Print() << "RELAX ym: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j-1,k) << ' '
                               << data_arr(i  ,j-1,k  ,n+icomp) << ' ' << arr_xlo(i  ,j-1,k,n) << "\n";
                amrex::Print() << "\n";
                //exit(0);
            }
            */
        }
    },
    bx_xhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // Corners with x boxes
        int j_lo  = std::min(j-dom_lo.y,width-1);
        int j_hi  = std::min(dom_hi.y-j,width-1);
        int jj    = std::min(j_lo,j_hi);
        int n_ind = std::min(dom_hi.x-i,jj) + 1;
        if (n_ind <= Spec_z) {
            rhs_arr(i,j,k,n) = 0.0;
        } else {
            amrex::Real Factor   = (num - amrex::Real(n_ind))/denom;
            amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
            amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
            amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
            amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
            amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
            amrex::Real delta    = arr_xhi(i  ,j  ,k,n) - d;
            amrex::Real delta_xp = arr_xhi(i+1,j  ,k,n) - d_ip1;
            amrex::Real delta_xm = arr_xhi(i-1,j  ,k,n) - d_im1;
            amrex::Real delta_yp = arr_xhi(i  ,j+1,k,n) - d_jp1;
            amrex::Real delta_ym = arr_xhi(i  ,j-1,k,n) - d_jm1;
            amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
            rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;

            /*
            if ((std::abs((arr_xhi(i  ,j  ,k,n)-data_arr(i  ,j  ,k  ,n+icomp))/data_arr(i  ,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xhi(i+1,j  ,k,n)-data_arr(i+1,j  ,k  ,n+icomp))/data_arr(i+1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xhi(i-1,j  ,k,n)-data_arr(i-1,j  ,k  ,n+icomp))/data_arr(i-1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xhi(i  ,j+1,k,n)-data_arr(i  ,j+1,k  ,n+icomp))/data_arr(i  ,j+1,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_xhi(i  ,j-1,k,n)-data_arr(i  ,j-1,k  ,n+icomp))/data_arr(i  ,j-1,k  ,n+icomp)) > 0.5) ) {
                amrex::Print() << "XHI error:\n";
                amrex::Print() << "RELAX x: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j,k) << ' '
                               << data_arr(i  ,j  ,k  ,n+icomp) << ' ' << arr_xhi(i  ,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i+1,j,k) << ' '
                               << data_arr(i+1,j  ,k  ,n+icomp) << ' ' << arr_xhi(i+1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xm: " << icomp << ' ' << n << ' ' << amrex::IntVect(i-1,j,k) << ' '
                               << data_arr(i-1,j  ,k  ,n+icomp) << ' ' << arr_xhi(i-1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX yp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j+1,k) << ' '
                               << data_arr(i  ,j+1,k  ,n+icomp) << ' ' << arr_xhi(i  ,j+1,k,n) << "\n";
                amrex::Print() << "RELAX ym: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j-1,k) << ' '
                               << data_arr(i  ,j-1,k  ,n+icomp) << ' ' << arr_xhi(i  ,j-1,k,n) << "\n";
                amrex::Print() << "\n";
                //exit(0);
            }
            */
        }
    });

    amrex::ParallelFor(bx_ylo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = j - dom_lo.y + 1;
        if (n_ind <= Spec_z) {
            rhs_arr(i,j,k,n) = 0.0;
        } else {
            amrex::Real Factor   = (num - amrex::Real(n_ind))/denom;
            amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
            amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
            amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
            amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
            amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
            amrex::Real delta    = arr_ylo(i  ,j  ,k,n) - d;
            amrex::Real delta_xp = arr_ylo(i+1,j  ,k,n) - d_ip1;
            amrex::Real delta_xm = arr_ylo(i-1,j  ,k,n) - d_im1;
            amrex::Real delta_yp = arr_ylo(i  ,j+1,k,n) - d_jp1;
            amrex::Real delta_ym = arr_ylo(i  ,j-1,k,n) - d_jm1;
            amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
            rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;

            /*
            if ((std::abs((arr_ylo(i  ,j  ,k,n)-data_arr(i  ,j  ,k  ,n+icomp))/data_arr(i  ,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_ylo(i+1,j  ,k,n)-data_arr(i+1,j  ,k  ,n+icomp))/data_arr(i+1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_ylo(i-1,j  ,k,n)-data_arr(i-1,j  ,k  ,n+icomp))/data_arr(i-1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_ylo(i  ,j+1,k,n)-data_arr(i  ,j+1,k  ,n+icomp))/data_arr(i  ,j+1,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_ylo(i  ,j-1,k,n)-data_arr(i  ,j-1,k  ,n+icomp))/data_arr(i  ,j-1,k  ,n+icomp)) > 0.5) ) {
                amrex::Print() << "YLO error:\n";
                amrex::Print() << "RELAX x: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j,k) << ' '
                               << data_arr(i  ,j  ,k  ,n+icomp) << ' ' << arr_ylo(i  ,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i+1,j,k) << ' '
                               << data_arr(i+1,j  ,k  ,n+icomp) << ' ' << arr_ylo(i+1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xm: " << icomp << ' ' << n << ' ' << amrex::IntVect(i-1,j,k) << ' '
                               << data_arr(i-1,j  ,k  ,n+icomp) << ' ' << arr_ylo(i-1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX yp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j+1,k) << ' '
                               << data_arr(i  ,j+1,k  ,n+icomp) << ' ' << arr_ylo(i  ,j+1,k,n) << "\n";
                amrex::Print() << "RELAX ym: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j-1,k) << ' '
                               << data_arr(i  ,j-1,k  ,n+icomp) << ' ' << arr_ylo(i  ,j-1,k,n) << "\n";
                amrex::Print() << "\n";
                //exit(0);
            }
            */
        }
    },
    bx_yhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = dom_hi.y - j + 1;
        if (n_ind <= Spec_z) {
            rhs_arr(i,j,k,n) = 0.0;
        } else {
            amrex::Real Factor   = (num - amrex::Real(n_ind))/denom;
            amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp) + delta_t*rhs_arr(i  , j  , k  ,n+icomp);
            amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i+1, j  , k  ,n+icomp);
            amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp) + delta_t*rhs_arr(i-1, j  , k  ,n+icomp);
            amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp) + delta_t*rhs_arr(i  , j+1, k  ,n+icomp);
            amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp) + delta_t*rhs_arr(i  , j-1, k  ,n+icomp);
            amrex::Real delta    = arr_yhi(i  ,j  ,k,n) - d;
            amrex::Real delta_xp = arr_yhi(i+1,j  ,k,n) - d_ip1;
            amrex::Real delta_xm = arr_yhi(i-1,j  ,k,n) - d_im1;
            amrex::Real delta_yp = arr_yhi(i  ,j+1,k,n) - d_jp1;
            amrex::Real delta_ym = arr_yhi(i  ,j-1,k,n) - d_jm1;
            amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
            rhs_arr(i,j,k,n) += (F1*delta - F2*Laplacian) * Factor;

            /*
            if ((std::abs((arr_yhi(i  ,j  ,k,n)-data_arr(i  ,j  ,k  ,n+icomp))/data_arr(i  ,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_yhi(i+1,j  ,k,n)-data_arr(i+1,j  ,k  ,n+icomp))/data_arr(i+1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_yhi(i-1,j  ,k,n)-data_arr(i-1,j  ,k  ,n+icomp))/data_arr(i-1,j  ,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_yhi(i  ,j+1,k,n)-data_arr(i  ,j+1,k  ,n+icomp))/data_arr(i  ,j+1,k  ,n+icomp)) > 0.5) ||
                (std::abs((arr_yhi(i  ,j-1,k,n)-data_arr(i  ,j-1,k  ,n+icomp))/data_arr(i  ,j-1,k  ,n+icomp)) > 0.5) ) {
                amrex::Print() << "YHI error:\n";
                amrex::Print() << "RELAX x: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j,k) << ' '
                               << data_arr(i  ,j  ,k  ,n+icomp) << ' ' << arr_yhi(i  ,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i+1,j,k) << ' '
                               << data_arr(i+1,j  ,k  ,n+icomp) << ' ' << arr_yhi(i+1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX xm: " << icomp << ' ' << n << ' ' << amrex::IntVect(i-1,j,k) << ' '
                               << data_arr(i-1,j  ,k  ,n+icomp) << ' ' << arr_yhi(i-1,j  ,k,n) << "\n";
                amrex::Print() << "RELAX yp: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j+1,k) << ' '
                               << data_arr(i  ,j+1,k  ,n+icomp) << ' ' << arr_yhi(i  ,j+1,k,n) << "\n";
                amrex::Print() << "RELAX ym: " << icomp << ' ' << n << ' ' << amrex::IntVect(i,j-1,k) << ' '
                               << data_arr(i  ,j-1,k  ,n+icomp) << ' ' << arr_yhi(i  ,j-1,k,n) << "\n";
                amrex::Print() << "\n";
                //exit(0);
            }
            */
        }
    });
}
#endif
