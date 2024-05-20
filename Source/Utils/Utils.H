#ifndef _UTILS_H_
#define _UTILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <DataStruct.H>
#include <IndexDefines.H>
#include <ABLMost.H>
#include <ERF_FillPatcher.H>

#include <Microphysics_Utils.H>
#include <ERF_Constants.H>
#include <EOS.H>

/*
 * Create the Jacobian for the metric transformation when use_terrain is true
 */
void make_J (const amrex::Geometry& geom,
             amrex::MultiFab& z_phys_nd,
             amrex::MultiFab& detJ_cc);

void make_areas (const amrex::Geometry& geom,
                 amrex::MultiFab& z_phys_nd,
                 amrex::MultiFab& ax,
                 amrex::MultiFab& ay,
                 amrex::MultiFab& az);

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
                         const amrex::MultiFab& zmom_in,
                         const amrex::Box& domain,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type_h);

/*
 * Convert velocity to momentum by multiplying by density averaged onto faces
 */
void VelocityToMomentum (const amrex::MultiFab& xvel_in,
                         const amrex::IntVect & xvel_ngrow,
                         const amrex::MultiFab& yvel_in,
                         const amrex::IntVect & yvel_ngrow,
                         const amrex::MultiFab& zvel_in,
                         const amrex::IntVect & zvel_ngrow,
                         const amrex::MultiFab& cons_in,
                         amrex::MultiFab& xmom_out,
                         amrex::MultiFab& ymom_out,
                         amrex::MultiFab& zmom_out,
                         const amrex::Box& domain,
                         const amrex::Vector<amrex::BCRec>& domain_bcs_type_h);

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
void realbdy_compute_interior_ghost_rhs (const std::string& init_type,
                                         const amrex::Real& bdy_time_interval,
                                         const amrex::Real& start_bdy_time,
                                         const amrex::Real& time,
                                         const amrex::Real& delta_t,
                                         int width,
                                         int set_width,
                                         const amrex::Geometry& geom,
                                         amrex::Vector<amrex::MultiFab>& S_rhs,
                                         amrex::Vector<amrex::MultiFab>& S_old_data,
                                         amrex::Vector<amrex::MultiFab>& S_cur_data,
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

/*
 * Accumulate time averaged velocity fields
 */
void
Time_Avg_Vel_atCC (const amrex::Real& dt,
                   amrex::Real& t_avg_cnt,
                   amrex::MultiFab* vel_t_avg,
                   amrex::MultiFab& xvel,
                   amrex::MultiFab& yvel,
                   amrex::MultiFab& zvel);

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
wrfbdy_set_rhs_in_spec_region (const amrex::Real& dt,
                               const int& icomp,
                               const int& num_var,
                               const int& width,
                               const int& set_width,
                               const amrex::Dim3& dom_lo,
                               const amrex::Dim3& dom_hi,
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
    int Spec_z  = set_width;
    amrex::ParallelFor(bx_xlo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // Corners with x boxes
        int j_lo  = std::min(j-dom_lo.y,width-1);
        int j_hi  = std::min(dom_hi.y-j,width-1);
        int jj    = std::min(j_lo,j_hi);
        int n_ind = std::min(i-dom_lo.x,jj) + 1;
        if (n_ind <= Spec_z) {
            amrex::Real tend = ( arr_xlo(i,j,k) - data_arr(i,j,k,n+icomp) ) / dt;
            rhs_arr(i,j,k,n+icomp) = tend;
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
            amrex::Real tend = ( arr_xhi(i,j,k) - data_arr(i,j,k,n+icomp) ) / dt;
            rhs_arr(i,j,k,n+icomp) = tend;
        }
    });

    amrex::ParallelFor(bx_ylo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = j - dom_lo.y + 1;
        if (n_ind <= Spec_z) {
            amrex::Real tend = ( arr_ylo(i,j,k) - data_arr(i,j,k,n+icomp) ) / dt;
            rhs_arr(i,j,k,n+icomp) = tend;
        }
    },
    bx_yhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = dom_hi.y - j + 1;
        if (n_ind <= Spec_z) {
            amrex::Real tend = ( arr_yhi(i,j,k) - data_arr(i,j,k,n+icomp) ) / dt;
            rhs_arr(i,j,k,n+icomp) = tend;
        }
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
wrfbdy_compute_laplacian_relaxation (const int& icomp,
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
    amrex::Real num     = amrex::Real(Spec_z + Relax_z);
    amrex::Real denom   = amrex::Real(Relax_z - 1);
    amrex::Real SpecExp = -std::log(0.01) / amrex::Real(width - Spec_z);
    amrex::ParallelFor(bx_xlo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // Corners with x boxes
        int j_lo  = std::min(j-dom_lo.y,width-1);
        int j_hi  = std::min(dom_hi.y-j,width-1);
        int jj    = std::min(j_lo,j_hi);
        int n_ind = std::min(i-dom_lo.x,jj) + 1;
        AMREX_ASSERT(n_ind > Spec_z);
        amrex::Real Factor   = (num - amrex::Real(n_ind))/denom
                              * std::exp(-SpecExp * amrex::Real(n_ind - Spec_z));
        amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp);
        amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp);
        amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp);
        amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp);
        amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp);
        amrex::Real delta    = arr_xlo(i  ,j  ,k,n) - d;
        amrex::Real delta_xp = arr_xlo(i+1,j  ,k,n) - d_ip1;
        amrex::Real delta_xm = arr_xlo(i-1,j  ,k,n) - d_im1;
        amrex::Real delta_yp = arr_xlo(i  ,j+1,k,n) - d_jp1;
        amrex::Real delta_ym = arr_xlo(i  ,j-1,k,n) - d_jm1;
        amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
        amrex::Real Temp = (F1*delta - F2*Laplacian) * Factor;
        rhs_arr(i,j,k,n+icomp) += Temp;
    },
    bx_xhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // Corners with x boxes
        int j_lo  = std::min(j-dom_lo.y,width-1);
        int j_hi  = std::min(dom_hi.y-j,width-1);
        int jj    = std::min(j_lo,j_hi);
        int n_ind = std::min(dom_hi.x-i,jj) + 1;
        AMREX_ASSERT(n_ind > Spec_z);
        amrex::Real Factor   = (num - amrex::Real(n_ind))/denom
                              * std::exp(-SpecExp * amrex::Real(n_ind - Spec_z));
        amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp);
        amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp);
        amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp);
        amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp);
        amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp);
        amrex::Real delta    = arr_xhi(i  ,j  ,k,n) - d;
        amrex::Real delta_xp = arr_xhi(i+1,j  ,k,n) - d_ip1;
        amrex::Real delta_xm = arr_xhi(i-1,j  ,k,n) - d_im1;
        amrex::Real delta_yp = arr_xhi(i  ,j+1,k,n) - d_jp1;
        amrex::Real delta_ym = arr_xhi(i  ,j-1,k,n) - d_jm1;
        amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
        amrex::Real Temp = (F1*delta - F2*Laplacian) * Factor;
        rhs_arr(i,j,k,n+icomp) += Temp;
    });

    amrex::ParallelFor(bx_ylo, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = j - dom_lo.y + 1;
        AMREX_ASSERT(n_ind > Spec_z);
        amrex::Real Factor   = (num - amrex::Real(n_ind))/denom
                              * std::exp(-SpecExp * amrex::Real(n_ind - Spec_z));
        amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp);
        amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp);
        amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp);
        amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp);
        amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp);
        amrex::Real delta    = arr_ylo(i  ,j  ,k,n) - d;
        amrex::Real delta_xp = arr_ylo(i+1,j  ,k,n) - d_ip1;
        amrex::Real delta_xm = arr_ylo(i-1,j  ,k,n) - d_im1;
        amrex::Real delta_yp = arr_ylo(i  ,j+1,k,n) - d_jp1;
        amrex::Real delta_ym = arr_ylo(i  ,j-1,k,n) - d_jm1;
        amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
        amrex::Real Temp = (F1*delta - F2*Laplacian) * Factor;
        rhs_arr(i,j,k,n+icomp) += Temp;
    },
    bx_yhi, num_var, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        // No corners for y boxes
        int n_ind = dom_hi.y - j + 1;
        AMREX_ASSERT(n_ind > Spec_z);
        amrex::Real Factor   = (num - amrex::Real(n_ind))/denom
                              * std::exp(-SpecExp * amrex::Real(n_ind - Spec_z));
        amrex::Real d        = data_arr(i  ,j  ,k  ,n+icomp);
        amrex::Real d_ip1    = data_arr(i+1,j  ,k  ,n+icomp);
        amrex::Real d_im1    = data_arr(i-1,j  ,k  ,n+icomp);
        amrex::Real d_jp1    = data_arr(i  ,j+1,k  ,n+icomp);
        amrex::Real d_jm1    = data_arr(i  ,j-1,k  ,n+icomp);
        amrex::Real delta    = arr_yhi(i  ,j  ,k,n) - d;
        amrex::Real delta_xp = arr_yhi(i+1,j  ,k,n) - d_ip1;
        amrex::Real delta_xm = arr_yhi(i-1,j  ,k,n) - d_im1;
        amrex::Real delta_yp = arr_yhi(i  ,j+1,k,n) - d_jp1;
        amrex::Real delta_ym = arr_yhi(i  ,j-1,k,n) - d_jm1;
        amrex::Real Laplacian = delta_xp + delta_xm + delta_yp + delta_ym - 4.0*delta;
        amrex::Real Temp = (F1*delta - F2*Laplacian) * Factor;
        rhs_arr(i,j,k,n+icomp) += Temp;
    });
}

/*
 * Effectively a Multiply for a MultiFab and an iMultiFab mask
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
ApplyMask (amrex::MultiFab& dst,
           const amrex::iMultiFab& imask,
           const int nghost = 0)
{
    for (amrex::MFIter mfi(dst,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto        dstFab = dst.array(mfi);
            const auto maskFab = imask.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dstFab(i,j,k) *= maskFab(i,j,k);
            });
        }
    }
}
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void
ApplyInvertedMask (amrex::MultiFab& dst,
                   const amrex::iMultiFab& imask,
                   const int nghost = 0)
{
    for (amrex::MFIter mfi(dst,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(nghost);
        if (bx.ok())
        {
            auto        dstFab = dst.array(mfi);
            const auto maskFab = imask.const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dstFab(i,j,k) *= (1-maskFab(i,j,k));
            });
        }
    }
}

/*
 * Compute the Brunt-Vaisala frequency for dry and moist environments
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
amrex::Real
Brunt_Vaisala_Freq (int& i,
                    int& j,
                    int& k,
                    amrex::Real& dzInv,
                    const amrex::Array4<const amrex::Real>& cell_data)
{
    // NOTE: The effect of the Brunt-Vaisala freq is included here. We account
    //       for moist and dry conditions similar to  DK 82
    //       (https://doi.org/10.1175/1520-0469(1982)039<2152:OTEOMO>2.0.CO;2)
    //       The exact implementation follows the WRF 2019 documentation and can
    //       be found in Sec 4.2.4
    amrex::Real N2;
    amrex::Real T, P;
    amrex::Real Th, Th_hi, Th_lo;
    bool use_moisture = (cell_data.ncomp >= RhoQ2_comp+1);
    Th = cell_data(i,j,k,RhoTheta_comp)/cell_data(i,j,k,Rho_comp);
    if (use_moisture) {
        amrex::Real qv, qv_hi, qv_lo;
        amrex::Real qc, qw_hi, qw_lo;
        amrex::Real qsat, qsat_hi, qsat_lo;

        // Qvs at cc
        qv = cell_data(i,j,k,RhoQ1_comp   )/cell_data(i,j,k,Rho_comp);
        qc = cell_data(i,j,k,RhoQ2_comp   )/cell_data(i,j,k,Rho_comp);
        T  = getTgivenRandRTh(cell_data(i,j,k,Rho_comp),
                              cell_data(i,j,k,RhoTheta_comp),
                              qv);
        P  = getPgivenRTh(cell_data(i,j,k,RhoTheta_comp), qv) * 0.01; // convert to mbar
        erf_qsatw(T, P, qsat);

        if ( (qv > qsat) || (qc > qc_crit) ) {
            // Compute A coefficient
            amrex::Real A = ( 1. + Rv_on_Rd * (mwwv/mwdair) * L_v * qv / (R_d * T) );
            A /= ( 1. + (mwwv/mwdair) * L_v * L_v * qv / (Cp_d * R_v * T * T) );
            A *= (1./Th);

            // Comput qw on hi and lo face
            qw_hi = cell_data(i,j,k+1,RhoQ1_comp)/cell_data(i,j,k+1,Rho_comp)
                  + cell_data(i,j,k+1,RhoQ2_comp)/cell_data(i,j,k+1,Rho_comp);
            qw_lo = cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp)
                  + cell_data(i,j,k-1,RhoQ2_comp)/cell_data(i,j,k-1,Rho_comp);

            // Compute Theta_e on hi face
            qv_hi = cell_data(i,j,k+1,RhoQ1_comp)/cell_data(i,j,k+1,Rho_comp);
            T  = getTgivenRandRTh(cell_data(i,j,k+1,Rho_comp),
                                  cell_data(i,j,k+1,RhoTheta_comp),
                                  qv_hi);
            P  = getPgivenRTh(cell_data(i,j,k+1,RhoTheta_comp), qv_hi) * 0.01; // convert to mbar
            erf_qsatw(T, P, qsat_hi);
            Th_hi  = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
            Th_hi *= ( 1. + L_v * qsat_hi / (Cp_d * T) );

            // Compute Theta_e on lo face
            qv_lo = cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp);
            T  = getTgivenRandRTh(cell_data(i,j,k-1,Rho_comp),
                                  cell_data(i,j,k-1,RhoTheta_comp),
                                  qv_lo);
            P  = getPgivenRTh(cell_data(i,j,k-1,RhoTheta_comp), qv_lo) * 0.01; // convert to mbar
            erf_qsatw(T, P, qsat_lo);
            Th_lo  = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);
            Th_lo *= ( 1. + L_v * qsat_lo / (Cp_d * T) );

            // Compute N2
            N2 = CONST_GRAV * ( A * 0.5 * ( Th_hi - Th_lo ) * dzInv
                                  - 0.5 * ( qw_hi - qw_lo ) * dzInv );
        } else {
            // Comput qw on hi and lo face
            qv_hi = cell_data(i,j,k+1,RhoQ1_comp)/cell_data(i,j,k+1,Rho_comp);
            qv_lo = cell_data(i,j,k-1,RhoQ1_comp)/cell_data(i,j,k-1,Rho_comp);

            // Comput qw on hi and lo face
            qw_hi = qv_hi + cell_data(i,j,k+1,RhoQ2_comp)/cell_data(i,j,k+1,Rho_comp);
            qw_lo = qv_lo + cell_data(i,j,k-1,RhoQ2_comp)/cell_data(i,j,k-1,Rho_comp);

            // Compute theta on hi and lo face
            Th_hi  = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
            Th_lo  = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);

            // Compute N2
            N2 = CONST_GRAV * ( (1./Th) * 0.5 * ( Th_hi - Th_lo ) * dzInv
                             + Rv_on_Rd * 0.5 * ( qv_hi - qv_lo ) * dzInv
                                        - 0.5 * ( qw_hi - qw_lo ) * dzInv );
        }
    // No moisture present
    } else {
        // Compute theta on hi and lo face
        Th_hi  = cell_data(i,j,k+1,RhoTheta_comp)/cell_data(i,j,k+1,Rho_comp);
        Th_lo  = cell_data(i,j,k-1,RhoTheta_comp)/cell_data(i,j,k-1,Rho_comp);

        // Compute N2
        N2 = CONST_GRAV * ( (1./Th) * 0.5 * ( Th_hi - Th_lo ) * dzInv );
    }

    return N2;
}
#endif
