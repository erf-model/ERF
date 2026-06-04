#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "AMReX_FArrayBox.H"
#include "AMReX_Print.H"

#include "ERF_NCWpsFile.H"
#include "ERF_IndexDefines.H"
#include "ERF_EOS.H"
#include "ERF_DataStruct.H"
#include "ERF_NCInterface.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

namespace WRFBdyTypes {
    enum {
        x_lo,
        x_hi,
        y_lo,
        y_hi
    };
}

Real
read_times_from_wrfbdy (const std::string& nc_bdy_file,
                        Vector<Vector<FArrayBox>>& bdy_data_xlo,
                        Vector<Vector<FArrayBox>>& bdy_data_xhi,
                        Vector<Vector<FArrayBox>>& bdy_data_ylo,
                        Vector<Vector<FArrayBox>>& bdy_data_yhi,
                        Real& start_bdy_time,
                        Real& final_bdy_time)
{
    Print() << "Loading boundary data from NetCDF file " << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // *******************************************************************************

    int ntimes;
    Real timeInterval;
    const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor())
    {
        // Read the time stamps
        using CharArray = NDArray<char>;
        Vector<CharArray> array_ts(1);
        Vector<int> success(1);
        ReadNetCDFFile(nc_bdy_file, {"Times"}, array_ts, success);

        ntimes = array_ts[0].get_vshape()[0];

        Vector<std::string> timeStamps;
        timeStamps.reserve(ntimes);

        const char* data = array_ts[0].get_data();
        auto dateStrLen  = array_ts[0].get_vshape()[1];

        for (int nt = 0; nt < ntimes; ++nt) {
            const char* begin = data + nt * dateStrLen;
            timeStamps.emplace_back(begin, begin + dateStrLen);
        }

        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string date = timeStamps[nt];
            auto epochTime = getEpochTime(date, dateTimeFormat);
            Print() << "  wrfbdy datetime " << nt << " : " << date << " " << epochTime << std::endl;
            epochTimes.push_back(epochTime);

            if (nt == 1) {
                timeInterval = static_cast<Real>(epochTimes[1] - epochTimes[0]);
            } else if (nt >= 1) {
                AMREX_ALWAYS_ASSERT(static_cast<Real>(epochTimes[nt] - epochTimes[nt-1]) == timeInterval);
            }
        }
        start_bdy_time = static_cast<Real>(epochTimes[0]);
        final_bdy_time = static_cast<Real>(epochTimes[ntimes-1] + timeInterval);
        Print() << "  start_bdy_time " << start_bdy_time << std::endl;
        Print() << "  final_bdy_time " << final_bdy_time << std::endl;
    }

    ParallelDescriptor::Bcast(&start_bdy_time,1,ioproc);
    ParallelDescriptor::Bcast(&final_bdy_time,1,ioproc);
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Our outermost loop is time
    bdy_data_xlo.resize(ntimes+1);
    bdy_data_xhi.resize(ntimes+1);
    bdy_data_ylo.resize(ntimes+1);
    bdy_data_yhi.resize(ntimes+1);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

void
convert_wrfbdy_data (const int itime,
                     const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     std::unique_ptr<MultiFab>& wrf_MUB,
                     std::unique_ptr<MultiFab>& wrf_C1H,
                     std::unique_ptr<MultiFab>& wrf_C2H,
                     std::unique_ptr<MultiFab>& wrf_PHB,
                     const iMultiFab* mask_u,
                     const iMultiFab* mask_v,
                     const iMultiFab* mask_c,
                     const bool& use_moist)
{
    // Temporary bdy data structures for global reductions
    int vsize = bdy_data[itime].size() - 3; // Don't do PH, MU, or PC
    amrex::Vector<amrex::FArrayBox> bdy_data_tmp; bdy_data_tmp.resize(vsize);
    for (int ivar(0); ivar < vsize; ++ivar) {
        bdy_data_tmp[ivar].resize(bdy_data[itime][ivar].box(),1,The_Managed_Arena());
        bdy_data_tmp[ivar].template setVal<RunOn::Device>(0);
    }

    // Temporary bdy data structures for interpolation
    amrex::Vector<amrex::FArrayBox> bdy_data_int; bdy_data_int.resize(vsize);
    for (int ivar(0); ivar < vsize; ++ivar) {
        bdy_data_int[ivar].resize(bdy_data[itime][ivar].box(),1,The_Managed_Arena());
        bdy_data_int[ivar].template setVal<RunOn::Device>(0);
    }

    // Temporary "NEW" heights (this is the source array)
    amrex::FArrayBox bdy_c_z_new, bdy_u_z_new, bdy_v_z_new;
    bdy_c_z_new.resize(bdy_data[itime][WRFBdyVars::T].box(),1,The_Managed_Arena());
    bdy_u_z_new.resize(bdy_data[itime][WRFBdyVars::U].box(),1,The_Managed_Arena());
    bdy_v_z_new.resize(bdy_data[itime][WRFBdyVars::V].box(),1,The_Managed_Arena());
    Array4<Real> bdy_c_z_src = bdy_c_z_new.array();
    Array4<Real> bdy_u_z_src = bdy_u_z_new.array();
    Array4<Real> bdy_v_z_src = bdy_v_z_new.array();

    // Temporary "OLD" heights (these are the heights to interpolate to at itime=0)
    amrex::FArrayBox bdy_c_z_old, bdy_u_z_old, bdy_v_z_old;
    bdy_c_z_old.resize(bdy_data[0][WRFBdyVars::T].box(),1,The_Managed_Arena());
    bdy_u_z_old.resize(bdy_data[0][WRFBdyVars::U].box(),1,The_Managed_Arena());
    bdy_v_z_old.resize(bdy_data[0][WRFBdyVars::V].box(),1,The_Managed_Arena());
    Array4<Real> bdy_c_z_dst = bdy_c_z_old.array();
    Array4<Real> bdy_u_z_dst = bdy_u_z_old.array();
    Array4<Real> bdy_v_z_dst = bdy_v_z_old.array();

    // BDY data
    Array4<Real> bdy_u_arr  = bdy_data[itime][WRFBdyVars::U].array();  // This is x-face-centered
    Array4<Real> bdy_v_arr  = bdy_data[itime][WRFBdyVars::V].array();  // This is y-face-centered
    Array4<Real> bdy_t_arr  = bdy_data[itime][WRFBdyVars::T].array();  // This is cell-centered
    Array4<Real> bdy_qv_arr = bdy_data[itime][WRFBdyVars::QV].array(); // This is cell-centered
    Array4<Real> mu_arr     = bdy_data[itime][WRFBdyVars::MU].array(); // This is cell-centered
    Array4<Real> bdy_ph_arr = bdy_data[itime][WRFBdyVars::PH].array(); // This is z-face-centered

    // For height interpolation (removes averaging error)
    Array4<Real> mu0_arr     = bdy_data[0][WRFBdyVars::MU].array(); // This is cell-centered
    Array4<Real> bdy_ph0_arr = bdy_data[0][WRFBdyVars::PH].array(); // This is z-face-centered

    // Bounds limiting
    int ilo  = domain.smallEnd()[0];
    int ihi  = domain.bigEnd()[0];
    int jlo  = domain.smallEnd()[1];
    int jhi  = domain.bigEnd()[1];
    int klo  = domain.smallEnd()[2];
    int khi  = domain.bigEnd()[2];

    // PH bounds limiting
    Box ph_bx  = bdy_data[itime][WRFBdyVars::PH].box();
    int ilo_ph = ph_bx.smallEnd()[0];
    int ihi_ph = ph_bx.bigEnd()[0];
    int jlo_ph = ph_bx.smallEnd()[1];
    int jhi_ph = ph_bx.bigEnd()[1];

    for ( MFIter mfi(*mask_c); mfi.isValid(); ++mfi )
    {
        Box tbx = mfi.tilebox();
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);

        const Box& bx_u  = (xbx & bdy_data[itime][WRFBdyVars::U].box());
        const Box& bx_v  = (ybx & bdy_data[itime][WRFBdyVars::V].box());
        const Box& bx_t  = (tbx & bdy_data[itime][WRFBdyVars::T].box());
        const Box& bx_qv = (tbx & bdy_data[itime][WRFBdyVars::QV].box());

        // TMP BDY data
        Array4<Real> bdy_u_tmp  = bdy_data_tmp[WRFBdyVars::U].array();  // This is x-face-centered
        Array4<Real> bdy_v_tmp  = bdy_data_tmp[WRFBdyVars::V].array();  // This is y-face-centered
        Array4<Real> bdy_t_tmp  = bdy_data_tmp[WRFBdyVars::T].array();  // This is cell-centered
        Array4<Real> bdy_qv_tmp = bdy_data_tmp[WRFBdyVars::QV].array(); // This is cell-centered

        // TMP INTERP BDY data
        Array4<Real> bdy_u_int  = bdy_data_int[WRFBdyVars::U].array();  // This is x-face-centered
        Array4<Real> bdy_v_int  = bdy_data_int[WRFBdyVars::V].array();  // This is y-face-centered
        Array4<Real> bdy_t_int  = bdy_data_int[WRFBdyVars::T].array();  // This is cell-centered
        Array4<Real> bdy_qv_int = bdy_data_int[WRFBdyVars::QV].array(); // This is cell-centered

        // Mask data
        const Array4<const int>& mask_c_arr = mask_c->const_array(mfi);
        const Array4<const int>& mask_u_arr = mask_u->const_array(mfi);
        const Array4<const int>& mask_v_arr = mask_v->const_array(mfi);

        // Populated from read wrfinput
        Array4<Real const> c1h_arr = wrf_C1H->const_array(mfi);
        Array4<Real const> c2h_arr = wrf_C2H->const_array(mfi);
        Array4<Real const> mub_arr = wrf_MUB->const_array(mfi);
        Array4<Real>       PHB_arr = wrf_PHB->array(mfi);

        // New z values
        ParallelFor(bx_t, bx_u, bx_v,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Mass coupling
            Real mu    = mu_arr(i ,j ,0)  + mub_arr(i ,j ,0);
            Real mu0   = mu0_arr(i ,j ,0) + mub_arr(i ,j ,0);

            // Pert and base geopotential
            Real P     = PHB_arr(i ,j ,k  ) + bdy_ph_arr(i ,j ,k  )/mu  ;
            Real P_kp  = PHB_arr(i ,j ,k+1) + bdy_ph_arr(i ,j ,k+1)/mu  ;
            Real P0    = PHB_arr(i ,j ,k  ) + bdy_ph0_arr(i ,j ,k  )/mu0;
            Real P0_kp = PHB_arr(i ,j ,k+1) + bdy_ph0_arr(i ,j ,k+1)/mu0;

            // New heights
            bdy_c_z_src(i,j,k) = Real(0.5  ) * ( P + P_kp ) / CONST_GRAV;

            // Original heights
            bdy_c_z_dst(i,j,k) = Real(0.5  ) * ( P0 + P0_kp ) / CONST_GRAV;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Prevent averaging outside domain and match init from WRF input
            int ii = std::max(std::min(i  ,ihi_ph),ilo_ph);
            int im = std::max(std::min(i-1,ihi_ph),ilo_ph);

            // Mass coupling
            Real mu      = mu_arr(ii,j ,0)  + mub_arr(ii,j ,0);
            Real mu_im   = mu_arr(im,j ,0)  + mub_arr(im,j ,0);
            Real mu0     = mu0_arr(ii,j ,0) + mub_arr(ii,j ,0);
            Real mu0_im  = mu0_arr(im,j ,0) + mub_arr(im,j ,0);

            // Pert and base geopotential
            Real P        = PHB_arr(ii,j ,k  ) + bdy_ph_arr(ii,j ,k  )/mu     ;
            Real P_im     = PHB_arr(im,j ,k  ) + bdy_ph_arr(im,j ,k  )/mu_im  ;
            Real P_kp     = PHB_arr(ii,j ,k+1) + bdy_ph_arr(ii,j ,k+1)/mu     ;
            Real P_im_kp  = PHB_arr(im,j ,k+1) + bdy_ph_arr(im,j ,k+1)/mu_im  ;
            Real P0       = PHB_arr(ii,j ,k  ) + bdy_ph0_arr(ii,j ,k  )/mu0   ;
            Real P0_im    = PHB_arr(im,j ,k  ) + bdy_ph0_arr(im,j ,k  )/mu0_im;
            Real P0_kp    = PHB_arr(ii,j ,k+1) + bdy_ph0_arr(ii,j ,k+1)/mu0   ;
            Real P0_im_kp = PHB_arr(im,j ,k+1) + bdy_ph0_arr(im,j ,k+1)/mu0_im;

            // New heights
            bdy_u_z_src(i,j,k) = Real(0.25) * ( P + P_kp + P_im + P_im_kp ) / CONST_GRAV;

            // Original heights
            bdy_u_z_dst(i,j,k) = Real(0.25) * ( P0 + P0_kp + P0_im + P0_im_kp ) / CONST_GRAV;
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Prevent averaging outside domain and match init from WRF input
            int jj = std::max(std::min(j  ,jhi_ph),jlo_ph);
            int jm = std::max(std::min(j-1,jhi_ph),jlo_ph);

            // Mass coupling
            Real mu      = mu_arr(i ,jj,0)  + mub_arr(i ,jj,0);
            Real mu_jm   = mu_arr(i ,jm,0)  + mub_arr(i ,jm,0);
            Real mu0     = mu0_arr(i ,jj,0) + mub_arr(i ,jj,0);
            Real mu0_jm  = mu0_arr(i ,jm,0) + mub_arr(i ,jm,0);

            // Pert and base geopotential
            Real P        = PHB_arr(i ,jj,k  ) + bdy_ph_arr(i ,jj,k  )/mu     ;
            Real P_jm     = PHB_arr(i ,jm,k  ) + bdy_ph_arr(i ,jm,k  )/mu_jm  ;
            Real P_kp     = PHB_arr(i ,jj,k+1) + bdy_ph_arr(i ,jj,k+1)/mu     ;
            Real P_jm_kp  = PHB_arr(i ,jm,k+1) + bdy_ph_arr(i ,jm,k+1)/mu_jm  ;
            Real P0       = PHB_arr(i ,jj,k  ) + bdy_ph0_arr(i ,jj,k  )/mu0   ;
            Real P0_jm    = PHB_arr(i ,jm,k  ) + bdy_ph0_arr(i ,jm,k  )/mu0_jm;
            Real P0_kp    = PHB_arr(i ,jj,k+1) + bdy_ph0_arr(i ,jj,k+1)/mu0   ;
            Real P0_jm_kp = PHB_arr(i ,jm,k+1) + bdy_ph0_arr(i ,jm,k+1)/mu0_jm;

            // New heights
            bdy_v_z_src(i,j,k) = Real(0.25) * ( P + P_kp + P_jm + P_jm_kp ) / CONST_GRAV;

            // Original heights
            bdy_v_z_dst(i,j,k) = Real(0.25) * ( P0 + P0_kp + P0_jm + P0_jm_kp ) / CONST_GRAV;
        });

        // Define u velocity
        ParallelFor(bx_u, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_u_arr(i,j,k)) {
                Real xmu;
                if (i == ilo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (i > ihi) {
                    xmu  = mu_arr(i-1,j,0) + mub_arr(i-1,j,0);
                } else {
                    xmu = (  mu_arr(i,j,0) +  mu_arr(i-1,j,0)
                          + mub_arr(i,j,0) + mub_arr(i-1,j,0)) * myhalf;
                }
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy     = bdy_u_arr(i,j,k) / xmu_mult;
                bdy_u_tmp(i,j,k) = new_bdy;
            }
        });

        // Define v velocity
        ParallelFor(bx_v, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_v_arr(i,j,k)) {
                Real xmu;
                if (j == jlo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (j > jhi) {
                    xmu  = mu_arr(i,j-1,0) + mub_arr(i,j-1,0);
                } else {
                    xmu =  (  mu_arr(i,j,0) +  mu_arr(i,j-1,0)
                           + mub_arr(i,j,0) + mub_arr(i,j-1,0) ) * myhalf;
                }
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy     = bdy_v_arr(i,j,k) / xmu_mult;
                bdy_v_tmp(i,j,k) = new_bdy;
            }
        });

        // Convert perturbational moist pot. temp. (Th_m) to dry pot. temp. (Th_d)
        const Real wrf_theta_ref = Real(300.);
        ParallelFor(bx_t, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_c_arr(i,j,k)) {
                Real xmu         = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_Th  = bdy_t_arr(i,j,k) / xmu_mult + wrf_theta_ref;
                Real qv_fac      = (one + (R_v/R_d) * bdy_qv_arr(i,j,k) / xmu_mult);
                new_bdy_Th      /= qv_fac;
                bdy_t_tmp(i,j,k) = new_bdy_Th;
            }
        });

        // Define Qv
        ParallelFor(bx_qv, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_c_arr(i,j,k)) {
                Real xmu          = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult     = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_QV   = bdy_qv_arr(i,j,k) / xmu_mult;
                bdy_qv_tmp(i,j,k) = (use_moist) ? new_bdy_QV : zero;
            }
        });

        // Interpolate in height
        ParallelFor(bx_t, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_c_arr(i,j,k)) {
                int kstart, kend;
                Real z_dst, z_hi_src, z_lo_src;

                kstart   = k - amrex::min(5,k);
                z_dst    = bdy_c_z_dst(i,j,k);
                z_lo_src = bdy_c_z_src(i,j,kstart);

                bool found = false;
                for (int lk(kstart+1); lk<khi; ++lk) {
                    z_hi_src = bdy_c_z_src(i,j,lk);
                    if (z_dst >= z_lo_src && z_dst <= z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    bdy_t_int(i,j,k) = (  bdy_t_tmp(i,j,kend) -  bdy_t_tmp(i,j,kstart) ) * dz_rat +  bdy_t_tmp(i,j,kstart);
                    bdy_qv_int(i,j,k) = ( bdy_qv_tmp(i,j,kend) - bdy_qv_tmp(i,j,kstart) ) * dz_rat + bdy_qv_tmp(i,j,kstart);
                } else {
                    bdy_t_int(i,j,k)  =  bdy_t_tmp(i,j,k);
                    bdy_qv_int(i,j,k) = bdy_qv_tmp(i,j,k);
                }
            }
        });

        ParallelFor(bx_u, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_u_arr(i,j,k)) {
                int kstart, kend;
                Real z_dst, z_hi_src, z_lo_src;

                kstart   = k - amrex::min(5,k);
                z_dst    = bdy_u_z_dst(i,j,k);
                z_lo_src = bdy_u_z_src(i,j,kstart);

                bool found = false;
                for (int lk(kstart+1); lk<khi; ++lk) {
                    z_hi_src = bdy_u_z_src(i,j,lk);
                    if (z_dst >= z_lo_src && z_dst <= z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    bdy_u_int(i,j,k) = ( bdy_u_tmp(i,j,kend) -  bdy_u_tmp(i,j,kstart) ) * dz_rat +  bdy_u_tmp(i,j,kstart);
                } else {
                    bdy_u_int(i,j,k) =  bdy_u_tmp(i,j,k);
                }
            }
        });

        ParallelFor(bx_v, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_v_arr(i,j,k)) {
                int kstart, kend;
                Real z_dst, z_hi_src, z_lo_src;

                kstart   = k - amrex::min(5,k);
                z_dst    = bdy_v_z_dst(i,j,k);
                z_lo_src = bdy_v_z_src(i,j,kstart);

                bool found = false;
                for (int lk(kstart+1); lk<khi; ++lk) {
                    z_hi_src = bdy_v_z_src(i,j,lk);
                    if (z_dst >= z_lo_src && z_dst <= z_hi_src) {
                        found = true;
                        kend  = lk;
                        break;
                    }
                    z_lo_src = z_hi_src;
                    kstart   = lk;
                }

                if (found) {
                    Real dz_rat = (z_dst - z_lo_src) / (z_hi_src - z_lo_src);
                    bdy_v_int(i,j,k) = ( bdy_v_tmp(i,j,kend) -  bdy_v_tmp(i,j,kstart) ) * dz_rat +  bdy_v_tmp(i,j,kstart);
                } else {
                    bdy_v_int(i,j,k) =  bdy_v_tmp(i,j,k);
                }
            }
        });
    } // mfi

    for (int ivar(0); ivar < vsize; ++ivar) {
        amrex::ParallelAllReduce::Sum(bdy_data_int[ivar].dataPtr(),
                                      bdy_data_int[ivar].size(),
                                      ParallelContext::CommunicatorAll());
        bdy_data[itime][ivar].template  copy<RunOn::Device>(bdy_data_int[ivar],0,0,1);
    }
}

void
read_and_convert_from_wrfbdy (const int itime, const std::string& nc_bdy_file,
                              Vector<Vector<FArrayBox>>& bdy_data_xlo,
                              Vector<Vector<FArrayBox>>& bdy_data_xhi,
                              Vector<Vector<FArrayBox>>& bdy_data_ylo,
                              Vector<Vector<FArrayBox>>& bdy_data_yhi,
                              std::unique_ptr<MultiFab>& wrf_MUB,
                              std::unique_ptr<MultiFab>& wrf_C1H,
                              std::unique_ptr<MultiFab>& wrf_C2H,
                              std::unique_ptr<MultiFab>& wrf_PHB,
                              const MultiFab& xvel, const MultiFab& yvel, const MultiFab& cons,
                              const Geometry& geom,
                              const bool& use_moist,
                              int real_width, Real bdy_time_interval,
                              bool do_conversion)
{
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // If we are trying to define the bdy data at the final time, we must do it by first reading the tendency from
    // the previous time, then adding the previous time value + dT * tendency.
    bool do_tendency = (itime == bdy_data_xlo.size()-1);

    // If we are going to extrapolate from the older time we need to re-grab it since it has been over-written already
    if (do_tendency)
    {
        read_and_convert_from_wrfbdy(itime-1,nc_bdy_file,
                                     bdy_data_xlo, bdy_data_xhi, bdy_data_ylo, bdy_data_yhi,
                                     wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB,
                                     xvel, yvel, cons, geom,
                                     use_moist, real_width, bdy_time_interval,
                                     false);
    }

    // Even though we may not read in all the variables, we need to make the arrays big enough for them (for now)
    int nvars = WRFBdyVars::NumTypes*4;

    const Box& domain = geom.Domain();
    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();
    const int khi = hi[2];

    IntVect plo(lo);
    IntVect phi(hi);

    // ******************************************************************
    // Read the netcdf file and fill these FABs
    // NOTE: the order and number of these must match the WRFBdyVars enum!
    // WRFBdyVars:  U, V, T, QV, PH, MU, PC
    //
    // These fields are at myhalf levels (unstaggered)
    // ******************************************************************
    Vector<std::string> nc_var_names;
    Vector<std::string> nc_var_prefix = {"U","V","T","QVAPOR","PH","MU","PC"};

    for (int ip = 0; ip < nc_var_prefix.size(); ++ip)
    {
        if (do_tendency) {
            nc_var_names.push_back(nc_var_prefix[ip] + "_BTXS");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BTXE");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BTYS");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BTYE");
        } else {
            nc_var_names.push_back(nc_var_prefix[ip] + "_BXS");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BXE");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BYS");
            nc_var_names.push_back(nc_var_prefix[ip] + "_BYE");
        }
    }

    using RARRAY = NDArray<float>;
    Vector<RARRAY> tslice(nc_var_names.size());

    int width; // size of bdy_width from wrfbdy

    if (ParallelDescriptor::IOProcessor())
    {
        Vector<int> success(nc_var_names.size());

        int itime_to_read = (do_tendency) ? itime-1 : itime;
        ReadTimeSliceFromNetCDFFile(nc_bdy_file, itime_to_read, nc_var_names, tslice, success);

        for (auto &istat:success) {
            AMREX_ALWAYS_ASSERT(istat==1);
        }

        // Width of the boundary region
        width = tslice[0].get_vshape()[1];

        if (width != real_width) {
            AMREX_ALWAYS_ASSERT(real_width < width);
            Print() << "Note: Requested boundary width is " << real_width
                << " < " << width << " (bdy_width size in file)" << std::endl;
        }
    }
    ParallelDescriptor::Bcast(&width,1,ioproc);

    // This loops over every variable on every face, so nvars should be 4 * number of "ivartype" below
    for (int iv = 0; iv < nvars; iv++)
    {
        // Print() << "Building FAB for the NetCDF variable : " << nc_var_names[iv] << std::endl;

        int bdyVarType;

        std::string first1 = nc_var_names[iv].substr(0,1);
        std::string first2 = nc_var_names[iv].substr(0,2);

        if        (first1 == "U") {
            bdyVarType = WRFBdyVars::U;
        } else if (first1 == "V") {
            bdyVarType = WRFBdyVars::V;
        } else if (first1 == "T") {
            bdyVarType = WRFBdyVars::T;
        } else if (first2 == "QV") {
            bdyVarType = WRFBdyVars::QV;
        } else if (first2 == "PH") {
            bdyVarType = WRFBdyVars::PH;
        } else if (first2 == "MU") {
            bdyVarType = WRFBdyVars::MU;
        } else if (first2 == "PC") {
            bdyVarType = WRFBdyVars::PC;
        } else {
            Print() << "Trying to read " << first1 << " or " << first2 << std::endl;
            Abort("dont know this variable");
        }

        std::string last3 = nc_var_names[iv].substr(nc_var_names[iv].size()-3, 3);
        std::string last4 = nc_var_names[iv].substr(nc_var_names[iv].size()-4, 4);
        int bdyType;

        if        (last3 == "BXS" || last4 == "BTXS") {
            bdyType = WRFBdyTypes::x_lo;
        } else if (last3 == "BXE" || last4 == "BTXE") {
            bdyType = WRFBdyTypes::x_hi;
        } else if (last3 == "BYS" || last4 == "BTYS") {
            bdyType = WRFBdyTypes::y_lo;
        } else if (last3 == "BYE" || last4 == "BTYE") {
            bdyType = WRFBdyTypes::y_hi;
        }

        Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
        Arena_Used = The_Pinned_Arena();
#endif

        if (bdyType == WRFBdyTypes::x_lo) {

            // *******************************************************************************
            // xlo bdy
            // *******************************************************************************
            plo[0] = lo[0]             ; plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = lo[0]+real_width-1; phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xlo(plo, phi);

            Box xlo_plane_no_stag(pbx_xlo);
            Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
            Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});
            Box xlo_plane_z_stag = convert(pbx_xlo, {0, 0, 1});

            Box xlo_line(IntVect(lo[0], lo[1], 0), IntVect(lo[0]+real_width-1, hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_x_stag, 1, Arena_Used));  // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_y_stag, 1, Arena_Used));  // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::PH) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_z_stag, 1, Arena_Used));  // PH
            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_line, 1, Arena_Used));          // MU/PC
            }

        } else if (bdyType == WRFBdyTypes::x_hi) {

            // *******************************************************************************
            // xhi bdy
            // *******************************************************************************
            plo[0] = hi[0]-real_width+1; plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = hi[0]             ; phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xhi(plo, phi);

            Box xhi_plane_no_stag(pbx_xhi);
            Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
            Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});
            Box xhi_plane_z_stag = convert(pbx_xhi, {0, 0, 1});

            Box xhi_line(IntVect(hi[0]-real_width+1, lo[1], 0), IntVect(hi[0], hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_x_stag, 1, Arena_Used));  // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_y_stag, 1, Arena_Used));  // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::PH) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_z_stag, 1, Arena_Used));  // PH
            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_line, 1, Arena_Used));          // MU/PC
            }

        } else if (bdyType == WRFBdyTypes::y_lo) {

            // *******************************************************************************
            // ylo bdy
            // *******************************************************************************
            plo[1] = lo[1]             ; plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = lo[1]+real_width-1; phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_ylo(plo, phi);

            Box ylo_plane_no_stag(pbx_ylo);
            Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
            Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);
            Box ylo_plane_z_stag = convert(pbx_ylo, {0, 0, 1});

            Box ylo_line(IntVect(lo[0], lo[1], 0), IntVect(hi[0], lo[1]+real_width-1, 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_x_stag, 1, Arena_Used));  // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_y_stag, 1, Arena_Used));  // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::PH) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_z_stag, 1, Arena_Used));  // PH
            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_line, 1, Arena_Used));          // MU/PC
            }

        } else if (bdyType == WRFBdyTypes::y_hi) {

            // *******************************************************************************
            // yhi bdy
            // *******************************************************************************
            plo[1] = hi[1]-real_width+1; plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = hi[1]             ; phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_yhi(plo, phi);

            Box yhi_plane_no_stag(pbx_yhi);
            Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
            Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);
            Box yhi_plane_z_stag = convert(pbx_yhi, {0, 0, 1});

            Box yhi_line(IntVect(lo[0], hi[1]-real_width+1, 0), IntVect(hi[0], hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_x_stag, 1, Arena_Used));  // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_y_stag, 1, Arena_Used));  // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::PH) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_z_stag, 1, Arena_Used));  // PH
            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_line, 1, Arena_Used));          // MU/PC
            }
        }

        long num_pts;

        // Now fill the data
        if (ParallelDescriptor::IOProcessor())
        {
            // Print() << "SHAPE0 " << tslice[iv].get_vshape()[0] << std::endl;
            // Print() << "SHAPE1 " << tslice[iv].get_vshape()[1] << std::endl;
            // Print() << "SHAPE2 " << tslice[iv].get_vshape()[2] << std::endl;
            // Print() << "SHAPE3 " << tslice[iv].get_vshape()[3] << std::endl;

            Array4<Real> fab_arr;

            if (bdyVarType == WRFBdyVars::U || bdyVarType == WRFBdyVars::V  ||
                bdyVarType == WRFBdyVars::T || bdyVarType == WRFBdyVars::QV ||
                bdyVarType == WRFBdyVars::PH)
            {
                // xlo,xhi dims: (Time, bdy_width, bottom_top, south_north)
                // ylo,yhi dims: (Time, bdy_width, bottom_top, west_east)

                int ns2 = tslice[iv].get_vshape()[2]; // vertical size
                int ns3 = tslice[iv].get_vshape()[3]; // lateral size, may be staggered

                int lkhi = (bdyVarType == WRFBdyVars::PH) ? khi+1 : khi;

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xlo[itime][bdyVarType].smallEnd()[0];
                    fab_arr  = bdy_data_xlo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / (ns2 * ns3);
                        if (i >= real_width) continue;
                        int k = (n - i * (ns2 * ns3)) / ns3;
                        if (k > lkhi) continue;
                        int j =  n - i * (ns2 * ns3) - k * ns3;
                        fab_arr(ioff+i, j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xhi[itime][bdyVarType].bigEnd()[0];
                    fab_arr  = bdy_data_xhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / (ns2 * ns3);
                        if (i >= real_width) continue;
                        int k = (n - i * (ns2 * ns3)) / ns3;
                        if (k > lkhi) continue;
                        int j =  n - i * (ns2 * ns3) - k * ns3;
                        fab_arr(ioff-i, j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_ylo[itime][bdyVarType].smallEnd()[1];
                    fab_arr  = bdy_data_ylo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / (ns2 * ns3);
                        if (j >= real_width) continue;
                        int k = (n - j * (ns2 * ns3)) / ns3;
                        if (k > lkhi) continue;
                        int i =  n - j * (ns2 * ns3) - k * ns3;
                        fab_arr(i, joff+j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_yhi[itime][bdyVarType].bigEnd()[1];
                    fab_arr  = bdy_data_yhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / (ns2 * ns3);
                        if (j >= real_width) continue;
                        int k = (n - j * (ns2 * ns3)) / ns3;
                        if (k > lkhi) continue;
                        int i =  n - j * (ns2 * ns3) - k * ns3;
                        fab_arr(i, joff-j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } // bdyType

            } else if (bdyVarType == WRFBdyVars::MU ||
                       bdyVarType == WRFBdyVars::PC) {
                // xlo,xhi dims: (Time, bdy_width, south_north)
                // ylo,yhi dims: (Time, bdy_width, west_east)

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xlo[itime][bdyVarType].smallEnd()[0];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_xlo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / ns2;
                        if (i >= real_width) continue;
                        int j = n - i * ns2;
                        fab_arr(ioff+i, j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xhi[itime][bdyVarType].bigEnd()[0];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_xhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / ns2;
                        if (i >= real_width) continue;
                        int j = n - i * ns2;
                        fab_arr(ioff-i, j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_ylo[itime][bdyVarType].smallEnd()[1];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_ylo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / ns2;
                        if (j >= real_width) continue;
                        int i = n - j * ns2;
                        fab_arr(i, joff+j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_yhi[itime][bdyVarType].bigEnd()[1];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_yhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / ns2;
                        if (j >= real_width) continue;
                        int i = n - j * ns2;
                        fab_arr(i, joff-j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                }
            } // bdyVarType
        } // if ParalleDescriptor::IOProcessor()
    } // nc_var_names

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int n_per_time = nc_var_prefix.size();
    for (int i = 0; i < n_per_time; i++)
    {
        ParallelDescriptor::Bcast(bdy_data_xlo[itime][i].dataPtr(),bdy_data_xlo[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_xhi[itime][i].dataPtr(),bdy_data_xhi[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_ylo[itime][i].dataPtr(),bdy_data_ylo[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_yhi[itime][i].dataPtr(),bdy_data_yhi[itime][i].box().numPts(),ioproc);
    }

    if (do_tendency) {
        for (int i = 0; i < n_per_time; i++)
        {
            // Multiply the tendency bdy_tend_prev (stored in bdy_data at itime) by dt to get difference between old and new
            bdy_data_xlo[itime][i].mult(bdy_time_interval,0,1);
            bdy_data_xhi[itime][i].mult(bdy_time_interval,0,1);
            bdy_data_ylo[itime][i].mult(bdy_time_interval,0,1);
            bdy_data_yhi[itime][i].mult(bdy_time_interval,0,1);

            // Add bdy_prev to dt*bdy_tend_prev to get bdy_current
            bdy_data_xlo[itime][i].plus(bdy_data_xlo[itime-1][i], 0, 0, 1);
            bdy_data_xhi[itime][i].plus(bdy_data_xhi[itime-1][i], 0, 0, 1);
            bdy_data_ylo[itime][i].plus(bdy_data_ylo[itime-1][i], 0, 0, 1);
            bdy_data_yhi[itime][i].plus(bdy_data_yhi[itime-1][i], 0, 0, 1);
        }
    }

    if (do_conversion)
    {
        // Owner masks for parallel reduce sum
        std::unique_ptr<iMultiFab> mask_c = OwnerMask(cons, geom.periodicity());
        std::unique_ptr<iMultiFab> mask_u = OwnerMask(xvel, geom.periodicity());
        std::unique_ptr<iMultiFab> mask_v = OwnerMask(yvel, geom.periodicity());

        if (do_tendency) {
            convert_wrfbdy_data(itime-1, domain, bdy_data_xlo, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
            convert_wrfbdy_data(itime-1, domain, bdy_data_xhi, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
            convert_wrfbdy_data(itime-1, domain, bdy_data_ylo, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
            convert_wrfbdy_data(itime-1, domain, bdy_data_yhi, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
        }

        convert_wrfbdy_data(itime, domain, bdy_data_xlo, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
        convert_wrfbdy_data(itime, domain, bdy_data_xhi, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
        convert_wrfbdy_data(itime, domain, bdy_data_ylo, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
        convert_wrfbdy_data(itime, domain, bdy_data_yhi, wrf_MUB, wrf_C1H, wrf_C2H, wrf_PHB, mask_u.get(), mask_v.get(), mask_c.get(), use_moist);
    }
}
#endif // ERF_USE_NETCDF
