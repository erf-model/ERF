/**
 * \file ERF_ReadFromWRFBdy.H
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_ProbCommon.H>
#include <ERF_DataStruct.H>
#include <ERF_SurfaceLayer.H>

#ifdef ERF_USE_NETCDF

amrex::Real
read_times_from_wrfbdy (const std::string& nc_bdy_file,
                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                        amrex::Real& start_bdy_time, amrex::Real& final_bdy_time);

void
read_and_convert_from_wrfbdy (const int itime, const std::string& nc_bdy_file,
                              amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                              amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                              amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                              amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                              std::unique_ptr<amrex::MultiFab>& wrf_MUB, std::unique_ptr<amrex::MultiFab>& wrf_C1H,
                              std::unique_ptr<amrex::MultiFab>& wrf_C2H, std::unique_ptr<amrex::MultiFab>& wrf_PHB,
                              const amrex::MultiFab& xvel, const amrex::MultiFab& yvel, const amrex::MultiFab& cons,
                              const amrex::Geometry& geom, const bool& use_moist,
                              int real_width, amrex::Real bdy_time_interval,
                              bool do_conversion = true);

amrex::Real
read_times_from_wrflow (const std::string& nc_low_file,
                        amrex::Vector<amrex::Vector<amrex::FArrayBox>>& low_data_zlo,
                        amrex::Real& start_low_time, amrex::Real& final_low_time);

void
read_from_wrflow (const int itime,
                  const std::string& nc_low_file,
                  const amrex::Box& domain,
                  amrex::Vector<amrex::Vector<amrex::FArrayBox>>& low_data_zlo);

void
update_sst_tsk (const int itime,
                const amrex::Geometry& geom,
                const amrex::BoxArray& ba2d_lev,
                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& sst_lev,
                amrex::Vector<std::unique_ptr<amrex::MultiFab>>& tsk_lev,
                std::unique_ptr<SurfaceLayer>& SurfLayer,
                const amrex::Vector<amrex::Vector<amrex::FArrayBox>>& low_data_zlo,
                const amrex::MultiFab& cons,
                const amrex::MultiFab& mf_PSFC_lev,
                const amrex::Real rdOcp,
                std::unique_ptr<amrex::iMultiFab>& lmask,
                const bool use_moist);

#endif // ERF_USE_NETCDF
