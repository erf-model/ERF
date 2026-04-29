#ifndef ERF_READ_FROM_ERFBDY_H_
#define ERF_READ_FROM_ERFBDY_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Vector.H>

// Read metadata and time values from erfbdy Header
// Returns bdy_time_interval
amrex::Real
read_times_from_erfbdy(const std::string& bdy_file_name,
                       int& ntimes,
                       int& nvars,
                       int& real_width,
                       amrex::Vector<amrex::Real>& bdy_times,
                       amrex::Real& start_bdy_time,
                       amrex::Real& final_bdy_time);

// Read a single time slice from erfbdy file
void
read_from_erfbdy(int itime,
                 const std::string& bdy_file_name,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                 int nvars, int real_width);

#endif
