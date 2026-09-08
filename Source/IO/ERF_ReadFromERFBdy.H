/**
 * \file ERF_ReadFromERFBdy.H
 */
#ifndef ERF_READ_FROM_ERFBDY_H_
#define ERF_READ_FROM_ERFBDY_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Vector.H>

// Read metadata and time values from erfbdy Header.
// Returns bdy_time_interval.
double
read_times_from_erfbdy(const std::string& bdy_file_name,
                       int& ntimes,
                       int& nvars,
                       int& real_width,
                       amrex::Vector<double>& bdy_times,
                       double& start_bdy_time,
                       double& final_bdy_time);

// Read a single time index from erfbdy file.
void
read_from_erfbdy(int itime,
                 const std::string& bdy_file_name,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                 amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi,
                 int nvars, int real_width);

#endif
