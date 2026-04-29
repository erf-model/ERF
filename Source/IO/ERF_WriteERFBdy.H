#ifndef ERF_WRITE_ERFBDY_H_
#define ERF_WRITE_ERFBDY_H_

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Vector.H>

// Initialize erfbdy file and write header (call once at start)
void InitERFBdyFile(const std::string& bdy_file_name,
                    int ntimes,
                    const amrex::Vector<amrex::Real>& bdy_times,
                    const amrex::Box& domain,
                    int nvars,
                    int real_width);

// Write a single time slice (call once per time)
void WriteERFBdyTimeSlice(const std::string& bdy_file_name,
                          int itime,
                          const amrex::Vector<amrex::FArrayBox>& bdy_data_xlo,
                          const amrex::Vector<amrex::FArrayBox>& bdy_data_xhi,
                          const amrex::Vector<amrex::FArrayBox>& bdy_data_ylo,
                          const amrex::Vector<amrex::FArrayBox>& bdy_data_yhi,
                          int nvars);

#endif
