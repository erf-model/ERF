#ifndef READ_PLOTFILE_H_
#define READ_PLOTFILE_H_

#include <string>
#include <AMReX.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

amrex::Vector<std::string> ReadVarNames(const std::string& filename);

// ------------------------------------------------------------
// ReadPlotFile
//
// var_filename  : text file containing variable names (one per line)
// plot_filename : AMReX plotfile directory
// mf            : MultiFab that will be defined and filled
// ------------------------------------------------------------
void ReadPlotFile(const std::string& var_filename,
                  amrex::PlotFileData& pf,
                  amrex::MultiFab& mf);

void
check_large_values(const amrex::MultiFab& mf);

void
ApplyNeumannBCs(const amrex::Geometry& geom,
                     amrex::MultiFab& mf_cc);

void WriteCustomDataFile(const amrex::Geometry& geom,
                         const amrex::MultiFab& mf_cc,
                         const std::string& filename_custom);

void AverageFineToCoarse(const amrex::MultiFab& mf_fine,
                         amrex::MultiFab&       mf_coarse,
                         const amrex::Geometry& geom_fine,
                         const amrex::Geometry& geom_coarse);
#endif
