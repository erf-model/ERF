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
GetCoarseMultiFabOnFineDMap(const amrex::Geometry& geom_coarse,
                            const amrex::Geometry& geom_fine,
                            const amrex::MultiFab& multifab_coarse,
                            const amrex::MultiFab& multifab_fine,
                            amrex::MultiFab& coarse_multifab_on_fine_dmap);

void
CreateNodalMultiFabFromCellCenteredMultiFab (amrex::MultiFab& mf_nc,      // output nodal MF
                    amrex::MultiFab& mf_cc,      // input cell-centered MF (coarse)
                    const amrex::Geometry& geom);

void
CreateCellCenteredMultiFabFromNodalMultiFab(amrex::MultiFab& tmp_cc,   // cell-centered output
                                            amrex::MultiFab& mf_nc);   // node-centered input

#endif
