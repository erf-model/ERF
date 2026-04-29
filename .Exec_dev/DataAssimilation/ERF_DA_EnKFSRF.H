#ifndef ERF_DA_EnKFSRF_H
#define ERF_DA_EnKFSRF_H

#include <ERF.H>
#include <ERF_DataStruct.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
void
read_plot_file(amrex::PlotFileData& pf,
               const std::vector<std::string> varnames,
               amrex::MultiFab& mf);

amrex::MultiFab
read_member_multifab(int n,
                     const std::string& pf_name,
                     const amrex::Vector<std::string>& varnames);

amrex::MultiFab
compute_ensemble_mean(const std::string& pf_name,
                      int Nens,
                      const amrex::Vector<std::string>& varnames);
#endif
