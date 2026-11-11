/**
 * \file ERF_NCPlotFile.H
 */
#ifndef ERF_NCPLOTFILE_H_
#define ERF_NCPLOTFILE_H_

#include <string>

#include <AMReX_Array.H>
#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include "ERF_DataStruct.H"

// Narrow declaration home for the ERF NetCDF plotfile writer used by the 2D
// plotfile assembly path. Keep this header small so the IO split remains easy
// to understand when NetCDF support is enabled or disabled.

#ifdef ERF_USE_NETCDF
/**
 * Git hashes of the code that produced the data, written as NetCDF global
 * attributes. Empty entries are skipped, so a build without Noah-MP simply
 * omits that attribute.
 *
 * The callers supply these rather than the writer reading build info itself:
 * PlotfileToNetCDF converts someone else's plotfile, so stamping the
 * converter's own build would misattribute the provenance. That tool leaves
 * all three empty.
 */
struct NCProvenance {
    std::string erf_git_hash;
    std::string amrex_git_hash;
    std::string noahmp_git_hash;
};

void
writeNCPlotFile (int lev, int which_subdomain, const std::string& dir,
                 const amrex::Vector<const amrex::MultiFab*> &plotMF,
                 const amrex::Vector<std::string> &plot_var_names,
                 const amrex::Vector<int>& level_steps,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> prob_lo,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> prob_hi,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> dx_in,
                 const amrex::Box& subdomain,
                 const double& time,
                 const double& start_bdy_time,
                 const SolverChoice& solverChoice,
                 const amrex::Vector<amrex::Real>& zlevels_stag,
                 const NCProvenance& provenance = {});
#endif

#endif
