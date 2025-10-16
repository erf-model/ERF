//
// Read in an amrex plotfile and write out in NetCDF format
//

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>
//#include <AMReX_BoxArray.H>
//#include <AMReX_FArrayBox.H>
//#include <AMReX_ParallelDescriptor.H>
//#include <AMReX_VisMF.H>

//#include <iomanip>
//#include <iostream>
//#include <string>
//#include <ctime>

#include "ERF_NCInterface.H"

using namespace amrex;

void
writeNCPlotFile (int lev, int which, const std::string& dir,
                 const amrex::Vector<const amrex::MultiFab*> &mf,
                 const amrex::Vector<std::string> &plot_var_names,
                 const amrex::Vector<int>& level_steps,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> prob_lo,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> prob_hi,
                 amrex::Array<amrex::Real,AMREX_SPACEDIM> dx,
                 const amrex::Box& bounding_region,
                 amrex::Real time, amrex::Real start_bdy_time);

static
void
PrintUsage ()
{
    amrex::Print()
        << "\n"
        << " Convert a multilevel AMReX-formatted plotfile into one netcdf file per level \n"
        << "\n"
        << " usage:\n"
        << "   main*exe [-v] plotfilename\n"
        << "\n"
        << " optional arguments:\n"
        << "    -v verbosity            : verbose if set, otherwise silent \n";
  exit(1);
}

void
main_main()
{
    const int narg = amrex::command_argument_count();

    if (narg == 0) {
      PrintUsage();
    }

    bool verbose = false;

    std::string iFile;
    int farg = 1;
    while (farg <= narg) {
        const std::string fname = amrex::get_command_argument(farg);
        if (fname == "-h" || fname == "--help") {
          PrintUsage();
        } else if (fname == "-v") {
          verbose = true;
        } else {
            iFile = fname;
        }
        ++farg;
    }

    if (iFile.empty()) {
        amrex::Print() << "No plotfilename specified " << std::endl;
        return;
    } else {
        // Remove trailing backslash if present
        if (iFile.back() == '/') {
            iFile.pop_back();
        }
        if (verbose) {
           amrex::Print() << "Reading " << iFile << std::endl;
        }
    }

    PlotFileData pf_data(iFile);

    int finest_level = pf_data.finestLevel();

    int ncomp = pf_data.nComp();

    if (verbose) {
        amrex::Print() << "Finished defining pf with finest level " << finest_level << " and ncomp " << ncomp << std::endl;
    }

    const Vector<std::string>& varnames = pf_data.varNames();

    Vector<int> istep{3};
    Real time = 0.;

    int max_grid_size = 64;

    Vector<MultiFab> mfvec(finest_level+1);
    Vector<Geometry> geom(finest_level+1);

    Real start_bdy_time = time;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        mfvec[lev].define(pf_data.boxArray(lev), pf_data.DistributionMap(lev), ncomp, 0);

        for (int icomp = 0; icomp < ncomp; icomp++) {
            MultiFab tmp_data = pf_data.get(lev, varnames[icomp]);
            MultiFab::Copy(mfvec[lev],tmp_data,0,icomp,1,0);
        }

        // We assume only one "subdomain" at each level that holds all the grids
        BoxArray ba(pf_data.boxArray(lev));
        Box bounding_region = (lev == 0) ? pf_data.probDomain(lev) : ba.minimalBox();

        // We assume only one "subdomain" at each level
        int which = 0;

        writeNCPlotFile(lev, which, iFile, GetVecOfConstPtrs(mfvec), varnames, istep,
                        pf_data.probLo(), pf_data.probHi(), pf_data.cellSize(lev), bounding_region, time, start_bdy_time);
    }
}

int
main (int   argc,
      char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    main_main();
    amrex::Finalize();
    return 1;
}
