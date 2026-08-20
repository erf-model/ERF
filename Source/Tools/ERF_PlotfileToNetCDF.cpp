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

// Do not redeclare writeNCPlotFile here -- pull in the real declaration so the
// tool cannot drift out of sync with the definition in Source/IO/ERF_NCPlotFile.cpp
#include "ERF_NCPlotFile.H"

using namespace amrex;

/**
 * Print the usage message and exit.
 */
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

/**
 * Convert a multilevel AMReX plotfile into NetCDF format.
 *
 * @return 0 on success, otherwise an error code.
 */
int
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
        return 1;
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
    double time = pf_data.time();

    int max_grid_size = 64;

    Vector<MultiFab> mfvec(finest_level+1);
    Vector<Geometry> geom(finest_level+1);

    double start_bdy_time = time;

    // A plotfile does not carry the solver configuration, and the only thing
    // writeNCPlotFile takes from it is the mesh type.  A plotfile also does not
    // record the staggered z levels needed for a stretched mesh, so we can only
    // write the (x,y,z) grid for a constant-dz mesh -- which is the default.
    SolverChoice solverChoice;
    Vector<Real> zlevels_stag;

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
                        pf_data.probLo(), pf_data.probHi(), pf_data.cellSize(lev), bounding_region,
                        time, start_bdy_time, solverChoice, zlevels_stag);
    }

    return 0;
}

/**
 * Application entry point.
 *
 * @param[in] argc Argument count.
 * @param[in] argv Argument vector.
 * @return Execution status.
 */
int
main (int   argc,
      char* argv[])
{
    amrex::SetVerbose(0);
    amrex::Initialize(argc, argv, false);
    int status = main_main();
    amrex::Finalize();
    return status;
}
