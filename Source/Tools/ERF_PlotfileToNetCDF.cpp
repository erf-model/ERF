//
// Read in an amrex plotfile and write out in NetCDF format
//

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

#include <AMReX_PlotFileUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include "AMReX_Utility.H"
#include <AMReX_VisMF.H>

#include <iomanip>
#include <iostream>
#include <string>
#include <ctime>

#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "ERF_NCInterface.H"

using namespace amrex;

void
writeNCPlotFile (int lev, int which, const std::string& dir,
                 const amrex::Vector<const amrex::MultiFab*> &mf,
                 const amrex::Vector<std::string> &plot_var_names,
                 const amrex::Vector<int>& level_steps,
                 const amrex::Vector<Geometry>& geom,
                 amrex::Vector<amrex::Vector<amrex::Box>> boxes_at_level,
                 amrex::Real time, amrex::Real start_bdy_time);

static
void
PrintUsage (char* progName)
{
  std::cout << "\nUsage:\n"
            << progName
            << "\n\tinfile = inputFileName"
            << "\n\t[-help]"
            << "\n\n";
  exit(1);
}

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);

  bool verbose = true;

  {
    if (argc == 1)
      PrintUsage(argv[0]);

    if (ParallelDescriptor::NProcs() > 1) {
      amrex::Error("This is an inherently serial program");
    }

    ParmParse pp;

    const std::string farg = amrex::get_command_argument(1);
    if (farg == "-h" || farg == "--help") {
      PrintUsage(argv[0]);
    }

    //
    // Scan the arguments.
    //
    std::string iFile;
    pp.query("infile", iFile);
    if (iFile.empty()) {
      amrex::Abort("You must specify `infile'");
    }

    PlotFileData pf_data(iFile);

    int finest_level = pf_data.finestLevel();

    int ncomp = pf_data.nComp();

    const Vector<std::string>& varnames = pf_data.varNames();

    int l_which = 0;
    Vector<int> istep{3};
    Real time = 0.;

    int max_grid_size = 64;

    Vector<MultiFab> mfvec(finest_level+1);
    Vector<Geometry> geom(finest_level+1);
    Vector<Vector<Box>> boxes_at_level(finest_level+1);

    Real start_bdy_time = time;

    for (int lev = 0; lev < finest_level; lev++)
    {
        mfvec[lev].define(pf_data.boxArray(lev), pf_data.DistributionMap(lev), ncomp, 0);

        for (int icomp = 0; icomp < ncomp; icomp++) {
            MultiFab tmp_data = pf_data.get(lev, varnames[icomp]);
            MultiFab::Copy(mfvec[lev],tmp_data,0,icomp,1,0);
        }

        std::string outfile = iFile+ "d01.nc";
        if (verbose) {
            std::cout << "Writing " << outfile << std::endl;
        }

       writeNCPlotFile(lev, l_which, outfile, GetVecOfConstPtrs(mfvec), varnames, istep, geom, boxes_at_level, time, start_bdy_time);
    }
  }

  amrex::Finalize();
}
