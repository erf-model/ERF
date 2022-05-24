#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include <ERF.H>

std::string inputs_name = "";

using namespace amrex;

// Set the refine_grid_layout flags to (1,1,0) by default
// since the ERF default is different from the amrex default (1,1,1)
// Also set max_grid_size to very large since the only reason for
// chopping grids is if Nprocs > Ngrids
void add_par () {
   ParmParse pp("amr");

   // Set the refine_grid_layout flags to (1,1,0) by default
   pp.add("refine_grid_layout_x",1);
   pp.add("refine_grid_layout_y",1);
   pp.add("refine_grid_layout_z",0);

   // n_proper is the minimum number of coarse cells between coarse-fine boundaries
   // between levels (ell and ell+1) and levels (ell-1 and ell).   We want this to be
   // greater than or equal to the stencil width (a function of spatial order) divided by
   // ref_ratio (which can be 2,3 or 4).  This ensures that fillpatch at level (ell)
   // does not need to reach beyond level (ell-1). Here to be conservative we set this to 2
   // (rather than the amrex default of 1).
   pp.add("n_proper",2);

   pp.add("max_grid_size",2048);
   pp.add("blocking_factor",1);
   pp.add("n_error_buf",0);
}

int main(int argc, char* argv[])
{
    // Check to see if the command line contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                ERF::writeBuildInfo(std::cout);
                return 0;
            }
        }
    }
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);

    // Save the inputs file name for later.
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

    // wallclock time
    const Real strt_total = amrex::second();

    {
        // constructor - reads in parameters from inputs file
        //             - sizes multilevel arrays and data structures
        ERF erf;

        // initialize AMR data
        erf.InitData();

        // advance solution to final time
        erf.Evolve();

        // wallclock time
        Real end_total = amrex::second() - strt_total;

        // print wallclock time
        ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
        if (erf.Verbose()) {
            amrex::Print() << "\nTotal Time: " << end_total << '\n';
        }
    }

    // destroy timer for profiling
    BL_PROFILE_VAR_STOP(pmain);

    amrex::Finalize();
}
