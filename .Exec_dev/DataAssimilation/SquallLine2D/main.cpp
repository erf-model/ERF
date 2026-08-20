#include <iostream>
#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include "ERF.H"
#include <filesystem>

std::string inputs_name;

using namespace amrex;
namespace fs = std::filesystem;

/**
 * Function to set the refine_grid_layout flags to (1,1,0) by default
 * since the ERF default is different from the amrex default (1,1,1)
 * Also set max_grid_size to very large since the only reason for
 * chopping grids is if Nprocs > Ngrids
*/
void add_par () {
   ParmParse pp("amr");

   pp.add("refine_grid_layout_x",1);
   pp.add("refine_grid_layout_y",1);
   pp.add("refine_grid_layout_z",0);

   pp.add("n_proper",2);

   int max_grid_size = 2048;
   pp.queryAdd("max_grid_size",max_grid_size);

   int blocking_factor = 1;
   pp.queryAdd("blocking_factor",blocking_factor);

   int n_error_buf = 0;
   pp.queryAdd("n_error_buf",n_error_buf);
}

/**
 * Main driver -- ensemble wrapper around ERF
*/
int main (int argc, char* argv[])
{

#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    if (argc < 2) {
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        const std::string param(argv[i]);
        if ((param == "--help") || (param == "-h") || (param == "--usage")) {
            ERF::print_banner(MPI_COMM_WORLD, std::cout);
            ERF::print_usage(MPI_COMM_WORLD, std::cout);
            return 0;
        }
    }

    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD,
            "Input file does not exist = " + std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "--describe") {
            ERF::writeBuildInfo(std::cout);
            return 0;
        }
    }

    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);

    if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
    }

    BL_PROFILE_VAR("main()", pmain);

    const Real strt_total = amrex::second();

    // ------------------------------------------------------------
    // Ensemble control
    // ------------------------------------------------------------

    ParmParse pp_erf("erf");
    bool is_init_for_ensemble = false;
    pp_erf.query("is_init_for_ensemble", is_init_for_ensemble);

    ParmParse pp_ens("ensemble");

    int n_ens = 1;
    pp_ens.query("n_members", n_ens);

    // Ensemble run loop
    for (int ens_no = 0; ens_no < n_ens; ++ens_no)
    {
        // --------------------------------------------------------
        // Fresh ERF instance per ensemble
        // --------------------------------------------------------
        ERF erf;

        if(is_init_for_ensemble) {
            erf.SetDirsForPlotfilesAndCheckpointsForDA(ens_no);
        }
        erf.InitData();
        erf.Evolve();

        Real end_total = amrex::second() - strt_total;
        ParallelDescriptor::ReduceRealMax(
            end_total, ParallelDescriptor::IOProcessorNumber());

        if (erf.Verbose()) {
            amrex::Print() << "Ensemble " << ens_no
                           << " wallclock time: "
                           << end_total << '\n';
        }

        // --------------------------------------------------------
        // MPI barrier to ensure all ranks finish Evolve
        // --------------------------------------------------------
        ParallelDescriptor::Barrier();

   } // Ensemble run loop complete

   ERF tmp_erf;
   // This is only a post-processing step for visualization
   //tmp_erf.ComputeAndWriteEnsemblePerturbations();

   // Perform data assimilation
   int da_iter = 0;
   tmp_erf.PerformDataAssimilation(da_iter);

   BL_PROFILE_VAR_STOP(pmain);

   amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
