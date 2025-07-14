#include <iostream>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

#include "ERF.H"

#ifdef ERF_USE_WW3_COUPLING
#include <mpi.h>
#include <AMReX_MPMD.H>
#endif

#ifdef ERF_USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

std::string inputs_name;

using namespace amrex;

#ifdef ERF_USE_KOKKOS
static void test_kokkos (MultiFab& mfa, MultiFab const& mfb)
{
    for (MFIter mfi(mfa); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto N0 = bx.length(0);
        auto N1 = bx.length(1);
        auto N2 = bx.length(2);
        auto* pa = mfa[mfi].dataPtr();
        auto const* pb = mfb[mfi].dataPtr();
        Kokkos::View<double***,Kokkos::LayoutLeft> a(pa,N0,N1,N2);
        Kokkos::View<double const***, Kokkos::LayoutLeft> b(pb,N0,N1,N2);
        Kokkos::parallel_for(Kokkos::MDRangePolicy({0,0,0},{N0,N1,N2}),
                             KOKKOS_LAMBDA (int i, int j, int k)
        {
            a(i,j,k) += 0.5*b(i,j,k);
        });
    }
}
static void test_amrex (MultiFab& mfa, MultiFab const& mfb)
{
    for (MFIter mfi(mfa); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        Array4<Real const> const& b = mfb.const_array(mfi);
        Array4<Real> const& a = mfa.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i,j,k) += 0.5*b(i,j,k);
        });
    }
}
#endif

/**
 * Function to set the refine_grid_layout flags to (1,1,0) by default
 * since the ERF default is different from the amrex default (1,1,1)
 * Also set max_grid_size to very large since the only reason for
 * chopping grids is if Nprocs > Ngrids
*/
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

   int max_grid_size = 2048;
   pp.queryAdd("max_grid_size",max_grid_size);

   // This will set the default value of blocking_factor to be 1, but will allow
   //     the user to override it in the inputs file or on command line
   int blocking_factor = 1;
   pp.queryAdd("blocking_factor",blocking_factor);

   int n_error_buf = 0;
   pp.queryAdd("n_error_buf",n_error_buf);
}

/**
 * Main driver -- creates the ERF object, calls ERF.InitData() and ERF.Evolve()
*/
int main (int argc, char* argv[])
{

#if defined(AMREX_MPI_THREAD_MULTIPLE)
    int requested = MPI_THREAD_MULTIPLE;
    int provided = -1;
    MPI_Init_thread(&argc, &argv, requested, &provided);
#elif defined(AMREX_USE_MPI)
    MPI_Init(&argc, &argv);
#endif

    if (argc < 2) {
        // Print usage and exit with error code if no input file was provided.
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD, "No input file provided. Exiting!!");
        return 1;
    }

    // Look for "-h" or "--help" flag and print usage
    for (auto i = 1; i < argc; i++) {
        const std::string param(argv[i]);
        if ((param == "--help") || (param == "-h") || (param == "--usage")) {
            ERF::print_banner(MPI_COMM_WORLD, std::cout);
            ERF::print_usage(MPI_COMM_WORLD, std::cout);
            return 0;
        }
    }

    if (!amrex::FileSystem::Exists(std::string(argv[1]))) {
        // Print usage and exit with error code if we cannot find the input file
        ERF::print_usage(MPI_COMM_WORLD, std::cout);
        ERF::print_error(
            MPI_COMM_WORLD, "Input file does not exist = " +
                                std::string(argv[1]) + ". Exiting!!");
        return 1;
    }

  //  print_banner(MPI_COMM_WORLD, std::cout);
    // Check to see if the command line contains --describe
    if (argc >= 2) {
        for (auto i = 1; i < argc; i++) {
            if (std::string(argv[i]) == "--describe") {
                ERF::writeBuildInfo(std::cout);
                return 0;
            }
        }
    }
#ifdef ERF_USE_WW3_COUPLING
    MPI_Comm comm = amrex::MPMD::Initialize(argc, argv);
    amrex::Initialize(argc,argv,true,comm,add_par);
#else
    amrex::Initialize(argc,argv,true,MPI_COMM_WORLD,add_par);
#endif

    // Save the inputs file name for later.
    if (!strchr(argv[1], '=')) {
      inputs_name = argv[1];
    }

    // timer for profiling
    BL_PROFILE_VAR("main()", pmain);

#ifdef ERF_USE_KOKKOS
    Kokkos::initialize(argc, argv);
    {
        BoxArray ba(Box(IntVect(0),IntVect(257)));
        DistributionMapping dm(ba);
        MultiFab mfa(ba, dm, 1, 0);
        MultiFab mfb(ba, dm, 1, 0);
        mfa.setVal(1.0);
        mfb.setVal(2.0);
        double tamrex = 1.e10;
        double tkokkos = 1.e10;
        for (int count = 0; count < 1; ++count) {
            amrex::Gpu::synchronize();
            double t0 = amrex::second();

            test_amrex(mfa, mfb);

            amrex::Gpu::synchronize();
            double t1 = amrex::second();

            test_kokkos(mfa, mfb);

            amrex::Gpu::synchronize();
            double t2 = amrex::second();

            tamrex = std::min(t1-t0, tamrex);
            tkokkos = std::min(t2-t1, tkokkos);
        }
        std::cout << "amrex  run time is " << tamrex << "\n"
                  << "kokkos run time is " << tkokkos << "\n";
    }
    Kokkos::finalize();
#endif

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
#ifdef ERF_USE_WW3_COUPLING
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    amrex::Finalize();
#ifdef AMREX_USE_MPI
#ifdef ERF_USE_WW3_COUPLING
    amrex::MPMD::Finalize();
#else
    MPI_Finalize();
#endif
#endif
}
