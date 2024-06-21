#ifdef ERF_USE_WW3_COUPLING

#include <ERF.H>
#include <Utils.H>
#include <mpi.h>
#include <AMReX_MPMD.H>

using namespace amrex;

void
ERF::read_waves (int lev)
{
    for ( MFIter mfi(*Hwave_onegrid[lev],false); mfi.isValid(); ++mfi)
    {
         //auto my_H_ptr = Hwave[lev]->array(mfi);
         //auto my_L_ptr = Lwave[lev]->array(mfi);
         const auto & bx = mfi.validbox();

        amrex::Print() <<  " HERE " << bx << std::endl;
         // How to declare my_H_ptr directly?
         amrex::Array4<Real> my_H_arr = Hwave_onegrid[lev]->array(mfi);
         amrex::Array4<Real> my_L_arr = Lwave_onegrid[lev]->array(mfi);

         Real* my_H_ptr = my_H_arr.dataPtr();
         Real* my_L_ptr = my_L_arr.dataPtr();

         int rank_offset = amrex::MPMD::MyProc() - amrex::ParallelDescriptor::MyProc();
         int this_root, other_root;
         if (rank_offset == 0) { // First program
             this_root = 0;
             other_root = amrex::ParallelDescriptor::NProcs();
         } else {
             this_root = rank_offset;
             other_root = 0;
         }


         int nx=2147483647; // sanity check
         int ny=2147483647; // sanity check

         //JUST RECV
         if (amrex::MPMD::MyProc() == this_root) {
             if (rank_offset == 0) // the first program
             {
                 MPI_Recv(&nx, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 MPI_Recv(&ny, 1, MPI_INT, other_root, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             }
             else // the second program
             {
                 MPI_Recv(&nx, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 MPI_Recv(&ny, 1, MPI_INT, other_root, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             }
             //This may not be necessary
             ParallelDescriptor::Bcast(&nx, 1);
             ParallelDescriptor::Bcast(&ny, 1);
         }
         if((nx)*(ny) > 0) {
             int nsealm = (nx)*ny;
             Print()<<nsealm<<std::endl;
             Print()<<" NX = " << nx<<std::endl;
             Print()<<" Ny = " << ny<<std::endl;
             Print()<<" BX  NPTS = " << bx.numPts() <<std::endl;
             Print()<<" BX       = " << bx <<std::endl;
             // Print()<<" FAB NPTS = " << Hwave_onegrid[0]->array(mfi).size() <<std::endl;
             Print()<<" FAB NPTS = " << (*Hwave_onegrid[0])[0].box().numPts()<<std::endl;
             Print()<<" FAB      = " << (*Hwave_onegrid[0])[0].box() <<std::endl;
             // AMREX_ALWAYS_ASSERT_WITH_MESSAGE((nx)*ny <= bx.numPts() || bx.numPts() == 0, "total number of points being filled exceeds the size of the current box\n");

             if (amrex::MPMD::MyProc() == this_root) {
                 if (rank_offset == 0) // the first program
                 {
                     MPI_Recv(my_H_ptr, nsealm, MPI_DOUBLE, other_root, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                     MPI_Recv(my_L_ptr, nsealm, MPI_DOUBLE, other_root, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 }
                 else // the second program
                 {
                     MPI_Recv(my_H_ptr, nsealm, MPI_DOUBLE, other_root, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                     MPI_Recv(my_L_ptr, nsealm, MPI_DOUBLE, other_root, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 }
             }
             amrex::Print()<<"Just recieved "<<nsealm<<"as a double*"<<std::endl;

             if(bx.contains(IntVect(192,2,0))) {
                 std::cout<<my_H_arr(192,92,0)<<std::endl;
                 std::cout<<my_L_arr(192,92,0)<<std::endl;
                 std::cout<<my_H_ptr[192-0+(92-0)*193]<<std::endl;
                 std::cout<<my_L_ptr[192-0+(92-0)*193]<<std::endl;
             }
             amrex::AllPrintToFile("output_HS_cpp.txt")<<FArrayBox(my_H_arr)<<std::endl;
             amrex::AllPrintToFile("output_L_cpp.txt")<<FArrayBox(my_L_arr)<<std::endl;
         } else {
             finished_wave = true;
         }
    }
    //May need to be Redistribute
    //    ParallelCopy(Hwave[lev],Hwave_onegrid[lev],0,0,1,0);
    Hwave[lev]->ParallelCopy(*Hwave_onegrid[lev]);
    Hwave[lev]->FillBoundary(geom[lev].periodicity());
    Lwave[lev]->ParallelCopy(*Lwave_onegrid[lev]);
    Lwave[lev]->FillBoundary(geom[lev].periodicity());
    amrex::Print() << "HWAVE BOX " << (*Hwave[lev])[0].box() << std::endl;

    print_state(*Hwave[lev],IntVect(103,-3,0),0,IntVect(3,3,0));
    print_state(*Hwave[lev],IntVect(103,88,0),0,IntVect(3,3,0));
}
#endif

