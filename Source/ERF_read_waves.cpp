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
         
         const auto & bx = mfi.validbox();

         amrex::Print() <<  " HERE " << bx << std::endl;
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


         int nx=2147483647; 
         int ny=2147483647; // sanity check

         //JUST RECEIVED
         if (amrex::MPMD::MyProc() == this_root) {
             if (rank_offset == 0) // First program
             {
                 MPI_Recv(&nx, 1, MPI_INT, other_root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 MPI_Recv(&ny, 1, MPI_INT, other_root, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             }
             else // Second program
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

         } 
    }


    //May need to be Redistribute
    //    ParallelCopy(Hwave[lev],Hwave_onegrid[lev],0,0,1,0);
    Hwave[lev]->ParallelCopy(*Hwave_onegrid[lev]);
    Hwave[lev]->FillBoundary(geom[lev].periodicity());
    Lwave[lev]->ParallelCopy(*Lwave_onegrid[lev]);
    Lwave[lev]->FillBoundary(geom[lev].periodicity());
    amrex::Print() << "HWAVE BOX " << (*Hwave[lev])[0].box() << std::endl;
    for (MFIter mfi(*Hwave[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        Box bx = mfi.tilebox();
        const Array4<Real const>& Hwave_arr = Hwave[lev]->const_array(mfi);
        const Array4<int>& Lmask_arr = lmask_lev[lev][0]->array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
            if (Hwave_arr(i,j,k)<0) {
                Lmask_arr(i,j,k) = 1;
             } else {
                Lmask_arr(i,j,k) = 0;
            }
        });
    }
    
}
#endif

