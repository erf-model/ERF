#ifdef ERF_USE_WW3_COUPLING

#include <ERF.H>
#include <Utils.H>
#include <mpi.h>
#include <AMReX_MPMD.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

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

             // amrex::AllPrintToFile("output_HS_cpp.txt")<<FArrayBox(my_H_arr)<<std::endl;
             // amrex::AllPrintToFile("output_L_cpp.txt")<<FArrayBox(my_L_arr)<<std::endl;

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

void
ERF::send_waves (int lev)
{
    int ncomp = 1; // number components
    auto& lev_new = vars_new[lev];
    const double PI = 3.1415926535897932384626433832795028841971693993751058209;

    // Access xvel, yvel from ABL
    MultiFab xvel_data(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_data(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());

    // Make local copy of xvel, yvel
    MultiFab::Copy (xvel_data, lev_new[Vars::xvel], 0, 0, 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab::Copy (yvel_data, lev_new[Vars::yvel], 0, 0, 1, lev_new[Vars::yvel].nGrowVect());


    MultiFab x_avg(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       ncomp, lev_new[Vars::cons].nGrow());
    MultiFab y_avg(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       ncomp, lev_new[Vars::cons].nGrow());

    MultiFab u_mag(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       ncomp, lev_new[Vars::cons].nGrow());

    MultiFab u_dir(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       ncomp, lev_new[Vars::cons].nGrow());
    x_avg.setVal(0.);
    y_avg.setVal(0.);
    u_mag.setVal(0.);
    u_dir.setVal(0.);

    for (MFIter mfi(x_avg,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        const Array4<Real>& u_vel = x_avg.array(mfi);
        const Array4<const Real>& velx_arr = xvel_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){
            u_vel(i,j,k)  = 0.5 *( velx_arr(i,j,k) + velx_arr(i+1,j,k) );

            // amrex::AllPrintToFile("uvel.txt") << amrex::IntVect(i,j,k) << " [" <<velx_arr(i,j,k) << "| avg:  " << u_vel(i,j,k)<< " | " << velx_arr(i+1,j,k) << "]" <<std::endl;
        });
    }

    for (MFIter mfi(y_avg,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        const Array4<Real>& v_vel = y_avg.array(mfi);
        const Array4<const Real>& vely_arr = yvel_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            v_vel(i,j,k)  = 0.5 *( vely_arr(i,j,k) + vely_arr(i,j+1,k) );

            // amrex::AllPrintToFile("vvel.txt") << amrex::IntVect(i,j,k) << " ["<<vely_arr(i,j,k)<<"| avg: " << v_vel(i,j,k)<< " | " << vely_arr(i,j+1,k) << "]" <<std::endl;
        });
    }

    for (MFIter mfi(u_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        const Array4<Real>& magnitude = u_mag.array(mfi);
        const Array4<const Real>& u = x_avg.array(mfi);
        const Array4<const Real>& v = y_avg.array(mfi);
        const Array4<Real>& theta = u_dir.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            magnitude(i,j,k)  = std::sqrt( pow(u(i,j,k), 2) + pow(v(i,j,k), 2) );

            double u_val = u(i, j, k);
            double v_val = v(i, j, k);
            if ( u_val == 0 ) {
                u_val = std::max( u_val, 1e-15 );  // Ensure u_val is non-zero
            }


            if ( u_val < 0 && v_val > 0 || u_val < 0 && v_val < 0 ) {

                theta(i,j,k) = PI + ( atan( v_val / u_val ) );

            } else {

                theta(i,j,k) = atan ( v_val / u_val );
            }


            // amrex::AllPrintToFile("mag_theta.txt") << amrex::IntVect(i,j,k) <<  " Magnitude: " << magnitude(i,j,k) << " Theta: " << theta(i,j,k) <<std::endl;
        });


            // MPI_Send to WW3
    // Calculate the number of elements in the current box
    int n_elements = bx.numPts();

    // Get the data pointers for the arrays
    Real* magnitude_ptr = &magnitude(bx.smallEnd());
    Real* theta_ptr = &theta(bx.smallEnd());

    // Initialize other_root as needed
    int other_root = 0; // Example initialization, replace with appropriate logic

    // Send magnitude array
//    MPI_Send(magnitude_ptr, n_elements, MPI_DOUBLE, other_root, 0, MPI_COMM_WORLD);
    // Send theta array
//    MPI_Send(theta_ptr, n_elements, MPI_DOUBLE, other_root, 1, MPI_COMM_WORLD);
// amrex::AllPrintToFile("debug_send.txt") << n_elements << std::endl;

}
/*
            for (MFIter mfi(u_mag); mfi.isValid(); ++mfi) {

                const Array4<Real>& magnitude = u_mag.array(mfi);
                const Array4<Real>& theta = u_dir.array(mfi);

                Real* magnitude_ptr = magnitude.dataPtr();
                Real* theta_ptr = theta.dataPtr();

                // Get number of elements in arrays
                int n_elements = mfi.validbox().numPts();
                int this_root = 0;
                int other_root = 1;

                // Send magnitude and theta
                MPI_Send(magnitude_ptr. n_elements, MPI_DOUBLE, this_root, 0, MPI_COMM_WORLD)
                MPI_Send(theta_ptr, n_elements, MPI_DOUBLE, other_root, 1, MPI_COMM_WORLD);
            }
*/
}
#endif

