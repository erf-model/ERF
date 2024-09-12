#ifdef ERF_USE_WW3_COUPLING

#include <ERF.H>
#include <ERF_Utils.H>
#include <mpi.h>
#include <AMReX_MPMD.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_Box.H>
#include <cmath>
#include <time.h>
#include <iostream>

using namespace amrex;


void
ERF::read_waves (int lev)
{
         double clkStart, timedif;
         clkStart = (double) clock() / CLOCKS_PER_SEC;

    for ( MFIter mfi(*Hwave_onegrid[lev],false); mfi.isValid(); ++mfi)
    {

         const auto & bx = mfi.validbox();

         amrex::Print() <<  " Just called ERF::read_waves to receive from WW3 " << bx << std::endl;
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

             amrex::AllPrintToFile("output_HS_cpp.txt")<<FArrayBox(my_H_arr)<<std::endl;
             amrex::AllPrintToFile("output_L_cpp.txt")<<FArrayBox(my_L_arr)<<std::endl;

         }
    }


         amrex::Print() <<  " Just called received HS and LM from WW3 "  << std::endl;
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
         //    amrex::Real myclock = ParallelDescriptor::second();

        //     amrex::AllPrintToFile("timer.txt") << "At " << myclock << " seconds, I reached the end of read_waves" << std::endl;

    timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - clkStart;

     amrex::AllPrintToFile("timer.txt") << "It took " << timedif << " seconds to reach the end of read_waves" << std::endl;
}

void
ERF::send_to_ww3 (int lev)
{
    int ncomp = 1; // number components
    auto& lev_new = vars_new[lev];
    const double PI = 3.1415926535897932384626433832795028841971693993751058209;

    int count_send = 0;
    int k_ref = 0;
    int nlevs_max = max_level + 1;
    const auto dz  = geom[lev].CellSize(2); //For 10m
    if (dz < 10){
        k_ref = std::floor( (10 / dz) - 0.5 );
    }
    double clkStart, timedif;
    clkStart = (double) clock() / CLOCKS_PER_SEC;

    // Access xvel, yvel from ABL
    MultiFab xvel_data(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());

    MultiFab yvel_data(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());

/*
    BoxArray ba_onegrid(geom[lev].Domain());
    BoxList bl2d = ba.boxList();
    Real* theta_ptr = &theta(domlo);
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));

    // create a new BoxArray and DistributionMapping for a MultiFab with 1 box
    BoxArray ba_onegrid(lev_new[Vars::cons].boxArray());
    BoxList bl2d_onegrid = ba_onegrid.boxList();
    for (auto& b : bl2d_onegrid) {
        b.setRange(2,0);
    }
    BoxArray ba2d_onegrid(std::move(bl2d_onegrid));
    Vector<int> pmap;
    pmap.resize(1);
    pmap[0]=0;
    DistributionMapping dm_onegrid(ba2d_onegrid);
    dm_onegrid.define(pmap);
*/

    // Make local copy of xvel, yvel
    MultiFab::Copy (xvel_data, lev_new[Vars::xvel], 0, 0, 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab::Copy (yvel_data, lev_new[Vars::yvel], 0, 0, 1, lev_new[Vars::yvel].nGrowVect());


    // Initialize Multifabs to store x_avg, y_avg, u_mag, u_dir at cell centers
    MultiFab x_avg(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, lev_new[Vars::cons].nGrowVect());

    MultiFab y_avg(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, lev_new[Vars::cons].nGrowVect());

    MultiFab u_mag(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, lev_new[Vars::cons].nGrowVect());


    MultiFab u_dir(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(), 1, lev_new[Vars::cons].nGrowVect());


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

            amrex::AllPrintToFile("uvel.txt") << amrex::IntVect(i,j,k) << " [" <<velx_arr(i,j,k) << "| avg:  " << u_vel(i,j,k)<< " | " << velx_arr(i+1,j,k) << "]" <<std::endl;
        });
    }

    for (MFIter mfi(y_avg,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        const Array4<Real>& v_vel = y_avg.array(mfi);
        const Array4<const Real>& vely_arr = yvel_data.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

            v_vel(i,j,k)  = 0.5 *( vely_arr(i,j,k) + vely_arr(i,j+1,k) );

            amrex::AllPrintToFile("vvel.txt") << "%ld" << amrex::IntVect(i,j,k) << " ["<<vely_arr(i,j,k)<<"| avg: " << v_vel(i,j,k)<< " | " << vely_arr(i,j+1,k) << "]" <<std::endl;
        });
    }


    for (MFIter mfi(u_mag, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box bx = mfi.tilebox();
        const Array4<Real>& magnitude = u_mag.array(mfi);
        const Array4<const Real>& u = x_avg.array(mfi);
        const Array4<const Real>& v = y_avg.array(mfi);
        const Array4<Real>& theta = u_dir.array(mfi);

        amrex::Vector<std::unique_ptr<amrex::MultiFab>> magnitude_onegrid;
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> theta_onegrid;


 // create a new BoxArray and DistributionMapping for a MultiFab with 1 box
    BoxArray ba_onegrid(geom[lev].Domain());
    BoxList bl2d_onegrid = ba_onegrid.boxList();
    for (auto& b : bl2d_onegrid) {
        b.setRange(2,0);
    }
    BoxArray ba2d_onegrid(std::move(bl2d_onegrid));
    Vector<int> pmap;
    pmap.resize(1);
    pmap[0]=0;
    DistributionMapping dm_onegrid(ba2d_onegrid);
    dm_onegrid.define(pmap);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k){

 //           magnitude(i,j,k)  = std::sqrt( pow(u(i,j,k), 2) + pow(v(i,j,k), 2) );

            double u_val = u(i, j, k);
            double v_val = v(i, j, k);

            if ( u_val == 0 ) {
                u_val = std::max( u_val, 1e-15 );  // Ensure u_val is non-zero
            }

            magnitude(i,j,k) = std::sqrt( pow(u_val, 2) + pow(v_val, 2) );

            if ( u_val < 0 && v_val > 0 || u_val < 0 && v_val < 0 ) {

                theta(i,j,k) = PI + ( atan( v_val / u_val ) );

            } else {

                theta(i,j,k) = atan ( v_val / u_val );
            }


            amrex::AllPrintToFile("mag_theta.txt") << amrex::IntVect(i,j,k) <<  " Magnitude: " << magnitude(i,j,k) << " Theta: " << theta(i,j,k) <<std::endl;
  });


    // Send the 2D slice at k_ref
        // Box slice_box = bx;
        amrex::IntVect boxSmall = bx.smallEnd();
        amrex::IntVect boxBig = bx.bigEnd();
        Box slice_box_ref = makeSlab(bx, 2, k_ref);

    // Calculate the number of elements in the current box
    int n_elements = slice_box_ref.numPts();

    // Initialize vectors to send to WW3
    std::vector<Real> magnitude_values(n_elements);
    std::vector<Real> theta_values(n_elements);
    std::vector<amrex::IntVect> indices(n_elements);
    // Copy values
    int counter = 0;
    for (BoxIterator bi(slice_box_ref); bi.ok(); ++bi) {
    IntVect iv = bi();
    magnitude_values[counter] = magnitude(iv);
    theta_values[counter] = theta(iv);
    indices[counter] = iv;
    ++counter;
    }
//    timedif2 = ( ((double) clock()) / CLOCKS_PER_SEC) - clkStart2;
//    amrex::AllPrintToFile("timer.txt") << "It took " << (double) timedif2 << " seconds to reach the part before sending" << std::endl;
//amrex::Print() << "It took " << (double) timedif2 << " seconds to reach the part before sending" << std::endl;

// Print magnitude values and corresponding IntVect indices
for (int j = 0; j < n_elements; ++j) {
    amrex::AllPrintToFile("debug_send.txt")
        << "dz, k_ref " << dz << ", " << k_ref << " "
        << "Index: " << j
        << ", IntVect: (" << indices[j][0] << ", " << indices[j][1] << ", " << indices[j][2] << ")"
        << ", Magnitude: " << magnitude_values[j]
        << ", Theta: " << theta_values[j]
        << std::endl;
}

         int rank_offset = amrex::MPMD::MyProc() - amrex::ParallelDescriptor::MyProc();
         int this_root, other_root;
         if (rank_offset == 0) { // First program
             this_root = 0;
             other_root = amrex::ParallelDescriptor::NProcs();
         } else {
             this_root = rank_offset;
             other_root = 0;
         }


         amrex::Print()<< "Sending " << n_elements << " from ERF::send_to_ww3 now" << std::endl;

         if (amrex::MPMD::MyProc() == this_root) {
             if (rank_offset == 0) // First program
             {
             MPI_Send(&n_elements, 1, MPI_INT, other_root, 11, MPI_COMM_WORLD);
MPI_Send(magnitude_values.data(), n_elements, MPI_DOUBLE, other_root, 13, MPI_COMM_WORLD);
MPI_Send(theta_values.data(), n_elements, MPI_DOUBLE, other_root, 15, MPI_COMM_WORLD);
             }
             else // Second program
             {
                 MPI_Send(&n_elements, 1, MPI_INT, other_root, 10, MPI_COMM_WORLD);
MPI_Send(magnitude_values.data(), n_elements, MPI_DOUBLE, other_root, 12, MPI_COMM_WORLD);
MPI_Send(theta_values.data(), n_elements, MPI_DOUBLE, other_root, 14, MPI_COMM_WORLD);
                 //MPI_Recv(&nx, 1, MPI_INT, other_root, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                 //MPI_Recv(&ny, 1, MPI_INT, other_root, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
             }
         }
    timedif = ( ((double) clock()) / CLOCKS_PER_SEC) - clkStart;

// amrex::Real myclock = ParallelDescriptor::second();
//    amrex::AllPrintToFile("timer.txt") << "At " << myclock << " seconds I reached the end of send_to_ww3" << std::endl;

     amrex::AllPrintToFile("timer.txt") << "It took " << timedif << " seconds to reach the end of send_to_WW3" << std::endl;


}

}
#endif

