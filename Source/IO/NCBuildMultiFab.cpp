#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "ERF.H"
#include "NCInterface.H"
#include "AMReX_Print.H"

//
// This is the example function that demonstrate how to build the MultiFab using the
// NetCDF variables.
//
void BuildMultiFabFromMetgridOutputFileDemo(const std::string &fname) {

   using RARRAY = NDArray<float>;
   amrex::Vector<RARRAY> arrays(1);
   amrex::Vector<std::string> nc_vars = {"PRES"};
   ReadMetgridOutputFile(fname, nc_vars, arrays);

   // build the box first using the shape information
   amrex::Print() << "Building MultiFab from the netcdf variable just read: " << nc_vars[0] << std::endl;

   //amrex::Print() << "Creating dimension, box types, etc. " << std::endl;
   amrex::IntVect smallEnd{0,0,0};
   amrex::IntVect boxtyp{0,0,0};
   amrex::IntVect bigEnd(3);

  // amrex::Print() << "Assigning shapes to the box based on dimensions read from NetCDF file. " << std::endl;
   bigEnd[0] = arrays[0].get_vshape()[2] - 1;
   bigEnd[1] = arrays[0].get_vshape()[3] - 1;
   bigEnd[2] = arrays[0].get_vshape()[1] - 1;
   int tot_size = arrays[0].get_vshape()[1]*arrays[0].get_vshape()[2]*arrays[0].get_vshape()[3];
   amrex::Print() << "Total points in the box: " << tot_size << std::endl;

   amrex::Box box = amrex::Box(smallEnd, bigEnd, boxtyp);

   BoxArray ba;
   ba.define(box);

   //amrex::Print() << "Creating distribution mapping. " << std::endl;
   // create a distribution mapping
   //DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
   DistributionMapping dm (ba);
   //amrex::Print() << "Successfully created distribution mapping. " << std::endl;

   MultiFab pres;
   pres.define(convert(ba,IntVect(0,0,0)), dm, 1, 0);
    MFIter::allowMultipleMFIters(true);

   //amrex::Print() << "Now filling the MultiFab. " << std::endl;
   // assign the data to multifab
   for (amrex::MFIter mfi(pres); mfi.isValid(); ++mfi) {
       auto ncomp   = pres.nComp();
       auto box     = pres.get(mfi).box();
       auto num_pts = pres.get(mfi).numPts();

       for (int k(0); k < ncomp; ++k) {
          auto dataPtr = pres.get(mfi).dataPtr(k);
          for (int n(0); n < num_pts; ++n) {
            *(dataPtr+n) = static_cast<amrex::Real>(*(arrays[0].get_data()+n));
//            if ((n >= tot_size - 100) && (n - tot_size)%10 == 0)
//            if (n % 200 == 0)
//                amrex::Print() << "n: " << n << ", data[n]: " << *(dataPtr+n) << std::endl;
          }
       }
  }
  //amrex::Print() << "Done filling he MultiFab. " << std::endl;
}

// Function to read a NetCDF variable and fill a corresponding MultiFab and Array4
Box BuildMultiFabFromIdealOutputFile(const std::string &fname, const std::string &nc_var_name, MultiFab &mf_var,
                                     Array4<Real> &array4_var, const enum NC_Data_Dims_Type &nc_data_dims_type) {

    amrex::Print() << "Reading the NetCDF File: " << fname << std::endl;
    //amrex::Print() << "NetCDF var name: " << nc_var_name << std::endl;

    using RARRAY = NDArray<float>;
    amrex::Vector<RARRAY> arrays(1);
    ReadIdealOutputFile(fname, {nc_var_name}, arrays);

    // build the box first using the shape information
    //amrex::Print() << "Building MultiFab from the netcdf variable just read: " << nc_var_name << std::endl;
    amrex::IntVect smallEnd{0,0,0};
    amrex::IntVect boxType{0, 0, 0};
    amrex::IntVect bigEnd{0, 0, 0};

    //amrex::Print() << "Assigning shapes to the box based on dimensions read from NetCDF file. " << std::endl;
//    bigEnd[0] = arrays[0].get_vshape()[1] - 1;
//    bigEnd[1] = arrays[0].get_vshape()[2] - 1;
//    bigEnd[2] = arrays[0].get_vshape()[3] - 1;

    switch (nc_data_dims_type) {
        case NC_Data_Dims_Type::Time_BT_SN_WE:
            bigEnd[0] = arrays[0].get_vshape()[3] - 1;
            bigEnd[1] = arrays[0].get_vshape()[2] - 1;
            bigEnd[2] = arrays[0].get_vshape()[1] - 1;
            break;
        case NC_Data_Dims_Type::Time_SN_WE:
            bigEnd[0] = arrays[0].get_vshape()[2] - 1;
            bigEnd[1] = arrays[0].get_vshape()[1] - 1;
            bigEnd[2] = 0;
            break;
        case NC_Data_Dims_Type::Time_BT:
            bigEnd[0] = 0;
            bigEnd[1] = 0;
            bigEnd[2] = arrays[0].get_vshape()[1] - 1;
            break;
        case NC_Data_Dims_Type::Time:
            bigEnd[0] = 0;
            bigEnd[1] = 0;
            bigEnd[2] = 0;
            break;
        default:
            amrex::Error("Unrecognized NetCDF data dimensions type");
    }
    int tot_size = (bigEnd[0] + 1) * (bigEnd[1] + 1) * (bigEnd[2] + 1);
    amrex::Print() << "Total points in the box constructed from netCDF variable: " << tot_size << std::endl;

    amrex::Box box = amrex::Box(smallEnd, bigEnd, boxType);

    BoxArray ba;
    ba.define(box);

    //amrex::Print() << "Creating distribution mapping. " << std::endl;
    // create a distribution mapping
    //DistributionMapping dm { ba, ParallelDescriptor::NProcs() };
    DistributionMapping dm { ba };
    //amrex::Print() << "Successfully created distribution mapping. " << std::endl;

    mf_var.define(convert(ba,IntVect(0,0,0)), dm, 1, 0);
    MFIter::allowMultipleMFIters(true);

    Box box_nc;
    // assign the data to multifab
    //amrex::Print() << "Now filling the MultiFab. " << std::endl;
    for (amrex::MFIter mfi(mf_var); mfi.isValid(); ++mfi) {
        // Retrieve some attributes of the MultiFab (Sanity Check)
        auto ncomp   = mf_var.nComp();
        box_nc     = mf_var.get(mfi).box();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE( (bigEnd == box_nc.bigEnd() && smallEnd == box_nc.smallEnd()),
                "Dimensions of box used to define the Multifab and the box retrieved don't match");
        auto num_pts = mf_var.get(mfi).numPts();
        //auto type = mf_var.get(mfi).getType();

        for (int k(0); k < ncomp; ++k) {
            auto dataPtr = mf_var.get(mfi).dataPtr(k);
            for (int n(0); n < num_pts; ++n) {
                *(dataPtr+n) = static_cast<amrex::Real>(*(arrays[0].get_data()+n));
//                if (n % 100 == 0)
//                    amrex::Print() << "n: " << n << ", data[n]: " << *(dataPtr+n) << std::endl;
            }
        }
        array4_var = mf_var.array(mfi);
    }
    amrex::Print() << "Done filling the MultiFab. " << std::endl << std::endl;

    return box_nc;
}