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
void BuildFABsFromMetgridOutputFileDemo(const std::string &fname) {

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
}

// Function to read a NetCDF variable and fill a corresponding MultiFab and Array4
void
BuildFABsFromIdealOutputFile(const std::string &fname,
                             Vector<std::string> nc_var_names,
                             Vector<amrex::FArrayBox*> fab_vars,
                             const enum NC_Data_Dims_Type &nc_data_dims_type)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        using RARRAY = NDArray<float>;
        amrex::Vector<RARRAY> arrays(nc_var_names.size());
        ReadIdealOutputFile(fname, nc_var_names, arrays);

        for (int i = 0; i < nc_var_names.size(); i++)
        {
            // Build the box first using the shape information
            amrex::IntVect smallEnd{0, 0, 0};
            amrex::IntVect bigEnd{0, 0, 0};

            switch (nc_data_dims_type) {
                case NC_Data_Dims_Type::Time_BT_SN_WE:
                    bigEnd[0] = arrays[i].get_vshape()[3] - 1;
                    bigEnd[1] = arrays[i].get_vshape()[2] - 1;
                    bigEnd[2] = arrays[i].get_vshape()[1] - 1;
                    break;
                case NC_Data_Dims_Type::Time_SN_WE:
                    bigEnd[0] = arrays[i].get_vshape()[2] - 1;
                    bigEnd[1] = arrays[i].get_vshape()[1] - 1;
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
            // amrex::Print() << "Total points in the box constructed from netCDF variable: " << tot_size << std::endl;

            amrex::IntVect boxType{0, 0, 0};
            if (nc_var_names[i] == "V") boxType = amrex::IntVect(0,1,0);
            if (nc_var_names[i] == "W") boxType = amrex::IntVect(0,0,1);
            if (nc_var_names[i] == "U") boxType = amrex::IntVect(1,0,0);
            amrex::Box bx = amrex::Box(smallEnd, bigEnd, boxType);

            AMREX_ALWAYS_ASSERT(bx == fab_vars[i]->box());

            int ncomp  = 1;
            auto num_pts = fab_vars[i]->box().numPts();
            for (int k(0); k < ncomp; ++k) {
                auto dataPtr = fab_vars[i]->dataPtr(k);
                for (int n(0); n < num_pts; ++n) {
                    *(dataPtr+n) = static_cast<amrex::Real>(*(arrays[i].get_data()+n));
                }
            }
        }
    } // if IOProcessor
}
