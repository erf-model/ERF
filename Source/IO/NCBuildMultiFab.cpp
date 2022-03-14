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
void BuildMultiFabFromNCFile(const std::string &fname) {

   using RARRAY = NDArray<float>;
   amrex::Vector<RARRAY> arrays(1);
   amrex::Vector<std::string> nc_vars = {"PRES"};
   ReadNCWpsFile(fname, nc_vars, arrays);

   //
   // build the box first using the shape information
   amrex::Print() << "Building MultiFab from the netcdf variable just read: " << nc_vars[0] << std::endl;

   amrex::Print() << "Creating dimension, box types, etc. " << std::endl;
   amrex::IntVect smallEnd{0,0,0};
   amrex::IntVect boxtyp{0,0,0};
   amrex::IntVect bigEnd(3);

   amrex::Print() << "Assigning shapes to the box based on dimensions read from NetCDF file. " << std::endl;
   bigEnd[0] = arrays[0].get_vshape()[1];
   bigEnd[1] = arrays[0].get_vshape()[2];
   bigEnd[2] = arrays[0].get_vshape()[3];

   amrex::Box box = amrex::Box(smallEnd, bigEnd, boxtyp);

   BoxArray ba;
   ba.define(box);

   amrex::Print() << "Creating distribution mapping. " << std::endl;
   // create a distribution mapping
   DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

   MultiFab pres;
   pres.define(convert(ba,IntVect(0,0,0)), dm, 1, 0);
    MFIter::allowMultipleMFIters(true);

   amrex::Print() << "Now filling the MultiFab. " << std::endl;
   // assign the data to multifab
   for (amrex::MFIter mfi(pres); mfi.isValid(); ++mfi) {
       auto ncomp   = pres.nComp();
       auto box     = pres.get(mfi).box();
       auto num_pts = pres.get(mfi).numPts();

       for (int k(0); k < ncomp; ++k) {
          auto dataPtr = pres.get(mfi).dataPtr(k);
          for (int n(0); n < num_pts; ++n) {
            *(dataPtr+n) = static_cast<amrex::Real>(*(arrays[0].get_data()+n));
          }
       }
  }
  amrex::Print() << "Done filling he MultiFab. " << std::endl;

}

// Function to read a NetCDF variable and fill a corresponding multifab
void BuildMultiFabFromNCFile(const std::string &fname, const std::string &nc_var_name, MultiFab& mf_var) {
    amrex::Print() << "Reading the NetCDF File: " << fname << std::endl;
    amrex::Print() << "NetCDF var name: " << nc_var_name << std::endl;

    using RARRAY = NDArray<float>;
    amrex::Vector<RARRAY> arrays(1);
    ReadNCWpsFile(fname, {nc_var_name}, arrays);

    //
    // build the box first using the shape information
    //
    amrex::IntVect smallEnd{0,0,0};
    amrex::IntVect boxtyp{0,0,0};
    amrex::IntVect bigEnd(3);

    bigEnd[0] = arrays[0].get_vshape()[1];
    bigEnd[1] = arrays[0].get_vshape()[2];
    bigEnd[2] = arrays[0].get_vshape()[3];

    amrex::Box box = amrex::Box(smallEnd, bigEnd, boxtyp);

    BoxArray ba;
    ba.define(box);

    // create a distribution mapping
    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

    mf_var.define(convert(ba,IntVect(0,0,0)), dm, 1, 0);

    // assign the data to multifab
    for (amrex::MFIter mfi(mf_var); mfi.isValid(); ++mfi) {
        auto ncomp   = mf_var.nComp();
        auto box     = mf_var.get(mfi).box();
        auto num_pts = mf_var.get(mfi).numPts();

        for (int k(0); k < ncomp; ++k) {
            auto dataPtr = mf_var.get(mfi).dataPtr(k);
            for (int n(0); n < num_pts; ++n) {
                *(dataPtr+n) = static_cast<amrex::Real>(*(arrays[0].get_data()+n));
            }
        }
    }

}