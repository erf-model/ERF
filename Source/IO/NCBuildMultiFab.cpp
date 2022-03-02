#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "ERF.H"
#include "NCInterface.H"

//
// This is the example function that demonstrate how to build the MultiFab using the
// NetCDF variables.
//
void BuildMultiFabFromNCFile(const std::string &fname) {

   using RARRAY = NDArray<float>;
   amrex::Vector<RARRAY> arrays(1);
   ReadNCWpsFile(fname, {"PRES"}, arrays);

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
   ba.set(0, box);

   // create a distribution mapping
   DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

   MultiFab pres;
   pres.define(convert(ba,IntVect(0,0,0)), dm, 1, 0);

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

}
