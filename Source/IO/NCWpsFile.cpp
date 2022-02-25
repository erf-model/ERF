#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>

#include "NCInterface.H"
#include "NCWpsFile.H"
#include <ERF.H>

template<typename DType>
void ReadNCWpsFile(const std::string  &fname, amrex::Vector<NDArray<DType> >& arrays, 
                  amrex::Vector<std::string> names ) {

    amrex::Print() << "Reading NetCDF WPS file: " << fname << "\n";
  
    AMREX_ASSERT(arrays.size() == names.size()); 

    if (amrex::ParallelDescriptor::IOProcessor())
    {
      auto ncf = ncutils::NCFile::create(fname, NC_CLOBBER | NC_NETCDF4);
      
      // get the dimension information
      int DateStrLen         = static_cast<int>(ncf.dim("DateStrLen").len());
      int Time               = static_cast<int>(ncf.dim("Time").len());
      int west_east          = static_cast<int>(ncf.dim("west_east").len());
      int south_north        = static_cast<int>(ncf.dim("south_north").len());
      int num_metgrid_levels = static_cast<int>(ncf.dim("num_metgrid_levels").len());
      int num_st_layers      = static_cast<int>(ncf.dim("num_st_layers").len());
      int num_sm_layers      = static_cast<int>(ncf.dim("num_sm_layers").len());
      int south_north_stag   = static_cast<int>(ncf.dim("south_north_stag").len());
      int west_east_stag     = static_cast<int>(ncf.dim("west_east_stag").len());
      int z_dimension0012    = static_cast<int>(ncf.dim("z-dimension0012").len());
      int z_dimension0016    = static_cast<int>(ncf.dim("z-dimension0016").len());
      int z_dimension0024    = static_cast<int>(ncf.dim("z-dimension0024").len());

      for (auto n=0; n<arrays.size(); ++n) {
         // get the pressure
         std::vector<size_t> shape = ncf.var(names[n]).shape();
         arrays[n] = NDArray<DType>(names[n],shape);
         DType* dataPtr = arrays[n].get_data();
         ncf.var(names[n]).get(dataPtr, {0}, {shape.size()});
      }
     ncf.close();
   }
}
