#ifndef _NCWPSFILE_H_
#define _NCWPSFILE_H_

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "NCInterface.H"

void BuildMultiFabFromNCFile(const std::string &fname);
void BuildMultiFabFromNCFile(const std::string &fname, const std::string &nc_var_name, MultiFab& mf_var);

//
// NDArray is the datatype designed to hold any data, including scalars, multidimensional
// arrays, that read from the NetCDF file.
//
// The data read from NetCDF file are stored in a continuous memory, and the data layout is descriped
// by using a vector (shape). AMRex Box can be constructed using the data shape information, and MultiFab
// data array can be setup using the data that stored in the NDArray.
//
template <typename DataType>
struct NDArray
{
   using DType = typename std::remove_const<DataType>::type;

   // constructor
   explicit NDArray (const std::string vname, const std::vector<size_t>& vshape)
                    : name{vname}, shape{vshape}, ref_counted{1}, owned{true} {
      data = new DType [this->ndim()];
   }

   // default constructor
   NDArray() : name{"null"}, data{nullptr}, ref_counted{1}, owned{false} {}

   // copy constructor
   NDArray(const NDArray& array) {
     name  = array.name;
     shape = array.shape;
     data  = array.data;
     owned = false;
     ref_counted.fetch_add(1, std::memory_order_relaxed);
   }

   // copy assignment
   NDArray& operator=(const NDArray& array) {
      name  = array.name;
      shape = array.shape;
      data  = array.data;
      owned = false;
      ref_counted.fetch_add(1, std::memory_order_relaxed);
      return *this;
    }

   // destructor
   ~NDArray() {
     ref_counted.fetch_sub(1, std::memory_order_acq_rel);
     if (ref_counted == 1 && owned) delete[] data;
   }

   // get the data pointer
   decltype(auto) get_data() {
     ref_counted.fetch_add(1, std::memory_order_relaxed);
     return data;
   }

   // get the variable name
   std::string get_vname() {
     return name;
   }

   // get the variable data shape
   std::vector<size_t> get_vshape() {
     return shape;
   }

   // return the total number of data
   size_t ndim() {
     size_t num = 1;
     for (auto i=0; i<shape.size(); ++i) num *= shape[i];
     return num;
   }

   // set the data shape information
   void set_vshape(std::vector<size_t> vshape) {
     shape = vshape;
    }

 private:
   std::string name;
   std::vector<size_t> shape;
   DType* data;
   bool owned;
   mutable std::atomic<size_t> ref_counted;
};

//
// This function read NetCDF file.
//
// inputs:
//   fname is the NetCDF filename
//   names is the variables names that read from the NetCDF file
//
// output:
//   arrays is the vector of NDArray that stores the data read from the NetCDF file.
//
template<typename DType>
void ReadNCWpsFile(const std::string  &fname, amrex::Vector<std::string> names,
                  amrex::Vector<NDArray<DType> >& arrays) {

    amrex::Print() << "Reading NetCDF WPS file: " << fname << "\n";

    AMREX_ASSERT(arrays.size() == names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
      auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

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

     // This following dimension defined based on different timestep, the dimensional name
     // varies at different NetCDF output.
     // TODO: should we keep reading this information from the NetCDF file? we actually don't use it.
     //
     // int z_dimension0012    = static_cast<int>(ncf.dim("z-dimension0012").len());
     // int z_dimension0016    = static_cast<int>(ncf.dim("z-dimension0016").len());
     // int z_dimension0024    = static_cast<int>(ncf.dim("z-dimension0024").len());

      for (auto n=0; n<arrays.size(); ++n) {
         // read the data from NetCDF file
         std::string vname         = names[n];
         std::vector<size_t> shape = ncf.var(vname).shape();
         arrays[n]                 = NDArray<DType>(vname,shape);
         auto numPts               = arrays[n].ndim();
         DType* dataPtr            = arrays[n].get_data();

         std::vector<size_t> start(shape.size(), 0);
         std::cout << "name= " << vname << " ;numpts= " << numPts  << std::endl;
         ncf.var(vname).get(dataPtr, start, shape);
      }
     ncf.close();
   }
    amrex::Print() << "Succesfully read the NetCDF File " << fname << std::endl;
}
#endif
