#ifndef _NCWPSFILE_H_
#define _NCWPSFILE_H_

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "AMReX_FArrayBox.H"
#include "AMReX_IArrayBox.H"
#include "NCInterface.H"

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

/*
 // Read from metgrid
 NetCDF variables of dimensions Time_BT_SN_WE: "UU", "VV", "TT", "RH", "PRES", "GHT"
 NetCDF variables of dimensions Time_SN_WE   : "HGT", "MAPFAC_U",  "MAPFAC_V",  "MAPFAC_M", "PSFC"
 NetCDF global attributes of type int        : "FLAG_PSFC", "FLAG_MAPFAC_U", "FLAG_MAPFAC_V", "FLAG_MAPFAC_M",
                                               "FLAG_HGT_M", "WEST-EAST_GRID_DIMENSION", "SOUTH-NORTH_GRID_DIMENSION"
 NetCDF global attributes of type string     : "SIMULATION_START_DATE"
 NetCDF global attributes of type real       : "DX", "DY"

 // Read from wrfbdy
 NetCDF variables of dimensions Time_BdyWidth_BT_SN         : "U_BXS", "U_BXE", "V_BXS", "V_BXE" etc.
 NetCDF variables of dimensions Time_BdyWidth_BT_WE         : "U_BYS", "U_BYE", "V_BYS", "V_BYE" etc.
 NetCDF variables of dimensions Time_BdyWidth_SN            : "MU_BXS", "MU_BXE", "PC_BXS", "PC_BXE", etc.
 NetCDF variables of dimensions Time_BdyWidth_WE            : "MU_BYS", "MU_BYE", "PC_BYS", "PC_BYE", etc.
*/
enum class NC_Data_Dims_Type {
    Time_BT_SN_WE,
    Time_SN_WE,
    Time_BT,
    Time,
    Time_BdyWidth_BT_SN,
    Time_BdyWidth_BT_WE,
    Time_BdyWidth_SN,
    Time_BdyWidth_WE
};

//
// NDArray is the datatype designed to hold any data, including scalars, multidimensional
// arrays, that read from the NetCDF file.
//
// The data read from NetCDF file are stored in a continuous memory, and the data layout is described
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
   NDArray () : name{"null"}, ref_counted{1}, owned{false}, data{nullptr} {}

   // copy constructor
   NDArray (const NDArray& array) {
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
   ~NDArray () {
     ref_counted.fetch_sub(1, std::memory_order_acq_rel);
     if (ref_counted == 1 && owned) delete[] data;
   }

   // get the data pointer
   decltype(auto) get_data () {
     ref_counted.fetch_add(1, std::memory_order_relaxed);
     return data;
   }

   // get the variable name
   std::string get_vname () {
     return name;
   }

   // get the variable data shape
   std::vector<size_t> get_vshape () {
     return shape;
   }

   // return the total number of data
   size_t ndim () {
     size_t num = 1;
     int isize = static_cast<int>(shape.size());
     for (auto i=0; i<isize; ++i) num *= shape[i];
     return num;
   }

   // set the data shape information
   void set_vshape (std::vector<size_t> vshape) {
     shape = vshape;
    }

 private:
   std::string name;
   std::vector<size_t> shape;
   mutable std::atomic<size_t> ref_counted;
   bool owned;
   DType* data;
};

int BuildFABsFromWRFBdyFile (const std::string &fname,
                             amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xlo,
                             amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_xhi,
                             amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_ylo,
                             amrex::Vector<amrex::Vector<amrex::FArrayBox>>& bdy_data_yhi);

AMREX_FORCE_INLINE
std::time_t
getEpochTime (const std::string& dateTime, const std::string& dateTimeFormat)
{
    // Create a stream which we will use to parse the string,
    // which we provide to constructor of stream to fill the buffer.
    std::istringstream ss{ dateTime };

    // Create a tm object to store the parsed date and time.
    std::tm tmTime;
    memset(&tmTime, 0, sizeof(tmTime));

    // Now we read from buffer using get_time manipulator
    // and formatting the input appropriately.
    strptime(dateTime.c_str(), dateTimeFormat.c_str(), &tmTime);

    // Convert the tm structure to time_t value and return.
    // Here we use timegm since the output should be relative to UTC.
    auto epoch = timegm(&tmTime);
    // Print() << "Time Stamp: "<< std::put_time(&tmTime, "%c")
    //         << " , Epoch: " << epoch << std::endl;

    return epoch;
}

template<typename DType>
void ReadNetCDFFile (const std::string& fname, amrex::Vector<std::string> names,
                     amrex::Vector<NDArray<DType> >& arrays)
{
    AMREX_ASSERT(arrays.size() == names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        auto ncf = ncutils::NCFile::open(fname, NC_CLOBBER | NC_NETCDF4);

        /*
        // get the dimension information
        int Time               = static_cast<int>(ncf.dim("Time").len());
        int DateStrLen         = static_cast<int>(ncf.dim("DateStrLen").len());
        int west_east          = static_cast<int>(ncf.dim("west_east").len());
        int south_north        = static_cast<int>(ncf.dim("south_north").len());
        int bottom_top         = static_cast<int>(ncf.dim("bottom_top").len());
        int bottom_top_stag    = static_cast<int>(ncf.dim("bottom_top_stag").len());
        int west_east_stag     = static_cast<int>(ncf.dim("west_east_stag").len());
        int south_north_stag   = static_cast<int>(ncf.dim("south_north_stag").len());
        int bdy_width          = static_cast<int>(ncf.dim("bdy_width").len());
        */

        // amrex::Print() << "Reading the dimensions from the netcdf file " << "\n";
        for (auto n=0; n<arrays.size(); ++n) {
            // read the data from NetCDF file
            std::string vname_to_write = names[n];
            std::string vname_to_read  = names[n];
            if (vname_to_read.substr(0,2) == "R_") {
                vname_to_read  = names[n+4]; // This allows us to read "T" instead -- we will over-write this later
            }

            /*
            amrex::AllPrint() << "About to read " << vname_to_read
                              << " while filling the array for " << vname_to_write << std::endl;
            */

            std::vector<size_t> shape = ncf.var(vname_to_read).shape();
            arrays[n]                 = NDArray<DType>(vname_to_read,shape);
            DType* dataPtr            = arrays[n].get_data();

            std::vector<size_t> start(shape.size(), 0);

            // auto numPts               = arrays[n].ndim();
            // amrex::Print() << "NetCDF Variable name = " << vname_to_read << std::endl;
            // amrex::Print() << "numPts read from NetCDF file/var = " << numPts << std::endl;
            // amrex::Print()  << "Dims of the variable = (";
            // for (auto &dim:shape)
            //     amrex::Print() << dim << ", " ;
            // amrex::Print()  << ")" << std::endl;

            ncf.var(vname_to_read).get(dataPtr, start, shape);
        }
        ncf.close();
    }
}

/**
 * Helper function for reading data from NetCDF file into a
 * provided FAB.
 *
 * @param iv Index for which variable we are going to fill
 * @param nc_arrays Arrays of data from NetCDF file
 * @param var_name Variable name
 * @param NC_dim_type Dimension type for the variable as stored in the NetCDF file
 * @param temp FAB where we store the variable data from the NetCDF Arrays
 */
template<class FAB,typename DType>
void
fill_fab_from_arrays (int iv,
                      amrex::Real& Latitude,
                      amrex::Real& Longitude,
                      std::string& Lat_var_name,
                      std::string& Lon_var_name,
                      amrex::Vector<NDArray<float>>& nc_arrays,
                      const std::string& var_name,
                      NC_Data_Dims_Type& NC_dim_type,
                      FAB& temp)
{
    int ns1, ns2, ns3;
    if (NC_dim_type == NC_Data_Dims_Type::Time_BT) {
        ns1 = nc_arrays[iv].get_vshape()[1];
        ns2 = 1;
        ns3 = 1;
        // amrex::Print() << "TYPE BT " << ns1 << std::endl;
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_SN_WE) {
        ns1 = 1;
        ns2 = nc_arrays[iv].get_vshape()[1];
        ns3 = nc_arrays[iv].get_vshape()[2];
        // amrex::Print() << "TYPE SN WE " << ns2 << " " << ns3 << std::endl;
    } else if (NC_dim_type == NC_Data_Dims_Type::Time_BT_SN_WE) {
        ns1 = nc_arrays[iv].get_vshape()[1];
        ns2 = nc_arrays[iv].get_vshape()[2];
        ns3 = nc_arrays[iv].get_vshape()[3];
        // amrex::Print() << "TYPE BT SN WE " << ns1 << " " << ns2 << " " << ns3 << std::endl;
    } else {
        amrex::Abort("Dont know this NC_Data_Dims_Type");
    }

    // TODO:  The box will only start at (0,0,0) at level 0 -- we need to generalize this
    amrex::Box my_box(amrex::IntVect(0,0,0), amrex::IntVect(ns3-1,ns2-1,ns1-1));
    // amrex::Print() <<" MY BOX " << my_box << std::endl;

    if (var_name == "U" || var_name == "UU" ||
        var_name == "MAPFAC_U" || var_name == "MAPFAC_UY") my_box.setType(amrex::IndexType(amrex::IntVect(1,0,0)));
    if (var_name == "V" || var_name == "VV" ||
        var_name == "MAPFAC_V" || var_name == "MAPFAC_VY") my_box.setType(amrex::IndexType(amrex::IntVect(0,1,0)));
    if (var_name == "W" || var_name == "WW") my_box.setType(amrex::IndexType(amrex::IntVect(0,0,1)));

    amrex::Arena* Arena_Used = amrex::The_Arena();
#ifdef AMREX_USE_GPU
    // Make sure temp lives on CPU since nc_arrays lives on CPU only
    Arena_Used = amrex::The_Pinned_Arena();
#endif
    temp.resize(my_box,1, Arena_Used);
    amrex::Array4<DType> fab_arr = temp.array();

    int ioff = temp.box().smallEnd()[0];
    int joff = temp.box().smallEnd()[1];

    auto num_pts = my_box.numPts();

    // amrex::Print() <<" ns1 * ns2 * ns3 " << ns1 * ns2 * ns3 << std::endl;
    // amrex::Print() <<" NUMPTS " << num_pts << std::endl;

    for (int n(0); n < num_pts; ++n) {
        int k  = n / (ns2*ns3);
        int j  = (n - k*(ns2*ns3)) / ns3 + joff;
        int i  = n - k*(ns2*ns3) - (j-joff) * ns3 + ioff;
        fab_arr(i,j,k,0) = static_cast<DType>(*(nc_arrays[iv].get_data()+n));
    }

    if (var_name == Lat_var_name) Latitude  = fab_arr(0,0,0);
    if (var_name == Lon_var_name) Longitude = fab_arr(0,0,0);
}

/**
 * Function to read NetCDF variables and fill the corresponding Array4's
 *
 * @param fname Name of the NetCDF file to be read
 * @param nc_var_names Variable names in the NetCDF file
 * @param NC_dim_types NetCDF data dimension types
 * @param fab_vars Fab data we are to fill
 */
template<class FAB,typename DType>
void
BuildFABsFromNetCDFFile (const amrex::Box& domain,
                         amrex::Real& Latitude,
                         amrex::Real& Longitude,
                         std::string& Lat_var_name,
                         std::string& Lon_var_name,
                         const std::string &fname,
                         amrex::Vector<std::string> nc_var_names,
                         amrex::Vector<enum NC_Data_Dims_Type> NC_dim_types,
                         amrex::Vector<FAB*> fab_vars)
{
    int ioproc = amrex::ParallelDescriptor::IOProcessorNumber();  // I/O rank

    amrex::Vector<NDArray<float>> nc_arrays(nc_var_names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        ReadNetCDFFile(fname, nc_var_names, nc_arrays);
    }

    for (int iv = 0; iv < nc_var_names.size(); iv++)
    {
        FAB tmp;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            fill_fab_from_arrays<FAB,DType>(iv, Latitude, Longitude,
                                            Lat_var_name, Lon_var_name,
                                            nc_arrays, nc_var_names[iv],
                                            NC_dim_types[iv], tmp);
        }

        int ncomp = tmp.nComp();
        amrex::Box box = tmp.box();

        amrex::ParallelDescriptor::Bcast(&box,   1, ioproc);
        amrex::ParallelDescriptor::Bcast(&ncomp, 1, ioproc);

        if (!amrex::ParallelDescriptor::IOProcessor()) {
#ifdef AMREX_USE_GPU
            tmp.resize(box,ncomp,amrex::The_Pinned_Arena());
#else
            tmp.resize(box,ncomp);
#endif
        }

        amrex::ParallelDescriptor::Bcast(tmp.dataPtr(), tmp.size(), ioproc);

        // Shift box by the domain lower corner
        amrex::Box  fab_bx = tmp.box();
        amrex::Dim3 dom_lb = lbound(domain);
        fab_bx += amrex::IntVect(dom_lb.x,dom_lb.y,dom_lb.z);
        // fab_vars points to data on device
        fab_vars[iv]->resize(fab_bx,1);
#ifdef AMREX_USE_GPU
        amrex::Gpu::copy(amrex::Gpu::hostToDevice,
                         tmp.dataPtr(), tmp.dataPtr() + tmp.size(),
                         fab_vars[iv]->dataPtr());
#else
        // Provided by BaseFab inheritance through FArrayBox
        fab_vars[iv]->copy(tmp,tmp.box(),0,fab_bx,0,1);
#endif
    }
}

#endif
