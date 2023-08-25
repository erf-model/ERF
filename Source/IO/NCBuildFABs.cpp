#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "DataStruct.H"
#include "NCInterface.H"
#include "NCWpsFile.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_IndexType.H"
#include "AMReX_Print.H"

using namespace amrex;

using RARRAY = NDArray<float>;

void
fill_fab_from_arrays(int iv, Vector<RARRAY>& nc_arrays,
                     const std::string& var_name,
                     NC_Data_Dims_Type& NC_dim_type,
                     FArrayBox& temp);

/**
 * Function to read NetCDF variables and fill the corresponding Array4's
 *
 * @param fname Name of the NetCDF file to be read
 * @param nc_var_names Variable names in the NetCDF file
 * @param NC_dim_types NetCDF data dimension types
 * @param fab_vars Fab data we are to fill
 */
void
BuildFABsFromNetCDFFile(const Box& domain,
                        const std::string &fname,
                        Vector<std::string> nc_var_names,
                        Vector<enum NC_Data_Dims_Type> NC_dim_types,
                        Vector<amrex::FArrayBox*> fab_vars)
{
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    amrex::Vector<RARRAY> nc_arrays(nc_var_names.size());

    if (amrex::ParallelDescriptor::IOProcessor())
    {
        ReadNetCDFFile(fname, nc_var_names, nc_arrays);
    }

    for (int iv = 0; iv < nc_var_names.size(); iv++)
    {
        amrex::FArrayBox tmp;
        if (amrex::ParallelDescriptor::IOProcessor()) {
            fill_fab_from_arrays(iv, nc_arrays, nc_var_names[iv], NC_dim_types[iv], tmp);
        }

        Box box = tmp.box();
        int ncomp = tmp.nComp();

        ParallelDescriptor::Bcast(&box, 1, ioproc);
        ParallelDescriptor::Bcast(&ncomp, 1, ioproc);

        if (!amrex::ParallelDescriptor::IOProcessor()) {
#ifdef AMREX_USE_GPU
            tmp.resize(box,ncomp, The_Pinned_Arena());
#else
            tmp.resize(box,ncomp);
#endif
        }

        ParallelDescriptor::Bcast(tmp.dataPtr(), tmp.size(), ioproc);

        // Shift box by the domain lower corner
        Box  fab_bx = tmp.box();
        Dim3 dom_lb = lbound(domain);
        fab_bx += IntVect(dom_lb.x,dom_lb.y,dom_lb.z);
        // fab_vars points to data on device
        fab_vars[iv]->resize(fab_bx,1);
#ifdef AMREX_USE_GPU
        Gpu::copy(Gpu::hostToDevice, tmp.dataPtr(), tmp.dataPtr() + tmp.size(), fab_vars[iv]->dataPtr());
#else
        // Provided by BaseFab inheritance through FArrayBox
        fab_vars[iv]->copy(tmp,tmp.box(),0,fab_bx,0,1);
#endif
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
void
fill_fab_from_arrays(int iv, Vector<RARRAY>& nc_arrays,
                     const std::string& var_name,
                     NC_Data_Dims_Type& NC_dim_type,
                     FArrayBox& temp)
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
    Box my_box(IntVect(0,0,0), IntVect(ns3-1,ns2-1,ns1-1));
    // amrex::Print() <<" MY BOX " << my_box << std::endl;

    if (var_name == "U" || var_name == "UU" ||
        var_name == "MAPFAC_U" || var_name == "MAPFAC_UY") my_box.setType(amrex::IndexType(IntVect(1,0,0)));
    if (var_name == "V" || var_name == "VV" ||
        var_name == "MAPFAC_V" || var_name == "MAPFAC_VY") my_box.setType(amrex::IndexType(IntVect(0,1,0)));
    if (var_name == "W" || var_name == "WW") my_box.setType(amrex::IndexType(IntVect(0,0,1)));

#ifdef AMREX_USE_GPU
    // Make sure temp lives on CPU since nc_arrays lives on CPU only
    temp.resize(my_box,1,The_Pinned_Arena());
#else
    temp.resize(my_box,1);
#endif
    Array4<Real> fab_arr = temp.array();

    int ioff = temp.box().smallEnd()[0];
    int joff = temp.box().smallEnd()[1];

    auto num_pts = my_box.numPts();

    // amrex::Print() <<" ns1 * ns2 * ns3 " << ns1 * ns2 * ns3 << std::endl;
    // amrex::Print() <<" NUMPTS " << num_pts << std::endl;

    for (int n(0); n < num_pts; ++n) {
        int k  = n / (ns2*ns3);
        int j  = (n - k*(ns2*ns3)) / ns3 + joff;
        int i  = n - k*(ns2*ns3) - (j-joff) * ns3 + ioff;
        fab_arr(i,j,k,0) = static_cast<Real>(*(nc_arrays[iv].get_data()+n));
    }
}
