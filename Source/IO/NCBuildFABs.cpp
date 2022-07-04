#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "DataStruct.H"
#include "NCInterface.H"
#include "NCWpsFile.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Print.H"

using namespace amrex;

using PlaneVector = amrex::Vector<amrex::FArrayBox>;

// Function to read a NetCDF variable and fill a corresponding MultiFab and Array4
// fname is the "wrfinput_d01" resulting from running ideal.exe or real.exe
void
BuildFABsFromWRFInputFile(const std::string &fname,
                          Vector<std::string> nc_var_names,
                          Vector<enum NC_Data_Dims_Type> NC_dim_types,
                          Vector<amrex::FArrayBox*> fab_vars)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        using RARRAY = NDArray<float>;
        amrex::Vector<RARRAY> arrays(nc_var_names.size());

        ReadWRFFile(fname, nc_var_names, arrays);

        for (int iv = 0; iv < nc_var_names.size(); iv++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[iv] << std::endl;
            // amrex::Print() << "  with box " << fab_vars[iv]->box() << std::endl;

            Array4<Real> fab_arr = fab_vars[iv]->array();

            int ioff = fab_vars[iv]->box().smallEnd()[0];
            int joff = fab_vars[iv]->box().smallEnd()[1];

            auto bx = fab_vars[iv]->box();
            auto num_pts = bx.numPts();

            int ns2, ns3;
            if (NC_dim_types[iv] == NC_Data_Dims_Type::Time_BT) {
                ns2 = 1;
                ns3 = 1;
            } else if (NC_dim_types[iv] == NC_Data_Dims_Type::Time_SN_WE) {
                ns2 = arrays[iv].get_vshape()[1];
                ns3 = arrays[iv].get_vshape()[2];
                AMREX_ALWAYS_ASSERT(num_pts == ns2 * ns3);
            } else if (NC_dim_types[iv] == NC_Data_Dims_Type::Time_BT_SN_WE) {
                ns2 = arrays[iv].get_vshape()[2];
                ns3 = arrays[iv].get_vshape()[3];
            } else {
                amrex::Abort("Dont know this NC_Data_Dims_Type");
            }

            for (int n(0); n < num_pts; ++n) {
                int k  = n / (ns2*ns3);
                int j  = (n - k*(ns2*ns3)) / ns3 + joff;
                int i  = n - k*(ns2*ns3) - (j-joff) * ns3 + ioff;
                fab_arr(i,j,k,0) = static_cast<Real>(*(arrays[iv].get_data()+n));
            }
        } // iv
    } // if IOProcessor
}
