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
                          Vector<amrex::FArrayBox*> fab_vars,
                          Vector<enum NC_Data_Dims_Type> NC_dim_types)
{
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        using RARRAY = NDArray<float>;
        amrex::Vector<RARRAY> arrays(nc_var_names.size());
        ReadWRFInputFile(fname, nc_var_names, arrays);

        for (int iv = 0; iv < nc_var_names.size(); iv++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[iv] << std::endl;

            /*
            switch (NC_dim_types[iv]) {
                case NC_Data_Dims_Type::Time_BT_SN_WE:
                    bigEnd[0] = arrays[iv].get_vshape()[3] - 1;
                    bigEnd[1] = arrays[iv].get_vshape()[2] - 1;
                    bigEnd[2] = arrays[iv].get_vshape()[1] - 1;
                    break;
                case NC_Data_Dims_Type::Time_SN_WE:
                    bigEnd[0] = arrays[iv].get_vshape()[2] - 1;
                    bigEnd[1] = arrays[iv].get_vshape()[1] - 1;
                    bigEnd[2] = 0;
                    break;
                case NC_Data_Dims_Type::Time_BT:
                    bigEnd[0] = 0;
                    bigEnd[1] = 0;
                    bigEnd[2] = arrays[iv].get_vshape()[1] - 1;
                    break;
                case NC_Data_Dims_Type::Time:
                    bigEnd[0] = 0;
                    bigEnd[1] = 0;
                    bigEnd[2] = 0;
                    break;
                default:
                    amrex::Error("Unrecognized NetCDF data dimensions type");
            }
            */

            // amrex::Print() << "FAB IN BUILD " << fab_vars[iv]->box() << std::endl;

            int ncomp  = 1;
            auto num_pts = fab_vars[iv]->box().numPts();

            // The two versions shown below are equivalent -- we keep both for clarity
#if 0
            for (int kcomp(0); kcomp < ncomp; ++kcomp)
            {
                auto dataPtr = fab_vars[iv]->dataPtr(kcomp);
                for (int n(0); n < num_pts; ++n) {
                    *(dataPtr+n) = static_cast<Real>(*(arrays[iv].get_data()+n));
                }
            }
#else
            Array4<Real> fab_arr = fab_vars[iv]->array();
            int ns2 = arrays[iv].get_vshape()[2];
            int ns3 = arrays[iv].get_vshape()[3];

            int ioff = fab_vars[iv]->box().smallEnd()[0];
            int joff = fab_vars[iv]->box().smallEnd()[1];

            for (int kcomp(0); kcomp < ncomp; ++kcomp)
            {
                for (int n(0); n < num_pts; ++n) {
                    int k  = n / (ns2*ns3);
                    int j  = (n - k*(ns2*ns3)) / ns3 + joff;
                    int i  = n - k*(ns2*ns3) - (j-joff) * ns3 + ioff;
                    fab_arr(i,j,k,kcomp) = static_cast<Real>(*(arrays[iv].get_data()+n));
                }
            }
#endif
        }
    } // if IOProcessor
}
