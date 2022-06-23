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

        for (int i = 0; i < nc_var_names.size(); i++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[i] << std::endl;

            // Build the box first using the shape information
            amrex::IntVect smallEnd{0, 0, 0};
            amrex::IntVect bigEnd{0, 0, 0};

            switch (NC_dim_types[i]) {
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
                    bigEnd[2] = arrays[i].get_vshape()[1] - 1;
                    break;
                case NC_Data_Dims_Type::Time:
                    bigEnd[0] = 0;
                    bigEnd[1] = 0;
                    bigEnd[2] = 0;
                    break;
                default:
                    amrex::Error("Unrecognized NetCDF data dimensions type");
            }

            amrex::IntVect boxType{0, 0, 0};
            if (nc_var_names[i] == "U") boxType = amrex::IntVect(1,0,0);
            if (nc_var_names[i] == "V") boxType = amrex::IntVect(0,1,0);
            if (nc_var_names[i] == "W" || nc_var_names[i] == "PH" || nc_var_names[i] == "PHB") {
                boxType = amrex::IntVect(0,0,1);
            }
            amrex::Box bx = amrex::Box(smallEnd, bigEnd, boxType);

            AMREX_ALWAYS_ASSERT(bx == fab_vars[i]->box());

            int ncomp  = 1;
            auto num_pts = fab_vars[i]->box().numPts();
            for (int k(0); k < ncomp; ++k) {
                auto dataPtr = fab_vars[i]->dataPtr(k);
                for (int n(0); n < num_pts; ++n) {
                    *(dataPtr+n) = static_cast<Real>(*(arrays[i].get_data()+n));
                }
            }
        }
    } // if IOProcessor
}
