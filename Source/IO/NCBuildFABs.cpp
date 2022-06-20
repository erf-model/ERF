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

int
BuildFABsFromWRFBdyFile(const std::string &fname,
                        amrex::Vector<amrex::Vector<FArrayBox>>& bdy_data_xlo,
                        amrex::Vector<amrex::Vector<FArrayBox>>& bdy_data_xhi,
                        amrex::Vector<amrex::Vector<FArrayBox>>& bdy_data_ylo,
                        amrex::Vector<amrex::Vector<FArrayBox>>& bdy_data_yhi)
{
    int ntimes;
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        // We only read in data at one resolution, which we assume to be that of level 0
        //int lev = 0;

        Vector<std::string> nc_var_names;
        Vector<enum NC_Data_Dims_Type> NC_dim_types;
        Vector<enum Stagger> staggerDir;

        nc_var_names.push_back("U_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("U_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("U_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::x);
        nc_var_names.push_back("U_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::x);

        nc_var_names.push_back("V_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::y);
        nc_var_names.push_back("V_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::y);
        nc_var_names.push_back("V_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("V_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::None);

        nc_var_names.push_back("W_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::z);
        nc_var_names.push_back("W_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::z);
        nc_var_names.push_back("W_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::z);
        nc_var_names.push_back("W_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::z);

        nc_var_names.push_back("T_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("T_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("T_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::None);
        nc_var_names.push_back("T_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);   staggerDir.push_back(Stagger::None);

        using RARRAY = NDArray<float>;
        amrex::Vector<RARRAY> arrays(nc_var_names.size());
        ReadWRFBdyFile(fname, nc_var_names, arrays);

        using CharArray = NDArray<char*>;
        amrex::Vector<CharArray> array_ts(1);
        //ReadWRFBdyFile(fname, {"Times"}, array_ts);

        ntimes = arrays[0].get_vshape()[0];

        // Our outermost loop is time
        bdy_data_xlo.resize(ntimes);
        bdy_data_xhi.resize(ntimes);
        bdy_data_ylo.resize(ntimes);
        bdy_data_yhi.resize(ntimes);

        for (int n = 1; n < ntimes; n++) {

            bdy_data_xlo[n].push_back(FArrayBox(bdy_data_xlo[0][0].box(),1)); // U
            bdy_data_xlo[n].push_back(FArrayBox(bdy_data_xlo[0][1].box(),1)); // V
            bdy_data_xlo[n].push_back(FArrayBox(bdy_data_xlo[0][2].box(),1)); // W
            bdy_data_xlo[n].push_back(FArrayBox(bdy_data_xlo[0][3].box(),1)); // T

            Vector<FArrayBox> xhi_plane_data;
            bdy_data_xhi[n].push_back(FArrayBox(bdy_data_xhi[0][0].box(),1)); // U
            bdy_data_xhi[n].push_back(FArrayBox(bdy_data_xhi[0][1].box(),1)); // V
            bdy_data_xhi[n].push_back(FArrayBox(bdy_data_xhi[0][2].box(),1)); // W
            bdy_data_xhi[n].push_back(FArrayBox(bdy_data_xhi[0][3].box(),1)); // T

            Vector<FArrayBox> ylo_plane_data;
            bdy_data_ylo[n].push_back(FArrayBox(bdy_data_ylo[0][0].box(),1)); // U
            bdy_data_ylo[n].push_back(FArrayBox(bdy_data_ylo[0][1].box(),1)); // V
            bdy_data_ylo[n].push_back(FArrayBox(bdy_data_ylo[0][2].box(),1)); // W
            bdy_data_ylo[n].push_back(FArrayBox(bdy_data_ylo[0][3].box(),1)); // T

            Vector<FArrayBox> yhi_plane_data;
            bdy_data_yhi[n].push_back(FArrayBox(bdy_data_yhi[0][0].box(),1)); // U
            bdy_data_yhi[n].push_back(FArrayBox(bdy_data_yhi[0][1].box(),1)); // V
            bdy_data_yhi[n].push_back(FArrayBox(bdy_data_yhi[0][2].box(),1)); // W
            bdy_data_yhi[n].push_back(FArrayBox(bdy_data_yhi[0][3].box(),1)); // T
        }

        for (int i = 0; i < nc_var_names.size(); i++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[i] << std::endl;

            auto shape_data = arrays[i].get_vshape();
            auto currentVar = nc_var_names[i];

            // Build the box first using the shape information
            amrex::IntVect smallEnd{0, 0, 0};
            amrex::IntVect bigEnd{0, 0, 0};

            // get_vshape()[0] : number of times
            // get_vshape()[1] : number of points in x-direction
            // get_vshape()[2] : number of points in z-direction
            // get_vshape()[3] : number of points in y-direction

            // Assert that all data has the same number of time snapshots
            AMREX_ALWAYS_ASSERT(arrays[i].get_vshape()[0] == ntimes);

            //
            // NOTE: we always define the box here to be on the low side, but this won't matter
            //       when we copy data into the FABs below.  It is only the shape that matters.
            //
            switch (NC_dim_types[i]) {
                case NC_Data_Dims_Type::Time_BdyWidth_BT_SN: // lo-i and hi-i
                    // amrex::Print() << "TYPE IS Time_BdyWidth_BT_SN " << std::endl;
                    // amrex::Print() << "SHAPE 0 " << arrays[i].get_vshape()[0] << std::endl;
                    // amrex::Print() << "SHAPE 1 " << arrays[i].get_vshape()[1] << std::endl;
                    // amrex::Print() << "SHAPE 2 " << arrays[i].get_vshape()[2] << std::endl;
                    // amrex::Print() << "SHAPE 3 " << arrays[i].get_vshape()[3] << std::endl;
                    smallEnd[0] = -1;
                    bigEnd[0]   = -1;
                    bigEnd[1] = arrays[i].get_vshape()[3] - 1;
                    bigEnd[2] = arrays[i].get_vshape()[2] - 1;
                    break;
                case NC_Data_Dims_Type::Time_BdyWidth_BT_WE: // lo-j and hi-j
                    // amrex::Print() << "TYPE IS Time_BdyWidth_BT_WE " << std::endl;
                    // amrex::Print() << "SHAPE 0 " << arrays[i].get_vshape()[0] << std::endl;
                    // amrex::Print() << "SHAPE 1 " << arrays[i].get_vshape()[1] << std::endl;
                    // amrex::Print() << "SHAPE 2 " << arrays[i].get_vshape()[2] << std::endl;
                    // amrex::Print() << "SHAPE 3 " << arrays[i].get_vshape()[3] << std::endl;
                    bigEnd[0] = arrays[i].get_vshape()[1] - 1;
                    bigEnd[0] = arrays[i].get_vshape()[3] - 1;
                    smallEnd[1] = -1;
                    bigEnd[1]   = -1;
                    bigEnd[2] = arrays[i].get_vshape()[2] - 1;
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

            int tot_size = shape_data[0]*shape_data[1]*shape_data[2]*shape_data[3];
            // amrex::Print() << "Total points in the box constructed from netCDF variable: " << tot_size << std::endl;

            //amrex::IntVect boxType{0, 0, 0};
            //amrex::Box bx = amrex::Box(smallEnd, bigEnd, boxType);
            // amrex::Print() << "BX OF INPUT DATA " << bx << std::endl;

            long num_pts;
            Real* data_ptr;

            std::string first1 = nc_var_names[i].substr(0,1);
            int ifab;
            if      (first1 == "U") ifab = 0;
            else if (first1 == "V") ifab = 1;
            else if (first1 == "W") ifab = 2;
            else if (first1 == "T") ifab = 3;

            std::string  last3 = nc_var_names[i].substr(nc_var_names[i].size()-3, 3);

            for (int nt(0); nt < ntimes; ++nt)
            {
                if (last3 == "BXS") {
                    num_pts  = bdy_data_xlo[nt][ifab].box().numPts();
                    data_ptr = bdy_data_xlo[nt][ifab].dataPtr();
                } else if (last3 == "BXE") {
                    num_pts  = bdy_data_xhi[nt][ifab].box().numPts();
                    data_ptr = bdy_data_xhi[nt][ifab].dataPtr();
                } else if (last3 == "BYS") {
                    num_pts  = bdy_data_ylo[nt][ifab].box().numPts();
                    data_ptr = bdy_data_ylo[nt][ifab].dataPtr();
                } else if (last3 == "BYE") {
                    num_pts  = bdy_data_yhi[nt][ifab].box().numPts();
                    data_ptr = bdy_data_yhi[nt][ifab].dataPtr();
                }

                //
                // Here we put the components into the array in the order they are read
                //
                for (int n(0); n < num_pts; ++n) {
                    *(data_ptr+n) = static_cast<Real>(*(arrays[i].get_data()+n));
                }
            } // nt
        } //nc_var_names
    } // if IOProcessor
    return ntimes;
}
