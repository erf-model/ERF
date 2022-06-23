#include "AMReX_FArrayBox.H"
#include "NCWpsFile.H"

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

//TODO: Move this function to a suitable file

// Converts UTC time string to a time_t value.
std::time_t getEpochTime(const std::wstring& dateTime)
{
    // Let's consider we are getting all the input in
    // this format: '2014-07-25_20:17:22'
    // A better approach would be to pass in the format as well.
    static const std::wstring dateTimeFormat{ L"%Y-%m-%d_%H:%M:%S" };

    // Create a stream which we will use to parse the string,
    // which we provide to constructor of stream to fill the buffer.
    std::wistringstream ss{ dateTime };

    // Create a tm object to store the parsed date and time.
    std::tm dt;

    // Now we read from buffer using get_time manipulator
    // and formatting the input appropriately.
    ss >> std::get_time(&dt, dateTimeFormat.c_str());

    // Convert the tm structure to time_t value and return.
    return std::mktime(&dt);
}

#ifdef ERF_USE_NETCDF
void
read_from_wrfbdy(std::string nc_bdy_file, const Box& domain,
                 Vector<Vector<FArrayBox>>& bdy_data_xlo,
                 Vector<Vector<FArrayBox>>& bdy_data_xhi,
                 Vector<Vector<FArrayBox>>& bdy_data_ylo,
                 Vector<Vector<FArrayBox>>& bdy_data_yhi)
{
    amrex::Print() << "Loading boundary data from NetCDF file " << std::endl;
    // *********************************************************
    // Allocate space for all of the boundary planes we may need
    // Here we make only one enough space for one time -- we will
    //    add space for the later time slices later
    // *********************************************************
    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

    amrex::IntVect plo(lo);
    amrex::IntVect phi(hi);

    /*
     bdy_data_xlo[time] (etc) contain 4 different variables U_BXS, V_BXS, W_BXS, T_BXS.
     The dimensions of these on xlo and xhi are:
     U_BXS -> dim -> 1     NY       NZ
     V_BXS -> dim -> 1    (NY+1)    NZ
     W_BXS -> dim -> 1     NY      (NZ+1)
     T_BXS -> dim -> 1     NY       NZ

     The dimensions of these on ylo and yhi are:
     U_BYS -> dim -> 1    (NX+1)    NZ
     V_BYS -> dim -> 1     NX       NZ
     W_BYS -> dim -> 1     NX      (NZ+1)
     T_BYS -> dim -> 1     NX       NZ
    */

    // *******************************************************************************
    // First allocate space for just one time
    // *******************************************************************************
    bdy_data_xlo.resize(1);
    bdy_data_xhi.resize(1);
    bdy_data_ylo.resize(1);
    bdy_data_yhi.resize(1);

    // *******************************************************************************
    // xlo bdy
    // *******************************************************************************
    plo[0] = -1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = -1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xlo(plo, phi);

    Box xlo_plane_no_stag(pbx_xlo);
    Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});
    Box xlo_plane_z_stag = convert(pbx_xlo, {0, 0, 1});

    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_no_stag, 1)); // U
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_y_stag , 1)); // V
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_z_stag , 1)); // W
    bdy_data_xlo[0].push_back(FArrayBox(xlo_plane_no_stag, 1)); // T

    // *******************************************************************************
    // xhi bdy
    // *******************************************************************************
    plo[0] = hi[0] + 1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0] + 1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xhi(plo, phi);

    Box xhi_plane_no_stag(pbx_xhi);
    Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});
    Box xhi_plane_z_stag = convert(pbx_xhi, {0, 0, 1});

    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_no_stag, 1)); // U
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_y_stag , 1)); // V
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_z_stag , 1)); // W
    bdy_data_xhi[0].push_back(FArrayBox(xhi_plane_no_stag, 1)); // T

    // *******************************************************************************
    // ylo bdy
    // *******************************************************************************
    plo[1] = -1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = -1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_ylo(plo, phi);

    Box ylo_plane_no_stag(pbx_ylo);
    Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
    Box ylo_plane_z_stag = convert(pbx_ylo, {0, 0, 1});

    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_x_stag , 1)); // U
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_no_stag, 1)); // V
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_z_stag , 1)); // W
    bdy_data_ylo[0].push_back(FArrayBox(ylo_plane_no_stag, 1)); // T

    // *******************************************************************************
    // yhi bdy
    // *******************************************************************************
    plo[1] = hi[1] + 1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = hi[1] + 1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_yhi(plo, phi);

    Box yhi_plane_no_stag(pbx_yhi);
    Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
    Box yhi_plane_z_stag = convert(pbx_yhi, {0, 0, 1});

    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_x_stag , 1)); // U
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_no_stag, 1)); // V
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_z_stag , 1)); // W
    bdy_data_yhi[0].push_back(FArrayBox(yhi_plane_no_stag, 1)); // T

    // *******************************************************************************

    int ntimes;
    if (ParallelDescriptor::IOProcessor())
    {
        // Read the netcdf file and fill these FABs

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
        ReadWRFBdyFile(nc_bdy_file, nc_var_names, arrays);

        // Read the time stamps
        using CharArray = NDArray<char>;
        amrex::Vector<CharArray> array_ts(1);
        ReadWRFBdyFile(nc_bdy_file, {"Times"}, array_ts);

        //ntimes = arrays[0].get_vshape()[0];
        ntimes = array_ts[0].get_vshape()[0];
        auto dateStrLen = array_ts[0].get_vshape()[1];
        auto numCharInTimes = ntimes*dateStrLen;
        char timeStamps[ntimes][dateStrLen];

        // Fill up the characters read
        for (int nt(0); nt < ntimes; nt++)
            for (int dateStrCt(0); dateStrCt < dateStrLen; dateStrCt++) {
                auto n = nt*dateStrLen + dateStrCt;
                timeStamps[nt][dateStrCt] = *(array_ts[0].get_data() + n);
            }

        Vector<std::string> timeStampsString;
        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string str(&timeStamps[nt][0], &timeStamps[nt][dateStrLen-1]+1);
            std::wstring wstr(str.begin(), str.end());
            auto epochTime = getEpochTime(wstr);
            timeStampsString.push_back(str);
            epochTimes.push_back(epochTime);
        }

        int nvars = nc_var_names.size();
        // Our outermost loop is time
        bdy_data_xlo.resize(ntimes);
        bdy_data_xhi.resize(ntimes);
        bdy_data_ylo.resize(ntimes);
        bdy_data_yhi.resize(ntimes);

        for (int i = 0; i < nvars; i++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[i] << std::endl;

            // arrays[i].get_vshape()[0] : number of times
            // arrays[i].get_vshape()[1] : number of points in x-direction
            // arrays[i].get_vshape()[2] : number of points in z-direction
            // arrays[i].get_vshape()[3] : number of points in y-direction

            // Assert that all data has the same number of time snapshots
            AMREX_ALWAYS_ASSERT(arrays[i].get_vshape()[0] == ntimes);

            // amrex::Print() << "SHAPE 0 " << arrays[i].get_vshape()[0] << std::endl;
            // amrex::Print() << "SHAPE 1 " << arrays[i].get_vshape()[1] << std::endl;
            // amrex::Print() << "SHAPE 2 " << arrays[i].get_vshape()[2] << std::endl;
            // amrex::Print() << "SHAPE 3 " << arrays[i].get_vshape()[3] << std::endl;

            long num_pts;
            Real* data_ptr;

            std::string first1 = nc_var_names[i].substr(0,1);
            int ivartype;
            if      (first1 == "U") ivartype = 0;
            else if (first1 == "V") ivartype = 1;
            else if (first1 == "W") ivartype = 2;
            else if (first1 == "T") ivartype = 3;

            std::string  last3 = nc_var_names[i].substr(nc_var_names[i].size()-3, 3);

            for (int nt(0); nt < ntimes; ++nt)
            {
                bdy_data_xlo[nt].push_back(FArrayBox(bdy_data_xlo[0][i].box(),1));
                bdy_data_xhi[nt].push_back(FArrayBox(bdy_data_xhi[0][i].box(),1));
                bdy_data_ylo[nt].push_back(FArrayBox(bdy_data_ylo[0][i].box(),1));
                bdy_data_yhi[nt].push_back(FArrayBox(bdy_data_yhi[0][i].box(),1));

                if (last3 == "BXS") {
                    num_pts  = bdy_data_xlo[nt][ivartype].box().numPts();
                    data_ptr = bdy_data_xlo[nt][ivartype].dataPtr();
                } else if (last3 == "BXE") {
                    num_pts  = bdy_data_xhi[nt][ivartype].box().numPts();
                    data_ptr = bdy_data_xhi[nt][ivartype].dataPtr();
                } else if (last3 == "BYS") {
                    num_pts  = bdy_data_ylo[nt][ivartype].box().numPts();
                    data_ptr = bdy_data_ylo[nt][ivartype].dataPtr();
                } else if (last3 == "BYE") {
                    num_pts  = bdy_data_yhi[nt][ivartype].box().numPts();
                    data_ptr = bdy_data_yhi[nt][ivartype].dataPtr();
                }

                //
                // Here we put the components into the array in the order they are read
                //
                for (int n(0); n < num_pts; ++n) {
                    *(data_ptr+n) = static_cast<Real>(*(arrays[i].get_data()+n));
                }
            } // nt
        } //nc_var_names
    } // if ParalleDescriptor::IOProcessor()

    amrex::Print() << "NTIMES " << ntimes << std::endl;

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int nvars = bdy_data_xlo[0].size();
    for (int nt = 0; nt < ntimes; nt++)
    {
        for (int i = 0; i < nvars; i++)
        {
            ParallelDescriptor::Bcast(bdy_data_xlo[nt][i].dataPtr(),bdy_data_xlo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_xhi[nt][i].dataPtr(),bdy_data_xhi[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_ylo[nt][i].dataPtr(),bdy_data_ylo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_yhi[nt][i].dataPtr(),bdy_data_yhi[nt][i].box().numPts(),ioproc);
        }

        // CONVERT (THETA - 300) to (RHO THETA)
        amrex::Real theta_ref = 300.;
        bdy_data_xlo[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_xhi[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_ylo[nt][3].template plus<RunOn::Device>(theta_ref);
        bdy_data_yhi[nt][3].template plus<RunOn::Device>(theta_ref);

        // Now multiply by rho to get (rho theta) instead of theta
        // bdy_data_xlo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_xhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_ylo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_yhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
    }

    amrex::Print() << "Successfully loaded data from the wrfbdy (output of 'real.exe') NetCDF file" << std::endl << std::endl;
}
#endif // ERF_USE_NETCDF
