#include "AMReX_FArrayBox.H"
#include "NCWpsFile.H"
#include "IndexDefines.H"

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

//#define _XOPEN_SOURCE
// Converts UTC time string to a time_t value.
std::time_t getEpochTime(const std::string& dateTime, const std::string& dateTimeFormat)
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
    auto epoch = std::mktime(&tmTime);
    Print() << "Time Stamp: "<< std::put_time(&tmTime, "%c")
            << " , Epoch: " << epoch << std::endl;

    return epoch;
}

#ifdef ERF_USE_NETCDF
Real
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

    int ntimes;
    Real timeInterval;
    const std::string dateTimeFormat ="%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor())
    {
        // Read the netcdf file and fill these FABs

        Vector<std::string> nc_var_names;
        Vector<enum NC_Data_Dims_Type> NC_dim_types;

        nc_var_names.push_back("U_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("U_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("U_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);
        nc_var_names.push_back("U_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);

        nc_var_names.push_back("V_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("V_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("V_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);
        nc_var_names.push_back("V_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);

        nc_var_names.push_back("W_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("W_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("W_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);
        nc_var_names.push_back("W_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);

        nc_var_names.push_back("T_BXS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("T_BXE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_SN);
        nc_var_names.push_back("T_BYS")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);
        nc_var_names.push_back("T_BYE")     ; NC_dim_types.push_back(NC_Data_Dims_Type::Time_BdyWidth_BT_WE);

        using RARRAY = NDArray<float>;
        amrex::Vector<RARRAY> arrays(nc_var_names.size());
        ReadWRFBdyFile(nc_bdy_file, nc_var_names, arrays);

        // Read the time stamps
        using CharArray = NDArray<char>;
        amrex::Vector<CharArray> array_ts(1);
        ReadWRFBdyFile(nc_bdy_file, {"Times"}, array_ts);

        ntimes = array_ts[0].get_vshape()[0];
        auto dateStrLen = array_ts[0].get_vshape()[1];
        // auto numCharInTimes = ntimes*dateStrLen;
        char timeStamps[ntimes][dateStrLen];

        // Fill up the characters read
        int str_len = static_cast<int>(dateStrLen);
        for (int nt(0); nt < ntimes; nt++)
            for (int dateStrCt(0); dateStrCt < str_len; dateStrCt++) {
                auto n = nt*dateStrLen + dateStrCt;
                timeStamps[nt][dateStrCt] = *(array_ts[0].get_data() + n);
            }

        Vector<std::string> timeStampsString;
        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string date(&timeStamps[nt][0], &timeStamps[nt][dateStrLen-1]+1);
            auto epochTime = getEpochTime(date, dateTimeFormat);
            timeStampsString.push_back(date);
            epochTimes.push_back(epochTime);

            if (nt == 1)
                timeInterval = epochTimes[1] - epochTimes[0];
            else if (nt >= 1)
                AMREX_ALWAYS_ASSERT(epochTimes[nt] - epochTimes[nt-1] == timeInterval);
        }

        // Our outermost loop is time
        bdy_data_xlo.resize(ntimes);
        bdy_data_xhi.resize(ntimes);
        bdy_data_ylo.resize(ntimes);
        bdy_data_yhi.resize(ntimes);

        amrex::IntVect plo(lo);
        amrex::IntVect phi(hi);

        int nvars = nc_var_names.size();

        // This loops over every variable on every face, so nvars should be 4 * number of "ivartype" below
        for (int iv = 0; iv < nvars; iv++)
        {
            amrex::Print() << "Building FAB for the the NetCDF variable : " << nc_var_names[iv] << std::endl;

            std::string first1 = nc_var_names[iv].substr(0,1);
            int ivartype;

            if        (first1 == "U") {
                ivartype = WRFBdyVars::U;
            } else if (first1 == "V") {
                ivartype = WRFBdyVars::V;
            } else if (first1 == "W") {
                ivartype = WRFBdyVars::W;
            } else if (first1 == "T") {
                ivartype = WRFBdyVars::T;
            }

            // arrays[i].get_vshape()[0] : number of times
            // arrays[i].get_vshape()[1] : number of points in x-direction
            // arrays[i].get_vshape()[2] : number of points in z-direction
            // arrays[i].get_vshape()[3] : number of points in y-direction

            // Assert that all data has the same number of time snapshots
            int itimes = static_cast<int>(arrays[iv].get_vshape()[0]);
            AMREX_ALWAYS_ASSERT(itimes == ntimes);

            // amrex::Print() << "SHAPE0 " << arrays[iv].get_vshape()[0] << std::endl;
            // amrex::Print() << "SHAPE1 " << arrays[iv].get_vshape()[1] << std::endl;
            // amrex::Print() << "SHAPE2 " << arrays[iv].get_vshape()[2] << std::endl;
            // amrex::Print() << "SHAPE3 " << arrays[iv].get_vshape()[3] << std::endl;

            std::string  last3 = nc_var_names[iv].substr(nc_var_names[iv].size()-3, 3);

            // Width of the boundary region (1 <= ng <= 5)
            int ng = arrays[iv].get_vshape()[1];

            if (last3 == "BXS") {

                // *******************************************************************************
                // xlo bdy
                // *******************************************************************************
                plo[0] = -ng; plo[1] = lo[1]; plo[2] = lo[2];
                phi[0] = -1 ; phi[1] = hi[1]; phi[2] = hi[2];
                const Box pbx_xlo(plo, phi);

                Box xlo_plane_no_stag(pbx_xlo);
                //Box xlo_plane_x_stag = Box(IntVect(-ng+1,lo[1],lo[2]),IntVect(0,hi[1],hi[2]),{1,0,0});
                Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});
                Box xlo_plane_z_stag = convert(pbx_xlo, {0, 0, 1});

                if        (first1 == "U") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_no_stag, 1)); // U
                    }
                } else if (first1 == "V") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_y_stag , 1)); // V
                    }
                } else if (first1 == "W") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_z_stag , 1)); // W
                    }
                } else if (first1 == "T") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_no_stag, 1)); // T
                    }
                }

            } else if (last3 == "BXE") {

                // *******************************************************************************
                // xhi bdy
                // *******************************************************************************
                plo[0] = hi[0] +  1; plo[1] = lo[1]; plo[2] = lo[2];
                phi[0] = hi[0] + ng; phi[1] = hi[1]; phi[2] = hi[2];
                const Box pbx_xhi(plo, phi);

                Box xhi_plane_no_stag(pbx_xhi);
                //Box xhi_plane_x_stag = Box(IntVect(plo[0],lo[1],lo[2]),IntVect(phi[0],hi[1],hi[2]),{1,0,0});
                Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});
                Box xhi_plane_z_stag = convert(pbx_xhi, {0, 0, 1});

                if        (first1 == "U") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_no_stag, 1)); // U
                    }
                } else if (first1 == "V") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_y_stag , 1)); // V
                    }
                } else if (first1 == "W") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_z_stag , 1)); // W
                    }
                } else if (first1 == "T") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_no_stag, 1)); // T
                    }
                }

            } else if (last3 == "BYS") {

                // *******************************************************************************
                // ylo bdy
                // *******************************************************************************
                plo[1] = -ng; plo[0] = lo[0]; plo[2] = lo[2];
                phi[1] =  -1; phi[0] = hi[0]; phi[2] = hi[2];
                const Box pbx_ylo(plo, phi);

                Box ylo_plane_no_stag(pbx_ylo);
                Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
                //Box ylo_plane_y_stag = Box(IntVect(lo[0],-ng+1,lo[2]),IntVect(hi[0],0,hi[2]),{0,1,0});
                Box ylo_plane_z_stag = convert(pbx_ylo, {0, 0, 1});

                if        (first1 == "U") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_x_stag , 1)); // U
                    }
                } else if (first1 == "V") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_no_stag, 1)); // V
                    }
                } else if (first1 == "W") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_z_stag , 1)); // W
                    }
                } else if (first1 == "T") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_no_stag, 1)); // T
                    }
                }

            } else if (last3 == "BYE") {

                // *******************************************************************************
                // yhi bdy
                // *******************************************************************************
                plo[1] = hi[1] +  1; plo[0] = lo[0]; plo[2] = lo[2];
                phi[1] = hi[1] + ng; phi[0] = hi[0]; phi[2] = hi[2];
                const Box pbx_yhi(plo, phi);

                Box yhi_plane_no_stag(pbx_yhi);
                Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
                //Box yhi_plane_y_stag = Box(IntVect(lo[0],plo[1],lo[2]),IntVect(hi[0],phi[1],hi[2]),{0,1,0});
                Box yhi_plane_z_stag = convert(pbx_yhi, {0, 0, 1});

                if        (first1 == "U") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_x_stag , 1)); // U
                    }
                } else if (first1 == "V") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_no_stag, 1)); // V
                    }
                } else if (first1 == "W") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_z_stag , 1)); // W
                    }
                } else if (first1 == "T") {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_no_stag, 1)); // T
                    }
                }
            }

            long num_pts;

            int ns2 = arrays[iv].get_vshape()[2];
            int ns3 = arrays[iv].get_vshape()[3];

            // Now fill the data
            for (int nt(0); nt < ntimes; ++nt)
            {
                Array4<Real> fab_arr;
                if (last3 == "BXS") {
                    num_pts  = bdy_data_xlo[nt][ivartype].box().numPts();
                    fab_arr  = bdy_data_xlo[nt][ivartype].array();
                    int ioff = bdy_data_xlo[nt][ivartype].smallEnd()[0];
                    for (int n(0); n < num_pts; ++n) {
                        int i  = n / (ns2*ns3) + ioff;
                        int k  = (n - (i-ioff)*(ns2*ns3)) / ns3;
                        int j  =  n - (i-ioff)*(ns2*ns3) - k * ns3;
                        fab_arr(i,j,k,0) = static_cast<Real>(*(arrays[iv].get_data()+n));
                        if (nt == 0 and first1 == "T" and i == -1 and j == 0 and k == 0) amrex::Print() << "T AT XLO " << IntVect(i,j,k) << " " << n << " " <<
                            fab_arr(i,j,k,0) << std::endl;
                    }
                } else if (last3 == "BXE") {
                    num_pts  = bdy_data_xhi[nt][ivartype].box().numPts();
                    fab_arr  = bdy_data_xhi[nt][ivartype].array();
                    int ioff = bdy_data_xhi[nt][ivartype].bigEnd()[0]-ng+1;
                    for (int n(0); n < num_pts; ++n) {
                        int i  = n / (ns2*ns3) + ioff;
                        int k  = (n - (i-ioff)*(ns2*ns3)) / ns3;
                        int j  =  n - (i-ioff)*(ns2*ns3) - k * ns3;
                        fab_arr(i,j,k,0) = static_cast<Real>(*(arrays[iv].get_data()+n));
                    }
                } else if (last3 == "BYS") {
                    num_pts = bdy_data_ylo[nt][ivartype].box().numPts();
                    fab_arr  = bdy_data_ylo[nt][ivartype].array();
                    int joff = bdy_data_ylo[nt][ivartype].smallEnd()[1];
                    for (int n(0); n < num_pts; ++n) {
                        int j  = n / (ns2*ns3) + joff;
                        int k  = (n - (j-joff)*(ns2*ns3)) / ns3;
                        int i  =  n - (j-joff)*(ns2*ns3) - k * ns3;
                        fab_arr(i,j,k,0) = static_cast<Real>(*(arrays[iv].get_data()+n));
                        if (nt == 0 and first1 == "T" and i == 0 and j == -1 and k == 0) amrex::Print() << "T AT YLO " << IntVect(i,j,k) << " " << n << " " <<
                            fab_arr(i,j,k,0) << std::endl;
                    }
                } else if (last3 == "BYE") {
                    num_pts  = bdy_data_yhi[nt][ivartype].box().numPts();
                    fab_arr  = bdy_data_yhi[nt][ivartype].array();
                    int joff = bdy_data_yhi[nt][ivartype].bigEnd()[1]-ng+1;
                    for (int n(0); n < num_pts; ++n) {
                        int j  = n / (ns2*ns3) + joff;
                        int k  = (n - (j-joff)*(ns2*ns3)) / ns3;
                        int i  =  n - (j-joff)*(ns2*ns3) - k * ns3;
                        fab_arr(i,j,k,0) = static_cast<Real>(*(arrays[iv].get_data()+n));
                    }
                }
            } // nt
        } //nc_var_names
    } // if ParalleDescriptor::IOProcessor()

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // Make sure all processors know how timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

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

        // TODO: WE NEED TO MULTIPLY THESE BY RHO TO GET (RHO THETA)
        // bdy_data_xlo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_xhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_ylo[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
        // bdy_data_yhi[nt][3].template mult<RunOn::Device>(NC_rho_fab[lev][idx],0,0,1);
    }

    amrex::Print() << "Successfully loaded data from the NetCDF file" << std::endl << std::endl;

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}
#endif // ERF_USE_NETCDF
