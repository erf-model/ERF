#include "AMReX_FArrayBox.H"
#include "NCWpsFile.H"
#include "IndexDefines.H"
#include "EOS.H"

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "DataStruct.H"
#include "NCInterface.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Print.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

namespace WRFBdyTypes {
    enum {
        x_lo,
        x_hi,
        y_lo,
        y_hi
    };
}

Real
read_from_wrfbdy (const std::string& nc_bdy_file, const Box& domain,
                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                  Vector<Vector<FArrayBox>>& bdy_data_yhi,
                  int& width, Real& start_bdy_time)
{
    amrex::Print() << "Loading boundary data from NetCDF file " << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

    // *******************************************************************************

    int ntimes;
    Real timeInterval;
    const std::string dateTimeFormat ="%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor())
    {
        // Read the time stamps
        using CharArray = NDArray<char>;
        amrex::Vector<CharArray> array_ts(1);
        ReadNetCDFFile(nc_bdy_file, {"Times"}, array_ts);

        ntimes = array_ts[0].get_vshape()[0];
        auto dateStrLen = array_ts[0].get_vshape()[1];
        char timeStamps[ntimes][dateStrLen];

        // Fill up the characters read
        int str_len = static_cast<int>(dateStrLen);
        for (int nt(0); nt < ntimes; nt++)
            for (int dateStrCt(0); dateStrCt < str_len; dateStrCt++) {
                auto n = nt*dateStrLen + dateStrCt;
                timeStamps[nt][dateStrCt] = *(array_ts[0].get_data() + n);
            }

        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string date(&timeStamps[nt][0], &timeStamps[nt][dateStrLen-1]+1);
            auto epochTime = getEpochTime(date, dateTimeFormat);
            epochTimes.push_back(epochTime);

            if (nt == 1) {
                timeInterval = epochTimes[1] - epochTimes[0];
            } else if (nt >= 1) {
                AMREX_ALWAYS_ASSERT(epochTimes[nt] - epochTimes[nt-1] == timeInterval);
            }
        }
        start_bdy_time = epochTimes[0];
    }

    ParallelDescriptor::Bcast(&start_bdy_time,1,ioproc);
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Even though we may not read in all the variables, we need to make the arrays big enough for them (for now)
    int nvars = WRFBdyVars::NumTypes*4;

    // Our outermost loop is time
    bdy_data_xlo.resize(ntimes);
    bdy_data_xhi.resize(ntimes);
    bdy_data_ylo.resize(ntimes);
    bdy_data_yhi.resize(ntimes);

    amrex::IntVect plo(lo);
    amrex::IntVect phi(hi);

    // ******************************************************************
    // Read the netcdf file and fill these FABs
    // NOTE: the order and number of these must match the WRFBdyVars enum!
    // WRFBdyVars:  U, V, R, T, QV, MU, PC
    // ******************************************************************
    Vector<std::string> nc_var_names;
    Vector<std::string> nc_var_prefix = {"U","V","R","T","QVAPOR","MU","PC"};

    for (int ip = 0; ip < nc_var_prefix.size(); ++ip)
    {
       nc_var_names.push_back(nc_var_prefix[ip] + "_BXS");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BXE");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BYS");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BYE");
    }

    using RARRAY = NDArray<float>;
    amrex::Vector<RARRAY> arrays(nc_var_names.size());

    if (ParallelDescriptor::IOProcessor())
    {
        ReadNetCDFFile(nc_bdy_file, nc_var_names, arrays);

        // Assert that the data has the same number of time snapshots
        int itimes = static_cast<int>(arrays[0].get_vshape()[0]);
        AMREX_ALWAYS_ASSERT(itimes == ntimes);

        // Width of the boundary region
        width = arrays[0].get_vshape()[1];

        AMREX_ALWAYS_ASSERT(1 <= width && width <= 5);
    }
    ParallelDescriptor::Bcast(&width,1,ioproc);

    // This loops over every variable on every face, so nvars should be 4 * number of "ivartype" below
    for (int iv = 0; iv < nvars; iv++)
    {
        amrex::Print() << "Building FAB for the NetCDF variable : " << nc_var_names[iv] << std::endl;

        int bdyVarType;

        std::string first1 = nc_var_names[iv].substr(0,1);
        std::string first2 = nc_var_names[iv].substr(0,2);

        if        (first1 == "U") {
            bdyVarType = WRFBdyVars::U;
        } else if (first1 == "V") {
            bdyVarType = WRFBdyVars::V;
        } else if (first1 == "R") {
            bdyVarType = WRFBdyVars::R;
        } else if (first1 == "T") {
            bdyVarType = WRFBdyVars::T;
        } else if (first2 == "QV") {
            bdyVarType = WRFBdyVars::QV;
        } else if (first2 == "MU") {
            bdyVarType = WRFBdyVars::MU;
        } else if (first2 == "PC") {
            bdyVarType = WRFBdyVars::PC;
        } else {
            amrex::Print() << "Trying to read " << first1 << " or " << first2 << std::endl;
            amrex::Abort("dont know this variable");
        }

        std::string  last3 = nc_var_names[iv].substr(nc_var_names[iv].size()-3, 3);
        int bdyType;

        if        (last3 == "BXS") {
            bdyType = WRFBdyTypes::x_lo;
        } else if (last3 == "BXE") {
            bdyType = WRFBdyTypes::x_hi;
        } else if (last3 == "BYS") {
            bdyType = WRFBdyTypes::y_lo;
        } else if (last3 == "BYE") {
            bdyType = WRFBdyTypes::y_hi;
        }

        Box y_line_no_stag(IntVect(lo[0], 0, 0), IntVect(hi[0], 0, 0));

        Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
        Arena_Used = The_Pinned_Arena();
#endif

        if (bdyType == WRFBdyTypes::x_lo) {

                // *******************************************************************************
                // xlo bdy
                // *******************************************************************************
                plo[0] = lo[0]        ; plo[1] = lo[1]; plo[2] = lo[2];
                phi[0] = lo[0]+width-1; phi[1] = hi[1]; phi[2] = hi[2];
                const Box pbx_xlo(plo, phi);

                Box xlo_plane_no_stag(pbx_xlo);
                Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
                Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});

                Box xlo_line(IntVect(lo[0], lo[1], 0), IntVect(lo[0]+width-1, hi[1], 0));

                if        (bdyVarType == WRFBdyVars::U) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_x_stag, 1, Arena_Used)); // U
                    }
                } else if (bdyVarType == WRFBdyVars::V) {
                    for (int nt(0); nt < ntimes; ++nt) {
                        bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_y_stag , 1, Arena_Used)); // V
                    }
                } else if (bdyVarType == WRFBdyVars::R) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // R
                    }
                } else if (bdyVarType == WRFBdyVars::T) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // T
                    }
                } else if (bdyVarType == WRFBdyVars::QV) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xlo[nt].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // QV
                    }
                } else if (bdyVarType == WRFBdyVars::MU ||
                           bdyVarType == WRFBdyVars::PC) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xlo[nt].push_back(FArrayBox(xlo_line, 1, Arena_Used));
                    }
                }

            } else if (bdyType == WRFBdyTypes::x_hi) {

                // *******************************************************************************
                // xhi bdy
                // *******************************************************************************
                plo[0] = hi[0]-width+1; plo[1] = lo[1]; plo[2] = lo[2];
                phi[0] = hi[0]        ; phi[1] = hi[1]; phi[2] = hi[2];
                const Box pbx_xhi(plo, phi);

                Box xhi_plane_no_stag(pbx_xhi);
                Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
                Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});

                Box xhi_line(IntVect(hi[0]-width+1, lo[1], 0), IntVect(hi[0], hi[1], 0));

                //amrex::Print() << "HI XBX NO STAG " << pbx_xhi << std::endl;
                //amrex::Print() << "HI XBX  X STAG " << xhi_plane_x_stag << std::endl;
                //amrex::Print() << "HI XBX  Y STAG " << xhi_plane_y_stag << std::endl;

                if        (bdyVarType == WRFBdyVars::U) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_x_stag, 1, Arena_Used)); // U
                    }
                } else if (bdyVarType == WRFBdyVars::V) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_y_stag , 1, Arena_Used)); // V
                    }
                } else if (bdyVarType == WRFBdyVars::R) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // R
                    }
                } else if (bdyVarType == WRFBdyVars::T) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // T
                    }
                } else if (bdyVarType == WRFBdyVars::QV) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // QV
                    }
                } else if (bdyVarType == WRFBdyVars::MU ||
                           bdyVarType == WRFBdyVars::PC) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_xhi[nt].push_back(FArrayBox(xhi_line, 1, Arena_Used)); // MU
                    }
                }

            } else if (bdyType == WRFBdyTypes::y_lo) {

                // *******************************************************************************
                // ylo bdy
                // *******************************************************************************
                plo[1] = lo[1]        ; plo[0] = lo[0]; plo[2] = lo[2];
                phi[1] = lo[1]+width-1; phi[0] = hi[0]; phi[2] = hi[2];
                const Box pbx_ylo(plo, phi);

                Box ylo_plane_no_stag(pbx_ylo);
                Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
                Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);

                Box ylo_line(IntVect(lo[0], lo[1], 0), IntVect(hi[0], lo[1]+width-1, 0));

                if        (bdyVarType == WRFBdyVars::U) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_x_stag , 1, Arena_Used)); // U
                    }
                } else if (bdyVarType == WRFBdyVars::V) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_y_stag, 1, Arena_Used)); // V
                    }
                } else if (bdyVarType == WRFBdyVars::R) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // R
                    }
                } else if (bdyVarType == WRFBdyVars::T) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // T
                    }
                } else if (bdyVarType == WRFBdyVars::QV) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // QV
                    }
                } else if (bdyVarType == WRFBdyVars::MU ||
                           bdyVarType == WRFBdyVars::PC) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_ylo[nt].push_back(FArrayBox(ylo_line, 1, Arena_Used)); // PC
                    }
                }

            } else if (bdyType == WRFBdyTypes::y_hi) {

                // *******************************************************************************
                // yhi bdy
                // *******************************************************************************
                plo[1] = hi[1]-width+1; plo[0] = lo[0]; plo[2] = lo[2];
                phi[1] = hi[1]        ; phi[0] = hi[0]; phi[2] = hi[2];
                const Box pbx_yhi(plo, phi);

                Box yhi_plane_no_stag(pbx_yhi);
                Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
                Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);

                Box yhi_line(IntVect(lo[0], hi[1]-width+1, 0), IntVect(hi[0], hi[1], 0));

                //amrex::Print() << "HI YBX NO STAG " << pbx_yhi << std::endl;
                //amrex::Print() << "HI YBX  X STAG " << yhi_plane_x_stag << std::endl;
                //amrex::Print() << "HI YBX  Y STAG " << yhi_plane_y_stag << std::endl;

                if        (bdyVarType == WRFBdyVars::U) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_x_stag , 1, Arena_Used)); // U
                    }
                } else if (bdyVarType == WRFBdyVars::V) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_y_stag, 1, Arena_Used)); // V
                    }
                } else if (bdyVarType == WRFBdyVars::R) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // R
                    }
                } else if (bdyVarType == WRFBdyVars::T) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // T
                    }
                } else if (bdyVarType == WRFBdyVars::QV) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // QV
                    }
                } else if (bdyVarType == WRFBdyVars::MU ||
                           bdyVarType == WRFBdyVars::PC) {
                    for (int nt(0); nt < ntimes; ++nt) {
                      bdy_data_yhi[nt].push_back(FArrayBox(yhi_line, 1, Arena_Used)); // PC
                    }
                }
        }

        long num_pts;

        // Now fill the data
        if (ParallelDescriptor::IOProcessor())
        {
            // amrex::Print() << "SHAPE0 " << arrays[iv].get_vshape()[0] << std::endl;
            // amrex::Print() << "SHAPE1 " << arrays[iv].get_vshape()[1] << std::endl;
            // amrex::Print() << "SHAPE2 " << arrays[iv].get_vshape()[2] << std::endl;
            // amrex::Print() << "SHAPE3 " << arrays[iv].get_vshape()[3] << std::endl;

            Array4<Real> fab_arr;
            if (bdyVarType == WRFBdyVars::U || bdyVarType == WRFBdyVars::V ||
                bdyVarType == WRFBdyVars::R || bdyVarType == WRFBdyVars::T ||
                bdyVarType == WRFBdyVars::QV)
            {
                int ns2 = arrays[iv].get_vshape()[2];
                int ns3 = arrays[iv].get_vshape()[3];

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = bdy_data_xlo[0][bdyVarType].box().numPts();
                    int ioff = bdy_data_xlo[0][bdyVarType].smallEnd()[0];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_xlo[nt][bdyVarType].array();
                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int i = n / (ns2 * ns3);
                            int k = (n - i * (ns2 * ns3)) / ns3;
                            int j =  n - i * (ns2 * ns3) - k * ns3;
                            fab_arr(ioff+i, j, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = bdy_data_xhi[0][bdyVarType].box().numPts();
                    int ioff = bdy_data_xhi[0][bdyVarType].bigEnd()[0];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_xhi[nt][bdyVarType].array();
                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int i = n / (ns2 * ns3);
                            int k = (n - i * (ns2 * ns3)) / ns3;
                            int j =  n - i * (ns2 * ns3) - k * ns3;
                            fab_arr(ioff-i, j, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = bdy_data_ylo[0][bdyVarType].box().numPts();
                    int joff = bdy_data_ylo[0][bdyVarType].smallEnd()[1];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_ylo[nt][bdyVarType].array();
                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int j = n / (ns2 * ns3);
                            int k = (n - j * (ns2 * ns3)) / ns3;
                            int i =  n - j * (ns2 * ns3) - k * ns3;
                            fab_arr(i, joff+j, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = bdy_data_yhi[0][bdyVarType].box().numPts();
                    int joff = bdy_data_yhi[0][bdyVarType].bigEnd()[1];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_yhi[nt][bdyVarType].array();
                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int j = n / (ns2 * ns3);
                            int k = (n - j * (ns2 * ns3)) / ns3;
                            int i =  n - j * (ns2 * ns3) - k * ns3;
                            fab_arr(i, joff-j, k, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } // bdyType

            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = bdy_data_xlo[0][bdyVarType].box().numPts();
                    int ioff = bdy_data_xlo[0][bdyVarType].smallEnd()[0];
                    int ns2 = arrays[iv].get_vshape()[2];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_xlo[nt][bdyVarType].array();

                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int i = n / ns2;
                            int j = n - i * ns2;
                            fab_arr(ioff+i, j, 0, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = bdy_data_xhi[0][bdyVarType].box().numPts();
                    int ioff = bdy_data_xhi[0][bdyVarType].bigEnd()[0];
                    int ns2 = arrays[iv].get_vshape()[2];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_xhi[nt][bdyVarType].array();

                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int i = n / ns2;
                            int j = n - i * ns2;
                            fab_arr(ioff-i, j, 0, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = bdy_data_ylo[0][bdyVarType].box().numPts();
                    int joff = bdy_data_ylo[0][bdyVarType].smallEnd()[1];
                    int ns2 = arrays[iv].get_vshape()[2];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_ylo[nt][bdyVarType].array();

                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int j = n / ns2;
                            int i = n - j * ns2;
                            fab_arr(i, joff+j, 0, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = bdy_data_yhi[0][bdyVarType].box().numPts();
                    int joff = bdy_data_yhi[0][bdyVarType].bigEnd()[1];
                    int ns2 = arrays[iv].get_vshape()[2];
                    for (int nt(0); nt < ntimes; ++nt)
                    {
                        fab_arr  = bdy_data_yhi[nt][bdyVarType].array();

                        int n_off = nt * num_pts;
                        for (int n(0); n < num_pts; ++n) {
                            int j = n / ns2;
                            int i = n - j * ns2;
                            fab_arr(i, joff-j, 0, 0) = static_cast<Real>(*(arrays[iv].get_data() + n + n_off));
                        }
                    }
                }
            } // bdyVarType
        } // if ParalleDescriptor::IOProcessor()
    } // nc_var_names

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    amrex::ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int n_per_time = nc_var_prefix.size();
    for (int nt = 0; nt < ntimes; nt++)
    {
        for (int i = 0; i < n_per_time; i++)
        {
            ParallelDescriptor::Bcast(bdy_data_xlo[nt][i].dataPtr(),bdy_data_xlo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_xhi[nt][i].dataPtr(),bdy_data_xhi[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_ylo[nt][i].dataPtr(),bdy_data_ylo[nt][i].box().numPts(),ioproc);
            ParallelDescriptor::Bcast(bdy_data_yhi[nt][i].dataPtr(),bdy_data_yhi[nt][i].box().numPts(),ioproc);
        }
    }

    // Make sure all processors know how timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

void
convert_wrfbdy_data (int which, const Box& domain, Vector<Vector<FArrayBox>>& bdy_data,
                     const FArrayBox& NC_MUB_fab,
                     const FArrayBox& NC_PH_fab, const FArrayBox& NC_PHB_fab,
                     const FArrayBox& NC_C1H_fab, const FArrayBox& NC_C2H_fab,
                     const FArrayBox& NC_RDNW_fab,
                     const FArrayBox& NC_xvel_fab, const FArrayBox& NC_yvel_fab,
                     const FArrayBox& NC_rho_fab, const FArrayBox& NC_rhotheta_fab,
                     const FArrayBox& NC_QVAPOR_fab)
{
    // These were filled from wrfinput
    Array4<Real const> c1h_arr  = NC_C1H_fab.const_array();
    Array4<Real const> c2h_arr  = NC_C2H_fab.const_array();
    Array4<Real const> rdnw_arr = NC_RDNW_fab.const_array();
    Array4<Real const> mub_arr  = NC_MUB_fab.const_array();

    Array4<Real const>  ph_arr  = NC_PH_fab.const_array();
    Array4<Real const> phb_arr  = NC_PHB_fab.const_array();

    int ntimes = bdy_data.size();
    for (int nt = 0; nt < ntimes; nt++)
    {
        Array4<Real> bdy_u_arr  = bdy_data[nt][WRFBdyVars::U].array();  // This is face-centered
        Array4<Real> bdy_v_arr  = bdy_data[nt][WRFBdyVars::V].array();
        Array4<Real> bdy_r_arr  = bdy_data[nt][WRFBdyVars::R].array();
        Array4<Real> bdy_t_arr  = bdy_data[nt][WRFBdyVars::T].array();
        Array4<Real> bdy_qv_arr = bdy_data[nt][WRFBdyVars::QV].array();
        Array4<Real> mu_arr     = bdy_data[nt][WRFBdyVars::MU].array(); // This is cell-centered

        int ilo  = domain.smallEnd()[0];
        int ihi  = domain.bigEnd()[0];
        int jlo  = domain.smallEnd()[1];
        int jhi  = domain.bigEnd()[1];

        if (nt==0) {
            bdy_data[0][WRFBdyVars::U].template copy<RunOn::Device>(NC_xvel_fab);
            bdy_data[0][WRFBdyVars::V].template copy<RunOn::Device>(NC_yvel_fab);
            bdy_data[0][WRFBdyVars::R].template copy<RunOn::Device>(NC_rho_fab);
            bdy_data[0][WRFBdyVars::T].template copy<RunOn::Device>(NC_rhotheta_fab);
            bdy_data[0][WRFBdyVars::QV].template copy<RunOn::Device>(NC_QVAPOR_fab);
            bdy_data[0][WRFBdyVars::QV].template mult<RunOn::Device>(NC_rho_fab);
        } else {
            // Define u velocity
            const auto & bx_u  = bdy_data[0][WRFBdyVars::U].box();
            ParallelFor(bx_u, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xmu;
                if (i == ilo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (i > ihi) {
                    xmu  = mu_arr(i-1,j,0) + mub_arr(i-1,j,0);
                } else {
                    xmu = ( mu_arr(i,j,0) +  mu_arr(i-1,j,0)
                           +mub_arr(i,j,0) + mub_arr(i-1,j,0)) * 0.5;
                }
                Real xmu_mult = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy = bdy_u_arr(i,j,k) / xmu_mult;
                bdy_u_arr(i,j,k) = new_bdy;
            });

            // Define v velocity
            const auto & bx_v  = bdy_data[0][WRFBdyVars::V].box();
            ParallelFor(bx_v, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xmu;
                if (j == jlo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (j > jhi) {
                    xmu  = mu_arr(i,j-1,0) + mub_arr(i,j-1,0);
                } else {
                    xmu =  ( mu_arr(i,j,0) +  mu_arr(i,j-1,0)
                            +mub_arr(i,j,0) + mub_arr(i,j-1,0) ) * 0.5;
                }
                Real xmu_mult = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy = bdy_v_arr(i,j,k) / xmu_mult;
                bdy_v_arr(i,j,k) = new_bdy;
            });

            // Define density
            const auto & bx_t = bdy_data[0][WRFBdyVars::T].box(); // Note this is currently "THM" aka the perturbational moist pot. temp.
            ParallelFor(bx_t, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xmu  = c1h_arr(0,0,k) * (mu_arr(i,j,0) + mub_arr(i,j,0)) + c2h_arr(0,0,k);
                Real dpht = (ph_arr(i,j,k+1) + phb_arr(i,j,k+1)) - (ph_arr(i,j,k) + phb_arr(i,j,k));
                bdy_r_arr(i,j,k) = -xmu / ( dpht * rdnw_arr(0,0,k) );
                //if (nt == 0 and std::abs(r_arr(i,j,k) - bdy_r_arr(i,j,k)) > 0.) {
                //    amrex::Print() << "INIT VS BDY DEN " << IntVect(i,j,k) << " " << r_arr(i,j,k) << " " << bdy_r_arr(i,j,k) <<
                //                    " " << std::abs(r_arr(i,j,k) - bdy_r_arr(i,j,k)) << std::endl;
                //}
            });

            // Define theta
            amrex::Real theta_ref = 300.;
            ParallelFor(bx_t, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xmu  = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_Th = bdy_t_arr(i,j,k) / xmu_mult + theta_ref;
                Real qv_fac = (1. + bdy_qv_arr(i,j,k) / 0.622 / xmu_mult);
                new_bdy_Th /= qv_fac;
                bdy_t_arr(i,j,k) = new_bdy_Th * bdy_r_arr(i,j,k);
                //if (nt == 0 and std::abs(rth_arr(i,j,k) - bdy_t_arr(i,j,k)) > 0.) {
                //    amrex::Print() << "INIT VS BDY TH " << IntVect(i,j,k) << " " << rth_arr(i,j,k) << " " << bdy_t_arr(i,j,k) <<
                //                    " " << std::abs(th_arr(i,j,k) - bdy_t_arr(i,j,k)) << std::endl;
                //}
            });

            // Define Qv
            const auto & bx_qv = bdy_data[0][WRFBdyVars::QV].box();
            ParallelFor(bx_qv, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real xmu  = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_QV = bdy_qv_arr(i,j,k) / xmu_mult;
                bdy_qv_arr(i,j,k) = new_bdy_QV * bdy_r_arr(i,j,k);
            });

        } // nt ==0
    } // ntimes
}
#endif // ERF_USE_NETCDF
