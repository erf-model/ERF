#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "AMReX_FArrayBox.H"
#include "AMReX_Print.H"

#include "ERF_NCWpsFile.H"
#include "ERF_IndexDefines.H"
#include "ERF_EOS.H"
#include "ERF_DataStruct.H"
#include "ERF_NCInterface.H"

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
read_times_from_wrfbdy (const std::string& nc_bdy_file,
                        Vector<Vector<FArrayBox>>& bdy_data_xlo,
                        Vector<Vector<FArrayBox>>& bdy_data_xhi,
                        Vector<Vector<FArrayBox>>& bdy_data_ylo,
                        Vector<Vector<FArrayBox>>& bdy_data_yhi,
                        Real& start_bdy_time)
{
    Print() << "Loading boundary data from NetCDF file " << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // *******************************************************************************

    int ntimes;
    Real timeInterval;
    const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor())
    {
        // Read the time stamps
        using CharArray = NDArray<char>;
        Vector<CharArray> array_ts(1);
        Vector<int> success(1);
        ReadNetCDFFile(nc_bdy_file, {"Times"}, array_ts, success);

        ntimes = array_ts[0].get_vshape()[0];

        auto dateStrLen = array_ts[0].get_vshape()[1];
        char timeStamps[ntimes][dateStrLen];

        // Fill up the characters read
        int str_len = static_cast<int>(dateStrLen);
        for (int nt(0); nt < ntimes; nt++) {
            for (int dateStrCt(0); dateStrCt < str_len; dateStrCt++) {
                auto n = nt*dateStrLen + dateStrCt;
                timeStamps[nt][dateStrCt] = *(array_ts[0].get_data() + n);
            }
        }

        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string date(&timeStamps[nt][0], &timeStamps[nt][dateStrLen-1]+1);
            auto epochTime = getEpochTime(date, dateTimeFormat);
            Print() << "  wrfbdy datetime " << nt << " : " << date << " " << epochTime << std::endl;
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

    // Our outermost loop is time
    bdy_data_xlo.resize(ntimes);
    bdy_data_xhi.resize(ntimes);
    bdy_data_ylo.resize(ntimes);
    bdy_data_yhi.resize(ntimes);

    // Make sure all processors know timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

void
read_from_wrfbdy (const int itime, const std::string& nc_bdy_file, const Box& domain,
                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                  Vector<Vector<FArrayBox>>& bdy_data_yhi,
                  int real_width)
{
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // Even though we may not read in all the variables, we need to make the arrays big enough for them (for now)
    int nvars = WRFBdyVars::NumTypes*4;

    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();
    const int khi = hi[2];

    IntVect plo(lo);
    IntVect phi(hi);

    // ******************************************************************
    // Read the netcdf file and fill these FABs
    // NOTE: the order and number of these must match the WRFBdyVars enum!
    // WRFBdyVars:  U, V, R, T, QV, MU, PC
    //
    // These fields are at half levels (unstaggered)
    // ******************************************************************
    Vector<std::string> nc_var_names;
    Vector<std::string> nc_var_prefix = {"U","V","T","QVAPOR","MU","PC"};

    for (int ip = 0; ip < nc_var_prefix.size(); ++ip)
    {
       nc_var_names.push_back(nc_var_prefix[ip] + "_BXS");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BXE");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BYS");
       nc_var_names.push_back(nc_var_prefix[ip] + "_BYE");
    }

    using RARRAY = NDArray<float>;
    Vector<RARRAY> tslice(nc_var_names.size());

    int width; // size of bdy_width from wrfbdy

    if (ParallelDescriptor::IOProcessor())
    {
        Vector<int> success(nc_var_names.size());

        ReadTimeSliceFromNetCDFFile(nc_bdy_file, itime, nc_var_names, tslice, success);

        for (auto &istat:success) {
            AMREX_ALWAYS_ASSERT(istat==1);
        }

        // Width of the boundary region
        width = tslice[0].get_vshape()[1];

        if (width != real_width) {
            AMREX_ALWAYS_ASSERT(real_width < width);
            Print() << "Note: Requested boundary width is " << real_width
                << " < " << width << " (bdy_width size in file)" << std::endl;
        }
    }
    ParallelDescriptor::Bcast(&width,1,ioproc);

    // This loops over every variable on every face, so nvars should be 4 * number of "ivartype" below
    for (int iv = 0; iv < nvars; iv++)
    {
        // Print() << "Building FAB for the NetCDF variable : " << nc_var_names[iv] << std::endl;

        int bdyVarType;

        std::string first1 = nc_var_names[iv].substr(0,1);
        std::string first2 = nc_var_names[iv].substr(0,2);

        if        (first1 == "U") {
            bdyVarType = WRFBdyVars::U;
        } else if (first1 == "V") {
            bdyVarType = WRFBdyVars::V;
        } else if (first1 == "T") {
            bdyVarType = WRFBdyVars::T;
        } else if (first2 == "QV") {
            bdyVarType = WRFBdyVars::QV;
        } else if (first2 == "MU") {
            bdyVarType = WRFBdyVars::MU;
        } else if (first2 == "PC") {
            bdyVarType = WRFBdyVars::PC;
        } else {
            Print() << "Trying to read " << first1 << " or " << first2 << std::endl;
            Abort("dont know this variable");
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

        Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
        Arena_Used = The_Pinned_Arena();
#endif

        if (bdyType == WRFBdyTypes::x_lo) {

            // *******************************************************************************
            // xlo bdy
            // *******************************************************************************
            plo[0] = lo[0]             ; plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = lo[0]+real_width-1; phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xlo(plo, phi);

            Box xlo_plane_no_stag(pbx_xlo);
            Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
            Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});

            Box xlo_line(IntVect(lo[0], lo[1], 0), IntVect(lo[0]+real_width-1, hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_x_stag, 1, Arena_Used)); // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_y_stag , 1, Arena_Used)); // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::MU ||
                       bdyVarType == WRFBdyVars::PC) {
                bdy_data_xlo[itime].push_back(FArrayBox(xlo_line, 1, Arena_Used));
            }

        } else if (bdyType == WRFBdyTypes::x_hi) {

            // *******************************************************************************
            // xhi bdy
            // *******************************************************************************
            plo[0] = hi[0]-real_width+1; plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = hi[0]             ; phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xhi(plo, phi);

            Box xhi_plane_no_stag(pbx_xhi);
            Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
            Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});

            Box xhi_line(IntVect(hi[0]-real_width+1, lo[1], 0), IntVect(hi[0], hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_x_stag, 1, Arena_Used)); // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_y_stag , 1, Arena_Used)); // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::MU ||
                       bdyVarType == WRFBdyVars::PC) {
                bdy_data_xhi[itime].push_back(FArrayBox(xhi_line, 1, Arena_Used)); // MU
            }

        } else if (bdyType == WRFBdyTypes::y_lo) {

            // *******************************************************************************
            // ylo bdy
            // *******************************************************************************
            plo[1] = lo[1]             ; plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = lo[1]+real_width-1; phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_ylo(plo, phi);

            Box ylo_plane_no_stag(pbx_ylo);
            Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
            Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);

            Box ylo_line(IntVect(lo[0], lo[1], 0), IntVect(hi[0], lo[1]+real_width-1, 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_x_stag , 1, Arena_Used)); // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_y_stag, 1, Arena_Used)); // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::MU ||
                       bdyVarType == WRFBdyVars::PC) {
                bdy_data_ylo[itime].push_back(FArrayBox(ylo_line, 1, Arena_Used)); // PC
            }

        } else if (bdyType == WRFBdyTypes::y_hi) {

            // *******************************************************************************
            // yhi bdy
            // *******************************************************************************
            plo[1] = hi[1]-real_width+1; plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = hi[1]             ; phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_yhi(plo, phi);

            Box yhi_plane_no_stag(pbx_yhi);
            Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
            Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);

            Box yhi_line(IntVect(lo[0], hi[1]-real_width+1, 0), IntVect(hi[0], hi[1], 0));

            if        (bdyVarType == WRFBdyVars::U) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_x_stag , 1, Arena_Used)); // U
            } else if (bdyVarType == WRFBdyVars::V) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_y_stag, 1, Arena_Used)); // V
            } else if (bdyVarType == WRFBdyVars::T) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // T
            } else if (bdyVarType == WRFBdyVars::QV) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_plane_no_stag, 1, Arena_Used)); // QV
            } else if (bdyVarType == WRFBdyVars::MU ||
                       bdyVarType == WRFBdyVars::PC) {
                bdy_data_yhi[itime].push_back(FArrayBox(yhi_line, 1, Arena_Used)); // PC
            }
        }

        long num_pts;

        // Now fill the data
        if (ParallelDescriptor::IOProcessor())
        {
            // Print() << "SHAPE0 " << tslice[iv].get_vshape()[0] << std::endl;
            // Print() << "SHAPE1 " << tslice[iv].get_vshape()[1] << std::endl;
            // Print() << "SHAPE2 " << tslice[iv].get_vshape()[2] << std::endl;
            // Print() << "SHAPE3 " << tslice[iv].get_vshape()[3] << std::endl;

            Array4<Real> fab_arr;

            if (bdyVarType == WRFBdyVars::U || bdyVarType == WRFBdyVars::V ||
                bdyVarType == WRFBdyVars::T || bdyVarType == WRFBdyVars::QV)
            {
                // xlo,xhi dims: (Time, bdy_width, bottom_top, south_north)
                // ylo,yhi dims: (Time, bdy_width, bottom_top, west_east)

                int ns2 = tslice[iv].get_vshape()[2]; // vertical size
                int ns3 = tslice[iv].get_vshape()[3]; // lateral size, may be staggered

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xlo[itime][bdyVarType].smallEnd()[0];
                    fab_arr  = bdy_data_xlo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / (ns2 * ns3);
                        if (i >= real_width) continue;
                        int k = (n - i * (ns2 * ns3)) / ns3;
                        if (k > khi) continue;
                        int j =  n - i * (ns2 * ns3) - k * ns3;
                        fab_arr(ioff+i, j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xhi[itime][bdyVarType].bigEnd()[0];
                    fab_arr  = bdy_data_xhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / (ns2 * ns3);
                        if (i >= real_width) continue;
                        int k = (n - i * (ns2 * ns3)) / ns3;
                        if (k > khi) continue;
                        int j =  n - i * (ns2 * ns3) - k * ns3;
                        fab_arr(ioff-i, j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_ylo[itime][bdyVarType].smallEnd()[1];
                    fab_arr  = bdy_data_ylo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / (ns2 * ns3);
                        if (j >= real_width) continue;
                        int k = (n - j * (ns2 * ns3)) / ns3;
                        if (k > khi) continue;
                        int i =  n - j * (ns2 * ns3) - k * ns3;
                        fab_arr(i, joff+j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_yhi[itime][bdyVarType].bigEnd()[1];
                    fab_arr  = bdy_data_yhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / (ns2 * ns3);
                        if (j >= real_width) continue;
                        int k = (n - j * (ns2 * ns3)) / ns3;
                        if (k > khi) continue;
                        int i =  n - j * (ns2 * ns3) - k * ns3;
                        fab_arr(i, joff-j, k, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } // bdyType

            } else if (bdyVarType == WRFBdyVars::MU || bdyVarType == WRFBdyVars::PC) {
                // xlo,xhi dims: (Time, bdy_width, south_north)
                // ylo,yhi dims: (Time, bdy_width, west_east)

                if (bdyType == WRFBdyTypes::x_lo) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xlo[itime][bdyVarType].smallEnd()[0];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_xlo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / ns2;
                        if (i >= real_width) continue;
                        int j = n - i * ns2;
                        fab_arr(ioff+i, j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::x_hi) {
                    num_pts  = tslice[iv].ndim();
                    int ioff = bdy_data_xhi[itime][bdyVarType].bigEnd()[0];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_xhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int i = n / ns2;
                        if (i >= real_width) continue;
                        int j = n - i * ns2;
                        fab_arr(ioff-i, j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_lo) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_ylo[itime][bdyVarType].smallEnd()[1];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_ylo[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / ns2;
                        if (j >= real_width) continue;
                        int i = n - j * ns2;
                        fab_arr(i, joff+j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                } else if (bdyType == WRFBdyTypes::y_hi) {
                    num_pts  = tslice[iv].ndim();
                    int joff = bdy_data_yhi[itime][bdyVarType].bigEnd()[1];
                    int ns2 = tslice[iv].get_vshape()[2];
                    fab_arr  = bdy_data_yhi[itime][bdyVarType].array();
                    for (int n(0); n < num_pts; ++n) {
                        int j = n / ns2;
                        if (j >= real_width) continue;
                        int i = n - j * ns2;
                        fab_arr(i, joff-j, 0, 0) = static_cast<Real>(*(tslice[iv].get_data() + n));
                    }
                }
            } // bdyVarType
        } // if ParalleDescriptor::IOProcessor()
    } // nc_var_names

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    int n_per_time = nc_var_prefix.size();
    for (int i = 0; i < n_per_time; i++)
    {
        ParallelDescriptor::Bcast(bdy_data_xlo[itime][i].dataPtr(),bdy_data_xlo[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_xhi[itime][i].dataPtr(),bdy_data_xhi[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_ylo[itime][i].dataPtr(),bdy_data_ylo[itime][i].box().numPts(),ioproc);
        ParallelDescriptor::Bcast(bdy_data_yhi[itime][i].dataPtr(),bdy_data_yhi[itime][i].box().numPts(),ioproc);
    }
}

void
convert_wrfbdy_data (const int itime,
                     const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     const MultiFab& mf_MUB,
                     const MultiFab& mf_C1H,
                     const MultiFab& mf_C2H,
                     const MultiFab& xvel,
                     const MultiFab& yvel,
                     const MultiFab& cons,
                     const Geometry& geom,
                     const bool& use_moist)
{
    // Owner masks for parallel reduce sum
    std::unique_ptr<iMultiFab> mask_c = OwnerMask(cons, geom.periodicity());
    std::unique_ptr<iMultiFab> mask_u = OwnerMask(xvel, geom.periodicity());
    std::unique_ptr<iMultiFab> mask_v = OwnerMask(yvel, geom.periodicity());

    // Temporary bdy data structures for global reductions
    int vsize = bdy_data[itime].size() - 2; // Don't do MU & PC
    amrex::Vector<amrex::FArrayBox> bdy_data_tmp; bdy_data_tmp.resize(vsize);
    for (int ivar(0); ivar < vsize; ++ivar) {
        bdy_data_tmp[ivar].resize(bdy_data[itime][ivar].box(),1,The_Managed_Arena());
        bdy_data_tmp[ivar].template setVal<RunOn::Device>(0.);
    }

    // BDY data
    Array4<Real> bdy_u_arr  = bdy_data[itime][WRFBdyVars::U].array();  // This is x-face-centered
    Array4<Real> bdy_v_arr  = bdy_data[itime][WRFBdyVars::V].array();  // This is y-face-centered
    Array4<Real> bdy_t_arr  = bdy_data[itime][WRFBdyVars::T].array();  // This is cell-centered
    Array4<Real> bdy_qv_arr = bdy_data[itime][WRFBdyVars::QV].array(); // This is cell-centered
    Array4<Real> mu_arr     = bdy_data[itime][WRFBdyVars::MU].array(); // This is cell-centered

    // Bounds limiting
    int ilo  = domain.smallEnd()[0];
    int ihi  = domain.bigEnd()[0];
    int jlo  = domain.smallEnd()[1];
    int jhi  = domain.bigEnd()[1];

    for ( MFIter mfi(cons); mfi.isValid(); ++mfi )
    {
        Box tbx = mfi.tilebox();
        Box xbx = mfi.nodaltilebox(0);
        Box ybx = mfi.nodaltilebox(1);

        const Box& bx_u  = (xbx & bdy_data[itime][WRFBdyVars::U].box());
        const Box& bx_v  = (ybx & bdy_data[itime][WRFBdyVars::V].box());
        const Box& bx_t  = (tbx & bdy_data[itime][WRFBdyVars::T].box());
        const Box& bx_qv = (tbx & bdy_data[itime][WRFBdyVars::QV].box());

        // TMP BDY data
        Array4<Real> bdy_u_tmp  = bdy_data_tmp[WRFBdyVars::U].array();  // This is x-face-centered
        Array4<Real> bdy_v_tmp  = bdy_data_tmp[WRFBdyVars::V].array();  // This is y-face-centered
        Array4<Real> bdy_t_tmp  = bdy_data_tmp[WRFBdyVars::T].array();  // This is cell-centered
        Array4<Real> bdy_qv_tmp = bdy_data_tmp[WRFBdyVars::QV].array(); // This is cell-centered

        // Mask data
        const Array4<const int>& mask_c_arr = mask_c->const_array(mfi);
        const Array4<const int>& mask_u_arr = mask_u->const_array(mfi);
        const Array4<const int>& mask_v_arr = mask_v->const_array(mfi);

        // Populated from read wrfinput
        Array4<Real const> c1h_arr  = mf_C1H.const_array(mfi);
        Array4<Real const> c2h_arr  = mf_C2H.const_array(mfi);
        Array4<Real const> mub_arr  = mf_MUB.const_array(mfi);

        // Define u velocity
        ParallelFor(bx_u, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_u_arr(i,j,k)) {
                Real xmu;
                if (i == ilo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (i > ihi) {
                    xmu  = mu_arr(i-1,j,0) + mub_arr(i-1,j,0);
                } else {
                    xmu = (  mu_arr(i,j,0) +  mu_arr(i-1,j,0)
                          + mub_arr(i,j,0) + mub_arr(i-1,j,0)) * 0.5;
                }
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy     = bdy_u_arr(i,j,k) / xmu_mult;
                bdy_u_tmp(i,j,k) = new_bdy;
            }
        });

        // Define v velocity
        ParallelFor(bx_v, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_v_arr(i,j,k)) {
                Real xmu;
                if (j == jlo) {
                    xmu  = mu_arr(i,j,0) + mub_arr(i,j,0);
                } else if (j > jhi) {
                    xmu  = mu_arr(i,j-1,0) + mub_arr(i,j-1,0);
                } else {
                    xmu =  (  mu_arr(i,j,0) +  mu_arr(i,j-1,0)
                           + mub_arr(i,j,0) + mub_arr(i,j-1,0) ) * 0.5;
                }
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy     = bdy_v_arr(i,j,k) / xmu_mult;
                bdy_v_tmp(i,j,k) = new_bdy;
            }
        });

        // Convert perturbational moist pot. temp. (Th_m) to dry pot. temp. (Th_d)
        const Real wrf_theta_ref = 300.;
        ParallelFor(bx_t, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_c_arr(i,j,k)) {
                Real xmu         = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult    = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_Th  = bdy_t_arr(i,j,k) / xmu_mult + wrf_theta_ref;
                Real qv_fac      = (1. + (R_v/R_d) * bdy_qv_arr(i,j,k) / xmu_mult);
                new_bdy_Th      /= qv_fac;
                bdy_t_tmp(i,j,k) = new_bdy_Th;
            }
        });

        // Define Qv
        ParallelFor(bx_qv, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (mask_c_arr(i,j,k)) {
                Real xmu          = (mu_arr(i,j,0) + mub_arr(i,j,0));
                Real xmu_mult     = c1h_arr(0,0,k) * xmu + c2h_arr(0,0,k);
                Real new_bdy_QV   = bdy_qv_arr(i,j,k) / xmu_mult;
                bdy_qv_tmp(i,j,k) = (use_moist) ? new_bdy_QV : 0.;
            }
        });
    } // mfi

    for (int ivar(0); ivar < vsize; ++ivar) {
        amrex::ParallelAllReduce::Sum(bdy_data_tmp[ivar].dataPtr(),
                                      bdy_data_tmp[ivar].size(),
                                      ParallelContext::CommunicatorAll());
        bdy_data[itime][ivar].template  copy<RunOn::Device>(bdy_data_tmp[ivar],0,0,1);
    }
}

void
convert_all_wrfbdy_data (const int itime,
                         const Box& domain,
                         Vector<Vector<FArrayBox>>& bdy_data_xlo,
                         Vector<Vector<FArrayBox>>& bdy_data_xhi,
                         Vector<Vector<FArrayBox>>& bdy_data_ylo,
                         Vector<Vector<FArrayBox>>& bdy_data_yhi,
                         const MultiFab& mf_MUB,
                         const MultiFab& mf_C1H,
                         const MultiFab& mf_C2H,
                         const MultiFab& xvel,
                         const MultiFab& yvel,
                         const MultiFab& cons,
                         const Geometry& geom,
                         const bool& use_moist)
{
    convert_wrfbdy_data(itime, domain, bdy_data_xlo,
                        mf_MUB, mf_C1H, mf_C2H,
                        xvel, yvel, cons, geom, use_moist);
    convert_wrfbdy_data(itime, domain, bdy_data_xhi,
                        mf_MUB, mf_C1H, mf_C2H,
                        xvel, yvel, cons, geom, use_moist);
    convert_wrfbdy_data(itime, domain, bdy_data_ylo,
                        mf_MUB, mf_C1H, mf_C2H,
                        xvel, yvel, cons, geom, use_moist);
    convert_wrfbdy_data(itime, domain, bdy_data_yhi,
                        mf_MUB, mf_C1H, mf_C2H,
                        xvel, yvel, cons, geom, use_moist);
}
#endif // ERF_USE_NETCDF
