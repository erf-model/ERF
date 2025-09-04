#include "AMReX_FArrayBox.H"
#include "ERF_NCWpsFile.H"
#include "ERF_IndexDefines.H"

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "ERF_DataStruct.H"
#include "ERF_NCInterface.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_Print.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

Real
read_from_wrflow (const std::string& nc_low_file, const Box& domain,
                  Vector<FArrayBox>& low_data_zlo,
                  Real& start_bdy_time)
{
    Print() << "Loading low boundary data from NetCDF file " << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

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
        ReadNetCDFFile(nc_low_file, {"Times"}, array_ts, success);

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

    // *******************************************************************************
    // Build low_data_zlo
    // *******************************************************************************
    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = The_Pinned_Arena();
#endif

    IntVect plo(lo);
    IntVect phi(hi);

    plo[0] = lo[0]; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0]; phi[1] = hi[1]; phi[2] = lo[2];
    const Box pbx_zlo(plo, phi);

    for (int nt(0); nt < ntimes; ++nt) {
        amrex::Print() << "MAKING FAB FOR LOW DATA AT TIME " << nt << " " << pbx_zlo << std::endl;
        low_data_zlo.push_back(FArrayBox(pbx_zlo, 1, Arena_Used));
    }

    // ******************************************************************
    // Read the NetCDF file and fill this FABs
    // ******************************************************************
    Vector<std::string> nc_var_names;
    nc_var_names.push_back("SST");

    using RARRAY = NDArray<float>;
    Vector<RARRAY> arrays(nc_var_names.size());

    if (ParallelDescriptor::IOProcessor())
    {
        Vector<int> success(nc_var_names.size());

        for (int itime=0; itime < ntimes; ++itime)
        {
            Vector<RARRAY> tslice(nc_var_names.size());

            ReadTimeSliceFromNetCDFFile(nc_low_file, itime, nc_var_names, tslice, success);

#if 1
            std::vector<size_t> vshape = tslice[0].get_vshape();
            size_t offset{1};
            for (auto &dim:vshape) offset *= dim;

            Print() << "Time " << itime << " " << nc_var_names[0] << "(";
            for (auto &dim:vshape) {
                amrex::Print() << dim << ",";
            }
            Print() << ") offset=" << offset << std::endl;

            if (itime==0) {
                // allocate full array
                vshape[0] = ntimes;
                arrays[0] = NDArray<float>(tslice[0].get_vname(), vshape);
            }

            float* fullPtr = arrays[0].get_data();
            float* slicePtr = tslice[0].get_data();
            std::copy(slicePtr, slicePtr+offset, fullPtr+itime*offset);
        }
#endif
    }


    // *******************************************************************************
    // Now fill the data
    // *******************************************************************************
    if (ParallelDescriptor::IOProcessor())
    {
        Array4<Real> fab_arr;
        int ns3 = arrays[0].get_vshape()[2];

        long num_pts  = low_data_zlo[0].box().numPts();
        int ioff      = low_data_zlo[0].smallEnd()[0];
        int joff      = low_data_zlo[0].smallEnd()[1];

        for (int nt(0); nt < ntimes; ++nt)
        {
            fab_arr  = low_data_zlo[nt].array();
            int n_off = nt * num_pts;
            for (int n(0); n < num_pts; ++n) {
                int j  = n / ns3;
                int i  = n - j*ns3;
                fab_arr(ioff+i,joff+j,0) = static_cast<Real>(*(arrays[0].get_data()+n+n_off));
            }
        }
    } // if ParalleDescriptor::IOProcessor()

    // *******************************************************************************
    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    // *******************************************************************************
    ParallelDescriptor::Barrier();

    // *******************************************************************************
    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    // *******************************************************************************
    for (int nt = 0; nt < ntimes; nt++)
    {
        ParallelDescriptor::Bcast(low_data_zlo[nt].dataPtr(),low_data_zlo[nt].box().numPts(),ioproc);
    }

    // Make sure all processors know how timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Make sure all processors know how many times are stored
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

#endif // ERF_USE_NETCDF
