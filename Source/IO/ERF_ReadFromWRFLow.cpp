/**
 * \file ERF_ReadFromWRFLow.cpp
 */
#include "AMReX_FArrayBox.H"
#include "ERF_NCWpsFile.H"
#include "ERF_NCInterface.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"

#include <sstream>
#include <string>
#include <ctime>
#include <limits>
#include <atomic>

#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "AMReX_Print.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

double
read_times_from_wrflow (const std::string& nc_low_file,
                        Vector<Vector<FArrayBox>>& low_data_zlo,
                        double& start_low_time, double& final_low_time)
{
    Print() << "Loading lower boundary data from NetCDF file " << std::endl;

    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank

    // *******************************************************************************

    int ntimes = 0;
    // With a single record there is no interval; the largest double makes
    // elapsed_time / interval select that record (index 0) at every time
    double timeInterval = std::numeric_limits<double>::max();
    const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor())
    {
        // Read the time stamps
        using CharArray = NDArray<char>;
        Vector<CharArray> array_ts(1);
        Vector<int> success(1);
        ReadNetCDFFile(nc_low_file, {"Times"}, array_ts, success);

        ntimes = array_ts[0].get_vshape()[0];

        Vector<std::string> timeStamps;
        timeStamps.reserve(ntimes);

        const char* data = array_ts[0].get_data();
        auto dateStrLen  = array_ts[0].get_vshape()[1];

        for (int nt = 0; nt < ntimes; ++nt) {
            const char* begin = data + nt * dateStrLen;
            timeStamps.emplace_back(begin, begin + dateStrLen);
        }


        Vector<std::time_t> epochTimes;
        for (int nt(0); nt < ntimes; nt++) {
            std::string date = timeStamps[nt];
            auto epochTime = getEpochTime(date, dateTimeFormat);
            Print() << "  wrflow datetime " << nt << " : " << date << " " << epochTime << std::endl;
            epochTimes.push_back(epochTime);

            if (nt == 1) {
                timeInterval = static_cast<double>(epochTimes[1] - epochTimes[0]);
            } else if (nt >= 1) {
                AMREX_ALWAYS_ASSERT(static_cast<double>(epochTimes[nt] - epochTimes[nt-1]) == timeInterval);
            }
        }
        start_low_time = static_cast<double>(epochTimes[0]);
        final_low_time = static_cast<double>(epochTimes[ntimes-1]);
    }

    ParallelDescriptor::Bcast(&start_low_time,1,ioproc);
    ParallelDescriptor::Bcast(&final_low_time,1,ioproc);
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Our outermost loop is time
    low_data_zlo.resize(ntimes);

    // Make sure all processors know timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

/**
 * Read time slice itime of the SST in a wrflow file into low_data_zlo[itime], on the 2D box
 * of cells the file covers: the domain for level 0, or the nest box for a finer level.
 *
 * @param[in]    itime        time slice to read
 * @param[in]    nc_low_file  name of the wrflow file
 * @param[in]    box          cells the file covers (only its lower z index is used)
 * @param[inout] low_data_zlo data read, by time slice and variable
 */
void
read_from_wrflow (const int itime, const std::string& nc_low_file, const Box& box,
                  Vector<Vector<FArrayBox>>& low_data_zlo)
{
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
                                                           //
    const auto& lo = box.loVect();
    const auto& hi = box.hiVect();

    IntVect plo(lo);
    IntVect phi(hi);

    plo[0] = lo[0]; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0]; phi[1] = hi[1]; phi[2] = lo[2];
    const Box pbx_zlo(plo, phi);

    // ******************************************************************
    // Read the NetCDF file and fill this FAB
    // ******************************************************************
    Vector<std::string> nc_var_names;
    nc_var_names.push_back("SST");

    using RARRAY = NDArray<float>;
    Vector<RARRAY> tslice(nc_var_names.size());

    if (ParallelDescriptor::IOProcessor())
    {
        Vector<int> success(nc_var_names.size());

        ReadTimeSliceFromNetCDFFile(nc_low_file, itime, nc_var_names, tslice, success);

        for (auto &istat:success) {
            AMREX_ALWAYS_ASSERT(istat==1);
        }
    }

    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    Arena_Used = The_Pinned_Arena();
#endif

    // *******************************************************************************
    // Now fill the data
    // *******************************************************************************
    int nvar = nc_var_names.size();
    for (int iv = 0; iv < nvar; iv++)
    {
        low_data_zlo[itime].push_back(FArrayBox(pbx_zlo, 1, Arena_Used));

        if (ParallelDescriptor::IOProcessor())
        {
            Array4<Real> fab_arr = low_data_zlo[itime][iv].array();

            // dims: (Time, south_north, west_east)
            int ns1 = tslice[iv].get_vshape()[1];
            int ns2 = tslice[iv].get_vshape()[2];
            if (ns2 != pbx_zlo.length(0) || ns1 != pbx_zlo.length(1)) {
                Abort("read_from_wrflow: " + nc_var_names[iv] + " in " + nc_low_file + " is " +
                      std::to_string(ns2) + " x " + std::to_string(ns1) + " but the cells it must fill (" +
                      std::to_string(pbx_zlo.length(0)) + " x " + std::to_string(pbx_zlo.length(1)) +
                      ") do not match");
            }

            long num_pts  = low_data_zlo[itime][iv].box().numPts();
            int ioff      = low_data_zlo[itime][iv].smallEnd()[0];
            int joff      = low_data_zlo[itime][iv].smallEnd()[1];

            for (int n(0); n < num_pts; ++n) {
                int j = n / ns2;
                int i = n - j*ns2;
                fab_arr(ioff+i,joff+j,0) = static_cast<Real>(*(tslice[iv].get_data() + n));
            }
        } // if ParalleDescriptor::IOProcessor()
    } // nc_var_names

    // We put a barrier here so the rest of the processors wait to do anything until they have the data
    ParallelDescriptor::Barrier();

    // When an FArrayBox is built, space is allocated on every rank.  However, we only
    //    filled the data in these FABs on the IOProcessor.  So here we broadcast
    //    the data to every rank.
    for (int iv = 0; iv < nvar; iv++)
    {
        ParallelDescriptor::Bcast(low_data_zlo[itime][iv].dataPtr(),low_data_zlo[itime][iv].box().numPts(),ioproc);
    }
}

void
update_sst_tsk (const int lev,
                const int itime,
                const Geometry& geom,
                const BoxArray& ba2d_lev,
                Vector<std::unique_ptr<MultiFab>>& sst_lev,
                Vector<std::unique_ptr<MultiFab>>& tsk_lev,
                std::unique_ptr<SurfaceLayer>& SurfLayer,
                const Vector<Vector<FArrayBox>>& low_data_zlo,
                const MultiFab& cons,
                const MultiFab& mf_PSFC_lev,
                const Real rdOcp,
                std::unique_ptr<iMultiFab>& lmask,
                const bool /*use_moist*/)
{
    // Temporary MFs for derived quantities
    auto& dm    = cons.DistributionMap();
    IntVect ng  = cons.nGrowVect();
    IntVect ngv = ng; ngv[2] = 0;

    // Bounds limiting: to the cells the file covers, which are the domain at level 0 and the
    // nest box at a finer level
    const Box& src_box = low_data_zlo[itime][0].box();
    int ilo = src_box.smallEnd(0);
    int ihi = src_box.bigEnd(0);
    int jlo = src_box.smallEnd(1);
    int jhi = src_box.bigEnd(1);

    if (itime > 0) {
        sst_lev[itime] = std::make_unique<MultiFab>(ba2d_lev,dm,1,ngv);
        tsk_lev[itime] = std::make_unique<MultiFab>(ba2d_lev,dm,1,ngv);

        // NOTE: we assume that TSK was in the wrfinput file so it is defined
        //       at the first time, but not in the wrflow file so not defined
        //       at later times.  Here we fill the later time arrays by copying
        //       from the initial time
        MultiFab::Copy(*tsk_lev[itime],*tsk_lev[0],0,0,1,tsk_lev[0]->nGrowVect());

        if (SurfLayer) {
            SurfLayer->update_sst_ptr(lev, itime, sst_lev[itime].get());
            SurfLayer->update_tsk_ptr(lev, itime, tsk_lev[itime].get());
        }
    }

    for ( MFIter mfi(*(sst_lev[itime]), false); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();

        const FArrayBox& src = low_data_zlo[itime][0];
        FArrayBox& sst_fab = (*(sst_lev[itime]))[mfi];
        FArrayBox& tsk_fab = (*(tsk_lev[itime]))[mfi];

        const Array4<      Real>& sst_arr   = sst_fab.array();
        const Array4<      Real>& tsk_arr   = tsk_fab.array();
        const Array4<const Real>& src_arr   = src.const_array();
        const Array4<const Real>& psfc_arr  = mf_PSFC_lev.const_array(mfi);
        const Array4<const  int>& lmask_arr = (lmask) ? lmask->const_array(mfi) : Array4<const int> {};

        // For simplicity, assume constant surface pressure = p0
        // ==> surface theta == SST.
        // We convert to theta for the surface layer scheme using the
        // the initial surface pressure since it's not available in the wrflowinp file

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = min(max(i, ilo), ihi);
            int lj = min(max(j, jlo), jhi);

            bool is_land = (lmask_arr) ? lmask_arr(li,lj,0) : true;
            sst_arr(i,j,0) = getThgivenTandP(src_arr(li,lj,0), psfc_arr(li,lj,0), rdOcp);
            if (!is_land && std::abs(sst_arr(i,j,0)) < Real(400.0)) {
                tsk_arr(i,j,0) = sst_arr(i,j,0);
            }
        });
    }

    sst_lev[itime]->FillBoundary(geom.periodicity());
    tsk_lev[itime]->FillBoundary(geom.periodicity());
}
#endif // ERF_USE_NETCDF
