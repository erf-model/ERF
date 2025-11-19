#include "AMReX_FArrayBox.H"
#include "ERF_NCWpsFile.H"
#include "ERF_NCInterface.H"
#include "ERF_IndexDefines.H"
#include "ERF_SurfaceLayer.H"

#include <sstream>
#include <string>
#include <ctime>
#include <atomic>

#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "AMReX_Print.H"

using namespace amrex;

#ifdef ERF_USE_NETCDF

Real
read_times_from_wrflow (const std::string& nc_low_file,
                        Vector<Vector<FArrayBox>>& low_data_zlo,
                        Real& start_low_time)
{
    Print() << "Loading lower boundary data from NetCDF file " << std::endl;

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
        ReadNetCDFFile(nc_low_file, {"Times"}, array_ts, success);

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
            Print() << "  wrflow datetime " << nt << " : " << date << " " << epochTime << std::endl;
            epochTimes.push_back(epochTime);

            if (nt == 1) {
                timeInterval = epochTimes[1] - epochTimes[0];
            } else if (nt >= 1) {
                AMREX_ALWAYS_ASSERT(epochTimes[nt] - epochTimes[nt-1] == timeInterval);
            }
        }
        start_low_time = epochTimes[0];
    }

    ParallelDescriptor::Bcast(&start_low_time,1,ioproc);
    ParallelDescriptor::Bcast(&ntimes,1,ioproc);
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Our outermost loop is time
    low_data_zlo.resize(ntimes);

    // Make sure all processors know timeInterval
    ParallelDescriptor::Bcast(&timeInterval,1,ioproc);

    // Return the number of seconds between the boundary plane data
    return timeInterval;
}

void
read_from_wrflow (const int itime, const std::string& nc_low_file, const Box& domain,
                  Vector<Vector<FArrayBox>>& low_data_zlo)
{
    int ioproc = ParallelDescriptor::IOProcessorNumber();  // I/O rank
                                                           //
    const auto& lo = domain.loVect();
    const auto& hi = domain.hiVect();

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
            int ns2 = tslice[iv].get_vshape()[2];

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
update_sst_tsk (const int itime,
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
    auto& domain = geom.Domain();

    // Temporary MFs for derived quantities
    auto& dm    = cons.DistributionMap();
    IntVect ng  = cons.nGrowVect();
    IntVect ngv = ng; ngv[2] = 0;

    // Bounds limiting
    int ilo = domain.smallEnd()[0];
    int ihi = domain.bigEnd()[0];
    int jlo = domain.smallEnd()[1];
    int jhi = domain.bigEnd()[1];

    if (itime > 0) {
        sst_lev[itime] = std::make_unique<MultiFab>(ba2d_lev,dm,1,ngv);
        tsk_lev[itime] = std::make_unique<MultiFab>(ba2d_lev,dm,1,ngv);
        if (SurfLayer) {
            SurfLayer->update_sst_ptr(0, itime, sst_lev[itime].get());
            SurfLayer->update_tsk_ptr(0, itime, tsk_lev[itime].get());
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
      //const Array4<const Real>& con_arr = cons.const_array(mfi);

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = min(max(i, ilo), ihi);
            int lj = min(max(j, jlo), jhi);

            // For simplicity, assume constant surface pressure = p0
            // ==> surface theta == SST
//            sst_arr(i,j,0) = src_arr(li,lj,0);

            // Use local density to convert to potential temperature -- but
            // this should be averaged over the lowinput time interval
//            Real rho = con_arr(li, lj, 0, Rho_comp);
//            Real qv = (use_moist) con_arr(li, lj, 0, RhoQ1_comp) / rho : 0.0;
//            sst_arr(i,j,0) = getThgivenRandT(rho, src_arr(li,lj,0), rdOcp, qv);

            // NOTE: we convert to potential temperature for the surface
            // layer scheme using the initial surface pressure since it's
            // not available in the wrflowinp file
            bool is_land = (lmask_arr) ? lmask_arr(li,lj,0) : true;
            sst_arr(i,j,0) = getThgivenTandP(src_arr(li,lj,0), psfc_arr(li,lj,0), rdOcp);
            if (!is_land && std::abs(sst_arr(i,j,0)) < 400.0) {
                tsk_arr(i,j,0) = sst_arr(i,j,0);
            }
        });
    }

    sst_lev[itime]->FillBoundary(geom.periodicity());
    tsk_lev[itime]->FillBoundary(geom.periodicity());
}
#endif // ERF_USE_NETCDF
