/**
 * \file ERF_InitFromMetgrid.cpp
 */
#include <ERF_Constants.H>
#include <ERF_MetgridUtils.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

#include "ERF_NCWpsFile.H"

/**
 * Reads start_time from the first metgrid file
 *
*/
Real
read_start_time_from_metgrid(int lev, const std::string& fname)
{
    Real NC_epochTime;
    const std::string dateTimeFormat = "%Y-%m-%d_%H:%M:%S";

    if (ParallelDescriptor::IOProcessor()) {
        // Read the time stamps
        using CharArray = NDArray<char>;
        Vector<CharArray> array_ts(1);
        Vector<int> success(1);
        ReadNetCDFFile(fname, {"Times"}, array_ts, success);

        int ntimes = array_ts[0].get_vshape()[0];
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

        // Extract the first time entry
        std::string date(&timeStamps[0][0], &timeStamps[0][dateStrLen-1]+1);
        auto epochTime = getEpochTime(date, dateTimeFormat);
        Print() << "  metgrid datetime 0 : " << date << " " << epochTime << std::endl;
        NC_epochTime = static_cast<Real>(epochTime);

        amrex::Print() << "Have read start_time string at level "<< lev << " is " << date << std::endl;
        amrex::Print() << "Have read start_time number at level "<< lev << " is " << NC_epochTime << std::endl;
    }

    amrex::ParallelDescriptor::Bcast(&NC_epochTime,1,amrex::ParallelDescriptor::IOProcessorNumber());

    return NC_epochTime;
}

/**
 * Initializes ERF data using metgrid data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_metgrid (int lev)
{
    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
    if (use_moisture) {
        Print() << "Init with met_em with valid moisture model." << std::endl;
    } else {
        Print() << "Init with met_em without moisture model." << std::endl;
    }

    int ntimes = num_files_at_level[lev];
    Print() << ntimes << " met_em.d0" << lev+1 << "*.nc files are listed" << std::endl;

    if (nc_init_file.empty())
        Error("NetCDF initialization file name must be provided via input");

    if (nc_init_file[lev].empty())
        Error("NetCDF initialization file name must be provided via input");

    // At least two met_em files are necessary to calculate tendency terms.
    if (lev == 0) {
        AMREX_ALWAYS_ASSERT(ntimes >= 2);
    } else {
        AMREX_ALWAYS_ASSERT(ntimes == 1);
    }

    // Ensure a reasonable value for the order of the vertical interpolation scheme.
    AMREX_ALWAYS_ASSERT(metgrid_order > 0 && metgrid_order <= 9);

    // Odd behavior can occur if not using the surface and also omitting some near-surface levels.
    if (metgrid_force_sfc_k > 0)
        AMREX_ALWAYS_ASSERT(metgrid_use_sfc);

    // Size the SST and LANDMASK
      sst_lev[lev].resize(ntimes);
      tsk_lev[lev].resize(ntimes);
    lmask_lev[lev].resize(ntimes);

    // *** FArrayBox's at this level for holding the metgrid data
    FArrayBox NC_xvel_fab;
    FArrayBox NC_yvel_fab;
    FArrayBox NC_temp_fab;
    FArrayBox NC_rhum_fab;
    FArrayBox NC_pres_fab;
    FArrayBox NC_ght_fab;
    FArrayBox NC_hgt_fab;
    FArrayBox NC_psfc_fab;
    FArrayBox NC_MSFU_fab;
    FArrayBox NC_MSFV_fab;
    FArrayBox NC_MSFM_fab;
    FArrayBox NC_sst_fab;
    FArrayBox NC_tsk_fab;
    FArrayBox NC_LAT_fab;
    FArrayBox NC_LON_fab;

    // *** IArrayBox's at this level for holding mask data
    IArrayBox NC_lmask_iab;

    // *** Variables at this level for holding metgrid file global attributes
    int flag_psfc = 0;
    int flag_msf = 0;
    int flag_sst = 0;
    int flag_tsk = 0;
    int flag_lmask = 0;
    int NC_nx;
    int NC_ny;
    Real NC_dx;
    Real NC_dy;
    Vector<std::string> NC_dateTime; NC_dateTime.resize( ntimes);
    Vector<Real> NC_epochTime;       NC_epochTime.resize(ntimes);

    // Define the arena to be used for data allocation
    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    // Make sure this lives on CPU and GPU
    Arena_Used = The_Pinned_Arena();
#endif

    if (metgrid_debug_quiescent)  Print() << "metgrid_debug_quiescent  = true" << std::endl;
    if (metgrid_debug_isothermal) Print() << "metgrid_debug_isothermal = true" << std::endl;
    if (metgrid_debug_dry)        Print() << "metgrid_debug_dry        = true" << std::endl;
    if (metgrid_debug_psfc)       Print() << "metgrid_debug_psfc       = true" << std::endl;
    if (metgrid_debug_msf)        Print() << "metgrid_debug_msf        = true" << std::endl;
    if (metgrid_interp_theta)     Print() << "metgrid_interp_theta     = true" << std::endl;
    if (metgrid_basic_linear)     Print() << "metgrid_basic_linear     = true" << std::endl;

    auto& lev_new = vars_new[lev];

    z_phys_nd[lev]->setVal(0);

    AMREX_ALWAYS_ASSERT(SolverChoice::terrain_type != TerrainType::None);

    auto& dm = lev_new[Vars::cons].DistributionMap();
    auto ngv = lev_new[Vars::cons].nGrowVect(); ngv[2] = 0;

    int i_lo = boxes_at_level[lev][0].smallEnd(0); int i_hi = boxes_at_level[lev][0].bigEnd(0);
    int j_lo = boxes_at_level[lev][0].smallEnd(1); int j_hi = boxes_at_level[lev][0].bigEnd(1);

    // Set up FABs to hold data that will be used to set lateral boundary conditions.
    int MetGridBdyEnd = MetGridBdyVars::NumTypes-1;
    if (use_moisture) MetGridBdyEnd = MetGridBdyVars::NumTypes;

    if (lev == 0) {
        // Zero out the bdy data to start off
        bdy_data_xlo.resize(ntimes);
        bdy_data_xhi.resize(ntimes);
        bdy_data_ylo.resize(ntimes);
        bdy_data_yhi.resize(ntimes);

        for (int itime(0); itime < ntimes; itime++) {
            bdy_data_xlo[itime].resize(MetGridBdyEnd);
            bdy_data_xhi[itime].resize(MetGridBdyEnd);
            bdy_data_ylo[itime].resize(MetGridBdyEnd);
            bdy_data_yhi[itime].resize(MetGridBdyEnd);

            // Build the boxes for each BDY FAB
            const auto& lo = geom[lev].Domain().loVect();
            const auto& hi = geom[lev].Domain().hiVect();
            IntVect plo(lo);
            IntVect phi(hi);

            plo[0] = lo[0];              plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = lo[0]+real_width-1; phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xlo(plo, phi);
            Box xlo_plane_no_stag(pbx_xlo);
            Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
            Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});

            plo[0] = hi[0]-real_width+1; plo[1] = lo[1]; plo[2] = lo[2];
            phi[0] = hi[0];              phi[1] = hi[1]; phi[2] = hi[2];
            const Box pbx_xhi(plo, phi);
            Box xhi_plane_no_stag(pbx_xhi);
            Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
            Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});

            plo[1] = lo[1];              plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = lo[1]+real_width-1; phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_ylo(plo, phi);
            Box ylo_plane_no_stag(pbx_ylo);
            Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
            Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);

            plo[1] = hi[1]-real_width+1; plo[0] = lo[0]; plo[2] = lo[2];
            phi[1] = hi[1];              phi[0] = hi[0]; phi[2] = hi[2];
            const Box pbx_yhi(plo, phi);
            Box yhi_plane_no_stag(pbx_yhi);
            Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
            Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);

            for (int nvar(0); nvar<MetGridBdyEnd; ++nvar) {
                if (nvar==MetGridBdyVars::U) {
                    bdy_data_xlo[itime][nvar].resize(xlo_plane_x_stag, 1, Arena_Used);
                    bdy_data_xhi[itime][nvar].resize(xhi_plane_x_stag, 1, Arena_Used);
                    bdy_data_ylo[itime][nvar].resize(ylo_plane_x_stag, 1, Arena_Used);
                    bdy_data_yhi[itime][nvar].resize(yhi_plane_x_stag, 1, Arena_Used);
                } else if (nvar==MetGridBdyVars::V) {
                    bdy_data_xlo[itime][nvar].resize(xlo_plane_y_stag, 1, Arena_Used);
                    bdy_data_xhi[itime][nvar].resize(xhi_plane_y_stag, 1, Arena_Used);
                    bdy_data_ylo[itime][nvar].resize(ylo_plane_y_stag, 1, Arena_Used);
                    bdy_data_yhi[itime][nvar].resize(yhi_plane_y_stag, 1, Arena_Used);
                } else {
                    bdy_data_xlo[itime][nvar].resize(xlo_plane_no_stag, 1, Arena_Used);
                    bdy_data_xhi[itime][nvar].resize(xhi_plane_no_stag, 1, Arena_Used);
                    bdy_data_ylo[itime][nvar].resize(ylo_plane_no_stag, 1, Arena_Used);
                    bdy_data_yhi[itime][nvar].resize(yhi_plane_no_stag, 1, Arena_Used);
                }
                bdy_data_xlo[itime][nvar].template setVal<RunOn::Device>(0);
                bdy_data_xhi[itime][nvar].template setVal<RunOn::Device>(0);
                bdy_data_ylo[itime][nvar].template setVal<RunOn::Device>(0);
                bdy_data_yhi[itime][nvar].template setVal<RunOn::Device>(0);
            }
        } // itime
    } // lev==0

    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    for (int itime(0); itime < ntimes; itime++) {
        Print() << " init_from_metgrid: reading nc_init_file[" << lev << "][" << itime << "]\t" << nc_init_file[lev][itime] << std::endl;
        read_from_metgrid(lev, itime,
                          boxes_at_level[lev][0], nc_init_file[lev][itime],
                          NC_dateTime[itime],  NC_epochTime[itime],
                          flag_psfc,    flag_msf,
                          flag_sst,     flag_tsk,    flag_lmask,
                          NC_nx,        NC_ny,       NC_dx,       NC_dy,
                          NC_xvel_fab,  NC_yvel_fab,
                          NC_temp_fab,  NC_rhum_fab, NC_pres_fab,
                          NC_ght_fab,   NC_hgt_fab,  NC_psfc_fab,
                          NC_MSFU_fab,  NC_MSFV_fab, NC_MSFM_fab,
                          NC_sst_fab,   NC_tsk_fab,  NC_LAT_fab,  NC_LON_fab,
                          NC_lmask_iab, geom[lev]);

        if (lev == 0) {
            if (itime == 0) {
                // Start at the earliest time in nc_init_file[lev].
                start_bdy_time = NC_epochTime[itime];

                //
                // Note that t_new and t_old carry *elapsed* time, not total time
                //
                Print() << "start_bdy_time is " << std::setprecision(timeprecision) << start_bdy_time
                        << " from metgrid file but note that time variable in simulation is elapsed time" << std::endl;
                t_new[lev] = zero;
                t_old[lev] = -Real(1.e200);
            } else {
                // Verify that files in nc_init_file[lev] are ordered from earliest to latest.
                AMREX_ALWAYS_ASSERT(NC_epochTime[itime] > NC_epochTime[itime-1]);

                // Determine the spacing between met_em files.
                bdy_time_interval = NC_epochTime[1]-NC_epochTime[0];

                // Verify that met_em files have even spacing in time.
                Real NC_dt = NC_epochTime[itime]-NC_epochTime[itime-1];
                Print() << " " << nc_init_file[lev][itime-1] << " / " << nc_init_file[lev][itime] << " are " << NC_dt << " seconds apart" << std::endl;
                if (NC_dt != bdy_time_interval) Error("Time interval between consecutive met_em files must be consistent.");
            } // itime==0

            if (itime == ntimes-1) {
                final_bdy_time = NC_epochTime[itime];
                Print() << "final_bdy_time is " << std::setprecision(timeprecision) << final_bdy_time
                        << " from metgrid file but note that time variable in simulation is elapsed time" << std::endl;
            }
        } // lev==0


        // Verify that the grid size and resolution from met_em file matches that in geom (from ERF inputs file).
        Real tol   = Real(1.0e-3);
        AMREX_ALWAYS_ASSERT(std::fabs(geom[lev].CellSizeArray()[0]-NC_dx) < tol);
        AMREX_ALWAYS_ASSERT(std::fabs(geom[lev].CellSizeArray()[1]-NC_dy) < tol);
        // NC_nx-2 because NC_nx is the number of staggered grid points indexed from one
        AMREX_ALWAYS_ASSERT(i_hi-i_lo == NC_nx-2);
        // NC_ny-2 because NC_ny is the number of staggered grid points indexed from one
        AMREX_ALWAYS_ASSERT(j_hi-j_lo == NC_ny-2);

        if (itime == 0) {

            for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                // This defines only the z(i,j,0) values given the FAB filled from the NetCDF input
                FArrayBox& z_phys_nd_fab = (*z_phys_nd[lev])[mfi];
                init_terrain_from_metgrid(z_phys_nd_fab, NC_hgt_fab);
            } // mf

            // This defines all the z(i,j,k) values given z(i,j,0) from above.
            make_terrain_fitted_coords(lev, geom[lev], *z_phys_nd[lev], zlevels_stag[lev], phys_bc_type);

            // This makes the Jacobian.
            make_J(geom[lev], *z_phys_nd[lev], *detJ_cc[lev]);
            make_areas(geom[lev], *z_phys_nd[lev], *ax[lev], *ay[lev], *az[lev]);

            // This defines z at w-cell faces.
            make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);

        } // itime==0

        // Copy LATITUDE, LONGITUDE, SST and LANDMASK data into MF and iMF data structures

        if (flag_sst) {
            sst_lev[lev][itime] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(sst_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(sst_lev[lev][itime]))[mfi];
                FArrayBox& src = NC_sst_fab;
                const Array4<      Real>& dst_arr = dst.array();
                const Array4<const Real>& src_arr = src.const_array();
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    int li = min(max(i, i_lo), i_hi);
                    int lj = min(max(j, j_lo), j_hi);
                    dst_arr(i,j,0) = src_arr(li,lj,0);
                });
            }
            sst_lev[lev][itime]->FillBoundary(geom[lev].periodicity());
        } else {
            sst_lev[lev][itime] = nullptr;
        }

        if (flag_tsk) {
            tsk_lev[lev][itime] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(tsk_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(tsk_lev[lev][itime]))[mfi];
                FArrayBox& src = NC_tsk_fab;
                const Array4<      Real>& dst_arr = dst.array();
                const Array4<const Real>& src_arr = src.const_array();
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    int li = min(max(i, i_lo), i_hi);
                    int lj = min(max(j, j_lo), j_hi);
                    dst_arr(i,j,0) = src_arr(li,lj,0);
                });
            }
            tsk_lev[lev][itime]->FillBoundary(geom[lev].periodicity());
        } else {
            tsk_lev[lev][itime] = nullptr;
        }

        if (flag_lmask) {
            lmask_lev[lev][itime] = std::make_unique<iMultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(lmask_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                IArrayBox& dst = (*(lmask_lev[lev][itime]))[mfi];
                IArrayBox& src = NC_lmask_iab;
                const Array4<      int>& dst_arr = dst.array();
                const Array4<const int>& src_arr = src.const_array();
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    int li = min(max(i, i_lo), i_hi);
                    int lj = min(max(j, j_lo), j_hi);
                    dst_arr(i,j,0) = src_arr(li,lj,0);
                });
            }
            lmask_lev[lev][itime]->FillBoundary(geom[lev].periodicity());
        }

        if (itime == 0) {
            lat_m[lev]    = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            sinPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            cosPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(lat_m[lev]))[mfi];
                FArrayBox& src = NC_LAT_fab;
                const Array4<      Real>& sin_arr = (sinPhi_m[lev])->array(mfi);
                const Array4<      Real>& cos_arr = (cosPhi_m[lev])->array(mfi);
                const Array4<      Real>& dst_arr = dst.array();
                const Array4<const Real>& src_arr = src.const_array();
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    int li = min(max(i, i_lo), i_hi);
                    int lj = min(max(j, j_lo), j_hi);
                    dst_arr(i,j,0) = src_arr(li,lj,0);

                    Real lat_rad = dst_arr(i,j,0) * (PI/Real(180.));
                    sin_arr(i,j,0) = std::sin(lat_rad);
                    cos_arr(i,j,0) = std::cos(lat_rad);
                });
            }

            lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(lon_m[lev]))[mfi];
                FArrayBox& src = NC_LON_fab;
                const Array4<      Real>& dst_arr = dst.array();
                const Array4<const Real>& src_arr = src.const_array();
                ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    int li = min(max(i, i_lo), i_hi);
                    int lj = min(max(j, j_lo), j_hi);
                    dst_arr(i,j,0) = src_arr(li,lj,0);
                });
            }
        } // itime==0

        // Set up a temporary MultiFab with  mixing ratio and potential
        // temperature on the origin model vertical grid. This is
        // necessary for input data that has relative humidity and
        // temperature, but does not include mixing ratio and potential
        // temperature. Also set up a temporary MultiFab with pressure
        // and temperature on the ERF vertical grid. This is necessary
        // for calculating potential temperature from interpolated
        // pressure and temperature.
        MultiFab tmp_src; // variables on origin model vertical levels
        MultiFab tmp_dst; // variables on ERF vertical levels
        {   // Start with the level's boxArray then modify dim-2 to span
            // from 0 to the number of levels in the origin model data.
            // boxArray->boxList->boxArray gets around a shallow copy problem.
            BoxArray ba_dst = lev_new[Vars::cons].boxArray();
            Box box_nc = NC_pres_fab.box();
            int kmax = ubound(box_nc).z;
            BoxList bl_src = ba_dst.boxList();
            for (auto& b : bl_src) {
                b.setSmall(2, 0);
                b.setBig(2, kmax);
            }
            BoxArray ba_src(std::move(bl_src));
            tmp_src.define(ba_src, dm, MetGridTmpSrcVars::NumTypes, 0);
            tmp_dst.define(ba_dst, dm, MetGridTmpDstVars::NumTypes, 0);
        }

        const Real l_rdOcp = solverChoice.rdOcp;
        std::unique_ptr<iMultiFab> mask_c = OwnerMask(lev_new[Vars::cons], geom[lev].periodicity());//, lev_new[Vars::cons].nGrowVect());
        std::unique_ptr<iMultiFab> mask_u = OwnerMask(lev_new[Vars::xvel], geom[lev].periodicity());//, lev_new[Vars::xvel].nGrowVect());
        std::unique_ptr<iMultiFab> mask_v = OwnerMask(lev_new[Vars::yvel], geom[lev].periodicity());//, lev_new[Vars::yvel].nGrowVect());
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for ( MFIter mfi(lev_new[Vars::cons], TileNoZ()); mfi.isValid(); ++mfi ) {
            Box tbxc = mfi.tilebox();
            Box tbxu = mfi.tilebox(IntVect(1,0,0));
            Box tbxv = mfi.tilebox(IntVect(0,1,0));

            // Define FABs for holding some of the initial data
            FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
            FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
            FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
            FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];
            FArrayBox& z_phys_nd_fab = (*z_phys_nd[lev])[mfi];

            // FABs holding temporary vars on origin and ERF grids, respectively.
            FArrayBox &tmp_src_fab = tmp_src[mfi];
            FArrayBox &tmp_dst_fab = tmp_dst[mfi];

            const Array4<const int>& mask_c_arr = mask_c->const_array(mfi);
            const Array4<const int>& mask_u_arr = mask_u->const_array(mfi);
            const Array4<const int>& mask_v_arr = mask_v->const_array(mfi);

            // Fill state data using origin data (initialization and BC arrays)
            //     x_vel   interpolated from origin levels
            //     y_vel   interpolated from origin levels
            //     z_vel   set to zero
            //     theta   (metgrid_interp_theta) calculate on origin levels then interpolate
            //             (!metgrid_interp_theta) interpolate P and T then calculate on ERF levels
            //     mxrat   convert RH -> Q on origin levels then interpolate
            init_state_from_metgrid(lev, itime, use_moisture, metgrid_interp_theta,
                                    metgrid_debug_quiescent, metgrid_debug_isothermal,
                                    metgrid_debug_dry, metgrid_basic_linear,
                                    metgrid_use_below_sfc, metgrid_use_sfc,
                                    metgrid_retain_sfc, metgrid_proximity,
                                    metgrid_order, metgrid_force_sfc_k, l_rdOcp,
                                    tbxc, tbxu, tbxv,
                                    cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                    z_phys_nd_fab,
                                    NC_ght_fab, NC_xvel_fab,
                                    NC_yvel_fab, NC_temp_fab, NC_rhum_fab,
                                    NC_pres_fab,
                                    tmp_src_fab, tmp_dst_fab,
                                    bdy_data_xlo, bdy_data_xhi,
                                    bdy_data_ylo, bdy_data_yhi,
                                    mask_c_arr, mask_u_arr, mask_v_arr);
        } // mf

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif

        if (itime == 0) {

            // Use map scale factors directly from the met_em files
            Print() << "[init_msfs_from_metgrid] lev = " << lev << ", itime = " << itime << std::endl;
            for ( MFIter mfi(*mapfac[lev][MapFacType::u_x], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                // Define fabs for holding the initial data
                FArrayBox &msfu_fab = (*mapfac[lev][MapFacType::u_x])[mfi];
                FArrayBox &msfv_fab = (*mapfac[lev][MapFacType::v_x])[mfi];
                FArrayBox &msfm_fab = (*mapfac[lev][MapFacType::m_x])[mfi];

                init_msfs_from_metgrid(metgrid_debug_msf,
                                       msfu_fab, msfv_fab, msfm_fab, flag_msf,
                                       NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
            } // mf

            Print() << "[init_base_state_from_metgrid] lev = " << lev << ", itime = " << itime << std::endl;
            for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                FArrayBox&     p_hse_fab = p_hse[mfi];
                FArrayBox&    pi_hse_fab = pi_hse[mfi];
                FArrayBox&    th_hse_fab = th_hse[mfi];
                FArrayBox&    qv_hse_fab = qv_hse[mfi];
                FArrayBox&     r_hse_fab = r_hse[mfi];
                FArrayBox&      cons_fab = lev_new[Vars::cons][mfi];
                FArrayBox& z_phys_cc_fab = (*z_phys_cc[lev])[mfi];
                FArrayBox& z_phys_nd_fab = (*z_phys_nd[lev])[mfi];

                // Fill base state data using origin data
                //     p_hse     calculate moist hydrostatic pressure
                //     r_hse     calculate moist hydrostatic density
                //     pi_hse    calculate Exner term given pressure
                //     th_hse    calculate potential temperature
                //     qv_hse    calculate qv
                const Box valid_bx = mfi.validbox();
                init_base_state_from_metgrid(use_moisture, metgrid_debug_psfc,
                                             l_rdOcp, valid_bx, flag_psfc, cons_fab,
                                             r_hse_fab, p_hse_fab, pi_hse_fab, th_hse_fab,
                                             qv_hse_fab, z_phys_nd_fab, z_phys_cc_fab, NC_psfc_fab);
            } // mf

        } // itime==0

    } // itime

    // FillBoundary to populate the internal halo cells
     r_hse.FillBoundary(geom[lev].periodicity());
     p_hse.FillBoundary(geom[lev].periodicity());
    pi_hse.FillBoundary(geom[lev].periodicity());
    th_hse.FillBoundary(geom[lev].periodicity());
    qv_hse.FillBoundary(geom[lev].periodicity());

    // NOTE: fabs_for_bcs is defined over the whole domain on each rank.
    //       However, the operations needed to define the data on the ERF
    //       grid are done over MultiFab boxes that are local to the rank.
    //       So when we save the data in fabs_for_bc, only regions owned
    //       by the rank are populated. Use an allreduce sum to make the
    //       complete data set; initialized to 0 above.
    if (lev == 0) {
        for (int itime(0); itime < ntimes; itime++) {
            for (int nvar(0); nvar<MetGridBdyEnd; ++nvar) {
                ParallelAllReduce::Sum(bdy_data_xlo[itime][nvar].dataPtr(),
                                       bdy_data_xlo[itime][nvar].size(),
                                       ParallelContext::CommunicatorAll());
                ParallelAllReduce::Sum(bdy_data_xhi[itime][nvar].dataPtr(),
                                       bdy_data_xhi[itime][nvar].size(),
                                       ParallelContext::CommunicatorAll());
                ParallelAllReduce::Sum(bdy_data_ylo[itime][nvar].dataPtr(),
                                       bdy_data_ylo[itime][nvar].size(),
                                       ParallelContext::CommunicatorAll());
                ParallelAllReduce::Sum(bdy_data_yhi[itime][nvar].dataPtr(),
                                       bdy_data_yhi[itime][nvar].size(),
                                       ParallelContext::CommunicatorAll());
            } // nvar
        } // itime
    } // lev==0

    Print() << "Running with relaxation width: " << real_width << std::endl;
}

/**
 * Helper function to initialize terrain nodal z coordinates given metgrid data.
 *
 * @param z_phys_nd_fab FArrayBox (Fab) holding the nodal z coordinates for terrain data we want to fill
 * @param NC_hgt_fab FArrayBox (Fab) holding height data read from the first NetCDF file of metgrid data
 */
void
init_terrain_from_metgrid (FArrayBox& z_phys_nd_fab,
                           FArrayBox& NC_hgt_fab)
{
   // This copies from NC_zphys on z-faces to z_phys_nd on nodes
   const Array4<Real      >&      z_arr = z_phys_nd_fab.array();
   const Array4<Real const>& nc_hgt_arr = NC_hgt_fab.const_array();

   const Box z_hgt_box = NC_hgt_fab.box();

   int ilo = z_hgt_box.smallEnd()[0];
   int ihi = z_hgt_box.bigEnd()[0];
   int jlo = z_hgt_box.smallEnd()[1];
   int jhi = z_hgt_box.bigEnd()[1];

   Box z_phys_box = z_phys_nd_fab.box();
   Box from_box = surroundingNodes(NC_hgt_fab.box());
   from_box.growHi(2,-1);

   Box bx = z_phys_box & from_box;

   ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
       int ii = std::max(std::min(i,ihi-1),ilo+1);
       int jj = std::max(std::min(j,jhi-1),jlo+1);
       z_arr(i,j,k) =  fourth * ( nc_hgt_arr (ii,jj  ,k) + nc_hgt_arr(ii-1,jj  ,k) +
                                nc_hgt_arr (ii,jj-1,k) + nc_hgt_arr(ii-1,jj-1,k) );
   });
}

/**
 * Helper function to initialize state and velocity data read from metgrid data.
 *
 * @param lev int
 * @param itime int
 * @param use_moisture bool True if solverChoice.moisture_type != MoistureType::None
 * @param metgrid_interp_theta bool calculate theta on origin levels, then interpolate
 * @param metgrid_debug_quiescent bool overwrite u and v with zero
 * @param metgrid_debug_isothermal bool overwrite theta with Real(300.0)
 * @param metgrid_debug_dry bool overwrite qv with zero
 * @param metgrid_basic_linear bool linear interpolation without quality control
 * @param metgrid_use_below_sfc bool quality control includes points below the surface
 * @param metgrid_use_sfc bool quality control includes the point at the surface
 * @param metgrid_retain_sfc bool set the lowest level of interpolated field to the surface value
 * @param metgrid_proximity Real pressure difference for quality control pruning
 * @param metgrid_order int interpolation order
 * @param metgrid_force_sfc_k int lower levels pruned by quality control
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param state_fab FArrayBox holding the state data to initialize
 * @param x_vel_fab FArrayBox holding the x-velocity data to initialize
 * @param y_vel_fab FArrayBox holding the y-velocity data to initialize
 * @param z_vel_fab FArrayBox holding the z-velocity data to initialize
 * @param z_phys_nd_fab FArrayBox holding nodal z coordinate data for terrain
 * @param NC_ght_fab  FArrayBox object holding metgrid data for height of cell centers
 * @param NC_xvel_fab FArrayBox object holding metgrid data for x-velocity
 * @param NC_yvel_fab FArrayBox object holding metgrid data for y-velocity
 * @param NC_zvel_fab FArrayBox object holding metgrid data for z-velocity
 * @param NC_temp_fab FArrayBox object holding metgrid data for temperature
 * @param NC_rhum_fab FArrayBox object holding metgrid data for relative humidity
 * @param NC_pres_fab FArrayBox object holding metgrid data for pressure
 * @param p_interp_fab FArrayBox object
 * @param t_interp_fab FArrayBox object
 * @param theta_fab FArrayBox object holding potential temperature calculated from temperature and pressure
 * @param mxrat_fab FArrayBox object holding vapor mixing ratio calculated from relative humidity
 * @param fabs_for_bcs Vector of Vector of FArrayBox objects holding MetGridBdyVars at each met_em time.
 * @param mask_c_arr
 * @param mask_u_arr
 * @param mask_v_arr
 */
void
init_state_from_metgrid (const int  lev,
                         const int  itime,
                         const bool use_moisture,
                         const bool metgrid_interp_theta,
                         const bool metgrid_debug_quiescent,
                         const bool metgrid_debug_isothermal,
                         const bool metgrid_debug_dry,
                         const bool metgrid_basic_linear,
                         const bool metgrid_use_below_sfc,
                         const bool metgrid_use_sfc,
                         const bool metgrid_retain_sfc,
                         const Real metgrid_proximity,
                         const int  metgrid_order,
                         const int  metgrid_force_sfc_k,
                         const Real l_rdOcp,
                         Box& tbxc,
                         Box& tbxu,
                         Box& tbxv,
                         FArrayBox& state_fab,
                         FArrayBox& x_vel_fab,
                         FArrayBox& y_vel_fab,
                         FArrayBox& z_vel_fab,
                         FArrayBox& z_phys_nd_fab,
                         const FArrayBox& NC_ght_fab,
                         const FArrayBox& NC_xvel_fab,
                         const FArrayBox& NC_yvel_fab,
                         const FArrayBox& NC_temp_fab,
                         const FArrayBox& NC_rhum_fab,
                         const FArrayBox& NC_pres_fab,
                         FArrayBox& tmp_src_fab,
                         FArrayBox& tmp_dst_fab,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_xlo,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_xhi,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_ylo,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_yhi,
                         const Array4<const int>& mask_c_arr,
                         const Array4<const int>& mask_u_arr,
                         const Array4<const int>& mask_v_arr)
{
    bool metgrid_exp_interp = false; // interpolate w.r.t. exp(z) for non-pressure variables.

    // Staging FABs for velocities. Interpolated velocities will be saved to these,
    // then copied into the initialization and boundaries as necessary.
    FArrayBox u_staging(x_vel_fab.box(), 1, The_Async_Arena());
    FArrayBox v_staging(y_vel_fab.box(), 1, The_Async_Arena());

    // ********************************************************
    // U
    // ********************************************************
    {
#ifndef AMREX_USE_GPU
    Print() << "[init_state_from_metgrid] vertical interpolation of u-velocity, lev = " << lev << ", itime = " << itime << std::endl;
#endif
    Box bx2d = NC_xvel_fab.box() & tbxu;
    bx2d.setRange(2,0);
    auto const orig_data = NC_xvel_fab.const_array();
    auto const orig_z    = NC_ght_fab.const_array();
    auto       new_data  = u_staging.array();
    auto const new_z     = z_phys_nd_fab.const_array();
    int kmax = ubound(tbxu).z;
    int src_indx = 0;
    int tmp_indx = 0;
    int dst_indx = 0;

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (metgrid_debug_quiescent) { // Debugging option to run quiescent.
            for (int k(0); k<=kmax; k++) {
                new_data(i,j,k,tmp_indx) = zero;
            }
        } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
            for (int k(0); k<=kmax; k++) {
                Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'X',src_indx,orig_z,orig_data,new_z);
                new_data(i,j,k,tmp_indx) = Interp_Val;
            }
        } else { // Vertical interpolation and quality control similar to that from WRF.
            interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                       metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                       metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                       'U', 'X', orig_z, orig_data, new_z, new_data);
        }
    });

    if (itime == 0) { // Copy results for initialization.
        auto       dst = x_vel_fab.array();
        auto const tmp = u_staging.const_array();
        ParallelFor(tbxu, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dst(i,j,k,dst_indx) = tmp(i,j,k,tmp_indx);
        });
    }
    }


    // ********************************************************
    // V
    // ********************************************************
    {
#ifndef AMREX_USE_GPU
    Print() << "[init_state_from_metgrid] vertical interpolation of v-velocity, lev = " << lev << ", itime = " << itime << std::endl;
#endif
    Box bx2d = NC_yvel_fab.box() & tbxv;
    bx2d.setRange(2,0);
    auto const orig_data = NC_yvel_fab.const_array();
    auto const orig_z    = NC_ght_fab.const_array();
    auto       new_data  = v_staging.array();
    auto const new_z     = z_phys_nd_fab.const_array();
    int kmax = ubound(tbxv).z;
    int src_indx = 0;
    int tmp_indx = 0;
    int dst_indx = 0;

    ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
    {
        if (metgrid_debug_quiescent) { // Debugging option to run quiescent.
            for (int k(0); k<=kmax; k++) {
                new_data(i,j,k,tmp_indx) = zero;
            }
        } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
            for (int k(0); k<=kmax; k++) {
                Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'Y',src_indx,orig_z,orig_data,new_z);
                new_data(i,j,k,tmp_indx) = Interp_Val;
            }
        } else {
            interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                       metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                       metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                       'V', 'Y', orig_z, orig_data, new_z, new_data);
        }
    });

    if (itime == 0) { // Copy results for initialization.
        auto       dst = y_vel_fab.array();
        auto const tmp = v_staging.const_array();
        ParallelFor(tbxv, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dst(i,j,k,dst_indx) = tmp(i,j,k,tmp_indx);
        });
    }
    }


    // ********************************************************
    // W
    // ********************************************************
    if (itime == 0) { // update at initialization
        z_vel_fab.template setVal<RunOn::Device>(0);
    }


    // ********************************************************
    // Initialize all state_fab variables to zero
    // ********************************************************
    if (itime == 0) { // update at initialization
        state_fab.template setVal<RunOn::Device>(0);
    }


    // ********************************************************
    // theta
    // ********************************************************
    if (metgrid_interp_theta) {
        // Calculate potential temperature on the origin model vertical levels
        // then interpolate that onto the ERF vertical levels.

        { // calculate potential temperature.
            Box bx = NC_temp_fab.box() & tbxc;
            auto const temp  = NC_temp_fab.const_array();
            auto const pres  = NC_pres_fab.const_array();
            auto       tmp_src_arr = tmp_src_fab.array();
            int src_indx = MetGridTmpSrcVars::Theta;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                tmp_src_arr(i,j,k,src_indx) = getThgivenTandP(temp(i,j,k),pres(i,j,k),l_rdOcp);
            });
        }

        { // vertical interpolation of potential temperature.
#ifndef AMREX_USE_GPU
        Print() << "[init_state_from_metgrid] vertical interpolation of potential temperature, lev = " << lev << ", itime = " << itime << std::endl;
#endif
        Box bx2d = NC_temp_fab.box() & tbxc;
        bx2d.setRange(2,0);
        auto const orig_data = tmp_src_fab.const_array();
        auto const orig_z    = NC_ght_fab.const_array();
        auto       new_data  = tmp_dst_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();
        int kmax = amrex::ubound(tbxc).z;
        int src_indx = MetGridTmpSrcVars::Theta;
        int tmp_indx = MetGridTmpDstVars::Theta;
        int dst_indx = RhoTheta_comp;

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (metgrid_debug_isothermal) { // Debugging option to run isothermal.
                for (int k(0); k<=kmax; k++) {
                    new_data(i,j,k,tmp_indx) = Real(300.0);
                }
            } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                for (int k(0); k<=kmax; k++) {
                    Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',src_indx,orig_z,orig_data,new_z);
                    new_data(i,j,k,tmp_indx) = Interp_Val;
                }
            } else { // Vertical interpolation and quality control similar to that from WRF.
                interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                           metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                           metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                           'T', 'M', orig_z, orig_data, new_z, new_data);
            }
        });

        if (itime == 0) { // Copy results for initialization.
            auto       dst = state_fab.array();
            auto const tmp = tmp_dst_fab.const_array();
            ParallelFor(tbxc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst(i,j,k,dst_indx) = tmp(i,j,k,tmp_indx);
            });
        }
        }

    } else { // metgrid_interp_theta == false

        { // vertical interpolation of pressure.
#ifndef AMREX_USE_GPU
        Print() << "[init_state_from_metgrid] vertical interpolation of pressure, lev = " << lev << ", itime = " << itime << std::endl;
#endif
        Box bx2d = tmp_dst_fab.box() & tbxc;
        bx2d.setRange(2,0);
        auto const orig_data = NC_pres_fab.const_array();
        auto const orig_z    = NC_ght_fab.const_array();
        auto       new_data  = tmp_dst_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();
        int kmax = ubound(tbxc).z;
        int src_indx = 0;
        int tmp_indx = MetGridTmpDstVars::P;

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                for (int k(0); k<=kmax; k++) {
                    Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',src_indx,orig_z,orig_data,new_z);
                    new_data(i,j,k,tmp_indx) = Interp_Val;
                }
            } else { // Vertical interpolation and quality control similar to that from WRF.
                // Interpolate pressure not w.r.t. z but rather p_0*exp(-CONST_GRAV*z/(t_0*R_d)).
                // This is akin to interpolating in pressure-space assuming a baroclinic atmosphere.
                interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, true,
                                           metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                           metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                           'T', 'M', orig_z, orig_data, new_z, new_data);
            }
        });
        }

        { // vertical interpolation of temperature.
#ifndef AMREX_USE_GPU
        Print() << "[init_state_from_metgrid] vertical interpolation of temperature, lev = " << lev << ", itime = " << itime << std::endl;
#endif
        Box bx2d = tmp_dst_fab.box() & tbxc;
        bx2d.setRange(2,0);
        auto const orig_data = NC_temp_fab.const_array();
        auto const orig_z    = NC_ght_fab.const_array();
        auto       new_data  = tmp_dst_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();
        int kmax = ubound(tbxc).z;
        int src_indx = 0;
        int tmp_indx = MetGridTmpDstVars::T;

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                for (int k(0); k<=kmax; k++) {
                    Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',src_indx,orig_z,orig_data,new_z);
                    new_data(i,j,k,tmp_indx) = Interp_Val;
                }
            } else { // Vertical interpolation and quality control similar to that from WRF.
                // According to WRF's code comments, "It is better to
                // interpolate temperature and potential temperature
                // in LOG(p), regardless of requested default."
                interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, false,
                                           metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                           metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                           'T', 'M', orig_z, orig_data, new_z, new_data);
            }
        });
        }

        { // calculate potential temperature on the ERF vertical levels.
        auto const tmp_dst_arr = tmp_dst_fab.array();
        auto       new_data = tmp_dst_fab.array();
        int T_indx = MetGridTmpDstVars::T;
        int P_indx = MetGridTmpDstVars::P;
        int tmp_indx = MetGridTmpDstVars::Theta;
        int dst_indx = RhoTheta_comp;

        ParallelFor(tbxc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real Calc_Val = getThgivenTandP(tmp_dst_arr(i,j,k,T_indx),tmp_dst_arr(i,j,k,P_indx),l_rdOcp);
            if (metgrid_debug_isothermal) Calc_Val = Real(300.0); // Debugging option to run isothermal.
            new_data(i,j,k,tmp_indx) = Calc_Val;
        });

        if (itime == 0) { // Copy results for initialization.
            auto       dst = state_fab.array();
            auto const tmp = tmp_dst_fab.const_array();
            ParallelFor(tbxc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                dst(i,j,k,dst_indx) = tmp(i,j,k,tmp_indx);
            });
        }
        }

    } // metgrid_interp_theta

    if (use_moisture) {
        // ********************************************************
        // specific humidity / relative humidity / mixing ratio
        // ********************************************************
        // TODO: we will need to check what input data we have for moisture
        // and then, if necessary, compute mixing ratio. For now, we will
        // focus on the case where we have relative humidity. Alternate cases
        // could be specific humidity or a mixing ratio.
        //
        { // calculate vapor mixing ratio from relative humidity.
            Box bx = NC_temp_fab.box() & tbxc;
            auto const rhum  = NC_rhum_fab.const_array();
            auto const temp  = NC_temp_fab.const_array();
            auto const pres  = NC_pres_fab.const_array();
            auto    tmp_src_arr = tmp_src_fab.array();
            int src_indx = MetGridTmpSrcVars::QV;

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                rh_to_mxrat(i, j, k, rhum, temp, pres, src_indx, tmp_src_arr);
            });
        }

        { // vertical interpolation of vapor mixing ratio.
#ifndef AMREX_USE_GPU
            Print() << "[init_state_from_metgrid] vertical interpolation of vapor mixing ratio, lev = " << lev << ", itime = " << itime << std::endl;
#endif
            Box bx2d = NC_temp_fab.box() & tbxc;
            bx2d.setRange(2,0);
            auto const orig_data = tmp_src_fab.const_array();
            auto const orig_z    = NC_ght_fab.const_array();
            auto       new_data  = tmp_dst_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();
            int kmax = ubound(tbxc).z;
            int src_indx = MetGridTmpSrcVars::QV;
            int tmp_indx = MetGridTmpDstVars::QV;
            int dst_indx = RhoQ1_comp;

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                if (metgrid_debug_dry) { // Debugging option to run dry.
                    for (int k(0); k<=kmax; k++) {
                        new_data(i,j,k,tmp_indx)   = zero;
                    }
                } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                    for (int k(0); k<=kmax; k++) {
                        Real Interp_Val  = interpolate_column_metgrid_linear(i,j,k,'M',src_indx,orig_z,orig_data,new_z);
                        new_data(i,j,k,tmp_indx)   = Interp_Val;
                    }
                } else { // Vertical interpolation and quality control similar to that from WRF.
                    interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                               metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                               metgrid_force_sfc_k, i, j, kmax, src_indx, tmp_indx,
                                               'Q', 'M', orig_z, orig_data, new_z, new_data);
                }
            });

            if (itime == 0) { // Copy results for initialization.
                auto       dst = state_fab.array();
                auto const tmp = tmp_dst_fab.const_array();
                ParallelFor(tbxc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    dst(i,j,k,dst_indx) = tmp(i,j,k,tmp_indx);
                });
            }
        }
    } // use_moisture

    if (lev == 0) { // Store boundary data if on level 0.

        { // U
            Box bx2d = NC_xvel_fab.box() & tbxu;
            bx2d.setRange(2,0);
            auto new_data = u_staging.const_array();
            Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::U].box();
            Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::U].box();
            Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::U].box();
            Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::U].box();
            auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::U].array();
            auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::U].array();
            auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::U].array();
            auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::U].array();
            int kmax = ubound(tbxu).z;
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k(0); k<=kmax; k++) {
                    if (mask_u_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_u_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_u_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_u_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = new_data(i,j,k,0);
                }
            });
        } // U

        { // V
            Box bx2d = NC_yvel_fab.box() & tbxv;
            bx2d.setRange(2,0);
            auto new_data = v_staging.const_array();
            Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::V].box();
            Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::V].box();
            Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::V].box();
            Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::V].box();
            auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::V].array();
            auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::V].array();
            auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::V].array();
            auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::V].array();
            int kmax = ubound(tbxv).z;
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k(0); k<=kmax; k++) {
                    if (mask_v_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_v_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_v_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = new_data(i,j,k,0);
                    if (mask_v_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = new_data(i,j,k,0);
                }
            });
        } // V

        { // theta
            Box bx2d = NC_temp_fab.box() & tbxc;
            bx2d.setRange(2,0);
            auto new_data = tmp_dst_fab.const_array();
            Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::T].box();
            Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::T].box();
            Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::T].box();
            Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::T].box();
            auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::T].array();
            auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::T].array();
            auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::T].array();
            auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::T].array();
            int kmax = amrex::ubound(tbxc).z;
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k(0); k<=kmax; k++) {
                    if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::Theta);
                    if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::Theta);
                    if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::Theta);
                    if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::Theta);
                }
            });
        } // theta

        if (use_moisture) {
            Box bx2d = NC_temp_fab.box() & tbxc;
            bx2d.setRange(2,0);
            auto       new_data  = tmp_dst_fab.const_array();
            Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::QV].box();
            Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::QV].box();
            Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::QV].box();
            Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::QV].box();
            auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::QV].array();
            auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::QV].array();
            auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::QV].array();
            auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::QV].array();
            int kmax = ubound(tbxc).z;
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                for (int k(0); k<=kmax; k++) {
                    if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::QV);
                    if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::QV);
                    if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::QV);
                    if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = new_data(i,j,k,MetGridTmpDstVars::QV);
                }
            });
        } // use_moisture

    } // lev==0

}


/**
 * Helper function for initializing hydrostatic base state data from metgrid data
 *
 * @param use_moisture bool True if solverChoice.moisture_type != MoistureType::None
 * @param metgrid_debug_psfc bool use 10**5 Pa as surface pressure when True
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param valid_bx Box specifying the index space we are to initialize
 * @param flag_psfc Int 1 if surface pressure is in metgrid data, 0 otherwise
 * @param state_fab FArrayBox object holding the state data to initialize
 * @param r_hse_fab FArrayBox object holding the hydrostatic base state density we are initializing
 * @param p_hse_fab FArrayBox object holding the hydrostatic base state pressure we are initializing
 * @param pi_hse_fab FArrayBox object holding the hydrostatic base Exner pressure we are initializing
 * @param th_hse_fab FArrayBox object holding the base state potential temperature we are initializing
 * @param qv_hse_fab FArrayBox object holding the base state qv we are initializing
 * @param z_phys_cc_fab FArrayBox object holding cell center z heights for terrain
 * @param NC_psfc_fab FArrayBox object holding metgrid data for surface pressure
 */
void
init_base_state_from_metgrid (const bool use_moisture,
                              const bool metgrid_debug_psfc,
                              const Real l_rdOcp,
                              const Box& valid_bx,
                              const int& flag_psfc,
                              FArrayBox& state_fab,
                              FArrayBox& r_hse_fab,
                              FArrayBox& p_hse_fab,
                              FArrayBox& pi_hse_fab,
                              FArrayBox& th_hse_fab,
                              FArrayBox& qv_hse_fab,
                              FArrayBox& z_phys_nd_fab,
                              FArrayBox& z_phys_cc_fab,
                              const FArrayBox& NC_psfc_fab)
{
    // NOTE: FOEXTRAP is utilized on the validbox but
    //       the FillBoundary call will populate the
    //       internal ghost cells and we are left with
    //       zero gradient at the domain boundaries.

    // Create halo boxes to populate the ghost cells of hse quantities
    Box gvbx_xlo(valid_bx); Box gvbx_xhi(valid_bx);
    Box gvbx_ylo(valid_bx); Box gvbx_yhi(valid_bx);
    Box gvbx_zlo(valid_bx); Box gvbx_zhi(valid_bx);
    gvbx_xlo.grow(IntVect(1,1,0)); gvbx_xhi.grow(IntVect(1,1,0));
    gvbx_ylo.grow(IntVect(1,1,0)); gvbx_yhi.grow(IntVect(1,1,0));
    gvbx_zlo.grow(1); gvbx_zhi.grow(1);
    gvbx_xlo.makeSlab(0,gvbx_xlo.smallEnd(0)); gvbx_xhi.makeSlab(0,gvbx_xhi.bigEnd(0));
    gvbx_ylo.makeSlab(1,gvbx_ylo.smallEnd(1)); gvbx_yhi.makeSlab(1,gvbx_yhi.bigEnd(1));
    gvbx_zlo.makeSlab(2,gvbx_zlo.smallEnd(2)); gvbx_zhi.makeSlab(2,gvbx_zhi.bigEnd(2));

    // Expose for GPU
    Real grav = CONST_GRAV;
    const int maxiter = 20;
    const Real tol    = Real(1.0e-10);

    //***********************************************************************************
    // Set the HSE base state only
    //***********************************************************************************
    {
        Box valid_bx2d = valid_bx;
        valid_bx2d.setRange(2,0);

        // Expose for copy to GPU
        int khi = ubound(valid_bx).z;
        int klo = lbound(valid_bx).z;

        // ARW V4 Constants (5.2.2 Reference State)
        const Real T00       = Real(290.0);
        const Real TLP       = Real(50.0);
        const Real TISO      = Real(200.0);
        const Real TLP_STRAT = Real(-11.0);
        const Real P_STRAT   = Real(0.);

        const Array4<Real>& r_hse_arr  = r_hse_fab.array();
        const Array4<Real>& p_hse_arr  = p_hse_fab.array();
        const Array4<Real>& pi_hse_arr = pi_hse_fab.array();
        const Array4<Real>& th_hse_arr = th_hse_fab.array();
        const Array4<Real>& qv_hse_arr = qv_hse_fab.array();
        auto const z_arr = z_phys_nd_fab.const_array();

        ParallelFor(valid_bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            // Surface values and constants
            Real dz, F, C;
            Real z_hi, Pd_hi, Td_hi, Rd_hi;
            Real z_lo  = Real(0.25)  * ( z_arr(i,j  ,klo  ) + z_arr(i+1,j  ,klo  )
                                       + z_arr(i,j+1,klo  ) + z_arr(i+1,j+1,klo  ) );
            Real Pd_lo = p_0 * std::exp( -T00/TLP + std::sqrt( (T00/TLP)*(T00/TLP) - two * grav * z_lo / (TLP * R_d) ) );
            Real Td_lo;
            if (P_STRAT > zero && Pd_lo < P_STRAT) {
                Td_lo = TISO + TLP_STRAT * std::log(Pd_lo/P_STRAT);
            } else {
                Td_lo = std::max(TISO, T00 + TLP * std::log(Pd_lo/p_0));
            }

            Real Rd_lo = getRhogivenTandPress(Td_lo, Pd_lo);
            for (int k(klo); k<=khi; ++k) {
                // Vertical grid spacing
                z_hi = Real(0.125) * ( z_arr(i,j,k  ) + z_arr(i+1,j,k  ) + z_arr(i,j+1,k  ) + z_arr(i+1,j+1,k  )
                                     + z_arr(i,j,k+1) + z_arr(i+1,j,k+1) + z_arr(i,j+1,k+1) + z_arr(i+1,j+1,k+1) );
                dz   = z_hi - z_lo;

                // Establish known constant
                C  = -Pd_lo + myhalf*Rd_lo*grav*dz;

                // Initial guess and residual
                Pd_hi = Pd_lo;
                Td_hi = Td_lo;
                Rd_hi = Rd_lo;
                F = Pd_hi + myhalf*Rd_hi*grav*dz + C;

                // Iterate to solution
                int niter = 0;
                while (std::fabs(F)>tol && niter<maxiter) {
                    Real dP      = amrex::max(Real(1.0e-3),Real(1.0e-3)*Pd_hi);
                    Real Pd_plus = Pd_hi + dP;
                    Real Td_plus;
                    if (P_STRAT > zero && Pd_plus < P_STRAT) {
                        Td_plus = TISO + TLP_STRAT * std::log(Pd_plus/P_STRAT);
                    } else {
                        Td_plus = std::max(TISO, T00 + TLP * std::log(Pd_plus/p_0));
                    }
                    Real Rd_plus = getRhogivenTandPress(Td_plus, Pd_plus);
                    Real F_plus  = Pd_plus + myhalf*Rd_plus*grav*dz + C;
                    Real dFdP    = (F_plus - F) / dP;

                    Pd_hi -= F / dFdP;
                    if (P_STRAT > zero && Pd_hi < P_STRAT) {
                        Td_hi = TISO + TLP_STRAT * std::log(Pd_hi/P_STRAT);
                    } else {
                        Td_hi = std::max(TISO, T00 + TLP * std::log(Pd_hi/p_0));
                    }
                    Rd_hi   = getRhogivenTandPress(Td_hi, Pd_hi);
                    F       = Pd_hi + myhalf*Rd_hi*grav*dz + C;
                    ++niter;
                }

                // Assign data
                r_hse_arr(i,j,k)  = Rd_hi;
                th_hse_arr(i,j,k) = getThgivenTandP(Td_hi, Pd_hi, l_rdOcp);
                qv_hse_arr(i,j,k) = Real(0.);
                p_hse_arr(i,j,k)  = Pd_hi;
                pi_hse_arr(i,j,k) = getExnergivenP(Pd_hi, l_rdOcp);

                // Transfer solution
                Pd_lo = Pd_hi;
                Td_lo = Td_hi;
                Rd_lo = Rd_hi;
                z_lo  = z_hi;
            }
        });

        // FOEXTRAP hse arrays to fill ghost cells. FillBoundary is
        // called later and will overwrite ghost cell values in the
        // interior of the domain, but the FOEXTRAP values will
        // remain along the lateral domain boundaries.
        ParallelFor(gvbx_xlo, gvbx_xhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int jj = max(j ,valid_bx.smallEnd(1));
                jj = min(jj,valid_bx.bigEnd(1));
            r_hse_arr(i,j,k) =  r_hse_arr(i+1,jj,k);
            p_hse_arr(i,j,k) =  p_hse_arr(i+1,jj,k);
            pi_hse_arr(i,j,k) = pi_hse_arr(i+1,jj,k);
            th_hse_arr(i,j,k) = th_hse_arr(i+1,jj,k);
            qv_hse_arr(i,j,k) = qv_hse_arr(i+1,jj,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int jj = max(j ,valid_bx.smallEnd(1));
                jj = min(jj,valid_bx.bigEnd(1));
             r_hse_arr(i,j,k) =  r_hse_arr(i-1,jj,k);
             p_hse_arr(i,j,k) =  p_hse_arr(i-1,jj,k);
            pi_hse_arr(i,j,k) = pi_hse_arr(i-1,jj,k);
            th_hse_arr(i,j,k) = th_hse_arr(i-1,jj,k);
            qv_hse_arr(i,j,k) = qv_hse_arr(i-1,jj,k);
        });
        ParallelFor(gvbx_ylo, gvbx_yhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r_hse_arr(i,j,k) =  r_hse_arr(i,j+1,k);
            p_hse_arr(i,j,k) =  p_hse_arr(i,j+1,k);
            pi_hse_arr(i,j,k) = pi_hse_arr(i,j+1,k);
            th_hse_arr(i,j,k) = th_hse_arr(i,j+1,k);
            qv_hse_arr(i,j,k) = qv_hse_arr(i,j+1,k);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r_hse_arr(i,j,k) =  r_hse_arr(i,j-1,k);
            p_hse_arr(i,j,k) =  p_hse_arr(i,j-1,k);
            pi_hse_arr(i,j,k) = pi_hse_arr(i,j-1,k);
            th_hse_arr(i,j,k) = th_hse_arr(i,j-1,k);
            qv_hse_arr(i,j,k) = qv_hse_arr(i,j-1,k);
        });
        ParallelFor(gvbx_zlo, gvbx_zhi,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r_hse_arr(i,j,k) =  r_hse_arr(i,j,k+1);
            p_hse_arr(i,j,k) =  p_hse_arr(i,j,k+1);
            pi_hse_arr(i,j,k) = pi_hse_arr(i,j,k+1);
            th_hse_arr(i,j,k) = th_hse_arr(i,j,k+1);
            qv_hse_arr(i,j,k) = qv_hse_arr(i,j,k+1);
        },
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r_hse_arr(i,j,k) =  r_hse_arr(i,j,k-1);
            p_hse_arr(i,j,k) =  p_hse_arr(i,j,k-1);
            pi_hse_arr(i,j,k) = pi_hse_arr(i,j,k-1);
            th_hse_arr(i,j,k) = th_hse_arr(i,j,k-1);
            qv_hse_arr(i,j,k) = qv_hse_arr(i,j,k-1);
        });
    }

    //***********************************************************************************
    // Set the state density only
    //***********************************************************************************
    {
        // Expose for GPU
        int RhoQ_comp = RhoQ1_comp;
        int kmax = ubound(valid_bx).z;

        Box valid_bx2d = valid_bx;
        valid_bx2d.setRange(2,0);
        auto const orig_psfc = NC_psfc_fab.const_array();
        auto       new_data  = state_fab.array();
        auto const new_z     = z_phys_cc_fab.const_array();

        ParallelFor(valid_bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            // Low and Hi column variables
            Real psurf;
            Real z_lo,   z_hi;
            Real p_lo,   p_hi;
            Real qv_lo, qv_hi;
            Real rd_lo, rd_hi;
            Real thetad_lo, thetad_hi;

            // Calculate or use pressure at the surface.
            if (metgrid_debug_psfc) {
                psurf = amrex::Math::powi<5>(10);
            } else if (flag_psfc == 1) {
                psurf = orig_psfc(i,j,0);
            } else {
                z_lo     = new_z(i,j,0);
                Real t_0 = Real(290.0); // WRF's model_config_rec%base_temp
                Real a   = Real(50.0);  // WRF's model_config_rec%base_lapse
                psurf = p_0*std::exp(-t_0/a + std::sqrt(std::pow(t_0/a, two)-two*grav*z_lo/(a*R_d)));
            }
            AMREX_ALWAYS_ASSERT(psurf > zero);
            AMREX_ALWAYS_ASSERT(new_data(i,j,0,RhoTheta_comp) > zero);

            // Iterations for the first CC point that is 1/2 dz off the surface
            {
                z_lo      = new_z(i,j,0);
                qv_lo     = (use_moisture) ? new_data(i,j,0,RhoQ_comp) : zero;
                rd_lo     = zero; // initial guess
                thetad_lo = new_data(i,j,0,RhoTheta_comp);
                // NOTE: The first iteration is from z=0 to z_cc(i,j,0) since the
                //       reference pressure (psurf) is at the ground.
                Real myhalf_dz = z_lo;
                Real qvf     = one+(R_v/R_d)*qv_lo;
                Real thetam  = thetad_lo*qvf;
                for (int it(0); it<maxiter; it++) {
                    p_lo = psurf-myhalf_dz*rd_lo*(one+qv_lo)*grav;
                    if (p_lo < zero) p_lo = zero;
                    rd_lo = (p_0/(R_d*thetam))*std::pow(p_lo/p_0, iGamma);
                } // it

                // Copy solution to state
                new_data(i,j,0,Rho_comp)       = rd_lo;
                new_data(i,j,0,RhoTheta_comp) *= rd_lo;
                if (use_moisture) {
                    new_data(i,j,0,RhoQ_comp) *= rd_lo;
                }
                for (int n(0); n < NSCALARS; n++) {
                    new_data(i,j,0,RhoScalar_comp+n) = zero;
                }
            }

            // Iterations for k \in [1 kmax]
            for (int k(1); k<=kmax; k++) {
                // Known hi data
                z_hi  = new_z(i,j,k);
                qv_hi = (use_moisture) ? new_data(i,j,k,RhoQ_comp) : zero;
                thetad_hi = new_data(i,j,k,RhoTheta_comp);

                // Initial guesses for hi data
                 p_hi = p_lo;
                rd_hi = getRhogivenThetaPress(thetad_hi,
                                              p_hi,
                                              R_d/Cp_d,
                                              qv_hi);

                // Vertical grid spacing
                Real dz = z_hi - z_lo;

                // Establish known constant
                Real rho_tot_lo = rd_lo * (one + qv_lo);
                Real C = -p_lo + myhalf*rho_tot_lo*grav*dz;

                // Initial residual
                Real rho_tot_hi = rd_hi * (one + qv_hi);
                Real F = p_hi + myhalf*rho_tot_hi*grav*dz + C;

                // Do iterations
                if (std::abs(F)>tol) HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                                  grav, C, thetad_hi,
                                                                  qv_hi, qv_hi, p_hi,
                                                                  rd_hi, F);

                // Copy solution to state
                new_data(i,j,k,Rho_comp)       = rd_hi;
                new_data(i,j,k,RhoTheta_comp) *= rd_hi;
                if (use_moisture) {
                    new_data(i,j,k,RhoQ_comp) *= rd_hi;
                }
                for (int n(0); n < NSCALARS; n++) {
                    new_data(i,j,k,RhoScalar_comp+n) = zero;
                }

                // Copy hi to lo
                z_lo  = z_hi;
                p_lo  = p_hi;
                qv_lo = qv_hi;
                rd_lo = rd_hi;
                thetad_lo = thetad_hi;
            }
        });
    }
}


/**
 * Helper function to initialize map factors from metgrid data
 *
 * @param msfu_fab FArrayBox specifying x-velocity map factors
 * @param msfv_fab FArrayBox specifying y-velocity map factors
 * @param msfm_fab FArrayBox specifying z-velocity map factors
 * @param flag_msf Integer 1 if map factors are in metgrid data, 0 otherwise
 * @param NC_MSFU_fab FArrayBox object holding metgrid data for x-velocity map factors
 * @param NC_MSFV_fab FArrayBox object holding metgrid data for y-velocity map factors
 * @param NC_MSFM_fab FArrayBox object holding metgrid data for z-velocity map factors
 */
void
init_msfs_from_metgrid (const bool metgrid_debug_msf,
                        FArrayBox& msfu_fab,
                        FArrayBox& msfv_fab,
                        FArrayBox& msfm_fab,
                        const int& flag_msf,
                        FArrayBox& NC_MSFU_fab,
                        FArrayBox& NC_MSFV_fab,
                        FArrayBox& NC_MSFM_fab)
{
    //
    // FArrayBox to FArrayBox copy does "copy on intersection"
    // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
    //

    // This copies or sets mapfac
    if ((flag_msf == 1) and (!metgrid_debug_msf)) {
        msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab);
        msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab);
        msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab);
    } else {
#ifndef AMREX_USE_GPU
        Print() << " map factors are not present in met_em files. Setting to one" << std::endl;
#endif
        msfm_fab.template setVal<RunOn::Device>(1);
        msfu_fab.template setVal<RunOn::Device>(1);
        msfv_fab.template setVal<RunOn::Device>(1);
    }
}
#endif // ERF_USE_NETCDF
