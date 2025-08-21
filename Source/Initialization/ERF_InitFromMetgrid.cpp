/**
 * \file ERF_InitFromMetgrid.cpp
 */
#include <ERF_Constants.H>
#include <ERF_MetgridUtils.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF
/**
 * Initializes ERF data using metgrid data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_metgrid (int lev)
{
    bool use_moisture = (solverChoice.moisture_type != MoistureType::None);
#ifndef AMREX_USE_GPU
    if (use_moisture) {
        Print() << "Init with met_em with valid moisture model." << std::endl;
    } else {
        Print() << "Init with met_em without moisture model." << std::endl;
    }
#endif

    int ntimes = num_files_at_level[lev];

    if (nc_init_file.empty())
        Error("NetCDF initialization file name must be provided via input");

    if (nc_init_file[lev].empty())
        Error("NetCDF initialization file name must be provided via input");

    // At least two met_em files are necessary to calculate tendency terms.
    AMREX_ALWAYS_ASSERT(ntimes >= 2);

    // At least three points are necessary if there is a relaxation zone.
    if (real_width > real_set_width)
        AMREX_ALWAYS_ASSERT(real_width-real_set_width >= 3);

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
    Vector<FArrayBox> NC_xvel_fab;   NC_xvel_fab.resize(ntimes);
    Vector<FArrayBox> NC_yvel_fab;   NC_yvel_fab.resize(ntimes);
    Vector<FArrayBox> NC_temp_fab;   NC_temp_fab.resize(ntimes);
    Vector<FArrayBox> NC_rhum_fab;   NC_rhum_fab.resize(ntimes);
    Vector<FArrayBox> NC_pres_fab;   NC_pres_fab.resize(ntimes);
    Vector<FArrayBox> NC_ght_fab;    NC_ght_fab.resize( ntimes);
    Vector<FArrayBox> NC_hgt_fab;    NC_hgt_fab.resize( ntimes);
    Vector<FArrayBox> NC_psfc_fab;   NC_psfc_fab.resize(ntimes);
    Vector<FArrayBox> NC_MSFU_fab;   NC_MSFU_fab.resize(ntimes);
    Vector<FArrayBox> NC_MSFV_fab;   NC_MSFV_fab.resize(ntimes);
    Vector<FArrayBox> NC_MSFM_fab;   NC_MSFM_fab.resize(ntimes);
    Vector<FArrayBox> NC_sst_fab;    NC_sst_fab.resize (ntimes);
    Vector<FArrayBox> NC_tsk_fab;    NC_tsk_fab.resize (ntimes);
    Vector<FArrayBox> NC_LAT_fab;    NC_LAT_fab.resize (ntimes);
    Vector<FArrayBox> NC_LON_fab;    NC_LON_fab.resize (ntimes);

    // *** IArrayBox's at this level for holding mask data
    Vector<IArrayBox> NC_lmask_iab; NC_lmask_iab.resize(ntimes);

    // *** Variables at this level for holding metgrid file global attributes
    Vector<int> flag_psfc;           flag_psfc.resize(   ntimes);
    Vector<int> flag_msf;            flag_msf.resize(    ntimes);
    Vector<int> flag_sst;            flag_sst.resize(    ntimes);
    Vector<int> flag_tsk;            flag_tsk.resize(    ntimes);
    Vector<int> flag_lmask;          flag_lmask.resize(  ntimes);
    Vector<int> NC_nx;               NC_nx.resize(       ntimes);
    Vector<int> NC_ny;               NC_ny.resize(       ntimes);
    Vector<std::string> NC_dateTime; NC_dateTime.resize( ntimes);
    Vector<Real> NC_epochTime;       NC_epochTime.resize(ntimes);
    Vector<Real> NC_dx;              NC_dx.resize(       ntimes);
    Vector<Real> NC_dy;              NC_dy.resize(       ntimes);

    // Define the arena to be used for data allocation
    Arena* Arena_Used = The_Arena();
#ifdef AMREX_USE_GPU
    // Make sure this lives on CPU and GPU
    Arena_Used = The_Pinned_Arena();
#endif

#ifndef AMREX_USE_GPU
    if (metgrid_debug_quiescent)  Print() << "metgrid_debug_quiescent  = true" << std::endl;
    if (metgrid_debug_isothermal) Print() << "metgrid_debug_isothermal = true" << std::endl;
    if (metgrid_debug_dry)        Print() << "metgrid_debug_dry        = true" << std::endl;
    if (metgrid_debug_psfc)       Print() << "metgrid_debug_psfc       = true" << std::endl;
    if (metgrid_debug_msf)        Print() << "metgrid_debug_msf        = true" << std::endl;
    if (metgrid_interp_theta)     Print() << "metgrid_interp_theta     = true" << std::endl;
    if (metgrid_basic_linear)     Print() << "metgrid_basic_linear     = true" << std::endl;
#endif

    for (int itime(0); itime < ntimes; itime++) {
#ifndef AMREX_USE_GPU
        Print() << " init_from_metgrid: reading nc_init_file[" << lev << "][" << itime << "]\t" << nc_init_file[lev][itime] << std::endl;
#endif
        read_from_metgrid(lev, boxes_at_level[lev][0], nc_init_file[lev][itime],
                          NC_dateTime[itime],  NC_epochTime[itime],
                          flag_psfc[itime],    flag_msf[itime],
                          flag_sst[itime],     flag_tsk[itime],    flag_lmask[itime],
                          NC_nx[itime],        NC_ny[itime],       NC_dx[itime],      NC_dy[itime],
                          NC_xvel_fab[itime],  NC_yvel_fab[itime],
                          NC_temp_fab[itime],  NC_rhum_fab[itime], NC_pres_fab[itime],
                          NC_ght_fab[itime],   NC_hgt_fab[itime],  NC_psfc_fab[itime],
                          NC_MSFU_fab[itime],  NC_MSFV_fab[itime], NC_MSFM_fab[itime],
                          NC_sst_fab[itime],   NC_tsk_fab[itime],  NC_LAT_fab[itime],  NC_LON_fab[itime],
                          NC_lmask_iab[itime], geom[lev]);
    } // itime

    // Verify that files in nc_init_file[lev] are ordered from earliest to latest.
    for (int itime(1); itime < ntimes; itime++) AMREX_ALWAYS_ASSERT(NC_epochTime[itime] > NC_epochTime[itime-1]);

    // Start at the earliest time in nc_init_file[lev].
    start_bdy_time = NC_epochTime[0];

    //
    // Note that t_new and t_old carry *elapsed* time, not total time
    //
    Print() << "start_bdy_time is " << std::setprecision(timeprecision) << start_bdy_time
            << " from metgrid file but note that time variable in simulation is elapsed time" << std::endl;
    t_new[lev] = 0.;
    t_old[lev] = -1.e200;

    // Determine the spacing between met_em files.
    bdy_time_interval = NC_epochTime[1]-NC_epochTime[0];

    // Verify that met_em files have even spacing in time.
    for (int itime(1); itime < ntimes; itime++) {
        Real NC_dt = NC_epochTime[itime]-NC_epochTime[itime-1];
#ifndef AMREX_USE_GPU
        Print() << " " << nc_init_file[lev][itime-1] << " / " << nc_init_file[lev][itime] << " are " << NC_dt << " seconds apart" << std::endl;
#endif
        if (NC_dt != bdy_time_interval) Error("Time interval between consecutive met_em files must be consistent.");
    }

    auto& lev_new = vars_new[lev];

    std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];

    AMREX_ALWAYS_ASSERT(SolverChoice::terrain_type != TerrainType::None);

    z_phys->setVal(0.0);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // This defines only the z(i,j,0) values given the FAB filled from the NetCDF input
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        init_terrain_from_metgrid(z_phys_nd_fab, NC_hgt_fab);
    } // mf

    // This defines all the z(i,j,k) values given z(i,j,0) from above.
    make_terrain_fitted_coords(lev, geom[lev], *z_phys, zlevels_stag[lev], phys_bc_type);

    // Copy LATITUDE, LONGITUDE, SST and LANDMASK data into MF and iMF data structures
    auto& dm = lev_new[Vars::cons].DistributionMap();
    auto ngv = lev_new[Vars::cons].nGrowVect(); ngv[2] = 0;

    int i_lo = geom[lev].Domain().smallEnd(0); int i_hi = geom[lev].Domain().bigEnd(0);
    int j_lo = geom[lev].Domain().smallEnd(1); int j_hi = geom[lev].Domain().bigEnd(1);

    if (flag_sst[0]) {
        for (int itime(0); itime < ntimes; ++itime) {
            sst_lev[lev][itime] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(sst_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(sst_lev[lev][itime]))[mfi];
                FArrayBox& src = NC_sst_fab[itime];
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
        }
    } else {
        for (int itime(0); itime < ntimes; ++itime) sst_lev[lev][itime] = nullptr;
    }

    if (flag_tsk[0]) {
        for (int itime(0); itime < ntimes; ++itime) {
            tsk_lev[lev][itime] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(tsk_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                FArrayBox& dst = (*(tsk_lev[lev][itime]))[mfi];
                FArrayBox& src = NC_tsk_fab[itime];
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
        }
    } else {
        for (int itime(0); itime < ntimes; ++itime) tsk_lev[lev][itime] = nullptr;
    }

    if (flag_lmask[0]) {
        for (int itime(0); itime < ntimes; ++itime) {
            lmask_lev[lev][itime] = std::make_unique<iMultiFab>(ba2d[lev],dm,1,ngv);
            for ( MFIter mfi(*(lmask_lev[lev][itime]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
                Box gtbx = mfi.growntilebox();
                IArrayBox& dst = (*(lmask_lev[lev][itime]))[mfi];
                IArrayBox& src = NC_lmask_iab[itime];
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
    }

    solverChoice.has_lat_lon = true;
    lat_m[lev]    = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
    sinPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
    cosPhi_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
    for ( MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();
        FArrayBox& dst = (*(lat_m[lev]))[mfi];
        FArrayBox& src = NC_LAT_fab[0];
        const Array4<      Real>& sin_arr = (sinPhi_m[lev])->array(mfi);
        const Array4<      Real>& cos_arr = (cosPhi_m[lev])->array(mfi);
        const Array4<      Real>& dst_arr = dst.array();
        const Array4<const Real>& src_arr = src.const_array();
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = min(max(i, i_lo), i_hi);
            int lj = min(max(j, j_lo), j_hi);
            dst_arr(i,j,0) = src_arr(li,lj,0);

            Real lat_rad = dst_arr(i,j,0) * (PI/180.);
            sin_arr(i,j,0) = std::sin(lat_rad);
            cos_arr(i,j,0) = std::cos(lat_rad);
        });
    }

    lon_m[lev] = std::make_unique<MultiFab>(ba2d[lev],dm,1,ngv);
    for ( MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();
        FArrayBox& dst = (*(lon_m[lev]))[mfi];
        FArrayBox& src = NC_LON_fab[0];
        const Array4<      Real>& dst_arr = dst.array();
        const Array4<const Real>& src_arr = src.const_array();
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = min(max(i, i_lo), i_hi);
            int lj = min(max(j, j_lo), j_hi);
            dst_arr(i,j,0) = src_arr(li,lj,0);
        });
    }

    for (int itime(0); itime < ntimes; itime++) {
        // Verify that the grid size and resolution from met_em file matches that in geom (from ERF inputs file).
        AMREX_ALWAYS_ASSERT(geom[lev].CellSizeArray()[0] == NC_dx[itime]);
        AMREX_ALWAYS_ASSERT(geom[lev].CellSizeArray()[1] == NC_dy[itime]);
        // NC_nx-2 because NC_nx is the number of staggered grid points indexed from 1.
        AMREX_ALWAYS_ASSERT(geom[lev].Domain().bigEnd(0) == NC_nx[itime]-2);
        // NC_ny-2 because NC_ny is the number of staggered grid points indexed from 1.
        AMREX_ALWAYS_ASSERT(geom[lev].Domain().bigEnd(1) == NC_ny[itime]-2);
    } // itime

    // This makes the Jacobian.
    make_J(geom[lev],*z_phys,  *detJ_cc[lev]);
    make_areas(geom[lev],*z_phys,*ax[lev],*ay[lev],*az[lev]);

    // This defines z at w-cell faces.
    make_zcc(geom[lev],*z_phys,*z_phys_cc[lev]);

    // Set up FABs to hold data that will be used to set lateral boundary conditions.
    int MetGridBdyEnd = MetGridBdyVars::NumTypes-1;
    if (use_moisture) MetGridBdyEnd = MetGridBdyVars::NumTypes;

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
            bdy_data_xlo[itime][nvar].template setVal<RunOn::Device>(0.0);
            bdy_data_xhi[itime][nvar].template setVal<RunOn::Device>(0.0);
            bdy_data_ylo[itime][nvar].template setVal<RunOn::Device>(0.0);
            bdy_data_yhi[itime][nvar].template setVal<RunOn::Device>(0.0);
        }
    }

    // Set up FABs for mixing ratio and potential temperature on the origin model vertical grid.
    // This is necessary for input data that has relative humidity and temperature, but does not
    // include mixing ratio and potential temperature.
    // TODO: add alternate pathways for other origin models where different combinations of variables may be present.
    Box NC_box_unstag = NC_rhum_fab[0].box();
    FArrayBox mxrat_fab, theta_fab;
    mxrat_fab.resize(NC_box_unstag, 1, Arena_Used);
    theta_fab.resize(NC_box_unstag, 1, Arena_Used);

    // Set up FABs to hold interpolated origin variables that are only needed during initialization.
    // Pressure will be iteratively calculated from z, qv, and theta, but the origin model pressure
    // is necessary if metgrid_interp_theta=True and we interpolate temperature instead of theta.
    FArrayBox p_interp_fab, t_interp_fab;
    if (!metgrid_interp_theta) {
        Box ldomain = geom[lev].Domain();
        p_interp_fab.resize(ldomain, 1, Arena_Used);
        p_interp_fab.template setVal<RunOn::Device>(0.0);
        t_interp_fab.resize(ldomain, 1, Arena_Used);
        t_interp_fab.template setVal<RunOn::Device>(0.0);
    }

    const Real l_rdOcp = solverChoice.rdOcp;
    std::unique_ptr<iMultiFab> mask_c = OwnerMask(lev_new[Vars::cons], geom[lev].periodicity());//, lev_new[Vars::cons].nGrowVect());
    std::unique_ptr<iMultiFab> mask_u = OwnerMask(lev_new[Vars::xvel], geom[lev].periodicity());//, lev_new[Vars::xvel].nGrowVect());
    std::unique_ptr<iMultiFab> mask_v = OwnerMask(lev_new[Vars::yvel], geom[lev].periodicity());//, lev_new[Vars::yvel].nGrowVect());
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box tbxc = mfi.tilebox();
        Box tbxu = mfi.tilebox(IntVect(1,0,0));
        Box tbxv = mfi.tilebox(IntVect(0,1,0));

        // Define FABs for hlding some of the initial data
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];
        FArrayBox &z_phys_cc_fab = (*z_phys_cc[lev])[mfi];

        const Array4<const int>& mask_c_arr = mask_c->const_array(mfi);
        const Array4<const int>& mask_u_arr = mask_u->const_array(mfi);
        const Array4<const int>& mask_v_arr = mask_v->const_array(mfi);

        // Fill state data using origin data (initialization and BC arrays)
        //     x_vel   interpolated from origin levels
        //     y_vel   interpolated from origin levels
        //     z_vel   set to 0.0
        //     theta   calculate on origin levels then interpolate
        //     mxrat   convert RH -> Q on origin levels then interpolate
        init_state_from_metgrid(use_moisture, metgrid_interp_theta, metgrid_debug_quiescent,
                                metgrid_debug_isothermal, metgrid_debug_dry,
                                metgrid_basic_linear,
                                metgrid_use_below_sfc, metgrid_use_sfc,
                                metgrid_retain_sfc, metgrid_proximity,
                                metgrid_order, metgrid_force_sfc_k, l_rdOcp,
                                tbxc, tbxu, tbxv,
                                cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                z_phys_cc_fab,
                                NC_hgt_fab, NC_ght_fab, NC_xvel_fab,
                                NC_yvel_fab, NC_temp_fab, NC_rhum_fab,
                                NC_pres_fab, p_interp_fab, t_interp_fab,
                                theta_fab, mxrat_fab,
                                bdy_data_xlo, bdy_data_xhi,
                                bdy_data_ylo, bdy_data_yhi,
                                mask_c_arr, mask_u_arr, mask_v_arr);
    } // mf


#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // Use map scale factors directly from the met_em files
    for ( MFIter mfi(*mapfac[lev][MapFacType::u_x], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Define fabs for holding the initial data
        FArrayBox &msfu_fab = (*mapfac[lev][MapFacType::u_x])[mfi];
        FArrayBox &msfv_fab = (*mapfac[lev][MapFacType::v_x])[mfi];
        FArrayBox &msfm_fab = (*mapfac[lev][MapFacType::m_x])[mfi];

        init_msfs_from_metgrid(metgrid_debug_msf,
                               msfu_fab, msfv_fab, msfm_fab, flag_msf[0],
                               NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf


    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        FArrayBox&     p_hse_fab = p_hse[mfi];
        FArrayBox&    pi_hse_fab = pi_hse[mfi];
        FArrayBox&    th_hse_fab = th_hse[mfi];
        FArrayBox&    qv_hse_fab = qv_hse[mfi];
        FArrayBox&     r_hse_fab = r_hse[mfi];
        FArrayBox&      cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];


        // Fill base state data using origin data (initialization and BC arrays)
        //     p_hse     calculate moist hydrostatic pressure
        //     r_hse     calculate moist hydrostatic density
        //     pi_hse    calculate Exner term given pressure
        //     th_hse    calculate potential temperature
        //     qv_hse    calculate qv
        const Box valid_bx = mfi.validbox();
        init_base_state_from_metgrid(use_moisture, metgrid_debug_psfc, l_rdOcp,
                                     valid_bx, flag_psfc,
                                     cons_fab, r_hse_fab, p_hse_fab, pi_hse_fab, th_hse_fab,
                                     qv_hse_fab, z_phys_nd_fab, NC_psfc_fab);
    } // mf

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
        }
    }

    // NOTE: We must guarantee one halo cell in the bdy file.
    //       Otherwise, we make the total width match the set width.
    if (real_width-1 <= real_set_width) real_width = real_set_width;
#ifndef AMREX_USE_GPU
    Print() << "Running with specification width: " << real_set_width
            << " and relaxation width: " << real_width - real_set_width << std::endl;
#endif
}

/**
 * Helper function to initialize terrain nodal z coordinates given metgrid data.
 *
 * @param z_phys_nd_fab FArrayBox (Fab) holding the nodal z coordinates for terrain data we want to fill
 * @param NC_hgt_fab Vector of FArrayBox objects holding height data read from NetCDF files for metgrid data
 */
void
init_terrain_from_metgrid (FArrayBox& z_phys_nd_fab,
                           const Vector<FArrayBox>& NC_hgt_fab)
{
   int ntimes = 1; // Use terrain from the first met_em file.

   for (int itime(0); itime < ntimes; itime++) {
        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys_nd_fab.array();
        const Array4<Real const>& nc_hgt_arr = NC_hgt_fab[itime].const_array();

        const Box z_hgt_box = NC_hgt_fab[itime].box();

        int ilo = z_hgt_box.smallEnd()[0];
        int ihi = z_hgt_box.bigEnd()[0];
        int jlo = z_hgt_box.smallEnd()[1];
        int jhi = z_hgt_box.bigEnd()[1];

        Box z_phys_box = z_phys_nd_fab.box();
        Box from_box = surroundingNodes(NC_hgt_fab[itime].box());
        from_box.growHi(2,-1);

        Box bx = z_phys_box & from_box;
        Box bxu = bx; bxu.growHi(0,1);
        Box bxv = bx; bxv.growHi(1,1);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            z_arr(i,j,k) =  0.25 * ( nc_hgt_arr (ii,jj  ,k) + nc_hgt_arr(ii-1,jj  ,k) +
                                     nc_hgt_arr (ii,jj-1,k) + nc_hgt_arr(ii-1,jj-1,k) );
        });
    } // itime
}

/**
 * Helper function to initialize state and velocity data read from metgrid data.
 *
 * @param use_moisture bool True if solverChoice.moisture_type != MoistureType::None
 * @param metgrid_interp_theta bool calculate theta on origin levels, then interpolate
 * @param metgrid_debug_quiescent bool overwrite u and v with 0.0
 * @param metgrid_debug_isothermal bool overwrite theta with 300.0
 * @param metgrid_debug_dry bool overwrite qv with 0.0
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
 * @param NC_hgt_fab  Vector of FArrayBox objects holding metgrid data for terrain height
 * @param NC_ght_fab  Vector of FArrayBox objects holding metgrid data for height of cell centers
 * @param NC_xvel_fab Vector of FArrayBox objects holding metgrid data for x-velocity
 * @param NC_yvel_fab Vector of FArrayBox objects holding metgrid data for y-velocity
 * @param NC_zvel_fab Vector of FArrayBox objects holding metgrid data for z-velocity
 * @param NC_temp_fab Vector of FArrayBox objects holding metgrid data for temperature
 * @param NC_rhum_fab Vector of FArrayBox objects holding metgrid data for relative humidity
 * @param NC_pres_fab Vector of FArrayBox objects holding metgrid data for pressure
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
init_state_from_metgrid (const bool use_moisture,
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
                         const Vector<FArrayBox>& NC_hgt_fab,
                         const Vector<FArrayBox>& NC_ght_fab,
                         const Vector<FArrayBox>& NC_xvel_fab,
                         const Vector<FArrayBox>& NC_yvel_fab,
                         const Vector<FArrayBox>& NC_temp_fab,
                         const Vector<FArrayBox>& NC_rhum_fab,
                         const Vector<FArrayBox>& NC_pres_fab,
                         FArrayBox& p_interp_fab,
                         FArrayBox& t_interp_fab,
                         FArrayBox& theta_fab,
                         FArrayBox& mxrat_fab,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_xlo,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_xhi,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_ylo,
                         Vector<Vector<FArrayBox>>& fabs_for_bcs_yhi,
                         const Array4<const int>& mask_c_arr,
                         const Array4<const int>& mask_u_arr,
                         const Array4<const int>& mask_v_arr)
{
    bool metgrid_exp_interp = false; // interpolate w.r.t. exp(z) for non-pressure variables.

    // Loop over each time in the origin data.
    int ntimes = NC_hgt_fab.size();
    for (int itime(0); itime < ntimes; itime++)
    {

        // ********************************************************
        // U
        // ********************************************************
        {
#ifndef AMREX_USE_GPU
        Print() << "[init_state_from_metgrid] vertical interpolation of u-velocity, itime = " << itime << std::endl;
#endif
        Box bx2d = NC_xvel_fab[itime].box() & tbxu;
        bx2d.setRange(2,0);
        auto const orig_data = NC_xvel_fab[itime].const_array();
        auto const orig_z    = NC_ght_fab[itime].const_array();
        auto       new_data  = x_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

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
            if (metgrid_debug_quiescent) { // Debugging option to run quiescent.
                for (int k(0); k<=kmax; k++) {
                  if (mask_u_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = 0.0;
                  if (mask_u_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = 0.0;
                  if (mask_u_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = 0.0;
                  if (mask_u_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = 0.0;
                    if (itime == 0) new_data(i,j,k,0) = 0.0;
                }
            } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                for (int k(0); k<=kmax; k++) {
                    Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'X',0,orig_z,orig_data,new_z);
                    if (mask_u_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = Interp_Val;
                    if (mask_u_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = Interp_Val;
                    if (mask_u_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = Interp_Val;
                    if (mask_u_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = Interp_Val;
                    if (itime == 0) new_data(i,j,k,0) = Interp_Val;
                }
            } else { // Vertical interpolation and quality control similar to that from WRF.
                interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                           metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                           metgrid_force_sfc_k, i, j, 0, itime, 'U', 'X',
                                           orig_z, orig_data, new_z, new_data, true,
                                           bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                           bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_u_arr);
            }
        });
        }


        // ********************************************************
        // V
        // ********************************************************
        {
#ifndef AMREX_USE_GPU
        Print() << "[init_state_from_metgrid] vertical interpolation of v-velocity, itime = " << itime << std::endl;
#endif
        Box bx2d = NC_yvel_fab[itime].box() & tbxv;
        bx2d.setRange(2,0);
        auto const orig_data = NC_yvel_fab[itime].const_array();
        auto const orig_z    = NC_ght_fab[itime].const_array();
        auto       new_data  = y_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

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
            if (metgrid_debug_quiescent) { // Debugging option to run quiescent.
                for (int k(0); k<=kmax; k++) {
                    if (mask_v_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = 0.0;
                    if (mask_v_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = 0.0;
                    if (mask_v_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = 0.0;
                    if (mask_v_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = 0.0;
                    if (itime == 0) new_data(i,j,k,0) = 0.0;
                }
            } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                for (int k(0); k<=kmax; k++) {
                    Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'Y',0,orig_z,orig_data,new_z);
                    if (mask_v_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = Interp_Val;
                    if (mask_v_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = Interp_Val;
                    if (mask_v_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = Interp_Val;
                    if (mask_v_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = Interp_Val;
                    if (itime == 0) new_data(i,j,k,0) = Interp_Val;
                }
            } else {
                interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                           metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                           metgrid_force_sfc_k, i, j, 0, itime, 'V', 'Y',
                                           orig_z, orig_data, new_z, new_data, true,
                                           bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                           bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_v_arr);
            }
        });
        }


        // ********************************************************
        // W
        // ********************************************************
        if (itime == 0) { // update at initialization
            z_vel_fab.template setVal<RunOn::Device>(0.0);
        }


        // ********************************************************
        // Initialize all state_fab variables to zero
        // ********************************************************
        if (itime == 0) { // update at initialization
            state_fab.template setVal<RunOn::Device>(0.0);
        }


        // ********************************************************
        // theta
        // ********************************************************
        if (metgrid_interp_theta) {
            // Calculate potential temperature on the origin model vertical levels
            // then interpolate that onto the ERF vertical levels.

            { // calculate potential temperature.
                Box bx = NC_rhum_fab[itime].box() & tbxc;
                auto const temp  = NC_temp_fab[itime].const_array();
                auto const pres  = NC_pres_fab[itime].const_array();
                auto       theta = theta_fab.array();

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    theta(i,j,k) = getThgivenTandP(temp(i,j,k),pres(i,j,k),l_rdOcp);
                });
            }

            { // vertical interpolation of potential temperature.
#ifndef AMREX_USE_GPU
            Print() << "[init_state_from_metgrid] vertical interpolation of potential temperature, itime = " << itime << std::endl;
#endif
            Box bx2d = NC_temp_fab[itime].box() & tbxc;
            bx2d.setRange(2,0);
            auto const orig_data = theta_fab.const_array();
            auto const orig_z    = NC_ght_fab[itime].const_array();
            auto       new_data  = state_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

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
                if (metgrid_debug_isothermal) { // Debugging option to run isothermal.
                    for (int k(0); k<=kmax; k++) {
                        if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = 300.0;
                        if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = 300.0;
                        if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = 300.0;
                        if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = 300.0;
                        if (itime == 0) new_data(i,j,k,RhoTheta_comp) = 300.0;
                    }
                } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                    for (int k(0); k<=kmax; k++) {
                        Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',0,orig_z,orig_data,new_z);
                        if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = Interp_Val;
                        if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = Interp_Val;
                        if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = Interp_Val;
                        if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = Interp_Val;
                        if (itime == 0) new_data(i,j,k,RhoTheta_comp) = Interp_Val;
                    }
                } else { // Vertical interpolation and quality control similar to that from WRF.
                    interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                               metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                               metgrid_force_sfc_k, i, j, RhoTheta_comp, itime, 'T', 'M',
                                               orig_z, orig_data, new_z, new_data, true,
                                               bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                               bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_c_arr);
                }
            });
            }

        } else { // metgrid_interp_theta == false

            { // vertical interpolation of pressure.
#ifndef AMREX_USE_GPU
            Print() << "[init_state_from_metgrid] vertical interpolation of pressure, itime = " << itime << std::endl;
#endif
            Box bx2d = p_interp_fab.box() & tbxc;
            bx2d.setRange(2,0);
            auto const orig_data = NC_pres_fab[itime].const_array();
            auto const orig_z    = NC_ght_fab[itime].const_array();
            auto       new_data  = p_interp_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            const amrex::Array4<amrex::Real> bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi;

            int kmax = ubound(tbxc).z;

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                    for (int k(0); k<=kmax; k++) {
                        Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',0,orig_z,orig_data,new_z);
                        new_data(i,j,k) = Interp_Val;
                    }
                } else { // Vertical interpolation and quality control similar to that from WRF.
                    // Interpolate pressure not w.r.t. z but rather p_0*exp(-CONST_GRAV*z/(t_0*R_d)).
                    // This is akin to interpolating in pressure-space assuming a baroclinic atmosphere.
                    interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, true,
                                               metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                               metgrid_force_sfc_k, i, j, 0, 0, 'T', 'M',
                                               orig_z, orig_data, new_z, new_data, false,
                                               bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                               bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_c_arr);
                }
            });
            }

            { // vertical interpolation of temperature.
#ifndef AMREX_USE_GPU
            Print() << "[init_state_from_metgrid] vertical interpolation of temperature, itime = " << itime << std::endl;
#endif
            Box bx2d = p_interp_fab.box() & tbxc;
            bx2d.setRange(2,0);
            auto const orig_data = NC_temp_fab[itime].const_array();
            auto const orig_z    = NC_ght_fab[itime].const_array();
            auto       new_data  = t_interp_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

            Box bx_xlo, bx_xhi, bx_ylo, bx_yhi;
            const amrex::Array4<amrex::Real> bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi;

            int kmax = ubound(tbxc).z;

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
            {
                if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                    for (int k(0); k<=kmax; k++) {
                        Real Interp_Val = interpolate_column_metgrid_linear(i,j,k,'M',0,orig_z,orig_data,new_z);
                        new_data(i,j,k) = Interp_Val;
                    }
                } else { // Vertical interpolation and quality control similar to that from WRF.
                    // According to WRF's code comments, "It is better to
                    // interpolate temperature and potential temperature
                    // in LOG(p), regardless of requested default."
                    interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, false,
                                               metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                               metgrid_force_sfc_k, i, j, 0, 0, 'T', 'M',
                                               orig_z, orig_data, new_z, new_data, false,
                                               bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                               bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_c_arr);
                }
            });
            }

            { // calculate potential temperature on the ERF vertical levels.
            auto const temp  = t_interp_fab.const_array();
            auto const pres  = p_interp_fab.const_array();
            auto       new_data = state_fab.array();

            Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::T].box();
            Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::T].box();
            Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::T].box();
            Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::T].box();
            auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::T].array();
            auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::T].array();
            auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::T].array();
            auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::T].array();

            ParallelFor(tbxc, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real Calc_Val = getThgivenTandP(temp(i,j,k),pres(i,j,k),l_rdOcp);
                if (metgrid_debug_isothermal) Calc_Val = 300.0; // Debugging option to run isothermal.
                if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = Calc_Val;
                if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = Calc_Val;
                if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = Calc_Val;
                if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = Calc_Val;
                if (itime == 0) new_data(i,j,k,RhoTheta_comp) = Calc_Val;
            });
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
                Box bx = NC_temp_fab[itime].box() & tbxc;
                auto const rhum  = NC_rhum_fab[itime].const_array();
                auto const temp  = NC_temp_fab[itime].const_array();
                auto const pres  = NC_pres_fab[itime].const_array();
                auto       mxrat = mxrat_fab.array();

                ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    rh_to_mxrat(i,j,k,rhum,temp,pres,mxrat);
                });
            }

            { // vertical interpolation of vapor mixing ratio.
#ifndef AMREX_USE_GPU
                Print() << "[init_state_from_metgrid] vertical interpolation of vapor mixing ratio, itime = " << itime << std::endl;
#endif
                Box bx2d = NC_temp_fab[itime].box() & tbxc;
                bx2d.setRange(2,0);
                auto const orig_data = mxrat_fab.const_array();
                auto const orig_z    = NC_ght_fab[itime].const_array();
                auto       new_data  = state_fab.array();
                auto const new_z     = z_phys_nd_fab.const_array();

                Box bx_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::QV].box();
                Box bx_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::QV].box();
                Box bx_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::QV].box();
                Box bx_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::QV].box();
                auto bc_data_xlo = fabs_for_bcs_xlo[itime][MetGridBdyVars::QV].array();
                auto bc_data_xhi = fabs_for_bcs_xhi[itime][MetGridBdyVars::QV].array();
                auto bc_data_ylo = fabs_for_bcs_ylo[itime][MetGridBdyVars::QV].array();
                auto bc_data_yhi = fabs_for_bcs_yhi[itime][MetGridBdyVars::QV].array();

                int kmax = ubound(tbxc).z;

                int state_indx = RhoQ1_comp;
                ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
                {
                    if (metgrid_debug_dry) { // Debugging option to run dry.
                        for (int k(0); k<=kmax; k++) {
                            if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = 0.0;
                            if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = 0.0;
                            if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = 0.0;
                            if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = 0.0;
                            if (itime == 0) new_data(i,j,k,state_indx)   = 0.0;
                        }
                    } else if (metgrid_basic_linear) { // Linear interpolation with no quality control.
                        for (int k(0); k<=kmax; k++) {
                            Real Interp_Val  = interpolate_column_metgrid_linear(i,j,k,'M',0,orig_z,orig_data,new_z);
                            if (mask_c_arr(i,j,k) && bx_xlo.contains(i,j,k)) bc_data_xlo(i,j,k,0) = Interp_Val;
                            if (mask_c_arr(i,j,k) && bx_xhi.contains(i,j,k)) bc_data_xhi(i,j,k,0) = Interp_Val;
                            if (mask_c_arr(i,j,k) && bx_ylo.contains(i,j,k)) bc_data_ylo(i,j,k,0) = Interp_Val;
                            if (mask_c_arr(i,j,k) && bx_yhi.contains(i,j,k)) bc_data_yhi(i,j,k,0) = Interp_Val;
                            if (itime == 0) new_data(i,j,k,state_indx)   = Interp_Val;
                        }
                    } else { // Vertical interpolation and quality control similar to that from WRF.
                        interpolate_column_metgrid(metgrid_use_below_sfc, metgrid_use_sfc, metgrid_exp_interp,
                                                   metgrid_retain_sfc, metgrid_proximity, metgrid_order,
                                                   metgrid_force_sfc_k, i, j, state_indx, itime, 'Q', 'M',
                                                   orig_z, orig_data, new_z, new_data, true,
                                                   bc_data_xlo, bc_data_xhi, bc_data_ylo, bc_data_yhi,
                                                   bx_xlo, bx_xhi, bx_ylo, bx_yhi, mask_c_arr);
                    }
                });
            }
        } // use_moisture

    } // itime
}


/**
 * Helper function for initializing hydrostatic base state data from metgrid data
 *
 * @param use_moisture bool True if solverChoice.moisture_type != MoistureType::None
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param valid_bx Box specifying the index space we are to initialize
 * @param flag_psfc Vector of Integer 1 if surface pressure is in metgrid data, 0 otherwise
 * @param state_fab FArrayBox holding the state data to initialize
 * @param r_hse_fab FArrayBox holding the hydrostatic base state density we are initializing
 * @param p_hse_fab FArrayBox holding the hydrostatic base state pressure we are initializing
 * @param pi_hse_fab FArrayBox holding the hydrostatic base Exner pressure we are initializing
 * @param th_hse_fab FArrayBox holding the base state potential temperature we are initializing
 * @param qv_hse_fab FArrayBox holding the base state qv we are initializing
 * @param z_phys_nd_fab FArrayBox holding nodal z coordinate data for terrain
 * @param NC_psfc_fab Vector of FArrayBox objects holding metgrid data for surface pressure
 * @param fabs_for_bcs Vector of Vector of FArrayBox objects holding MetGridBdyVars at each met_em time.
 * @param mask_c_arr
 */
void
init_base_state_from_metgrid (const bool use_moisture,
                              const bool metgrid_debug_psfc,
                              const Real l_rdOcp,
                              const Box& valid_bx,
                              const Vector<int>& flag_psfc,
                              FArrayBox& state_fab,
                              FArrayBox& r_hse_fab,
                              FArrayBox& p_hse_fab,
                              FArrayBox& pi_hse_fab,
                              FArrayBox& th_hse_fab,
                              FArrayBox& qv_hse_fab,
                              FArrayBox& z_phys_cc_fab,
                              const Vector<FArrayBox>& NC_psfc_fab)
{
    int RhoQ_comp = RhoQ1_comp;
    int kmax = ubound(valid_bx).z;

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

    // Device vectors for psfc flags
    Gpu::DeviceVector<int>flag_psfc_d(flag_psfc.size());
    Gpu::copy(Gpu::hostToDevice, flag_psfc.begin(), flag_psfc.end(), flag_psfc_d.begin());
    int* flag_psfc_vec = flag_psfc_d.data();

    // Expose for copy to GPU
    Real grav = CONST_GRAV;

    { // set pressure and density at initialization.
        const Array4<Real>& r_hse_arr  = r_hse_fab.array();
        const Array4<Real>& p_hse_arr  = p_hse_fab.array();
        const Array4<Real>& pi_hse_arr = pi_hse_fab.array();
        const Array4<Real>& th_hse_arr = th_hse_fab.array();
        const Array4<Real>& qv_hse_arr = qv_hse_fab.array();
        auto psfc_flag = flag_psfc_vec[0];

        // ********************************************************
        // calculate density and pressure for initial conditions.
        // ********************************************************
        Box valid_bx2d = valid_bx;
        valid_bx2d.setRange(2,0);
        auto const orig_psfc = NC_psfc_fab[0].const_array();
        auto       new_data  = state_fab.array();
        auto const new_z     = z_phys_cc_fab.const_array();

        ParallelFor(valid_bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            const int maxiter = 10;
            const amrex::Real tol = 1.0e-10;

            // Low and Hi column variables
            Real psurf;
            Real z_lo,   z_hi;
            Real p_lo,   p_hi;
            Real qv_lo, qv_hi;
            Real rd_lo, rd_hi;
            Real thetad_lo, thetad_hi;

            // Calculate or use pressure at the surface.
            if (metgrid_debug_psfc) {
                psurf = std::pow(10, 5);
            } else if (psfc_flag == 1) {
                psurf = orig_psfc(i,j,0);
            } else {
                z_lo     = new_z(i,j,0);
                Real t_0 = 290.0; // WRF's model_config_rec%base_temp
                Real a   = 50.0;  // WRF's model_config_rec%base_lapse
                psurf = p_0*exp(-t_0/a+std::pow((std::pow(t_0/a, 2.)-2.0*grav*z_lo/(a*R_d)), 0.5));
            }

            // Iterations for the first CC point that is 1/2 dz off the surface
            {
                z_lo      = new_z(i,j,0);
                qv_lo     = (use_moisture) ? new_data(i,j,0,RhoQ_comp) : 0.0;
                rd_lo     = 0.0; // initial guess
                thetad_lo = new_data(i,j,0,RhoTheta_comp);
                Real half_dz = z_lo;
                Real qvf     = 1.0+(R_v/R_d)*qv_lo;
                Real thetam  = thetad_lo*qvf;
                for (int it(0); it<maxiter; it++) {
                    p_lo = psurf-half_dz*rd_lo*(1.0+qv_lo)*grav;
                    if (p_lo < 0.0) p_lo = 0.0;
                    rd_lo = (p_0/(R_d*thetam))*std::pow(p_lo/p_0, iGamma);
                } // it
                 p_hse_arr(i,j,0) =  p_lo;
                 r_hse_arr(i,j,0) = rd_lo;
                qv_hse_arr(i,j,0) = qv_lo;
            }

            // Iterations for k \in [1 kmax]
            for (int k(1); k<=kmax; k++) {
                // Known hi data
                z_hi  = new_z(i,j,k);
                qv_hi = (use_moisture) ? new_data(i,j,k,RhoQ_comp) : 0.0;
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
                Real rho_tot_lo = rd_lo * (1. + qv_lo);
                Real C = -p_lo + 0.5*rho_tot_lo*grav*dz;

                // Initial residual
                Real rho_tot_hi = rd_hi * (1. + qv_hi);
                Real F = p_hi + 0.5*rho_tot_hi*grav*dz + C;

                // Do iterations
                if (std::abs(F)>tol) HSEutils::Newton_Raphson_hse(tol, R_d/Cp_d, dz,
                                                                  grav, C, thetad_hi,
                                                                  qv_hi, qv_hi, p_hi,
                                                                  rd_hi, F);

                // Copy solution to base state
                p_hse_arr(i,j,k) =  p_hi;
                r_hse_arr(i,j,k) = rd_hi;

                // Copy hi to lo
                z_lo  = z_hi;
                p_lo  = p_hi;
                qv_lo = qv_hi;
                rd_lo = rd_hi;
                thetad_lo = thetad_hi;
            }
        });

        ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            // Multiply by Rho to get conserved vars
            Real Qv = 0.0;
            new_data(i,j,k,Rho_comp)       = r_hse_arr(i,j,k);
            new_data(i,j,k,RhoTheta_comp) *= r_hse_arr(i,j,k);
            if (use_moisture) {
                Qv = new_data(i,j,k,RhoQ_comp);
                new_data(i,j,k,RhoQ_comp) *= r_hse_arr(i,j,k);
            }
            for (int n(0); n < NSCALARS; n++) {
                new_data(i,j,k,RhoScalar_comp+n) = 0.0;
            }

            // r_hse needs to include the moisture (account for that here)
            r_hse_arr(i,j,k) *= (1.0 + Qv);

            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
            th_hse_arr(i,j,k) = getRhoThetagivenP(p_hse_arr(i,j,k), Qv) / new_data(i,j,k,Rho_comp);
            qv_hse_arr(i,j,k) = Qv;
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
            qv_hse_arr(i,j,k) = qv_hse_arr(i,j,k+1);
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
 * @param NC_MSFU_fab Vector of FArrayBox objects holding metgrid data for x-velocity map factors
 * @param NC_MSFV_fab Vector of FArrayBox objects holding metgrid data for y-velocity map factors
 * @param NC_MSFM_fab Vector of FArrayBox objects holding metgrid data for z-velocity map factors
 */
void
init_msfs_from_metgrid (const bool metgrid_debug_msf,
                        FArrayBox& msfu_fab,
                        FArrayBox& msfv_fab,
                        FArrayBox& msfm_fab,
                        const int& flag_msf,
                        const Vector<FArrayBox>& NC_MSFU_fab,
                        const Vector<FArrayBox>& NC_MSFV_fab,
                        const Vector<FArrayBox>& NC_MSFM_fab)
{
//    int ntimes = NC_MSFU_fab.size();
    int ntimes = 1;
    for (int itime(0); itime < ntimes; itime++) {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //

        // This copies or sets mapfac
        if ((flag_msf == 1) and (!metgrid_debug_msf)) {
            msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab[itime]);
            msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab[itime]);
            msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab[itime]);
        } else {
#ifndef AMREX_USE_GPU
            Print() << " map factors are not present in met_em files. Setting to 1.0" << std::endl;
#endif
            msfm_fab.template setVal<RunOn::Device>(1.0);
            msfu_fab.template setVal<RunOn::Device>(1.0);
            msfv_fab.template setVal<RunOn::Device>(1.0);
        }
    } // itime
}
#endif // ERF_USE_NETCDF
