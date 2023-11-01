/**
 * \file ERF_init_from_metgrid.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>

using namespace amrex;

void
read_from_metgrid (int lev, const Box& domain, const std::string& fname,
                   std::string& NC_dateTime, Real& NC_epochTime,
                   int& flag_psfc, int& flag_msfu, int& flag_msfv, int& flag_msfm,
                   int& flag_hgt,  int& flag_sst,  int& flag_lmask,
                   int& NC_nx,     int& NC_ny,
                   Real& NC_dx,    Real& NC_dy,
                   FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                   FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                   FArrayBox& NC_pres_fab, FArrayBox& NC_ght_fab,
                   FArrayBox& NC_hgt_fab,  FArrayBox& NC_psfc_fab,
                   FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                   FArrayBox& NC_msfm_fab, FArrayBox& NC_sst_fab,
                   IArrayBox& NC_lmask_iab);

void
interpolate_column (int i, int j, int src_comp, int dest_comp,
                    const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                    const Array4<Real const>&  new_z, const Array4<Real>&  new_data);

void
interpolate_column_metgrid (int i, int j, char stag, int src_comp, int dest_comp,
                            const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                            const Array4<Real const>&  new_z, const Array4<Real>&  new_data);

void
init_terrain_from_metgrid (FArrayBox& z_phys_nd_fab,
                           const Vector<FArrayBox>& NC_hgt_fab);

void
init_state_from_metgrid (const Real l_rdOcp,
                         FArrayBox& state_fab,
                         FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                         FArrayBox& z_vel_fab, FArrayBox& z_phys_nd_fab,
                         const Vector<FArrayBox>& NC_hgt_fab,
                         const Vector<FArrayBox>& NC_ght_fab,
                         const Vector<FArrayBox>& NC_xvel_fab,
                         const Vector<FArrayBox>& NC_yvel_fab,
                         const Vector<FArrayBox>& NC_zvel_fab,
                         const Vector<FArrayBox>& NC_temp_fab,
                         const Vector<FArrayBox>& NC_rhum_fab,
                         Vector<FArrayBox>& theta_fab,
                         Vector<FArrayBox>& mxrat_fab,
                         amrex::Vector<amrex::Vector<FArrayBox>>& fabs_for_bcs);

void
init_msfs_from_metgrid (FArrayBox& msfu_fab,
                        FArrayBox& msfv_fab,
                        FArrayBox& msfm_fab,
                        const int& flag_msfu,
                        const int& flag_msfv,
                        const int& flag_msfm,
                        const Vector<FArrayBox>& NC_MSFU_fab,
                        const Vector<FArrayBox>& NC_MSFV_fab,
                        const Vector<FArrayBox>& NC_MSFM_fab);

void
init_base_state_from_metgrid (const Real l_rdOcp,
                              const Box& valid_bx,
                              const Vector<int>& flag_psfc,
                              FArrayBox& state,
                              FArrayBox& r_hse_fab,
                              FArrayBox& p_hse_fab,
                              FArrayBox& pi_hse_fab,
                              FArrayBox& z_phys_nd_fab,
                              const Vector<FArrayBox>& NC_ght_fab,
                              const Vector<FArrayBox>& NC_psfc_fab,
                              amrex::Vector<amrex::Vector<FArrayBox>>& fabs_for_bcs);

void
rh_to_mxrat (int i, int j, int k,
             const Array4<Real const>& rhum,
             const Array4<Real const>& temp,
             const Array4<Real const>& pres,
             const Array4<Real>& mxrat);

void
calc_rho_p (const int& kmax,
            const int& flag_psfc,
            const Real& psfc,
            const Vector<Real>& Theta_vec,
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
            const Vector<Real>& Q_vec,
#endif
            const Vector<Real>& z_vec,
            Vector<Real>& Rhod_vec,
            Vector<Real>& Pd_vec);

#ifdef ERF_USE_NETCDF
/**
 * Initializes ERF data using metgrid data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_metgrid (int lev)
{
#ifndef AMREX_USE_GPU
#if defined(ERF_USE_MOISTURE)
    amrex::Print() << "Init with met_em with ERF_USE_MOISTURE" << std::endl;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    amrex::Print() << "Init with met_em with ERF_USE_WARM_NO_PRECIP" << std::endl;
#else
    amrex::Print() << "Init with met_em without moisture" << std::endl;
#endif
#endif

    int nboxes = num_boxes_at_level[lev];
    int ntimes = num_files_at_level[lev];

    if (nc_init_file.empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    if (nc_init_file[lev].empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    // At least two met_em files are necessary to calculate tendency terms.
    AMREX_ALWAYS_ASSERT(ntimes >= 2);

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

    // *** IArrayBox's at this level for holding mask data
    Vector<IArrayBox> NC_lmask_iab; NC_lmask_iab.resize(ntimes);

    // *** Variables at this level for holding metgrid file global attributes
    Vector<int> flag_psfc;           flag_psfc.resize(   ntimes);
    Vector<int> flag_msfu;           flag_msfu.resize(   ntimes);
    Vector<int> flag_msfv;           flag_msfv.resize(   ntimes);
    Vector<int> flag_msfm;           flag_msfm.resize(   ntimes);
    Vector<int> flag_hgt;            flag_hgt.resize(    ntimes);
    Vector<int> flag_sst;            flag_sst.resize(    ntimes);
    Vector<int> flag_lmask;          flag_lmask.resize(  ntimes);
    Vector<int> NC_nx;               NC_nx.resize(       ntimes);
    Vector<int> NC_ny;               NC_ny.resize(       ntimes);
    Vector<std::string> NC_dateTime; NC_dateTime.resize( ntimes);
    Vector<Real> NC_epochTime;       NC_epochTime.resize(ntimes);
    Vector<Real> NC_dx;              NC_dx.resize(       ntimes);
    Vector<Real> NC_dy;              NC_dy.resize(       ntimes);

    for (int it = 0; it < ntimes; it++) {
#ifndef AMREX_USE_GPU
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << ": nc_init_file[" << lev << "][" << it << "]\t" << nc_init_file[lev][it] << std::endl;
#endif
        read_from_metgrid(lev, boxes_at_level[lev][0], nc_init_file[lev][it],
                          NC_dateTime[it], NC_epochTime[it],
                          flag_psfc[it],   flag_msfu[it],   flag_msfv[it],  flag_msfm[it],
                          flag_hgt[it],    flag_sst[it],    flag_lmask[it],
                          NC_nx[it],       NC_ny[it],       NC_dx[it],      NC_dy[it],
                          NC_xvel_fab[it], NC_yvel_fab[it],
                          NC_temp_fab[it], NC_rhum_fab[it], NC_pres_fab[it],
                          NC_ght_fab[it],  NC_hgt_fab[it],  NC_psfc_fab[it],
                          NC_MSFU_fab[it], NC_MSFV_fab[it], NC_MSFM_fab[it],
                          NC_sst_fab[it],  NC_lmask_iab[it]);
#ifndef AMREX_USE_GPU
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": flag_psfc   \t" << flag_psfc[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": flag_msfu   \t" << flag_msfu[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": flag_msfv   \t" << flag_msfv[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": flag_msfm   \t" << flag_msfm[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": flag_hgt    \t" << flag_hgt[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": NC_nx       \t" << NC_nx[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": NC_ny       \t" << NC_ny[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": NC_dx       \t" << NC_dx[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": NC_dy       \t" << NC_dy[it] << std::endl;
        amrex::AllPrint() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << " it-" << it << ": NC_epochTime\t" << NC_epochTime[it] << std::endl;
        // NC_dateTime is only on the IOProcessor.
        amrex::Print() << " DJW init_from_metgrid proc-" << ParallelDescriptor::MyProc() << ": NC_dateTime \t" << NC_dateTime[it] << std::endl;
#endif
    } // it

    // Verify that files in nc_init_file[lev] are ordered from earliest to latest.
    for (int it = 1; it < ntimes; it++) AMREX_ALWAYS_ASSERT(NC_epochTime[it] > NC_epochTime[it-1]);

    // Start at the earliest time in nc_init_file[lev].
    start_bdy_time = NC_epochTime[0];

    // Determine the spacing between met_em files.
    bdy_time_interval = NC_epochTime[1]-NC_epochTime[0];

    // Verify that met_em files have even spacing in time.
    for (int it = 1; it < ntimes; it++) {
        Real NC_dt = NC_epochTime[it]-NC_epochTime[it-1];
#ifndef AMREX_USE_GPU
        amrex::Print() << " " << nc_init_file[lev][it-1] << " / " << nc_init_file[lev][it] << " are " << NC_dt << " seconds apart" << std::endl;
#endif
        if (NC_dt != bdy_time_interval) amrex::Error("Time interval between consecutive met_em files must be consistent.");
    }

    // Set up a FAB for mixing ratio and another for potential temperature.
    // Necessary because the input data has relative humidity and temperature, not mixing ratio and potential temperature.
    // TODO: add alternate pathways for other origin models where different combinations of variables may be present.
    Vector<FArrayBox> mxrat_fab; mxrat_fab.resize(ntimes);
    Vector<FArrayBox> theta_fab; theta_fab.resize(ntimes);
    for (int it = 0; it < ntimes; it++) {
        Box NC_box_unstag = NC_rhum_fab[it].box();
#ifdef AMREX_USE_GPU
        mxrat_fab[it].resize(NC_box_unstag, 1, The_Pinned_Arena());
        theta_fab[it].resize(NC_box_unstag, 1, The_Pinned_Arena());
#else
        mxrat_fab[it].resize(NC_box_unstag, 1);
        theta_fab[it].resize(NC_box_unstag, 1);
#endif
    }

    auto& lev_new = vars_new[lev];

    std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];

    AMREX_ALWAYS_ASSERT(solverChoice.use_terrain);

    // Verify that the terrain height (HGT_M) was in each met_em file.
    for (int it = 0; it < ntimes; it++) AMREX_ALWAYS_ASSERT(flag_hgt[it] == 1);

    z_phys->setVal(0.0);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // This defines only the z(i,j,0) values given the FAB filled from the NetCDF input
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        init_terrain_from_metgrid(z_phys_nd_fab, NC_hgt_fab);
    } // mf

    // This defines all the z(i,j,k) values given z(i,j,0) from above.
    init_terrain_grid(geom[lev], *z_phys);

    // Copy SST and LANDMASK data into MF and iMF data structures
    auto& ba = lev_new[Vars::cons].boxArray();
    auto& dm = lev_new[Vars::cons].DistributionMap();
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));
    if (flag_sst[0]) {
        sst_lev[lev] = std::make_unique<MultiFab>(ba2d,dm,1,0);
        for ( MFIter mfi(*(sst_lev[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Box tbx = mfi.tilebox();
            FArrayBox& dst = (*(sst_lev[lev]))[mfi];
            FArrayBox& src = NC_sst_fab[0];
            dst.copy(src, tbx, 0, tbx, 0, 1);
        }
    } else {
        sst_lev[lev] = nullptr;
    }
    if (flag_lmask[0]) {
        lmask_lev[lev] = std::make_unique<iMultiFab>(ba2d,dm,1,0);
        for ( MFIter mfi(*(lmask_lev[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
            Box tbx = mfi.tilebox();
            IArrayBox& dst = (*(lmask_lev[lev]))[mfi];
            IArrayBox& src = NC_lmask_iab[0];
            dst.copy(src, tbx, 0, tbx, 0, 1);
        }
    } else {
        lmask_lev[lev] = nullptr;
    }
    amrex::AllPrint() << "Copied SST and LANDMASK data. Dying now...\n";
    exit(0);

    for (int it = 0; it < ntimes; it++) {
        // Verify that the grid size and resolution from met_em file matches that in geom (from ERF inputs file).
        AMREX_ALWAYS_ASSERT(geom[lev].CellSizeArray()[0] == NC_dx[it]);
        AMREX_ALWAYS_ASSERT(geom[lev].CellSizeArray()[1] == NC_dy[it]);
        // NC_nx-2 because NC_nx is the number of staggered grid points indexed from 1.
        AMREX_ALWAYS_ASSERT(geom[lev].Domain().bigEnd(0) == NC_nx[it]-2);
        // NC_ny-2 because NC_ny is the number of staggered grid points indexed from 1.
        AMREX_ALWAYS_ASSERT(geom[lev].Domain().bigEnd(1) == NC_ny[it]-2);
    } // it

    // This makes the Jacobian.
    make_J(geom[lev],*z_phys,  *detJ_cc[lev]);

    // This defines z at w-cell faces.
    make_zcc(geom[lev],*z_phys,*z_phys_cc[lev]);

    // Set up FABs to hold data that will be used to set lateral boundary conditions.
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
    int MetGridBdyEnd = MetGridBdyVars::NumTypes;
#else
    int MetGridBdyEnd = MetGridBdyVars::NumTypes-1;
#endif
    //amrex::Vector<amrex::Vector<amrex::MultiFab> > fabs_for_bcs;
    amrex::Vector<amrex::Vector<FArrayBox>> fabs_for_bcs;
    fabs_for_bcs.resize(ntimes);
    for (int it(0); it < ntimes; it++) {
        fabs_for_bcs[it].resize(MetGridBdyEnd);
    }

    const Real l_rdOcp = solverChoice.rdOcp;

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Define FABs for holding some of the initial data
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        // Fill state data using origin data (initialization and BC arrays)
        //     x_vel   interpolated from origin levels
        //     y_vel   interpolated from origin levels
        //     z_vel   set to 0.0
        //     theta   calculate on origin levels then interpolate
        //     mxrat   convert RH -> Q on origin levels then interpolate
        init_state_from_metgrid(l_rdOcp,
                                cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                z_phys_nd_fab,
                                NC_hgt_fab, NC_ght_fab, NC_xvel_fab,
                                NC_yvel_fab, NC_temp_fab, NC_rhum_fab,
                                NC_pres_fab, theta_fab, mxrat_fab,
                                fabs_for_bcs);
    } // mf

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // Use map scale factors directly from the met_em files
    for ( MFIter mfi(*mapfac_u[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        // Define fabs for holding the initial data
        FArrayBox &msfu_fab = (*mapfac_u[lev])[mfi];
        FArrayBox &msfv_fab = (*mapfac_v[lev])[mfi];
        FArrayBox &msfm_fab = (*mapfac_m[lev])[mfi];

        init_msfs_from_metgrid(msfu_fab, msfv_fab, msfm_fab,
                               flag_msfu[0], flag_msfv[0], flag_msfm[0],
                               NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        FArrayBox&     p_hse_fab = p_hse[mfi];
        FArrayBox&    pi_hse_fab = pi_hse[mfi];
        FArrayBox&     r_hse_fab = r_hse[mfi];
        FArrayBox&      cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];

        // Fill base state data using origin data (initialization and BC arrays)
        //     p_hse     calculate dry pressure
        //     r_hse     calculate dry density
        //     pi_hse    calculate Exner term given pressure
        const Box valid_bx = mfi.validbox();
        init_base_state_from_metgrid(l_rdOcp,
                                     valid_bx,
                                     flag_psfc,
                                     cons_fab, r_hse_fab, p_hse_fab, pi_hse_fab,
                                     z_phys_nd_fab, NC_ght_fab, NC_psfc_fab,
                                     fabs_for_bcs);
    } // mf

    // Set up boxes for lateral boundary arrays.
    bdy_data_xlo.resize(ntimes);
    bdy_data_xhi.resize(ntimes);
    bdy_data_ylo.resize(ntimes);
    bdy_data_yhi.resize(ntimes);

    const auto& lo = geom[lev].Domain().loVect();
    const auto& hi = geom[lev].Domain().hiVect();
    amrex::IntVect plo(lo);
    amrex::IntVect phi(hi);

    plo[0] = lo[0];                     plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = lo[0]+metgrid_bdy_width-1; phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xlo(plo, phi);
    Box xlo_plane_no_stag(pbx_xlo);
    Box xlo_plane_x_stag = pbx_xlo; xlo_plane_x_stag.shiftHalf(0,-1);
    Box xlo_plane_y_stag = convert(pbx_xlo, {0, 1, 0});

    plo[0] = hi[0]-metgrid_bdy_width+1; plo[1] = lo[1]; plo[2] = lo[2];
    phi[0] = hi[0];                     phi[1] = hi[1]; phi[2] = hi[2];
    const Box pbx_xhi(plo, phi);
    Box xhi_plane_no_stag(pbx_xhi);
    Box xhi_plane_x_stag = pbx_xhi; xhi_plane_x_stag.shiftHalf(0,1);
    Box xhi_plane_y_stag = convert(pbx_xhi, {0, 1, 0});

    plo[1] = lo[1];                     plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = lo[1]+metgrid_bdy_width-1; phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_ylo(plo, phi);
    Box ylo_plane_no_stag(pbx_ylo);
    Box ylo_plane_x_stag = convert(pbx_ylo, {1, 0, 0});
    Box ylo_plane_y_stag = pbx_ylo; ylo_plane_y_stag.shiftHalf(1,-1);

    plo[1] = hi[1]-metgrid_bdy_width+1; plo[0] = lo[0]; plo[2] = lo[2];
    phi[1] = hi[1];                     phi[0] = hi[0]; phi[2] = hi[2];
    const Box pbx_yhi(plo, phi);
    Box yhi_plane_no_stag(pbx_yhi);
    Box yhi_plane_x_stag = convert(pbx_yhi, {1, 0, 0});
    Box yhi_plane_y_stag = pbx_yhi; yhi_plane_y_stag.shiftHalf(1,1);

#ifndef AMREX_USE_GPU
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tplo               \t" << plo << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tphi               \t" << phi << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txlo_plane_no_stag \t" << xlo_plane_no_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txlo_plane_x_stag  \t" << xlo_plane_x_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txlo_plane_y_stag  \t" << xlo_plane_y_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tplo               \t" << plo << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tphi               \t" << phi << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txhi_plane_no_stag \t" << xhi_plane_no_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txhi_plane_x_stag  \t" << xhi_plane_x_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \txhi_plane_y_stag  \t" << xhi_plane_y_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tplo               \t" << plo << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tphi               \t" << phi << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tylo_plane_no_stag \t" << ylo_plane_no_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tylo_plane_x_stag  \t" << ylo_plane_x_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tylo_plane_y_stag  \t" << ylo_plane_y_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tplo               \t" << plo << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tphi               \t" << phi << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tyhi_plane_no_stag \t" << yhi_plane_no_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tyhi_plane_x_stag  \t" << yhi_plane_x_stag << std::endl;
    amrex::AllPrint() << " proc-" << ParallelDescriptor::MyProc() << " \tyhi_plane_y_stag  \t" << yhi_plane_y_stag << std::endl;
#endif

    for (int ivar(MetGridBdyVars::U); ivar < MetGridBdyEnd; ivar++) {
        for (int it(0); it < ntimes; it++) {
            if (ivar == MetGridBdyVars::U) {
                bdy_data_xlo[it].push_back(FArrayBox(xlo_plane_x_stag, 1));
                bdy_data_xhi[it].push_back(FArrayBox(xhi_plane_x_stag, 1));
                bdy_data_ylo[it].push_back(FArrayBox(ylo_plane_x_stag, 1));
                bdy_data_yhi[it].push_back(FArrayBox(yhi_plane_x_stag, 1));
            } else if (ivar == MetGridBdyVars::V) {
                bdy_data_xlo[it].push_back(FArrayBox(xlo_plane_y_stag, 1));
                bdy_data_xhi[it].push_back(FArrayBox(xhi_plane_y_stag, 1));
                bdy_data_ylo[it].push_back(FArrayBox(ylo_plane_y_stag, 1));
                bdy_data_yhi[it].push_back(FArrayBox(yhi_plane_y_stag, 1));
            } else if (ivar == MetGridBdyVars::R) {
                bdy_data_xlo[it].push_back(FArrayBox(xlo_plane_no_stag, 1));
                bdy_data_xhi[it].push_back(FArrayBox(xhi_plane_no_stag, 1));
                bdy_data_ylo[it].push_back(FArrayBox(ylo_plane_no_stag, 1));
                bdy_data_yhi[it].push_back(FArrayBox(yhi_plane_no_stag, 1));
            } else if (ivar == MetGridBdyVars::T) {
                bdy_data_xlo[it].push_back(FArrayBox(xlo_plane_no_stag, 1));
                bdy_data_xhi[it].push_back(FArrayBox(xhi_plane_no_stag, 1));
                bdy_data_ylo[it].push_back(FArrayBox(ylo_plane_no_stag, 1));
                bdy_data_yhi[it].push_back(FArrayBox(yhi_plane_no_stag, 1));
            } else if (ivar == MetGridBdyVars::QV) {
                bdy_data_xlo[it].push_back(FArrayBox(xlo_plane_no_stag, 1));
                bdy_data_xhi[it].push_back(FArrayBox(xhi_plane_no_stag, 1));
                bdy_data_ylo[it].push_back(FArrayBox(ylo_plane_no_stag, 1));
                bdy_data_yhi[it].push_back(FArrayBox(yhi_plane_no_stag, 1));
            } else {
#ifndef AMREX_USE_GPU
                amrex::Print() << "Unexpected ivar " << ivar << std::endl;
#endif
                amrex::Abort("See Initialization/ERF_init_from_metgrid.cpp");
            }
        } // it
    } // ivar

    // Earlier we processed the entire domain at each time from the met_em files, even though
    // we only need the whole domain processed at initialization and the lateral boundaries
    // at subsequent times. We can optimize this later if needed. For now, we need to fill
    // the lateral boundary arrays using the info set aside earlier.
    for (int it(0); it < ntimes; it++) {

        const Array4<Real const>& R_bcs_arr = fabs_for_bcs[it][MetGridBdyVars::R].const_array();

        for (int ivar(MetGridBdyVars::U); ivar < MetGridBdyEnd; ivar++) {

            auto xlo_arr = bdy_data_xlo[it][ivar].array();
            auto xhi_arr = bdy_data_xhi[it][ivar].array();
            auto ylo_arr = bdy_data_ylo[it][ivar].array();
            auto yhi_arr = bdy_data_yhi[it][ivar].array();
            const Array4<Real const>& fabs_for_bcs_arr = fabs_for_bcs[it][ivar].const_array();

            if (ivar == MetGridBdyVars::U) {
                // xvel at west boundary
                amrex::ParallelFor(xlo_plane_x_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xlo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // xvel at east boundary
                amrex::ParallelFor(xhi_plane_x_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // xvel at south boundary
                amrex::ParallelFor(ylo_plane_x_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ylo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // xvel at north boundary
                amrex::ParallelFor(yhi_plane_x_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    yhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
            } // MetGridBdyVars::U

            if (ivar == MetGridBdyVars::V) {
                // yvel at west boundary
                amrex::ParallelFor(xlo_plane_y_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xlo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // yvel at east boundary
                amrex::ParallelFor(xhi_plane_y_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // yvel at south boundary
                amrex::ParallelFor(ylo_plane_y_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ylo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // yvel at north boundary
                amrex::ParallelFor(yhi_plane_y_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    yhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
            } // MetGridBdyVars::V

            if (ivar == MetGridBdyVars::R) {
                // density at west boundary
                amrex::ParallelFor(xlo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xlo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // density at east boundary
                amrex::ParallelFor(xhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // density at south boundary
                amrex::ParallelFor(ylo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ylo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // density at north boundary
                amrex::ParallelFor(yhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    yhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
            } // MetGridBdyVars::R

            if (ivar == MetGridBdyVars::T) {
                // potential temperature at west boundary
                amrex::ParallelFor(xlo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xlo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // potential temperature at east boundary
                amrex::ParallelFor(xhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // potential temperature at south boundary
                amrex::ParallelFor(ylo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ylo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
                // potential temperature at north boundary
                amrex::ParallelFor(yhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    yhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                });
            } // MetGridBdyVars::T

            if (ivar == MetGridBdyVars::QV) {
                // mixing ratio at west boundary
                amrex::ParallelFor(xlo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xlo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                    xlo_arr(i,j,k,0) *= R_bcs_arr(i,j,k);
                });
                // mixing ratio at east boundary
                amrex::ParallelFor(xhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    xhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                    xhi_arr(i,j,k,0) *= R_bcs_arr(i,j,k);
                });
                // mixing ratio at south boundary
                amrex::ParallelFor(ylo_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    ylo_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                    ylo_arr(i,j,k,0) *= R_bcs_arr(i,j,k);
                });
                // mixing ratio at north boundary
                amrex::ParallelFor(yhi_plane_no_stag, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    yhi_arr(i,j,k,0) = fabs_for_bcs_arr(i,j,k);
                    yhi_arr(i,j,k,0) *= R_bcs_arr(i,j,k);
                });
            } // MetGridBdyVars::QV

        } // ivar

    } // it

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

   for (int it = 0; it < ntimes; it++) {
#ifndef AMREX_USE_GPU
        amrex::Print() << " SIZE OF HGT FAB " << NC_hgt_fab[it].box() << std::endl;
        amrex::Print() << " SIZE OF ZP FAB "  << z_phys_nd_fab.box() << std::endl;
#endif

        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys_nd_fab.array();
        const Array4<Real const>& nc_hgt_arr = NC_hgt_fab[it].const_array();

        const Box z_hgt_box = NC_hgt_fab[it].box();

        int ilo = z_hgt_box.smallEnd()[0];
        int ihi = z_hgt_box.bigEnd()[0];
        int jlo = z_hgt_box.smallEnd()[1];
        int jhi = z_hgt_box.bigEnd()[1];

        Box z_phys_box = z_phys_nd_fab.box();
        Box from_box = surroundingNodes(NC_hgt_fab[it].box());
        from_box.growHi(2,-1);

        Box bx = z_phys_box & from_box;
        Box bxu = bx; bxu.growHi(0,1);
        Box bxv = bx; bxv.growHi(1,1);

#ifndef AMREX_USE_GPU
        amrex::Print() << "FROM BOX " << from_box << std::endl;
        amrex::Print() << "BX " << bx << std::endl;
#endif

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            z_arr(i,j,k) =  0.25 * ( nc_hgt_arr (ii,jj  ,k) + nc_hgt_arr(ii-1,jj  ,k) +
                                     nc_hgt_arr (ii,jj-1,k) + nc_hgt_arr(ii-1,jj-1,k) );
        });
    } // it
}

/**
 * Helper function to initialize state and velocity data read from metgrid data.
 *
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param state_fab FArrayBox holding the state data to initialize
 * @param x_vel_fab FArrayBox holding the x-velocity data to initialize
 * @param y_vel_fab FArrayBox holding the y-velocity data to initialize
 * @param z_vel_fab FArrayBox holding the z-velocity data to initialize
 * @param z_phys_nd_fab FArrayBox holding nodal z coordinate data for terrain
 * @param NC_hgt_fab Vector of FArrayBox obects holding metgrid data for terrain height
 * @param NC_ght_fab Vector of FArrayBox objects holding metgrid data for height of cell centers
 * @param NC_xvel_fab Vector of FArrayBox obects holding metgrid data for x-velocity
 * @param NC_yvel_fab Vector of FArrayBox obects holding metgrid data for y-velocity
 * @param NC_zvel_fab Vector of FArrayBox obects holding metgrid data for z-velocity
 * @param NC_temp_fab Vector of FArrayBox obects holding metgrid data for temperature
 * @param NC_rhum_fab Vector of FArrayBox obects holding metgrid data for relative humidity
 * @param NC_pres_fab Vector of FArrayBox obects holding metgrid data for pressure
 * @param theta_fab Vector of FArrayBox obects holding potential temperature calculated from temperature and pressure
 * @param mxrat_fab Vector of FArrayBox obects holding vapor mixing ratio calculated from relative humidity
 * @param fabs_for_bcs Vector of Vector of FArrayBox objects holding MetGridBdyVars at each met_em time.
 */
void
init_state_from_metgrid (const Real l_rdOcp,
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
                         Vector<FArrayBox>& theta_fab,
                         Vector<FArrayBox>& mxrat_fab,
                         amrex::Vector<amrex::Vector<FArrayBox>>& fabs_for_bcs)
{
    int ntimes = NC_hgt_fab.size();
    for (int it = 0; it < ntimes; it++)
    {
        // ********************************************************
        // U
        // ********************************************************
        if (it == 0) { // update at initialization
            Box bx2d = NC_xvel_fab[it].box() & x_vel_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = NC_xvel_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = x_vel_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'X',0,0,orig_z,orig_data,new_z,new_data);
            });
        }

        { // save for setting boundary conditions
#ifdef AMREX_USE_GPU
            fabs_for_bcs[it][MetGridBdyVars::U].resize(x_vel_fab.box(), 1, The_Pinned_Arena());
#else
            fabs_for_bcs[it][MetGridBdyVars::U].resize(x_vel_fab.box(), 1);
#endif
            Box bx2d = NC_xvel_fab[it].box() & x_vel_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = NC_xvel_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = fabs_for_bcs[it][MetGridBdyVars::U].array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'X',0,0,orig_z,orig_data,new_z,new_data);
            });
        }

        // ********************************************************
        // V
        // ********************************************************
        if (it == 0) { // update at initialization
            Box bx2d = NC_yvel_fab[it].box() & y_vel_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = NC_yvel_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = y_vel_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'Y',0,0,orig_z,orig_data,new_z,new_data);
            });
        }

        { // save for setting boundary conditions
#ifdef AMREX_USE_GPU
            fabs_for_bcs[it][MetGridBdyVars::V].resize(y_vel_fab.box(), 1, The_Pinned_Arena());
#else
            fabs_for_bcs[it][MetGridBdyVars::V].resize(y_vel_fab.box(), 1);
#endif
            Box bx2d = NC_yvel_fab[it].box() & y_vel_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = NC_yvel_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = fabs_for_bcs[it][MetGridBdyVars::V].array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'Y',0,0,orig_z,orig_data,new_z,new_data);
            });
        }

        // ********************************************************
        // W
        // ********************************************************
        if (it == 0) { // update at initialization
            z_vel_fab.template setVal<RunOn::Device>(0.0);
        }

        // ********************************************************
        // Initialize all state_fab variables to zero
        // ********************************************************
        if (it == 0) { // update at initialization
            state_fab.template setVal<RunOn::Device>(0.0);
        }

        // ********************************************************
        // theta
        // ********************************************************
        { // calculate potential temperature.
            Box bx = NC_rhum_fab[it].box() & state_fab.box();
            auto const temp  = NC_temp_fab[it].const_array();
            auto const pres  = NC_pres_fab[it].const_array();
            auto       theta = theta_fab[it].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                theta(i,j,k) = getThgivenPandT(temp(i,j,k),pres(i,j,k),l_rdOcp);
//                theta(i,j,k) = 300.0; // TODO: Remove when not needed. Force an isothermal atmosphere for debugging.
            });
        }

        // vertical interpolation of potential temperature.
        if (it == 0) { // update at initialization
            Box bx2d = NC_temp_fab[it].box() & state_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = theta_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = state_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'M',0,RhoTheta_comp,orig_z,orig_data,new_z,new_data);
            });
        }

        { // save for setting boundary conditions
#ifdef AMREX_USE_GPU
            fabs_for_bcs[it][MetGridBdyVars::T].resize(state_fab.box(), 1, The_Pinned_Arena());
#else
            fabs_for_bcs[it][MetGridBdyVars::T].resize(state_fab.box(), 1);
#endif
            Box bx2d = NC_temp_fab[it].box() & state_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = theta_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = fabs_for_bcs[it][MetGridBdyVars::T].array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'Y',0,0,orig_z,orig_data,new_z,new_data);
            });
        }

        // ********************************************************
        // specific humidity / relative humidity / mixing ratio
        // ********************************************************
        // TODO: we will need to check what input data we have for moisture
        // and then, if necessary, compute mixing ratio. For now, we will
        // focus on the case where we have relative humidity. Alternate cases
        // could be specific humidity or a mixing ratio.
        //
        { // calculate vapor mixing ratio from relative humidity.
            Box bx = NC_rhum_fab[it].box() & state_fab.box();
            auto const rhum  = NC_rhum_fab[it].const_array();
            auto const temp  = NC_temp_fab[it].const_array();
            auto const pres  = NC_pres_fab[it].const_array();
            auto       mxrat = mxrat_fab[it].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                rh_to_mxrat(i,j,k,rhum,temp,pres,mxrat);
            });
        }

        // vertical interpolation of vapor mixing ratio.
        if (it == 0) { // update at initialization
            Box bx2d = NC_rhum_fab[it].box() & state_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = mxrat_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = state_fab.array();
            auto const new_z     = z_phys_nd_fab.const_array();

#if defined(ERF_USE_MOISTURE)
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'M',0,RhoQt_comp,orig_z,orig_data,new_z,new_data);
            });
#elif defined(ERF_USE_WARM_NO_PRECIP)
            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'M',0,RhoQv_comp,orig_z,orig_data,new_z,new_data);
            });
#endif
        }

        { // save for setting boundary conditions
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
#ifdef AMREX_USE_GPU
            fabs_for_bcs[it][MetGridBdyVars::QV].resize(state_fab.box(), 1, The_Pinned_Arena());
#else
            fabs_for_bcs[it][MetGridBdyVars::QV].resize(state_fab.box(), 1);
#endif
            Box bx2d = NC_rhum_fab[it].box() & state_fab.box();
            bx2d.setRange(2,0);
            auto const orig_data = mxrat_fab[it].const_array();
            auto const orig_z    = NC_ght_fab[it].const_array();
            auto       new_data  = fabs_for_bcs[it][MetGridBdyVars::QV].array();
            auto const new_z     = z_phys_nd_fab.const_array();

            ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
                interpolate_column_metgrid(i,j,'M',0,0,orig_z,orig_data,new_z,new_data);
            });
#endif
        }

        // TODO: TEMPORARY CODE TO RUN QUIESCENT, REMOVE WHEN NOT NEEDED.
//        if (it == 0) {
//            x_vel_fab.template setVal<RunOn::Device>(0.0); // TODO: temporary code to initialize with quiescent atmosphere.
//            y_vel_fab.template setVal<RunOn::Device>(0.0); // TODO: temporary code to initialize with quiescent atmosphere.
//        }
//        fabs_for_bcs[it][MetGridBdyVars::U].template setVal<RunOn::Device>(0.0); // TODO: temporary code to force with quiescent atmosphere.
//        fabs_for_bcs[it][MetGridBdyVars::V].template setVal<RunOn::Device>(0.0); // TODO: temporary code to force with quiescent atmosphere.

    } // it
}

/**
 * Helper function to calculate vapor mixing ratio from relative humidity, pressure and temperature.
 *
 * @param i Integer specifying the x-dimension index for the column
 * @param j Integer specifying the y-dimension index for the column
 * @param k Integer specifying the z-dimension index for the column
 * @param rhum Array4 object holding relative humidty from the metgrid data
 * @param temp Array4 object holding temperature from the metgrid data
 * @param pres Array4 object holding pressure from the metgrid data
 * @param mxrat Array4 object to hold the calculated vapor mixing ratio
 */
void
rh_to_mxrat (int i, int j, int k,
             const Array4<Real const>& rhum,
             const Array4<Real const>& temp,
             const Array4<Real const>& pres,
             const Array4<Real>& mxrat)
{
    Real qv_max_p_safe = 10000.0; // WRF default value
    Real qv_max_flag = std::pow(10.0, -5); // WRF default value
    Real qv_max_value = 3.0*std::pow(10.0, -6); // WRF default value
    Real qv_min_p_safe = 110000.0; // WRF default value
    Real qv_min_flag = std::pow(10.0, -6); // WRF default value
    Real qv_min_value = 1.0*std::pow(10.0, -6); // WRF default value
    Real eps = 0.622;
    Real svp1 = 0.6112;
    Real svp2 = 17.67;
    Real svp3 = 29.65;
    Real svpt0 = 273.15;
    // WRF's method when model_config_rec%rh2qv_wrt_liquid=.true. (default behavior)
    if (temp(i,j,k) != 0.0) {
        Real es=0.01*rhum(i,j,k)*svp1*10.0*exp(svp2*(temp(i,j,k)-svpt0)/(temp(i,j,k)-svp3));
        if (es >= pres(i,j,k)/100.0) {
            // vapor pressure exceeds total pressure
            mxrat(i,j,k) = std::pow(10.0, -6);
        }
        else {
            mxrat(i,j,k) = max(eps*es/(pres(i,j,k)/100.0-es), std::pow(10.0, -6));
        }
    }
    else {
        // I don't know why there's a fringe case handled in WRF where T is absolute zero...
        // Let's just deal with it here in case we also end up needing it.
        mxrat(i,j,k) = std::pow(10.0, -6);
    }
    // See the below comment from WRF dyn_em/module_initialize_real.F rh_to_mxrat1.
    // For pressures above a defined level, reasonable Qv values should be
    // a certain value or smaller. If they are larger than this, the input data
    // probably had "missing" RH, and we filled in some values. This is an
    // attempt to catch those. Also, set the minimum value for the entire
    // domain that is above the selected pressure level.
    if (pres(i,j,k) < qv_max_p_safe) {
        if (mxrat(i,j,k) > qv_max_flag) {
            mxrat(i,j,k) = qv_max_value;
        }
    }
    if (pres(i,j,k) < qv_min_p_safe) {
        if (mxrat(i,j,k) < qv_min_flag) {
            mxrat(i,j,k) = qv_min_value;
        }
    }
}

/**
 * Helper function for calculating density and pressure in hydrostatic balance given metgrid data.
 *
 * @param kmax Integer specifying the number of nodal points in the z-dimension
 * @param flag_psfc Integer 1 if surface pressure is in metgrid data, 0 otherwise
 * @param psfc Real holding surface pressure from the metgrid data (if available)
 * @param Theta_vec Vector of Real holding potential temperature calculated from the metgrid data
 * @param Q_vec Vector of Real holding mixing ratio calcualted from the metgrid data
 * @param z_vec Vector of Real holding the z-coordinate of the ERF grid
 * @param Rhod_vec Vector of Real to hold dry density calculated within
 * @param Pd_vec Vector of Real to hold dry pressure calculated within
 */
void
calc_rho_p (const int& kmax,
            const int& flag_psfc,
            const Real& psfc,
            const Vector<Real>& Theta_vec,
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
            const Vector<Real>& Q_vec,
#endif
            const Vector<Real>& z_vec,
            Vector<Real>& Rhod_vec,
            Vector<Real>& Pd_vec)
{

    const int maxiter = 10;

    // Set up vectors for some components that are calculated as we integrate up and down the column.
    std::vector<Real> theta_v(kmax);
    std::vector<Real> Pm_vec(kmax);
    std::vector<Real> Rhom_vec(kmax);

    // Calculate or use moist pressure at the surface.
    if (flag_psfc == 1) {
        Pm_vec[0] = psfc;
    } else {
        Real t_0 = 290.0; // WRF's model_config_rec%base_temp
        Real a = 50.0; // WRF's model_config_rec%base_lapse
        Pm_vec[0] = p_0*exp(-t_0/a+std::pow((std::pow(t_0/a, 2)-2.0*CONST_GRAV*z_vec[0]/(a*R_d)), 0.5));
    }

    // calculate virtual potential temperature at the surface.
    {
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
        Real qvf = 1.0+(R_v/R_d+1.0)*Q_vec[0];
#else
        Real qvf = 1.0;
#endif
        theta_v[0] = Theta_vec[0]*qvf;
    }

    // calculate moist density at the surface.
    Rhom_vec[0] = 1.0/(R_d/p_0*theta_v[0]*std::pow(Pm_vec[0]/p_0, -iGamma));

    // integrate from the surface to the top boundary.
    for (int k=1; k < kmax; k++) {
        Real dz = z_vec[k]-z_vec[k-1];
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
        Real qvf = 1.0+(R_v/R_d+1.0)*Q_vec[k];
#else
        Real qvf = 1.0;
#endif
        theta_v[k] = Theta_vec[k]*qvf;
        Rhom_vec[k] = Rhom_vec[k-1]; // an initial guess.
        for (int it=0; it < maxiter; it++) {
            Pm_vec[k] = Pm_vec[k-1]-0.5*dz*(Rhom_vec[k]+Rhom_vec[k-1])*CONST_GRAV;
            if (Pm_vec[k] < 0.0) Pm_vec[k] = 0.0;
            Rhom_vec[k] = 1.0/(R_d/p_0*theta_v[k]*std::pow(Pm_vec[k]/p_0, -iGamma));
        } // it
    } // k

    // integrate from the top back down to get dry pressure and density.
    Pd_vec[kmax-1] = Pm_vec[kmax-1];
    Rhod_vec[kmax-1] = 1.0/(R_d/p_0*theta_v[kmax-1]*std::pow(Pd_vec[kmax-1]/p_0, -iGamma));
    for (int k=kmax-2; k >= 0; k--) {
        Real dz = z_vec[k+1]-z_vec[k];
        Rhod_vec[k] = Rhod_vec[k+1]; // an initial guess.
        for (int it=0; it < maxiter; it++) {
            Pd_vec[k] = Pd_vec[k+1]+0.5*dz*(Rhod_vec[k]+Rhod_vec[k+1])*CONST_GRAV;
            if (Pd_vec[k] < 0.0) Pd_vec[k] = 0.0;
            Rhod_vec[k] = 1.0/(R_d/p_0*theta_v[k]*std::pow(Pd_vec[k]/p_0, -iGamma));
        } // it
    } // k
}

/**
 * Helper function to initialize map factors from metgrid data
 *
 * @param msfu_fab FArrayBox specifying x-velocity map factors
 * @param msfv_fab FArrayBox specifying y-velocity map factors
 * @param msfm_fab FArrayBox specifying z-velocity map factors
 * @param flag_msfu Integer 1 if u-staggered map factor is in metgrid data, 0 otherwise
 * @param flag_msfv Integer 1 if v-staggered map factor is in metgrid data, 0 otherwise
 * @param flag_msfm Integer 1 if cell center map factor is in metgrid data, 0 otherwise
 * @param NC_MSFU_fab Vector of FArrayBox objects holding metgrid data for x-velocity map factors
 * @param NC_MSFV_fab Vector of FArrayBox objects holding metgrid data for y-velocity map factors
 * @param NC_MSFM_fab Vector of FArrayBox objects holding metgrid data for z-velocity map factors
 */
void
init_msfs_from_metgrid (FArrayBox& msfu_fab,
                        FArrayBox& msfv_fab,
                        FArrayBox& msfm_fab,
                        const int& flag_msfu,
                        const int& flag_msfv,
                        const int& flag_msfm,
                        const Vector<FArrayBox>& NC_MSFU_fab,
                        const Vector<FArrayBox>& NC_MSFV_fab,
                        const Vector<FArrayBox>& NC_MSFM_fab)
{
//    int ntimes = NC_MSFU_fab.size();
    int ntimes = 1;
    for (int it = 0; it < ntimes; it++) {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //

        // This copies or sets mapfac_m
        if (flag_msfm == 1) {
            msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab[it]);
        } else {
#ifndef AMREX_USE_GPU
            amrex::Print() << " MAPFAC_M not present in met_em files. Setting to 1.0" << std::endl;
#endif
            msfm_fab.template setVal<RunOn::Device>(1.0);
        }

        // This copies or sets mapfac_u
        if (flag_msfu == 1) {
            msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab[it]);
        } else {
#ifndef AMREX_USE_GPU
            amrex::Print() << " MAPFAC_U not present in met_em files. Setting to 1.0" << std::endl;
#endif
            msfu_fab.template setVal<RunOn::Device>(1.0);
        }

        // This copies or sets mapfac_v
        if (flag_msfv == 1) {
            msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab[it]);
        } else {
#ifndef AMREX_USE_GPU
            amrex::Print() << " MAPFAC_V not present in met_em files. Setting to 1.0" << std::endl;
#endif
            msfv_fab.template setVal<RunOn::Device>(1.0);
        }

    } // it
}

/**
 * Helper function for initializing hydrostatic base state data from metgrid data
 *
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param valid_bx Box specifying the index space we are to initialize
 * @param flag_psfc Vector of Integer 1 if surface pressure is in metgrid data, 0 otherwise
 * @param state_fab FArrayBox holding the state data to initialize
 * @param r_hse_fab FArrayBox holding the hydrostatic base state density we are initializing
 * @param p_hse_fab FArrayBox holding the hydrostatic base state pressure we are initializing
 * @param pi_hse_fab FArrayBox holding the hydrostatic base Exner pressure we are initializing
 * @param z_phys_nd_fab FArrayBox holding nodal z coordinate data for terrain
 * @param NC_ght_fab Vector of FArrayBox objects holding metgrid data for height of cell centers
 * @param NC_psfc_fab Vector of FArrayBox objects holding metgrid data for surface pressure
 * @param fabs_for_bcs Vector of Vector of FArrayBox objects holding MetGridBdyVars at each met_em time.
 */
void
init_base_state_from_metgrid (const Real l_rdOcp,
                              const Box& valid_bx,
                              const Vector<int>& flag_psfc,
                              FArrayBox& state_fab,
                              FArrayBox& r_hse_fab,
                              FArrayBox& p_hse_fab,
                              FArrayBox& pi_hse_fab,
                              FArrayBox& z_phys_nd_fab,
                              const Vector<FArrayBox>& NC_ght_fab,
                              const Vector<FArrayBox>& NC_psfc_fab,
                              amrex::Vector<amrex::Vector<FArrayBox>>& fabs_for_bcs)
{
#if defined(ERF_USE_MOISTURE)
    int RhoQ_comp = RhoQt_comp;
#elif defined(ERF_USE_WARM_NO_PRECIP)
    int RhoQ_comp = RhoQv_comp;
#endif
    int kmax = amrex::ubound(state_fab.box()).z-1;
    { // set pressure and density at initialization.
        const Array4<Real>& r_hse_arr = r_hse_fab.array();
        const Array4<Real>& p_hse_arr = p_hse_fab.array();
        const Array4<Real>& pi_hse_arr = pi_hse_fab.array();

        // ********************************************************
        // calculate dry density and dry pressure
        // ********************************************************
        // calculate density and dry pressure on the new grid.
        Box valid_bx2d = valid_bx;
        valid_bx2d.setRange(2,0);
        auto const orig_psfc = NC_psfc_fab[0].const_array();
        auto       new_data = state_fab.array();
        auto const new_z = z_phys_nd_fab.const_array();
        amrex::ParallelFor(valid_bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
            Vector<Real> z_vec;     z_vec.resize(kmax+1);
            Vector<Real> Theta_vec; Theta_vec.resize(kmax);
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
            Vector<Real> Q_vec;     Q_vec.resize(kmax);
#endif
            Vector<Real> Rhod_vec;  Rhod_vec.resize(kmax);
            Vector<Real> Pd_vec;    Pd_vec.resize(kmax);
            for (int k=0; k < kmax+1; k++) z_vec[k] = new_z(i,j,k);
            for (int k=0; k < kmax; k++) {
                Theta_vec[k] = new_data(i,j,k,RhoTheta_comp);
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
                Q_vec[k] = new_data(i,j,k,RhoQ_comp);
#endif
            }

            calc_rho_p(kmax,flag_psfc[0],orig_psfc(i,j,0),Theta_vec,
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
                       Q_vec,
#endif
                       z_vec,Rhod_vec,Pd_vec);
            for (int k=0; k < kmax; k++) {
                p_hse_arr(i,j,k) = Pd_vec[k];
                r_hse_arr(i,j,k) = Rhod_vec[k];
            }
        });

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            new_data(i,j,k,Rho_comp) = r_hse_arr(i,j,k);
            new_data(i,j,k,RhoScalar_comp) = 0.0;
            // RhoTheta and RhoQt or RhoQv currently hold Theta and Qt or Qv. Multiply by Rho.
            Real RhoTheta = r_hse_arr(i,j,k)*new_data(i,j,k,RhoTheta_comp);
            new_data(i,j,k,RhoTheta_comp) = RhoTheta;
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
            Real RhoQ = r_hse_arr(i,j,k)*new_data(i,j,k,RhoQ_comp);
            new_data(i,j,k,RhoQ_comp) = RhoQ;
#endif
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
//            pi_hse_arr(i,j,k) = getExnergivenRTh(RhoTheta, l_rdOcp);
        });
    }

    int ntimes = NC_psfc_fab.size();
    for (int it=0; it < ntimes; it++) {
        FArrayBox p_hse_bcs_fab;
        FArrayBox pi_hse_bcs_fab;
#ifdef AMREX_USE_GPU
        fabs_for_bcs[it][MetGridBdyVars::R].resize(state_fab.box(), 1, The_Pinned_Arena());
        p_hse_bcs_fab.resize(state_fab.box(), 1, The_Pinned_Arena());
#else
        fabs_for_bcs[it][MetGridBdyVars::R].resize(state_fab.box(), 1);
        p_hse_bcs_fab.resize(state_fab.box(), 1);
#endif

        // ********************************************************
        // calculate dry density and dry pressure
        // ********************************************************
        // calculate density and dry pressure on the new grid.
        Box valid_bx2d = valid_bx;
        valid_bx2d.setRange(2,0);
        auto const orig_psfc = NC_psfc_fab[it].const_array();
        auto const new_z = z_phys_nd_fab.const_array();
        auto       Theta_arr = fabs_for_bcs[it][MetGridBdyVars::T].array();
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
        auto       Q_arr = fabs_for_bcs[it][MetGridBdyVars::QV].array();
#endif
        auto r_hse_arr = fabs_for_bcs[it][MetGridBdyVars::R].array();
        auto p_hse_arr = p_hse_bcs_fab.array();
        amrex::ParallelFor(valid_bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
            Vector<Real> z_vec;     z_vec.resize(kmax+1);
            Vector<Real> Theta_vec; Theta_vec.resize(kmax);
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
            Vector<Real> Q_vec;     Q_vec.resize(kmax);
#endif
            Vector<Real> Rhod_vec;  Rhod_vec.resize(kmax);
            Vector<Real> Pd_vec;    Pd_vec.resize(kmax);
            for (int k=0; k < kmax+1; k++) z_vec[k] = new_z(i,j,k);
            for (int k=0; k < kmax; k++) {
                Theta_vec[k] = Theta_arr(i,j,k);
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
                Q_vec[k] = Q_arr(i,j,k);
#endif
            } // k
            calc_rho_p(kmax,flag_psfc[it],orig_psfc(i,j,0),Theta_vec,
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
                       Q_vec,
#endif
                       z_vec,Rhod_vec,Pd_vec);
            for (int k=0; k < kmax; k++) {
                r_hse_arr(i,j,k) = Rhod_vec[k];
                p_hse_arr(i,j,k) = Pd_vec[k];
#if defined(ERF_USE_MOISTURE) || defined(ERF_USE_WARM_NO_PRECIP)
                Q_arr(i,j,k) = Rhod_vec[k]*Q_vec[k];
#endif
                Theta_arr(i,j,k) = Rhod_vec[k]*Theta_vec[k];
            } // k
        });
    } // it
}

/**
 * Helper function to interpolate data and its associated z-coordinates to a new set of z-coordinates.
 *
 * Operates on a column of unsorted data (with fixed x,y coordinates) at a time.
 *
 * @param i Integer specifying the x-dimension index for the column
 * @param j Integer specifying the y-dimension index for the column
 * @param stag Char specifying the variable staggering ("X", "Y", or "M")
 * @param src_comp Integer specifying the source component of the data to start intepolation from
 * @param dest_comp Integer specifying the destination component of the data to interpolate into
 * @param orig_z Array4 object containing the original z-coordinates to interpolate from
 * @param orig_data Array4 object containing the data at the original z-coordinates to interpolate
 * @param new_z Array4 object containing the new z-coordinates to interpolate into
 * @param new_data Array4 object into which we interpolate the data at the new z-coordinates
 */
AMREX_GPU_DEVICE
void
interpolate_column_metgrid (int i, int j, char stag, int src_comp, int dest_comp,
                            const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                            const Array4<Real const>&  new_z, const Array4<Real>&  new_data)
{
    // This subroutine is a bit ham-handed and can be cleaned up later.
    int kmax = amrex::ubound(Box(new_data)).z;
    int imax_orig = amrex::ubound(Box(orig_data)).x;
    int jmax_orig = amrex::ubound(Box(orig_data)).y;
    int kmax_orig = amrex::ubound(Box(orig_data)).z;
    for (int k = 0; k < kmax-1; k++) {
        Real z;
        if (stag == 'X') {
            z = 0.25*(new_z(i,j,k)+new_z(i,j+1,k)+new_z(i,j,k+1)+new_z(i,j+1,k+1));
        }
        else if (stag == 'Y') {
            z = 0.25*(new_z(i,j,k)+new_z(i+1,j,k)+new_z(i,j,k+1)+new_z(i+1,j,k+1));
        }
        else if (stag == 'M') {
            z = 0.125*(new_z(i,j,k  )+new_z(i,j+1,k  )+new_z(i+1,j,k  )+new_z(i+1,j+1,k  )+
                       new_z(i,j,k+1)+new_z(i,j+1,k+1)+new_z(i+1,j,k+1)+new_z(i+1,j+1,k+1));
        }
        Real z0, z1;
        int klow = -1;
        int khi0 = -1;
        Real dzlow = 10000.0;
        Real dzhi0 = -10000.0;
        for (int kk = 0; kk < kmax_orig; kk++) {
            Real orig_z_stag = 0.0;
            if (stag == 'M') {
                orig_z_stag = orig_z(i,j,kk);
            }
            if (stag == 'X') {
                if (i == 0) {
                    orig_z_stag = orig_z(i,j,kk);
                }
                else if (i == imax_orig) {
                    orig_z_stag = orig_z(imax_orig-1,j,kk);
                }
                else {
                    orig_z_stag = (orig_z(i,j,kk)+orig_z(i-1,j,kk))/2.0;
                }
            }
            else if (stag == 'Y') {
                if (j == 0) {
                    orig_z_stag = orig_z(i,j,kk);
                }
                else if (j == jmax_orig) {
                    orig_z_stag = orig_z(i,jmax_orig-1,kk);
                }
                else {
                    orig_z_stag = (orig_z(i,j,kk)+orig_z(i,j-1,kk))/2.0;
                }
            }
            Real dz = z - orig_z_stag;
            if ((dz < 0.0) && (dz > dzhi0)) {
                dzhi0 = dz;
                khi0 = kk;
                z1 = orig_z_stag;
            }
            if ((dz >= 0.0) && (dz < dzlow)) {
                dzlow = dz;
                klow = kk;
                z0 = orig_z_stag;
            }
        }
        if (klow == -1) {
            // extrapolate
            int khi1 = -1;
            Real dzhi1 = -10000.0;
            for (int kk = 0; kk < kmax_orig; k++) {
                Real orig_z_stag = 0.0;
                if (stag == 'M') {
                    orig_z_stag = orig_z(i,j,kk);
                }
                else if (stag == 'X') {
                    if (i == 0) {
                        orig_z_stag = orig_z(i,j,kk);
                    }
                    else if (i == imax_orig) {
                        orig_z_stag = orig_z(imax_orig-1,j,kk);
                    }
                    else {
                        orig_z_stag = (orig_z(i,j,kk)+orig_z(i-1,j,kk))/2.0;
                    }
                }
                else if (stag == 'Y') {
                    if (j == 0) {
                        orig_z_stag = orig_z(i,j,kk);
                    }
                    else if (j == jmax_orig) {
                        orig_z_stag = orig_z(i,jmax_orig-1,kk);
                    }
                    else {
                        orig_z_stag = (orig_z(i,j,kk)+orig_z(i,j-1,kk))/2.0;
                    }
                }
                Real dz = z - orig_z_stag;
                if ((dz < 0.0) && (dz > dzhi1) && (kk != khi0)) {
                    dzhi1 = dz;
                    khi1 = kk;
                    z1 = orig_z_stag;
                }
            }
            Real y0 = orig_data(i,j,khi0,src_comp);
            Real y1 = orig_data(i,j,khi1,src_comp);
            new_data(i,j,k,dest_comp) = y0-(y1-y0)/(z1-z0)*(z0-z);
        } else {
            // interpolate
            Real y0 = orig_data(i,j,klow,src_comp);
            Real y1 = orig_data(i,j,khi0,src_comp);
            new_data(i,j,k,dest_comp) = y0+(y1-y0)/(z1-z0)*(z-z0);
        }
    }
}

/**
 * Helper function to interpolate data and its associated z-coordinates to a new set of z-coordinates.
 *
 * Operates on a column of data (with fixed x,y coordinates) at a time.
 *
 * @param i Integer specifying the x-dimension index for the column
 * @param j Integer specifying the y-dimension index for the column
 * @param src_comp Integer specifying the source component of the data to start intepolation from
 * @param dest_comp Integer specifying the destination component of the data to interpolate into
 * @param orig_z Array4 object containing the original z-coordinates to interpolate from
 * @param orig_data Array4 object containing the data at the original z-coordinates to interpolate
 * @param new_z Array4 object containing the new z-coordinates to interpolate into
 * @param new_data Array4 object into which we interpolate the data at the new z-coordinates
 */
AMREX_GPU_DEVICE
void
interpolate_column (int i, int j, int src_comp, int dest_comp,
                    const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                    const Array4<Real const>&  new_z, const Array4<Real>&  new_data)
{
    // CAVEAT: we only consider interpolation here - if we go past end of array this won't work for now

    int kmax     = amrex::ubound(Box(new_data)).z;
    int kmax_old = amrex::ubound(Box(orig_data)).z;
#ifndef AMREX_USE_GPU
    amrex::Print() << "KMAX     IN INTERP " << kmax << std::endl;
    amrex::Print() << "KMAX_OLD IN INTERP " << kmax_old << std::endl;
#endif
    int klast = 0;
    for (int k = 0; k < kmax; k++) {

        Real z = new_z(i,j,k);

        bool kfound = false;

#ifndef AMREX_USE_GPU
        amrex::Print() << "AT THIS K WE ARE USING KLAST " << k << " " << klast << std::endl;
        amrex::Print() << "z = " << z << std::endl;
#endif
        for (int kk = klast; kk < kmax_old; kk++) {
           if (z >= orig_z(i,j,k)) {
              Real y0 = orig_data(i,j,kk  ,src_comp);
              Real y1 = orig_data(i,j,kk+1,src_comp);
              Real z0 = orig_z(i,j,kk  );
              Real z1 = orig_z(i,j,kk+1);
              new_data(i,j,k,dest_comp) = y0 + (y1 - y0)*(z - z0) / (z1 - z0);
              klast = kk;
              return;
           }
        }
        if (!kfound) amrex::Error("Wasnt able to interpolate");
    }
}

#endif // ERF_USE_NETCDF
