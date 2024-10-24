/**
 * \file ERF_init_from_wrfinput.cpp
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_prob_common.H>
#include <ERF_DataStruct.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

void
read_from_wrfinput (int lev, const Box& domain, const std::string& fname,
                    FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                    FArrayBox& NC_zvel_fab, FArrayBox& NC_rho_fab,
                    FArrayBox& NC_rhop_fab, FArrayBox& NC_theta_fab,
                    FArrayBox& NC_MUB_fab ,
                    FArrayBox& NC_MSFU_fab, FArrayBox& NC_MSFV_fab,
                    FArrayBox& NC_MSFM_fab, FArrayBox& NC_SST_fab,
                    FArrayBox& NC_LANDMSK_fab,
                    FArrayBox& NC_C1H_fab , FArrayBox& NC_C2H_fab,
                    FArrayBox& NC_RDNW_fab,
                    FArrayBox& NC_QVAPOR_fab,
                    FArrayBox& NC_QCLOUD_fab,
                    FArrayBox& NC_QRAIN_fab,
                    FArrayBox& NC_PH_fab,
                    FArrayBox& NC_P_fab,
                    FArrayBox& NC_PHB_fab,
                    FArrayBox& NC_ALB_fab,
                    FArrayBox& NC_PB_fab,
                    FArrayBox& NC_LAT_fab,
                    FArrayBox& NC_LON_fab,
                    MoistureType moisture_type,
                    Real& Latitude,
                    Real& Longitude,
                    Geometry& geom);

Real
read_from_wrfbdy (const std::string& nc_bdy_file, const Box& domain,
                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                  Vector<Vector<FArrayBox>>& bdy_data_yhi,
                  int& width, Real& start_bdy_time);

void
convert_wrfbdy_data (const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     const FArrayBox& NC_MUB_fab,
                     const FArrayBox& NC_C1H_fab,
                     const FArrayBox& NC_C2H_fab,
                     const FArrayBox& NC_xvel_fab,
                     const FArrayBox& NC_yvel_fab,
                     const FArrayBox& NC_theta_fab,
                     const FArrayBox& NC_QVAPOR_fab);

void
init_state_from_wrfinput (int lev,
                          FArrayBox& cons_fab,
                          FArrayBox& x_vel_fab,
                          FArrayBox& y_vel_fab,
                          FArrayBox& z_vel_fab,
                          const Vector<FArrayBox>& NC_QVAPOR_fab,
                          const Vector<FArrayBox>& NC_QCLOUD_fab,
                          const Vector<FArrayBox>& NC_QRAIN_fab,
                          const Vector<FArrayBox>& NC_xvel_fab,
                          const Vector<FArrayBox>& NC_yvel_fab,
                          const Vector<FArrayBox>& NC_zvel_fab,
                          const Vector<FArrayBox>& NC_rho_fab,
                          const Vector<FArrayBox>& NC_theta_fab,
                          const int n_qstate,
                          const int RhoQr_comp);

void
init_msfs_from_wrfinput (int lev, FArrayBox& msfu_fab,
                         FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                         const Vector<FArrayBox>& NC_MSFU_fab,
                         const Vector<FArrayBox>& NC_MSFV_fab,
                         const Vector<FArrayBox>& NC_MSFM_fab);

void
verify_terrain_top_boundary (const Real& z_top,
                             const Vector<FArrayBox>& NC_PH_fab,
                             const Vector<FArrayBox>& NC_PHB_fab);

void
init_terrain_from_wrfinput (int lev, const Real& z_top,
                            const Box& domain, FArrayBox& z_phys,
                            const Vector<FArrayBox>& NC_PH_fab,
                            const Vector<FArrayBox>& NC_PHB_fab);

void
init_base_state_from_wrfinput (int lev,
                               const Box& gtbx,
                               const Box& domain,
                               Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate,
                               FArrayBox& cons_fab,
                               FArrayBox& p_hse,
                               FArrayBox& pi_hse,
                               FArrayBox& r_hse,
                               const Vector<FArrayBox>& NC_ALB_fab,
                               const Vector<FArrayBox>& NC_PB_fab,
                               const Vector<FArrayBox>& NC_P_fab);

/**
 * ERF function that initializes data from a WRF dataset
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_wrfinput (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_xvel_fab;   NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab;   NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_zvel_fab;   NC_zvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rho_fab;    NC_rho_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhop_fab;   NC_rhop_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_theta_fab;  NC_theta_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MUB_fab;    NC_MUB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFU_fab;   NC_MSFU_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFV_fab;   NC_MSFV_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFM_fab;   NC_MSFM_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_SST_fab;    NC_SST_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_LANDMSK_fab;NC_LANDMSK_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C1H_fab;    NC_C1H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C2H_fab;    NC_C2H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_RDNW_fab;   NC_RDNW_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PH_fab;     NC_PH_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_P_fab;      NC_P_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PHB_fab;    NC_PHB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ALB_fab;    NC_ALB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PB_fab;     NC_PB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QVAPOR_fab; NC_QVAPOR_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QCLOUD_fab; NC_QCLOUD_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QRAIN_fab ; NC_QRAIN_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_LAT_fab;    NC_LAT_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_LON_fab;    NC_LON_fab.resize(num_boxes_at_level[lev]);

    // Print() << "Building initial FABS from file " << nc_init_file[lev][idx] << std::endl;
    if (nc_init_file.empty()) {
        amrex::Error("NetCDF initialization file name must be provided via input");
    }

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_from_wrfinput(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                           NC_xvel_fab[idx]  , NC_yvel_fab[idx]   , NC_zvel_fab[idx] , NC_rho_fab[idx],
                           NC_rhop_fab[idx]  , NC_theta_fab[idx]  , NC_MUB_fab[idx]  ,
                           NC_MSFU_fab[idx]  , NC_MSFV_fab[idx]   , NC_MSFM_fab[idx] ,
                           NC_SST_fab[idx]   , NC_LANDMSK_fab[idx], NC_C1H_fab[idx]  ,
                           NC_C2H_fab[idx]   , NC_RDNW_fab[idx]   ,
                           NC_QVAPOR_fab[idx], NC_QCLOUD_fab[idx] , NC_QRAIN_fab[idx],
                           NC_PH_fab[idx]    , NC_P_fab[idx]      , NC_PHB_fab[idx]  ,
                           NC_ALB_fab[idx]   , NC_PB_fab[idx]     ,
                           NC_LAT_fab[idx]   , NC_LON_fab[idx]    ,
                           solverChoice.moisture_type, Latitude, Longitude, geom[lev]);
    }

    auto& lev_new = vars_new[lev];
    int n_qstate = micro->Get_Qstate_Size();
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // INITIAL DATA common for "ideal" as well as "real" simulation
    // Don't tile this since we are operating on full FABs in this routine
    for ( MFIter mfi(lev_new[Vars::cons], false); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];

        init_state_from_wrfinput(lev, cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                 NC_QVAPOR_fab, NC_QCLOUD_fab, NC_QRAIN_fab,
                                 NC_xvel_fab, NC_yvel_fab, NC_zvel_fab,
                                 NC_rho_fab, NC_theta_fab,
                                 n_qstate, solverChoice.RhoQr_comp);
    } // mf

    auto& ba = lev_new[Vars::cons].boxArray();
    auto& dm = lev_new[Vars::cons].DistributionMap();
    auto ngv = lev_new[Vars::cons].nGrowVect(); ngv[2] = 0;
    BoxList bl2d = ba.boxList();
    for (auto& b : bl2d) {
        b.setRange(2,0);
    }
    BoxArray ba2d(std::move(bl2d));
    int i_lo = geom[lev].Domain().smallEnd(0); int i_hi = geom[lev].Domain().bigEnd(0);
    int j_lo = geom[lev].Domain().smallEnd(1); int j_hi = geom[lev].Domain().bigEnd(1);

    lat_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,ngv);
    for ( MFIter mfi(*(lat_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();
        const Array4<      Real>& dst_arr = (lat_m[lev])->array(mfi);
        const Array4<const Real>& src_arr = NC_LAT_fab[0].const_array();
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);
            dst_arr(i,j,0) = src_arr(li,lj,0);
        });
    }
    lon_m[lev] = std::make_unique<MultiFab>(ba2d,dm,1,ngv);
    for ( MFIter mfi(*(lon_m[lev]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();
        const Array4<      Real>& dst_arr = (lon_m[lev])->array(mfi);
        const Array4<const Real>& src_arr = NC_LON_fab[0].const_array();
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);
            dst_arr(i,j,0) = src_arr(li,lj,0);
        });
    }

    for ( MFIter mfi(*(lmask_lev[lev][0]), TilingIfNotGPU()); mfi.isValid(); ++mfi ) {
        Box gtbx = mfi.growntilebox();
        const Array4<       int>& dst_arr = lmask_lev[lev][0]->array(mfi);
        const Array4<const Real>& src_arr = NC_LANDMSK_fab[0].const_array();
        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept
        {
            int li = amrex::min(amrex::max(i, i_lo), i_hi);
            int lj = amrex::min(amrex::max(j, j_lo), j_hi);
            dst_arr(i,j,0) = static_cast<int>(src_arr(li,lj,0));
        });
    }
    (lmask_lev[lev])[0]->FillBoundary(geom[lev].periodicity());

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // Map scale factors common for "ideal" as well as "real" simulation
    for ( MFIter mfi(*mapfac_u[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &msfu_fab = (*mapfac_u[lev])[mfi];
        FArrayBox &msfv_fab = (*mapfac_v[lev])[mfi];
        FArrayBox &msfm_fab = (*mapfac_m[lev])[mfi];

        init_msfs_from_wrfinput(lev, msfu_fab, msfv_fab, msfm_fab,
                                NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf

    const Box& domain = geom[lev].Domain();
    const Real& z_top = geom[lev].ProbHi(2);
    if (solverChoice.use_terrain)
    {
        if (ParallelDescriptor::IOProcessor()) {
            verify_terrain_top_boundary(z_top, NC_PH_fab, NC_PHB_fab);
        }

        std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
            init_terrain_from_wrfinput(lev, z_top, domain, z_phys_nd_fab, NC_PH_fab, NC_PHB_fab);
        } // mf

        make_J  (geom[lev],*z_phys_nd[lev],*  detJ_cc[lev]);
        make_areas(geom[lev],*z_phys_nd[lev],*ax[lev],*ay[lev],*az[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);

    } // use_terrain

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    IntVect ng = p_hse.nGrowVect();
    const Real l_rdOcp = solverChoice.rdOcp;

    if (init_type == InitType::Real) {
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox&   cons_fab = lev_new[Vars::cons][mfi];
            FArrayBox&  p_hse_fab = p_hse[mfi];
            FArrayBox& pi_hse_fab = pi_hse[mfi];
            FArrayBox&  r_hse_fab = r_hse[mfi];

            const Box gtbx = mfi.tilebox(IntVect(0), ng);
            init_base_state_from_wrfinput(lev, gtbx, domain, l_rdOcp, solverChoice.moisture_type, n_qstate,
                                          cons_fab, p_hse_fab, pi_hse_fab, r_hse_fab,
                                          NC_ALB_fab, NC_PB_fab, NC_P_fab);
        }

        // FillBoundary to populate the internal halo cells
         r_hse.FillBoundary(geom[lev].periodicity());
         p_hse.FillBoundary(geom[lev].periodicity());
        pi_hse.FillBoundary(geom[lev].periodicity());
    }

    if (init_type == InitType::Real && (lev == 0)) {
        if (nc_bdy_file.empty()) {
            amrex::Error("NetCDF boundary file name must be provided via input");
        }
        bdy_time_interval = read_from_wrfbdy(nc_bdy_file,geom[0].Domain(),
                                             bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                             real_width, start_bdy_time);

        Print() << "Read in boundary data with width "  << real_width << std::endl;
        Print() << "Running with specification width: " << real_set_width
                << " and relaxation width: " << real_width - real_set_width << std::endl;

        convert_wrfbdy_data(domain,bdy_data_xlo,
                            NC_MUB_fab[0] , NC_C1H_fab[0] , NC_C2H_fab[0],
                            NC_xvel_fab[0], NC_yvel_fab[0], NC_theta_fab[0], NC_QVAPOR_fab[0]);
        convert_wrfbdy_data(domain,bdy_data_xhi,
                            NC_MUB_fab[0] , NC_C1H_fab[0] , NC_C2H_fab[0],
                            NC_xvel_fab[0], NC_yvel_fab[0], NC_theta_fab[0], NC_QVAPOR_fab[0]);
        convert_wrfbdy_data(domain,bdy_data_ylo,
                            NC_MUB_fab[0] , NC_C1H_fab[0] , NC_C2H_fab[0],
                            NC_xvel_fab[0], NC_yvel_fab[0], NC_theta_fab[0], NC_QVAPOR_fab[0]);
        convert_wrfbdy_data(domain,bdy_data_yhi,
                            NC_MUB_fab[0] , NC_C1H_fab[0] , NC_C2H_fab[0] ,
                            NC_xvel_fab[0], NC_yvel_fab[0], NC_theta_fab[0], NC_QVAPOR_fab[0]);
    }

    // Start at the earliest time (read_from_wrfbdy)
    t_new[lev] = start_bdy_time;
    t_old[lev] = start_bdy_time - 1.e200;
}

/**
 * Helper function to initialize state and velocity data in a Fab from a WRF dataset.
 *
 * @param lev Integer specifying current level
 * @param state_fab FArrayBox object holding the state data we initialize
 * @param x_vel_fab FArrayBox object holding the x-velocity data we initialize
 * @param y_vel_fab FArrayBox object holding the y-velocity data we initialize
 * @param z_vel_fab FArrayBox object holding the z-velocity data we initialize
 * @param NC_xvel_fab Vector of FArrayBox objects with the WRF dataset specifying x-velocity
 * @param NC_yvel_fab Vector of FArrayBox objects with the WRF dataset specifying y-velocity
 * @param NC_zvel_fab Vector of FArrayBox objects with the WRF dataset specifying z-velocity
 * @param NC_rho_fab Vector of FArrayBox objects with the WRF dataset specifying density
 * @param NC_theta_fab Vector of FArrayBox objects with the WRF dataset specifying density*(potential temperature)
 */
void
init_state_from_wrfinput (int /*lev*/,
                          FArrayBox& state_fab,
                          FArrayBox& x_vel_fab,
                          FArrayBox& y_vel_fab,
                          FArrayBox& z_vel_fab,
                          const Vector<FArrayBox>& NC_QVAPOR_fab,
                          const Vector<FArrayBox>& NC_QCLOUD_fab,
                          const Vector<FArrayBox>& NC_QRAIN_fab,
                          const Vector<FArrayBox>& NC_xvel_fab,
                          const Vector<FArrayBox>& NC_yvel_fab,
                          const Vector<FArrayBox>& NC_zvel_fab,
                          const Vector<FArrayBox>& NC_rho_fab,
                          const Vector<FArrayBox>& NC_theta_fab,
                          const int n_qstate,
                          const int RhoQr_comp)
{
    int nboxes = NC_xvel_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        // This copies x-vel
        x_vel_fab.template copy<RunOn::Device>(NC_xvel_fab[idx]);

        // This copies y-vel
        y_vel_fab.template copy<RunOn::Device>(NC_yvel_fab[idx]);

        // This copies z-vel
        z_vel_fab.template copy<RunOn::Device>(NC_zvel_fab[idx]);

        // We first initialize all state_fab variables to zero
        state_fab.template setVal<RunOn::Device>(0.);

        // This copies the density
        state_fab.template copy<RunOn::Device>(NC_rho_fab[idx], 0, Rho_comp, 1);

        // This copies (rho*theta)
        state_fab.template copy<RunOn::Device>(NC_theta_fab[idx], 0, RhoTheta_comp, 1);
        state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]  , 0, RhoTheta_comp, 1);

        if (n_qstate >= 1) {
          state_fab.template copy<RunOn::Device>(NC_QVAPOR_fab[idx], 0, RhoQ1_comp, 1);
          state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]   , 0, RhoQ1_comp, 1);
        }

        if (n_qstate >= 2) {
          state_fab.template copy<RunOn::Device>(NC_QCLOUD_fab[idx], 0, RhoQ2_comp, 1);
          state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]   , 0, RhoQ2_comp, 1);
        }

        if (n_qstate >= 3) {
            state_fab.template copy<RunOn::Device>(NC_QRAIN_fab[idx], 0, RhoQr_comp, 1);
            state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]  , 0, RhoQr_comp, 1);
        }
    } // idx
}

/**
 * Helper function initializing velocity map factors from WRF data
 *
 * @param lev Integer specifying the current level
 * @param msfu_fab FArrayBox specifying the x-velocity map factors we initialize
 * @param msfv_fab FArrayBox specifying the y-velocity map factors we initialize
 * @param msfm_fab FArrayBox specifying the z-velocity map factors we initialize
 * @param NC_MSFU_fab Vector of FArrayBoxes holding WRF data specifying x-velocity map factors
 * @param NC_MSFV_fab Vector of FArrayBoxes holding WRF data specifying y-velocity map factors
 * @param NC_MSFM_fab Vector of FArrayBoxes holding WRF data specifying z-velocity map factors
 */
void
init_msfs_from_wrfinput (int /*lev*/, FArrayBox& msfu_fab,
                         FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                         const Vector<FArrayBox>& NC_MSFU_fab,
                         const Vector<FArrayBox>& NC_MSFV_fab,
                         const Vector<FArrayBox>& NC_MSFM_fab)
{
    int nboxes = NC_MSFU_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        // This copies mapfac_u
        msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab[idx]);

        // This copies mapfac_v
        msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab[idx]);

        // This copies mapfac_m
        msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab[idx]);
    } // idx
}

/**
 * Helper function to initialize hydrostatic base state data from WRF dataset
 *
 * @param lev Integer specifying current level
 * @param valid_bx Box specifying the index space we are to initialize
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param p_hse FArrayBox specifying the hydrostatic base state pressure we initialize
 * @param pi_hse FArrayBox specifying the hydrostatic base state Exner pressure we initialize
 * @param r_hse FArrayBox specifying the hydrostatic base state density we initialize
 * @param NC_ALB_fab Vector of FArrayBox objects containing WRF data specifying 1/density
 * @param NC_PB_fab Vector of FArrayBox objects containing WRF data specifying pressure
 */
void
init_base_state_from_wrfinput (int /*lev*/,
                               const Box& gtbx,
                               const Box& domain,
                               const Real l_rdOcp,
                               MoistureType moisture_type,
                               const int& n_qstate,
                               FArrayBox& cons_fab,
                               FArrayBox& p_hse,
                               FArrayBox& pi_hse,
                               FArrayBox& r_hse,
                               const Vector<FArrayBox>& NC_ALB_fab,
                               const Vector<FArrayBox>& NC_PB_fab,
                               const Vector<FArrayBox>& NC_P_fab)
{
    const auto& dom_lo = lbound(domain);
    const auto& dom_hi = ubound(domain);

    int nboxes = NC_ALB_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        const Array4<Real      >&   cons_arr = cons_fab.array();
        const Array4<Real      >&  p_hse_arr = p_hse.array();
        const Array4<Real      >& pi_hse_arr = pi_hse.array();
        const Array4<Real      >&  r_hse_arr = r_hse.array();
        //const Array4<Real const>&  alpha_arr = NC_ALB_fab[idx].const_array();
        const Array4<Real const>&  nc_pb_arr = NC_PB_fab[idx].const_array();
        const Array4<Real const>&   nc_p_arr = NC_P_fab[idx].const_array();

        ParallelFor(gtbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            // Base state needs ghost cells filled, protect FAB access
            int ii = std::max(i , dom_lo.x);
                ii = std::min(ii, dom_hi.x);
            int jj = std::max(j , dom_lo.y);
                jj = std::min(jj, dom_hi.y);
            int kk = std::max(k , dom_lo.z);
                kk = std::min(kk, dom_hi.z);

            // Base plus perturbational pressure
            Real Ptot = nc_pb_arr(ii,jj,kk) + nc_p_arr(ii,jj,kk);

            // Compute pressure from EOS
            Real Qv    = (moisture_type != MoistureType::None) ?
                         cons_arr(ii,jj,kk,RhoQ1_comp) / cons_arr(ii,jj,kk,Rho_comp) : 0.0;
            Real RT    = cons_arr(ii,jj,kk,RhoTheta_comp);
            Real P_eos = getPgivenRTh(RT, Qv);
            Real DelP  = std::fabs(Ptot - P_eos);
            AMREX_ASSERT_WITH_MESSAGE((DelP < 1.0), "Initial state is inconsistent with EOS!");

            // Compute rhse
            Real Rhse_Sum = cons_arr(ii,jj,kk,Rho_comp);
            for (int q_offset(0); q_offset<n_qstate; ++q_offset) Rhse_Sum += cons_arr(ii,jj,kk,RhoQ1_comp+q_offset);

            p_hse_arr(i,j,k)  = Ptot;
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
            r_hse_arr(i,j,k)  = Rhse_Sum;
        });
    } // idx
}

/**
 * Helper function for verifying the top boundary is valid.
 *
 * @param z_top Real user specified top boundary
 * @param NC_PH_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PH)
 * @param NC_PHB_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PHB)
 */
void
verify_terrain_top_boundary (const Real& z_top,
                             const Vector<FArrayBox>& NC_PH_fab,
                             const Vector<FArrayBox>& NC_PHB_fab)
{
    int nboxes = NC_PH_fab.size();
    for (int idx = 0; idx < nboxes; idx++) {
        Gpu::HostVector  <Real> MaxMax_h(2,-1.0e16);
        Gpu::DeviceVector<Real> MaxMax_d(2);
        Gpu::copy(Gpu::hostToDevice, MaxMax_h.begin(), MaxMax_h.end(), MaxMax_d.begin());

        Real* mm_d = MaxMax_d.data();

        Box Fab2dBox_hi (NC_PHB_fab[idx].box()); Fab2dBox_hi.makeSlab(2,Fab2dBox_hi.bigEnd(2));
        Box Fab2dBox_lo (NC_PHB_fab[idx].box()); Fab2dBox_lo.makeSlab(2,Fab2dBox_lo.bigEnd(2)-1);

        Box nodal_box = amrex::surroundingNodes(NC_PHB_fab[idx].box());
        int ilo = nodal_box.smallEnd()[0];
        int ihi = nodal_box.bigEnd()[0];
        int jlo = nodal_box.smallEnd()[1];
        int jhi = nodal_box.bigEnd()[1];

        auto const& phb = NC_PHB_fab[idx].const_array();
        auto const& ph  = NC_PH_fab[idx].const_array();

        ParallelFor(Fab2dBox_hi, Fab2dBox_lo,
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            Real z_calc = 0.25 * ( ph (ii,jj  ,k) + ph (ii-1,jj  ,k) +
                                   ph (ii,jj-1,k) + ph (ii-1,jj-1,k) +
                                   phb(ii,jj  ,k) + phb(ii-1,jj  ,k) +
                                   phb(ii,jj-1,k) + phb(ii-1,jj-1,k) ) / CONST_GRAV;
            amrex::Gpu::Atomic::Max(&(mm_d[0]),z_calc);
        },
        [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            Real z_calc = 0.25 * ( ph (ii,jj  ,k) + ph (ii-1,jj  ,k) +
                                   ph (ii,jj-1,k) + ph (ii-1,jj-1,k) +
                                   phb(ii,jj  ,k) + phb(ii-1,jj  ,k) +
                                   phb(ii,jj-1,k) + phb(ii-1,jj-1,k) ) / CONST_GRAV;
            amrex::Gpu::Atomic::Max(&(mm_d[1]),z_calc);
        });

        Gpu::copy(Gpu::deviceToHost, MaxMax_d.begin(), MaxMax_d.end(), MaxMax_h.begin());
        if ((z_top > MaxMax_h[0]) || (z_top < MaxMax_h[1])) {
            Print() << "Z problem extent " << z_top << " does not match NETCDF file min "
                    << MaxMax_h[1] << " and max " << MaxMax_h[0] << "!\n";
            Print() << "To run you must set the z component of prob_hi or prob_extent to lie between the netcdf bounds" << std::endl;
            Abort("Domain specification error");
        }
    }
}

/**
 * Helper function for initializing terrain coordinates from a WRF dataset.
 *
 * @param lev Integer specifying the current level
 * @param z_phys FArrayBox specifying the node-centered z coordinates of the terrain
 * @param NC_PH_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PH)
 * @param NC_PHB_fab Vector of FArrayBox objects storing WRF terrain coordinate data (PHB)
 */
void
init_terrain_from_wrfinput (int /*lev*/, const Real& z_top,
                            const Box& domain, FArrayBox& z_phys,
                            const Vector<FArrayBox>& NC_PH_fab,
                            const Vector<FArrayBox>& NC_PHB_fab)
{
    int nboxes = NC_PH_fab.size();
    for (int idx = 0; idx < nboxes; idx++) {
        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys.array();
        const Array4<Real const>& nc_phb_arr = NC_PHB_fab[idx].const_array();
        const Array4<Real const>& nc_ph_arr  = NC_PH_fab[idx].const_array();

        const Box& z_phys_box(z_phys.box());

        Box nodal_box = amrex::surroundingNodes(NC_PHB_fab[idx].box());

        int ilo = nodal_box.smallEnd()[0];
        int ihi = nodal_box.bigEnd()[0];
        int jlo = nodal_box.smallEnd()[1];
        int jhi = nodal_box.bigEnd()[1];
        int klo = nodal_box.smallEnd()[2];
        int khi = nodal_box.bigEnd()[2]-1;

        const auto domlo = amrex::lbound(domain);
        const auto domhi = amrex::ubound(domain);

        //
        // We must be careful not to read out of bounds of the WPS data
        //
        ParallelFor(z_phys_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            if (k < 0) {
                Real z_klo   =  0.25 * ( nc_ph_arr (ii,jj  ,klo  ) +  nc_ph_arr(ii-1,jj  ,klo  ) +
                                         nc_ph_arr (ii,jj-1,klo  ) + nc_ph_arr (ii-1,jj-1,klo) +
                                         nc_phb_arr(ii,jj  ,klo  ) + nc_phb_arr(ii-1,jj  ,klo  ) +
                                         nc_phb_arr(ii,jj-1,klo  ) + nc_phb_arr(ii-1,jj-1,klo) ) / CONST_GRAV;
                Real z_klop1 =  0.25 * ( nc_ph_arr (ii,jj  ,klo+1) +  nc_ph_arr(ii-1,jj  ,klo+1) +
                                         nc_ph_arr (ii,jj-1,klo+1) + nc_ph_arr (ii-1,jj-1,klo+1) +
                                         nc_phb_arr(ii,jj  ,klo+1) + nc_phb_arr(ii-1,jj  ,klo+1) +
                                         nc_phb_arr(ii,jj-1,klo+1) + nc_phb_arr(ii-1,jj-1,klo+1) ) / CONST_GRAV;
                z_arr(i, j, k) = 2.0 * z_klo - z_klop1;
            } else if (k > khi) {
                Real z_khim1 =  0.25 * ( nc_ph_arr (ii,jj  ,khi-1) + nc_ph_arr (ii-1,jj  ,khi-1) +
                                         nc_ph_arr (ii,jj-1,khi-1) + nc_ph_arr (ii-1,jj-1,khi-1) +
                                         nc_phb_arr(ii,jj  ,khi-1) + nc_phb_arr(ii-1,jj  ,khi-1) +
                                         nc_phb_arr(ii,jj-1,khi-1) + nc_phb_arr(ii-1,jj-1,khi-1) ) / CONST_GRAV;
                z_arr(i, j, k) = 2.0 * z_top - z_khim1;
            } else if (k == khi) {
                z_arr(i, j, k) = z_top;
            } else {
                z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj  ,k) + nc_ph_arr (ii-1,jj  ,k) +
                                          nc_ph_arr (ii,jj-1,k) + nc_ph_arr (ii-1,jj-1,k) +
                                          nc_phb_arr(ii,jj  ,k) + nc_phb_arr(ii-1,jj  ,k) +
                                          nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;
            } // k

            //
            // Fill values outside the domain on the fly (we will need these to make detJ in ghost cells)
            //
            if (i == domlo.x) z_arr(i-1,j,k) = z_arr(i,j,k);
            if (i == domhi.x) z_arr(i+1,j,k) = z_arr(i,j,k);
            if (j == domlo.y) z_arr(i,j-1,k) = z_arr(i,j,k);
            if (j == domhi.y) z_arr(i,j+1,k) = z_arr(i,j,k);

            if (i == domlo.x && j == domlo.y) z_arr(i-1,j-1,k) = z_arr(i,j,k);
            if (i == domlo.x && j == domhi.y) z_arr(i-1,j+1,k) = z_arr(i,j,k);
            if (i == domhi.x && j == domlo.y) z_arr(i+1,j-1,k) = z_arr(i,j,k);
            if (i == domhi.x && j == domhi.y) z_arr(i+1,j+1,k) = z_arr(i,j,k);
        });
    } // idx
}
#endif // ERF_USE_NETCDF
