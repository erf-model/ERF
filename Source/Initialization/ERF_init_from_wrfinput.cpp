/**
 * \file ERF_init_from_wrfinput.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>
#include <DataStruct.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

void
read_from_wrfinput (int lev, const Box& domain, const std::string& fname,
                    FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                    FArrayBox& NC_zvel_fab, FArrayBox& NC_rho_fab,
                    FArrayBox& NC_rhop_fab, FArrayBox& NC_rhotheta_fab,
                    FArrayBox& NC_MUB_fab ,
                    FArrayBox& NC_MSFU_fab, FArrayBox& NC_MSFV_fab,
                    FArrayBox& NC_MSFM_fab, FArrayBox& NC_SST_fab,
                    FArrayBox& NC_C1H_fab , FArrayBox& NC_C2H_fab,
                    FArrayBox& NC_RDNW_fab,
                    FArrayBox& NC_QVAPOR_fab,
                    FArrayBox& NC_QCLOUD_fab,
                    FArrayBox& NC_QRAIN_fab,
                    FArrayBox& NC_PH_fab  , FArrayBox& NC_PHB_fab,
                    FArrayBox& NC_ALB_fab , FArrayBox& NC_PB_fab,
                    MoistureType moisture_type);

Real
read_from_wrfbdy (const std::string& nc_bdy_file, const Box& domain,
                  Vector<Vector<FArrayBox>>& bdy_data_xlo,
                  Vector<Vector<FArrayBox>>& bdy_data_xhi,
                  Vector<Vector<FArrayBox>>& bdy_data_ylo,
                  Vector<Vector<FArrayBox>>& bdy_data_yhi,
                  int& width, amrex::Real& start_bdy_time);

void
convert_wrfbdy_data (int which, const Box& domain,
                     Vector<Vector<FArrayBox>>& bdy_data,
                     const FArrayBox& NC_MUB_fab,
                     const FArrayBox& NC_MSFU_fab,
                     const FArrayBox& NC_MSFV_fab,
                     const FArrayBox& NC_MSFM_fab,
                     const FArrayBox& NC_PH_fab,
                     const FArrayBox& NC_PHB_fab,
                     const FArrayBox& NC_C1H_fab,
                     const FArrayBox& NC_C2H_fab,
                     const FArrayBox& NC_RDNW_fab,
                     const FArrayBox& NC_xvel_fab,
                     const FArrayBox& NC_yvel_fab,
                     const FArrayBox& NC_rho_fab,
                     const FArrayBox& NC_rhoth_fab);

void
init_state_from_wrfinput (int lev, FArrayBox& state_fab,
                          FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                          FArrayBox& z_vel_fab,
                          const Vector<FArrayBox>& NC_QVAPOR_fab,
                          const Vector<FArrayBox>& NC_QCLOUD_fab,
                          const Vector<FArrayBox>& NC_QRAIN_fab,
                          const Vector<FArrayBox>& NC_xvel_fab,
                          const Vector<FArrayBox>& NC_yvel_fab,
                          const Vector<FArrayBox>& NC_zvel_fab,
                          const Vector<FArrayBox>& NC_rho_fab,
                          const Vector<FArrayBox>& NC_rhotheta_fab,
                          MoistureType moisture_type);

void
init_msfs_from_wrfinput (int lev, FArrayBox& msfu_fab,
                         FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                         const Vector<FArrayBox>& NC_MSFU_fab,
                         const Vector<FArrayBox>& NC_MSFV_fab,
                         const Vector<FArrayBox>& NC_MSFM_fab);
void
init_terrain_from_wrfinput (int lev, const Box& domain, FArrayBox& z_phys,
                            const Vector<FArrayBox>& NC_PH_fab,
                            const Vector<FArrayBox>& NC_PHB_fab);

void
init_base_state_from_wrfinput (int lev, const Box& bx, Real l_rdOcp,
                               FArrayBox& p_hse, FArrayBox& pi_hse,
                               FArrayBox& r_hse,
                               const Vector<FArrayBox>& NC_ALB_fab,
                               const Vector<FArrayBox>& NC_PB_fab);

/**
 * ERF function that initializes data from a WRF dataset
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_wrfinput (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_zvel_fab ; NC_zvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rho_fab  ; NC_rho_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhop_fab ; NC_rhop_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhoth_fab; NC_rhoth_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MUB_fab  ; NC_MUB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFU_fab ; NC_MSFU_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFV_fab ; NC_MSFV_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFM_fab ; NC_MSFM_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_SST_fab  ; NC_SST_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C1H_fab  ; NC_C1H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_C2H_fab  ; NC_C2H_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_RDNW_fab ; NC_RDNW_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PH_fab   ; NC_PH_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PHB_fab  ; NC_PHB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ALB_fab  ; NC_ALB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_PB_fab   ; NC_PB_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QVAPOR_fab; NC_QVAPOR_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QCLOUD_fab; NC_QCLOUD_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_QRAIN_fab ; NC_QRAIN_fab.resize(num_boxes_at_level[lev]);

    // amrex::Print() << "Building initial FABS from file " << nc_init_file[lev][idx] << std::endl;
    if (nc_init_file.empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_from_wrfinput(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                           NC_xvel_fab[idx], NC_yvel_fab[idx],  NC_zvel_fab[idx], NC_rho_fab[idx],
                           NC_rhop_fab[idx], NC_rhoth_fab[idx], NC_MUB_fab[idx],
                           NC_MSFU_fab[idx], NC_MSFV_fab[idx],  NC_MSFM_fab[idx],
                           NC_SST_fab[idx],  NC_C1H_fab[idx],   NC_C2H_fab[idx],  NC_RDNW_fab[idx],
                           NC_QVAPOR_fab[idx], NC_QCLOUD_fab[idx], NC_QRAIN_fab[idx],
                           NC_PH_fab[idx],NC_PHB_fab[idx],NC_ALB_fab[idx],NC_PB_fab[idx],
                           solverChoice.moisture_type);
    }

    auto& lev_new = vars_new[lev];

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
                                 NC_rho_fab, NC_rhoth_fab, solverChoice.moisture_type);
    } // mf

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

    if (solverChoice.use_terrain)
    {
        std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
            init_terrain_from_wrfinput(lev, domain, z_phys_nd_fab, NC_PH_fab, NC_PHB_fab);
        } // mf

        make_J  (geom[lev],*z_phys_nd[lev],*  detJ_cc[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);

    } // use_terrain

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    const Real l_rdOcp = solverChoice.rdOcp;

    if (init_type == "real") {
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox&  p_hse_fab = p_hse[mfi];
            FArrayBox& pi_hse_fab = pi_hse[mfi];
            FArrayBox&  r_hse_fab = r_hse[mfi];

            const Box valid_bx = mfi.validbox();
            init_base_state_from_wrfinput(lev, valid_bx, l_rdOcp,
                                          p_hse_fab, pi_hse_fab, r_hse_fab,
                                          NC_ALB_fab, NC_PB_fab);
        }
    }

    if (init_type == "real" && (lev == 0)) {
        if (nc_bdy_file.empty())
            amrex::Error("NetCDF boundary file name must be provided via input");
        bdy_time_interval = read_from_wrfbdy(nc_bdy_file,geom[0].Domain(),
                                             bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi,
                                             wrfbdy_width, start_bdy_time);

        if (wrfbdy_width-1 <= wrfbdy_set_width) wrfbdy_set_width = wrfbdy_width;
        amrex::Print() << "Read in boundary data with width "  << wrfbdy_width << std::endl;
        amrex::Print() << "Running with specification width: " << wrfbdy_set_width
                       << " and relaxation width: " << wrfbdy_width - wrfbdy_set_width << std::endl;

        // NOTE: Last WRF BDY cell is a ghost cell for Laplacian relaxation.
        //       Without relaxation zones, we must augment this value by 1.
        if (wrfbdy_width == wrfbdy_set_width) wrfbdy_width += 1;

        convert_wrfbdy_data(0,domain,bdy_data_xlo,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(1,domain,bdy_data_xhi,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(2,domain,bdy_data_ylo,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
        convert_wrfbdy_data(3,domain,bdy_data_yhi,
                            NC_MUB_fab[0], NC_MSFU_fab[0], NC_MSFV_fab[0], NC_MSFM_fab[0],
                            NC_PH_fab[0] , NC_PHB_fab[0],
                            NC_C1H_fab[0], NC_C2H_fab[0], NC_RDNW_fab[0],
                            NC_xvel_fab[0],NC_yvel_fab[0],NC_rho_fab[0],NC_rhoth_fab[0]);
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
 * @param NC_rhotheta_fab Vector of FArrayBox objects with the WRF dataset specifying density*(potential temperature)
 */
void
init_state_from_wrfinput (int lev,
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
                          const Vector<FArrayBox>& NC_rhotheta_fab,
                          MoistureType moisture_type)
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
        state_fab.template copy<RunOn::Device>(NC_rhotheta_fab[idx], 0, RhoTheta_comp, 1);

        if (moisture_type != MoistureType::None)
        {
            state_fab.template copy<RunOn::Device>(NC_QVAPOR_fab[idx], 0, RhoQ1_comp, 1);
            state_fab.template plus<RunOn::Device>(NC_QCLOUD_fab[idx], 0, RhoQ1_comp, 1);
            state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]   , 0, RhoQ1_comp, 1);

            state_fab.template copy<RunOn::Device>(NC_QRAIN_fab[idx], 0, RhoQ2_comp, 1);
            state_fab.template mult<RunOn::Device>(NC_rho_fab[idx]  , 0, RhoQ2_comp, 1);
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
init_msfs_from_wrfinput (int lev, FArrayBox& msfu_fab,
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
init_base_state_from_wrfinput (int lev, const Box& valid_bx, const Real l_rdOcp,
                               FArrayBox& p_hse, FArrayBox& pi_hse, FArrayBox& r_hse,
                               const Vector<FArrayBox>& NC_ALB_fab,
                               const Vector<FArrayBox>& NC_PB_fab)
{
    int nboxes = NC_ALB_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //
        const Array4<Real      >&  p_hse_arr =  p_hse.array();
        const Array4<Real      >& pi_hse_arr = pi_hse.array(); const Array4<Real      >&  r_hse_arr =  r_hse.array();
        const Array4<Real const>& alpha_arr = NC_ALB_fab[idx].const_array();
        const Array4<Real const>& nc_pb_arr = NC_PB_fab[idx].const_array();

        amrex::ParallelFor(valid_bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            p_hse_arr(i,j,k)  = nc_pb_arr(i,j,k);
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
            r_hse_arr(i,j,k)  = 1.0 / alpha_arr(i,j,k);

        });
    } // idx
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
init_terrain_from_wrfinput (int lev, const Box& domain, FArrayBox& z_phys,
                            const Vector<FArrayBox>& NC_PH_fab,
                            const Vector<FArrayBox>& NC_PHB_fab)
{
    int nboxes = NC_PH_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
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
        amrex::ParallelFor(z_phys_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
                Real z_khi   =  0.25 * ( nc_ph_arr (ii,jj  ,khi  ) + nc_ph_arr (ii-1,jj  ,khi  ) +
                                         nc_ph_arr (ii,jj-1,khi  ) + nc_ph_arr (ii-1,jj-1,khi) +
                                         nc_phb_arr(ii,jj  ,khi  ) + nc_phb_arr(ii-1,jj  ,khi  ) +
                                         nc_phb_arr(ii,jj-1,khi  ) + nc_phb_arr(ii-1,jj-1,khi) ) / CONST_GRAV;
                Real z_khim1 =  0.25 * ( nc_ph_arr (ii,jj  ,khi-1) + nc_ph_arr (ii-1,jj  ,khi-1) +
                                         nc_ph_arr (ii,jj-1,khi-1) + nc_ph_arr (ii-1,jj-1,khi-1) +
                                         nc_phb_arr(ii,jj  ,khi-1) + nc_phb_arr(ii-1,jj  ,khi-1) +
                                         nc_phb_arr(ii,jj-1,khi-1) + nc_phb_arr(ii-1,jj-1,khi-1) ) / CONST_GRAV;
                z_arr(i, j, k) = 2.0 * z_khi - z_khim1;
              } else {
                z_arr(i, j, k) = 0.25 * ( nc_ph_arr (ii,jj  ,k) + nc_ph_arr (ii-1,jj  ,k) +
                                          nc_ph_arr (ii,jj-1,k) + nc_ph_arr (ii-1,jj-1,k) +
                                          nc_phb_arr(ii,jj  ,k) + nc_phb_arr(ii-1,jj  ,k) +
                                          nc_phb_arr(ii,jj-1,k) + nc_phb_arr(ii-1,jj-1,k) ) / CONST_GRAV;

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

            } // k
        });
    } // idx
}
#endif // ERF_USE_NETCDF
