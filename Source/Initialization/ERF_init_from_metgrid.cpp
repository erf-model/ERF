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
                   FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                   FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                   FArrayBox& NC_pres_fab, FArrayBox& NC_hgt_fab,
                   FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                   FArrayBox& NC_msfm_fab);
void
interpolate_column (int i, int j, int src_comp, int dest_comp,
                    const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                    const Array4<Real const>&  new_z, const Array4<Real>&  new_data);

void
init_terrain_from_metgrid (int lev, FArrayBox& z_phys_nd_fab,
                           const Vector<FArrayBox>& NC_hgt_fab);

void
init_state_from_metgrid (int lev, FArrayBox& state_fab,
                         FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                         FArrayBox& z_vel_fab, FArrayBox& z_phys_nd_fab,
                         const Vector<FArrayBox>& NC_hgt_fab,
                         const Vector<FArrayBox>& NC_xvel_fab,
                         const Vector<FArrayBox>& NC_yvel_fab,
                         const Vector<FArrayBox>& NC_zvel_fab,
                         const Vector<FArrayBox>& NC_rho_fab,
                         const Vector<FArrayBox>& NC_rhotheta_fab);
void
init_msfs_from_metgrid (int lev, FArrayBox& msfu_fab,
                        FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                        const Vector<FArrayBox>& NC_MSFU_fab,
                        const Vector<FArrayBox>& NC_MSFV_fab,
                        const Vector<FArrayBox>& NC_MSFM_fab);
void
init_base_state_from_metgrid (int lev, const Box& valid_bx, Real l_rdOcp,
                              FArrayBox& p_hse, FArrayBox& pi_hse, FArrayBox& r_hse,
                              const Vector<FArrayBox>& NC_ALB_fab,
                              const Vector<FArrayBox>& NC_PB_fab);

#ifdef ERF_USE_NETCDF
/**
 * Initializes ERF data using metgrid data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_metgrid (int lev)
{
    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_xvel_fab ; NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab ; NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_temp_fab ; NC_temp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhum_fab ; NC_rhum_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pres_fab ; NC_pres_fab.resize(num_boxes_at_level[lev]);

    Vector<FArrayBox> NC_hgt_fab ; NC_hgt_fab.resize(num_boxes_at_level[lev]);

    Vector<FArrayBox> NC_MSFU_fab ; NC_MSFU_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFV_fab ; NC_MSFV_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFM_fab ; NC_MSFM_fab.resize(num_boxes_at_level[lev]);

    int nboxes = num_boxes_at_level[lev];

    if (nc_init_file.empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    if (nc_init_file[lev].empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    for (int idx = 0; idx < nboxes; idx++)
    {
        read_from_metgrid(lev, boxes_at_level[lev][idx], nc_init_file[lev][idx],
                          NC_xvel_fab[idx], NC_yvel_fab[idx],
                          NC_temp_fab[idx], NC_rhum_fab[idx],
                          NC_pres_fab[idx], NC_hgt_fab[idx],
                          NC_MSFU_fab[idx], NC_MSFV_fab[idx], NC_MSFM_fab[idx] );
    }

    auto& lev_new = vars_new[lev];

    std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];

    AMREX_ALWAYS_ASSERT(solverChoice.use_terrain);

    z_phys->setVal(0.);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // This defines only the z(i,j,0) values given the FAB filled from the NetCDF input
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        init_terrain_from_metgrid(lev, z_phys_nd_fab, NC_hgt_fab);
    } // mf

    // This defines all the z(i,j,k) values given z(i,j,0) from above.
    init_terrain_grid(geom[lev], *z_phys, zlevels_stag);

    // This makes the Jacobian.
    make_J  (geom[lev],*z_phys,  *detJ_cc[lev]);

    // This defines z at w-cell faces.
    make_zcc(geom[lev],*z_phys,*z_phys_cc[lev]);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];

        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        init_state_from_metgrid(lev, cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                z_phys_nd_fab,
                                NC_hgt_fab, NC_xvel_fab, NC_yvel_fab, NC_temp_fab,
                                NC_rhum_fab, NC_pres_fab);
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

        init_msfs_from_metgrid(lev, msfu_fab, msfv_fab, msfm_fab,
                               NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf

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

            const Box& bx = mfi.validbox();
            //init_base_state_from_metgrid(lev, bx, l_rdOcp, p_hse_fab, pi_hse_fab, r_hse_fab,
            //                             NC_ALB_fab, NC_PB_fab);
        }
    }
    exit(0);
}


/**
 * Helper function to initialize terrain nodal z coordinates for a Fab
 * given metgrid data.
 *
 * @param lev Integer specifying the current level
 * @param z_phys_nd_fab FArrayBox (Fab) holding the nodal z coordinates for terrain data we want to fill
 * @param NC_hgt_fab Vector of FArrayBox objects holding height data read from NetCDF files for metgrid data
 */
void
init_terrain_from_metgrid (int lev, FArrayBox& z_phys_nd_fab,
                           const Vector<FArrayBox>& NC_hgt_fab)
{
   int nboxes = NC_hgt_fab.size();

   // NOTE NOTE NOTE -- this routine currently only fills the k=0 value
   // TODO: we need to fill the rest of the values from the pressure (?) variable...

   for (int idx = 0; idx < nboxes; idx++)
   {
#ifndef AMREX_USE_GPU
        amrex::Print() << " SIZE OF HGT FAB " << NC_hgt_fab[idx].box() << std::endl;
        amrex::Print() << " SIZE OF ZP FAB "  << z_phys_nd_fab.box() << std::endl;
#endif

        // This copies from NC_zphys on z-faces to z_phys_nd on nodes
        const Array4<Real      >&      z_arr = z_phys_nd_fab.array();
        const Array4<Real const>& nc_hgt_arr = NC_hgt_fab[idx].const_array();

        const Box z_hgt_box = NC_hgt_fab[idx].box();

        int ilo = z_hgt_box.smallEnd()[0];
        int ihi = z_hgt_box.bigEnd()[0];
        int jlo = z_hgt_box.smallEnd()[1];
        int jhi = z_hgt_box.bigEnd()[1];

        Box z_phys_box = z_phys_nd_fab.box();
        Box   from_box = surroundingNodes(NC_hgt_fab[idx].box()); from_box.growHi(2,-1);
        Box bx = z_phys_box & from_box;

#ifndef AMREX_USE_GPU
        amrex::Print() << "FROM BOX " << from_box << std::endl;
        amrex::Print() << "BX " << bx << std::endl;
#endif

        //
        // We must be careful not to read out of bounds of the WPS data
        //
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
        {
            int ii = std::max(std::min(i,ihi-1),ilo+1);
            int jj = std::max(std::min(j,jhi-1),jlo+1);
            z_arr(i,j,k) =  0.25 * ( nc_hgt_arr (ii,jj  ,k) + nc_hgt_arr(ii-1,jj  ,k) +
                                     nc_hgt_arr (ii,jj-1,k) + nc_hgt_arr(ii-1,jj-1,k) );
        });
    } // idx
}

/**
 * Helper function to initialize state and velocity data
 * read from metgrid data.
 *
 * @param lev Integer specifying the current level
 * @param state_fab FArrayBox holding the state data to initialize
 * @param x_vel_fab FArrayBox holding the x-velocity data to initialize
 * @param y_vel_fab FArrayBox holding the y-velocity data to initialize
 * @param z_vel_fab FArrayBox holding the z-velocity data to initialize
 * @param z_phys_nd_fab FArrayBox holding nodal z coordinate data for terrain
 * @param NC_hgt_fab Vector of FArrayBox obects holding metgrid data for height
 * @param NC_xvel_fab Vector of FArrayBox obects holding metgrid data for x-velocity
 * @param NC_yvel_fab Vector of FArrayBox obects holding metgrid data for y-velocity
 * @param NC_zvel_fab Vector of FArrayBox obects holding metgrid data for z-velocity
 * @param NC_rho_fab Vector of FArrayBox obects holding metgrid data for density
 * @param NC_rhotheta_fab Vector of FArrayBox obects holding metgrid data for (density * potential temperature)
 */
void
init_state_from_metgrid (int lev, FArrayBox& state_fab,
                         FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                         FArrayBox& z_vel_fab, FArrayBox& z_phys_nd_fab,
                         const Vector<FArrayBox>& NC_hgt_fab,
                         const Vector<FArrayBox>& NC_xvel_fab,
                         const Vector<FArrayBox>& NC_yvel_fab,
                         const Vector<FArrayBox>& NC_zvel_fab,
                         const Vector<FArrayBox>& NC_rho_fab,
                         const Vector<FArrayBox>& NC_rhotheta_fab)
{
    int nboxes = NC_hgt_fab.size();

#ifndef AMREX_USE_GPU
    amrex::Print() << " U FROM NC " << NC_xvel_fab[0].box() << std::endl;
    amrex::Print() << " U INTO FAB " << x_vel_fab.box() << std::endl;
    exit(0);
#endif
    for (int idx = 0; idx < nboxes; idx++)
    {
        // ********************************************************
        // U
        // ********************************************************
        {
        Box bx2d = NC_xvel_fab[idx].box() & x_vel_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = NC_xvel_fab[idx].const_array();
        auto const orig_z    = NC_hgt_fab[idx].const_array();

        auto       new_data  = x_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column(i,j,0,0,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // V
        // ********************************************************
        {
        Box bx2d = NC_yvel_fab[idx].box() & y_vel_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = NC_yvel_fab[idx].const_array();
        auto const orig_z    = NC_hgt_fab[idx].const_array();

        auto       new_data  = y_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column(i,j,0,0,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // W
        // ********************************************************
        z_vel_fab.template setVal<RunOn::Device>(0.);

        // ********************************************************
        // rho
        // ********************************************************
        {
        Box bx2d = NC_rho_fab[idx].box() & state_fab.box();
        bx2d.setRange(2,0);

        auto const orig_data = NC_rho_fab[idx].const_array();
        auto const orig_z    = NC_hgt_fab[idx].const_array();
        auto        new_data  = state_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column(i,j,0,Rho_comp,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // rho_theta
        // ********************************************************
        {
        Box bx2d = NC_rhotheta_fab[idx].box() & state_fab.box();
        bx2d.setRange(2,0);

        auto const orig_data = NC_rhotheta_fab[idx].const_array();
        auto const orig_z    = NC_hgt_fab[idx].const_array();
        auto       new_data  = state_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column(i,j,0,RhoTheta_comp,orig_z,orig_data,new_z,new_data);
        });
        }
    } // idx
}

/**
 * Helper function to initialize map factors from metgrid data
 *
 * @param lev Integer specifying current level
 * @param msfu_fab FArrayBox specifying x-velocity map factors
 * @param msfv_fab FArrayBox specifying y-velocity map factors
 * @param msfm_fab FArrayBox specifying z-velocity map factors
 * @param NC_MSFU_fab Vector of FArrayBox objects holding metgrid data for x-velocity map factors
 * @param NC_MSFV_fab Vector of FArrayBox objects holding metgrid data for y-velocity map factors
 * @param NC_MSFM_fab Vector of FArrayBox objects holding metgrid data for z-velocity map factors
 */
void
init_msfs_from_metgrid (int lev, FArrayBox& msfu_fab,
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
 * Helper function for initializing hydrostatic base state data from metgrid data
 *
 * @param lev Integer specifying the current level
 * @param valid_bx Box specifying the index space we are initializing
 * @param l_rdOcp Real constant specifying Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param p_hse FArrayBox holding the hydrostatic base state pressure we are initializing
 * @param pi_hse FArrayBox holding the hydrostatic base Exner pressure we are initializing
 * @param r_hse FArrayBox holding the hydrostatic base state density we are initializing
 * @param NC_ALB_fab Vector of FArrayBox objects holding metgrid data specifying 1/density
 * @param NC_PB_fab Vector of FArrayBox objects holding metgrid data specifying pressure
 */
void
init_base_state_from_metgrid (int lev, const Box& valid_bx, const Real l_rdOcp,
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
        const Array4<Real      >& pi_hse_arr = pi_hse.array();
        const Array4<Real      >&  r_hse_arr =  r_hse.array();
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
