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
read_from_metgrid(int lev, const std::string& fname,
                  FArrayBox& NC_xvel_fab, FArrayBox& NC_yvel_fab,
                  FArrayBox& NC_temp_fab, FArrayBox& NC_rhum_fab,
                  FArrayBox& NC_pres_fab, FArrayBox& NC_ght_fab,
                  FArrayBox& NC_hgt_fab,  FArrayBox& NC_psfc_fab, 
                  FArrayBox& NC_msfu_fab, FArrayBox& NC_msfv_fab,
                  FArrayBox& NC_msfm_fab);

void
interpolate_column(int i, int j, int src_comp, int dest_comp,
                   const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                   const Array4<Real const>&  new_z, const Array4<Real>&  new_data);

void
interpolate_column_metgrid(int i, int j, char stag, int src_comp, int dest_comp,
                           const Array4<Real const>& orig_z, const Array4<Real const>& orig_data,
                           const Array4<Real const>&  new_z, const Array4<Real>&  new_data);

void
init_terrain_from_metgrid(FArrayBox& z_phys_nd_fab,
                          const Vector<FArrayBox>& NC_hgt_fab);

void
init_state_from_metgrid(FArrayBox& state_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox& z_vel_fab, FArrayBox& z_phys_nd_fab,
                        const Vector<FArrayBox>& NC_hgt_fab,
                        const Vector<FArrayBox>& NC_ght_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_zvel_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_rhum_fab,
                              Vector<FArrayBox>& mxrat_fab);

void
init_msfs_from_metgrid(FArrayBox& msfu_fab,
                       FArrayBox& msfv_fab,
                       FArrayBox& msfm_fab,
                       const Vector<FArrayBox>& NC_MSFU_fab,
                       const Vector<FArrayBox>& NC_MSFV_fab,
                       const Vector<FArrayBox>& NC_MSFM_fab);

void
init_base_state_from_metgrid(const Real l_rdOcp,
                             FArrayBox& state,
                             FArrayBox& r_hse_fab,
                             FArrayBox& p_hse_fab,
                             FArrayBox& pi_hse_fab,
                             FArrayBox& z_phys_nd_fab,
                             const Vector<FArrayBox>& NC_ght_fab,
                             const Vector<FArrayBox>& NC_psfc_fab);

void
rh_to_mxrat(int i, int j, int k,
            const Array4<Real const>& rhum,
            const Array4<Real const>& temp,
            const Array4<Real const>& pres,
            const Array4<Real>& mxrat);

void
calc_rho_p(int i, int j,
           const Array4<Real const>& psfc,
           const Array4<Real const>& new_z,
           const Array4<Real>& new_data,
           amrex::Array4<amrex::Real> const &p_hse_arr,
           amrex::Array4<amrex::Real> const &r_hse_arr);

#ifdef ERF_USE_NETCDF
/**
 * Initializes ERF data using metgrid data supplied by an external NetCDF file.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_metgrid(int lev)
{
#if defined(ERF_USE_MOISTURE)
#elif defined(ERF_USE_WARM_NO_PRECIP)
#else
    amrex::Error("Init with met_em requires either ERF_USE_MOISTURE or ERF_USE_WARM_NO_PRECIP");
#endif

    // *** FArrayBox's at this level for holding the INITIAL data
    Vector<FArrayBox> NC_xvel_fab;  NC_xvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_yvel_fab;  NC_yvel_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_temp_fab;  NC_temp_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_rhum_fab;  NC_rhum_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_pres_fab;  NC_pres_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_ght_fab;   NC_ght_fab.resize(num_boxes_at_level[lev]);

    Vector<FArrayBox> NC_hgt_fab;   NC_hgt_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_psfc_fab;  NC_psfc_fab.resize(num_boxes_at_level[lev]);

    Vector<FArrayBox> NC_MSFU_fab;  NC_MSFU_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFV_fab;  NC_MSFV_fab.resize(num_boxes_at_level[lev]);
    Vector<FArrayBox> NC_MSFM_fab;  NC_MSFM_fab.resize(num_boxes_at_level[lev]);

    int nboxes = num_boxes_at_level[lev];

    if (nc_init_file.empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    if (nc_init_file[lev].empty())
        amrex::Error("NetCDF initialization file name must be provided via input");

    for (int idx = 0; idx < nboxes; idx++)
    {
        read_from_metgrid(lev,nc_init_file[lev][idx],NC_xvel_fab[idx],NC_yvel_fab[idx],
                          NC_temp_fab[idx], NC_rhum_fab[idx], NC_pres_fab[idx],
                          NC_ght_fab[idx], NC_hgt_fab[idx], NC_psfc_fab[idx],
                          NC_MSFU_fab[idx], NC_MSFV_fab[idx], NC_MSFM_fab[idx] );
    }

    // Set up a FAB for mixing ratio, since our input data only has relative humidity.
    // Is there an easier way to set up an empty variable that we can calculate later?
    Vector<FArrayBox> mxrat_fab; mxrat_fab.resize(num_boxes_at_level[lev]);
    for (int idx = 0; idx < nboxes; idx++)
    {
        Box my_box = NC_rhum_fab[idx].box();
#ifdef AMREX_USE_GPU
        mxrat_fab[idx].resize(my_box, 1, The_Pinned_Arena());
#else
        mxrat_fab[idx].resize(my_box, 1);
#endif
    }

    auto& lev_new = vars_new[lev];

    std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];

    AMREX_ALWAYS_ASSERT(solverChoice.use_terrain);

    z_phys->setVal(0.);

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // This defines only the z(i,j,0) values given the FAB filled from the NetCDF input
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
        init_terrain_from_metgrid(z_phys_nd_fab, NC_hgt_fab);
    } // mf

    // This defines all the z(i,j,k) values given z(i,j,0) from above.
    init_terrain_grid(geom[lev], *z_phys);

#ifndef AMREX_USE_GPU
    const amrex::Box& domain = geom[lev].Domain();
    auto dx = geom[lev].CellSizeArray();
    for (int i = 0; i < dx.size(); i++)
    {
        amrex::Print() << " geom.CellSizeArray()[" << i << "]: " << dx[i] << std::endl;
    }
    amrex::Print() << " geom[" << lev << "].Domain.smallEnd(0): " << domain.smallEnd(0) << "    geom[" << lev << "].Domain.bigEnd(0): " << domain.bigEnd(0) << std::endl;
    amrex::Print() << " geom[" << lev << "].Domain.smallEnd(1): " << domain.smallEnd(1) << "    geom[" << lev << "].Domain.bigEnd(1): " << domain.bigEnd(1) << std::endl;
    amrex::Print() << " geom[" << lev << "].Domain.smallEnd(2): " << domain.smallEnd(2) << "    geom[" << lev << "].Domain.bigEnd(2): " << domain.bigEnd(2) << std::endl;
    amrex::Print() << " geom[" << lev << "].Domain.length(0): " << domain.length(0) << std::endl;
    amrex::Print() << " geom[" << lev << "].Domain.length(1): " << domain.length(1) << std::endl;
    amrex::Print() << " geom[" << lev << "].Domain.length(2): " << domain.length(2) << std::endl;
#endif

    // TODO: Verify that DX and DY from met_em file matches that in geom (from ERF inputs file).

    // This makes the Jacobian.
    make_J(geom[lev],*z_phys,  *detJ_cc[lev]);

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
        init_state_from_metgrid(cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                z_phys_nd_fab,
                                NC_hgt_fab, NC_ght_fab, NC_xvel_fab,
                                NC_yvel_fab, NC_temp_fab, NC_rhum_fab,
                                NC_pres_fab, mxrat_fab);
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

        // Zap the unnecessary ghost cells added by default when the MultiFab was created.
//        Box bxu = msfu_fab.box();
//        bxu.grow(2,-3);
//        msfu_fab.resize(bxu);
//        Box bxv = msfv_fab.box();
//        bxv.grow(2,-3);
//        msfv_fab.resize(bxv);
//        Box bxm = msfm_fab.box();
//        bxm.grow(2,-3);
//        msfm_fab.resize(bxm);
        init_msfs_from_metgrid(msfu_fab, msfv_fab, msfm_fab,
                               NC_MSFU_fab, NC_MSFV_fab, NC_MSFM_fab);
    } // mf

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    const Real l_rdOcp = solverChoice.rdOcp;

    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        FArrayBox&     p_hse_fab = p_hse[mfi];
        FArrayBox&    pi_hse_fab = pi_hse[mfi];
        FArrayBox&     r_hse_fab = r_hse[mfi];
        FArrayBox&      cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];

        init_base_state_from_metgrid(l_rdOcp,
                                     cons_fab, r_hse_fab, p_hse_fab, pi_hse_fab,
                                     z_phys_nd_fab, NC_ght_fab, NC_psfc_fab);
    }
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
init_terrain_from_metgrid(FArrayBox& z_phys_nd_fab,
                          const Vector<FArrayBox>& NC_hgt_fab)
{
   int nboxes = NC_hgt_fab.size();

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
        Box from_box = surroundingNodes(NC_hgt_fab[idx].box());
        from_box.growHi(2,-1);

        Box bx = z_phys_box & from_box;
        Box bxu = bx; bxu.growHi(0,1);
        Box bxv = bx; bxv.growHi(1,1);

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
init_state_from_metgrid(FArrayBox& state_fab,
                        FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                        FArrayBox& z_vel_fab, FArrayBox& z_phys_nd_fab,
                        const Vector<FArrayBox>& NC_hgt_fab,
                        const Vector<FArrayBox>& NC_ght_fab,
                        const Vector<FArrayBox>& NC_xvel_fab,
                        const Vector<FArrayBox>& NC_yvel_fab,
                        const Vector<FArrayBox>& NC_temp_fab,
                        const Vector<FArrayBox>& NC_rhum_fab,
                        const Vector<FArrayBox>& NC_pres_fab,
                              Vector<FArrayBox>& mxrat_fab)
{
    int nboxes = NC_hgt_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        // ********************************************************
        // U
        // ********************************************************
        {
        Box bx2d = NC_xvel_fab[idx].box() & x_vel_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = NC_xvel_fab[idx].const_array();
        auto const orig_z    = NC_ght_fab[idx].const_array();

        auto       new_data  = x_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column_metgrid(i,j,'X',0,0,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // V
        // ********************************************************
        {
        Box bx2d = NC_yvel_fab[idx].box() & y_vel_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = NC_yvel_fab[idx].const_array();
        auto const orig_z    = NC_ght_fab[idx].const_array();

        auto       new_data  = y_vel_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column_metgrid(i,j,'Y',0,0,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // W
        // ********************************************************
        z_vel_fab.template setVal<RunOn::Device>(0.0);

        // ********************************************************
        // Initialize all state_fab variables to zero
        // ********************************************************
        state_fab.template setVal<RunOn::Device>(0.0);

        // ********************************************************
        // theta
        // ********************************************************
        {
        Box bx2d = NC_temp_fab[idx].box() & state_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = NC_temp_fab[idx].const_array();
        auto const orig_z    = NC_ght_fab[idx].const_array();

        auto       new_data  = state_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();

        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column_metgrid(i,j,'M',0,RhoTheta_comp,orig_z,orig_data,new_z,new_data);
        });
        }

        // ********************************************************
        // specific humidity / relative humidity / mixing ratio
        // ********************************************************
        // In the future we will need to check what input data we have for moisture
        // and then, if necessary, compute mixing ratio. For now, we will focus on
        // the case where we have relative humidity.
        //
        { // calculate vapor mixing ratio from relative humidity.
        Box bx = NC_rhum_fab[idx].box() & state_fab.box();
        auto const rhum  = NC_rhum_fab[idx].const_array();
        auto const temp  = NC_temp_fab[idx].const_array();
        auto const pres  = NC_pres_fab[idx].const_array();
        auto       mxrat = mxrat_fab[idx].array();

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            rh_to_mxrat(i,j,k,rhum,temp,pres,mxrat);
        });
        }

        { // vertical interpolation of vapor mixing ratio.
        Box bx2d = NC_rhum_fab[idx].box() & state_fab.box();
        bx2d.setRange(2,0);
        auto const orig_data = mxrat_fab[idx].const_array();
        auto const orig_z    = NC_ght_fab[idx].const_array();
        auto       new_data  = state_fab.array();
        auto const new_z     = z_phys_nd_fab.const_array();
#if defined(ERF_USE_MOISTURE)
        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column_metgrid(i,j,'M',0,RhoQt_comp,orig_z,orig_data,new_z,new_data);
        });
#elif defined(ERF_USE_WARM_NO_PRECIP)
        ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            interpolate_column_metgrid(i,j,'M',0,RhoQv_comp,orig_z,orig_data,new_z,new_data);
        });
#endif

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
rh_to_mxrat(int i, int j, int k, 
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

void
calc_rho_p(int i, int j,
           const Array4<Real const>& psfc,
           const Array4<Real const>& new_z,
           const Array4<Real>& new_data,
           amrex::Array4<amrex::Real> const &p_hse_arr,
           amrex::Array4<amrex::Real> const &r_hse_arr)
{

    const int maxiter = 10;
    int kmax = amrex::ubound(Box(new_data)).z-1;

    // Set up vectors for some components that are calculated as we integrate up and down the column.
    // Is this a clean-enough / proper way of handling that?
    std::vector<Real> theta_v(kmax);
    std::vector<Real> pm_integ(kmax);
    std::vector<Real> rhom_integ(kmax);
    std::vector<Real> pd_integ(kmax);

    // TODO: check if we have pressure at the surface. If not, estimate it!
//    Real t_0 = 290.0; // WRF's model_config_rec%base_temp
//    Real a = 50.0; // WRF's model_config_rec%base_lapse
//    Real p_surf = p_0*exp(-t_0/a+std::pow((std::pow(t_0/a, 2)-2*CONST_GRAV*hgt(i,j)/(a*R_d)), 0.5);

    // calculate virtual potential temperature at the surface.
#if defined(ERF_USE_MOISTURE)
    theta_v[0] = new_data(i,j,0,RhoTheta_comp)*(1.0+(R_v/R_d-1.0)*new_data(i,j,0,RhoQt_comp));
#elif defined(ERF_USE_WARM_NO_PRECIP)
    theta_v[0] = new_data(i,j,0,RhoTheta_comp)*(1.0+(R_v/R_d-1.0)*new_data(i,j,0,RhoQv_comp));
#endif

    // calculate or use moist pressure at the surface.
    pm_integ[0] = psfc(i,j,0);

    // calculate moist density at the surface.
    rhom_integ[0] = 1.0/(R_d/p_0*theta_v[0]*std::pow(pm_integ[0]/p_0, -iGamma));

    // integrate from the surface to the top boundary.
    for (int k=1; k < kmax; ++k) {
        Real dz = new_z(i,j,k)-new_z(i,j,k-1);
#if defined(ERF_USE_MOISTURE)
        Real qvf = 1.0+(R_v/R_d-1.0)*new_data(i,j,k,RhoQt_comp);
        theta_v[k] = new_data(i,j,k,RhoTheta_comp)*(1.0+(R_v/R_d-1.0)*new_data(i,j,k,RhoQt_comp));
#elif defined(ERF_USE_WARM_NO_PRECIP)
        Real qvf = 1.0+(R_v/R_d-1.0)*new_data(i,j,k,RhoQv_comp);
        theta_v[k] = new_data(i,j,k,RhoTheta_comp)*(1.0+(R_v/R_d-1.0)*new_data(i,j,k,RhoQv_comp));
#else
        Real qvf = 0.0;
#endif
        rhom_integ[k] = rhom_integ[k-1]; // an initial guess.
        for (int it=0; it < maxiter; ++it) {
            pm_integ[k] = pm_integ[k-1]-0.5*dz*(rhom_integ[k]+rhom_integ[k-1])*CONST_GRAV;
            if (pm_integ[k] <= 0.0) pm_integ[k] = 0.0;
            rhom_integ[k] = 1.0/(R_d/p_0*theta_v[k]*qvf*std::pow(pm_integ[k]/p_0, -iGamma));
        }
    }

    // integrate from the top back down to get dry pressure and density.
    pd_integ[kmax-1] = pm_integ[kmax-1];
    new_data(i,j,kmax-1,Rho_comp) = 1.0/(R_d/p_0*new_data(i,j,kmax-1,RhoTheta_comp)*std::pow(pd_integ[kmax-1]/p_0, -iGamma));
    for (int k=kmax-2; k >= 0; --k) {
        Real dz = new_z(i,j,k+1)-new_z(i,j,k);
        new_data(i,j,k,Rho_comp) = new_data(i,j,k+1,Rho_comp); // an initial guess.
        for (int it=0; it < maxiter; ++it) {
            pd_integ[k] = pd_integ[k+1]+0.5*dz*(new_data(i,j,k,Rho_comp)+new_data(i,j,k+1,Rho_comp))*CONST_GRAV;
            if (pd_integ[k] <= 0.0) pd_integ[k] = 0.0;
            new_data(i,j,k,Rho_comp) = 1.0/(R_d/p_0*new_data(i,j,k,RhoTheta_comp)*std::pow(pd_integ[k]/p_0, -iGamma));
        }
        p_hse_arr(i,j,k)  = pd_integ[k];
        r_hse_arr(i,j,k)  = new_data(i,j,k,Rho_comp);
    }
}

void
init_msfs_from_metgrid(FArrayBox& msfu_fab, FArrayBox& msfv_fab, FArrayBox& msfm_fab,
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

        // This copies mapfac_m
        msfm_fab.template copy<RunOn::Device>(NC_MSFM_fab[idx]);

        // This copies mapfac_u
        msfu_fab.template copy<RunOn::Device>(NC_MSFU_fab[idx]);

        // This copies mapfac_v
        msfv_fab.template copy<RunOn::Device>(NC_MSFV_fab[idx]);

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
init_base_state_from_metgrid(const Real l_rdOcp,
                             FArrayBox& state_fab,
                             FArrayBox& r_hse_fab,
                             FArrayBox& p_hse_fab,
                             FArrayBox& pi_hse_fab,
                             FArrayBox& z_phys_nd_fab,
                             const Vector<FArrayBox>& NC_ght_fab,
                             const Vector<FArrayBox>& NC_psfc_fab)
{
    int nboxes = NC_psfc_fab.size();
    for (int idx = 0; idx < nboxes; idx++)
    {
        const Array4<Real>& r_hse_arr = r_hse_fab.array();
        const Array4<Real>& p_hse_arr = p_hse_fab.array();
        const Array4<Real>& pi_hse_arr = pi_hse_fab.array();

        // ********************************************************
        // calculate dry density and dry pressure
        // ********************************************************
        // calculate density and dry pressure on the new grid.
        {
        Box bx2d = state_fab.box();
        bx2d.setRange(2,0);
        bx2d.grow(0,-3); // TODO: is there a cleaner way to avoid the boundary / ghost cells?
        bx2d.grow(1,-3);
        auto const orig_psfc = NC_psfc_fab[idx].const_array();
        auto       new_data = state_fab.array();
        auto const new_z = z_phys_nd_fab.const_array();

        amrex::ParallelFor(bx2d, [=] AMREX_GPU_DEVICE (int i, int j, int) noexcept {
            calc_rho_p(i,j,orig_psfc,new_z,new_data,p_hse_arr,r_hse_arr);
        });

        Box bx3d = state_fab.box();
        bx3d.grow(0,-3); // TODO: is there a cleaner way to avoid the boundary / ghost cells?
        bx3d.grow(1,-3);
        bx3d.grow(2,-3);
        amrex::ParallelFor(bx3d, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            new_data(i,j,k,Rho_comp) = r_hse_arr(i,j,k);
            new_data(i,j,k,RhoScalar_comp) = 0.0;
            // RhoTheta and RhoQt are currently just Theta and Qt. Multiply by Rho.
            Real RhoTheta = r_hse_arr(i,j,k)*new_data(i,j,k,RhoTheta_comp);
            new_data(i,j,k,RhoTheta_comp) = RhoTheta;
#if defined(ERF_USE_MOISTURE)
            Real RhoQ = r_hse_arr(i,j,k)*new_data(i,j,k,RhoQt_comp);
            new_data(i,j,k,RhoQt_comp) = RhoQ;
#elif defined(ERF_USE_WARM_NO_PRECIP)
            Real RhoQ = r_hse_arr(i,j,k)*new_data(i,j,k,RhoQv_comp);
            new_data(i,j,k,RhoQv_comp) = RhoQ;
#endif
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), l_rdOcp);
//            pi_hse_arr(i,j,k) = getExnergivenRTh(RhoTheta, l_rdOcp);
        });
        }
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
interpolate_column_metgrid(int i, int j, char stag, int src_comp, int dest_comp,
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

AMREX_GPU_DEVICE
void
interpolate_column(int i, int j, int src_comp, int dest_comp,
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
