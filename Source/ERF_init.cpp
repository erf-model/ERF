/**
 * \file ERF_init.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>

using namespace amrex;

#ifdef ERF_USE_NETCDF

Real
read_from_wrfbdy(std::string nc_bdy_file, const Box& domain,
                 Vector<Vector<FArrayBox>>& bdy_data_xlo,
                 Vector<Vector<FArrayBox>>& bdy_data_xhi,
                 Vector<Vector<FArrayBox>>& bdy_data_ylo,
                 Vector<Vector<FArrayBox>>& bdy_data_yhi);

void
convert_wrfbdy_data(int which, const Box& domain,
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
ERF::init_from_wrfinput(int lev)
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

    if (nc_init_file.size() == 0)
        amrex::Error("NetCDF initialization file name must be provided via input");

    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        read_from_wrfinput(lev,idx,NC_xvel_fab,NC_yvel_fab,NC_zvel_fab,NC_rho_fab,
                           NC_rhop_fab,NC_rhoth_fab,NC_MUB_fab,
                           NC_MSFU_fab,NC_MSFV_fab,NC_MSFM_fab,
                           NC_SST_fab,
                           NC_C1H_fab,NC_C2H_fab,NC_RDNW_fab,
                           NC_PH_fab,NC_PHB_fab,NC_ALB_fab,NC_PB_fab);
    }

    auto& lev_new = vars_new[lev];

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    // INITIAL DATA common for "ideal" as well as "real" simulation
    for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
    {
        // Define fabs for holding the initial data
        FArrayBox &cons_fab = lev_new[Vars::cons][mfi];
        FArrayBox &xvel_fab = lev_new[Vars::xvel][mfi];
        FArrayBox &yvel_fab = lev_new[Vars::yvel][mfi];
        FArrayBox &zvel_fab = lev_new[Vars::zvel][mfi];

        init_state_from_wrfinput(lev, cons_fab, xvel_fab, yvel_fab, zvel_fab,
                                 NC_xvel_fab, NC_yvel_fab, NC_zvel_fab,
                                 NC_rho_fab, NC_rhoth_fab);
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

    if (solverChoice.use_terrain)
    {
        std::unique_ptr<MultiFab>& z_phys = z_phys_nd[lev];
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox& z_phys_nd_fab = (*z_phys)[mfi];
            init_terrain_from_wrfinput(lev, z_phys_nd_fab, NC_PH_fab, NC_PHB_fab);
        } // mf

        make_J  (geom[lev],*z_phys_nd[lev],*  detJ_cc[lev]);
        make_zcc(geom[lev],*z_phys_nd[lev],*z_phys_cc[lev]);

    } // use_terrain

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    if (init_type == "real") {
        for ( MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi )
        {
            FArrayBox&  p_hse_fab = p_hse[mfi];
            FArrayBox& pi_hse_fab = pi_hse[mfi];
            FArrayBox&  r_hse_fab = r_hse[mfi];

            const Box& bx = mfi.validbox();
            init_base_state_from_wrfinput(lev, bx, p_hse_fab, pi_hse_fab, r_hse_fab,
                                         NC_ALB_fab, NC_PB_fab);
        }
    }

    if (init_type == "real" && (lev == 0)) {
        if (nc_bdy_file.empty())
            amrex::Error("NetCDF boundary file name must be provided via input");
        bdy_time_interval = read_from_wrfbdy(nc_bdy_file,geom[0].Domain(),bdy_data_xlo,bdy_data_xhi,bdy_data_ylo,bdy_data_yhi);

        const Box& domain = geom[lev].Domain();

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
}

void
ERF::init_state_from_wrfinput(int lev, FArrayBox& state_fab,
                              FArrayBox& x_vel_fab, FArrayBox& y_vel_fab,
                              FArrayBox& z_vel_fab,
                              const Vector<FArrayBox>& NC_xvel_fab,
                              const Vector<FArrayBox>& NC_yvel_fab,
                              const Vector<FArrayBox>& NC_zvel_fab,
                              const Vector<FArrayBox>& NC_rho_fab,
                              const Vector<FArrayBox>& NC_rhotheta_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
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
    } // idx
}

void
ERF::init_msfs_from_wrfinput(int lev, FArrayBox& msfu_fab,
                             FArrayBox& msfv_fab, FArrayBox& msfm_fab,
                             const Vector<FArrayBox>& NC_MSFU_fab,
                             const Vector<FArrayBox>& NC_MSFV_fab,
                             const Vector<FArrayBox>& NC_MSFM_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
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

void
ERF::init_base_state_from_wrfinput(int lev, const Box& bx, FArrayBox& p_hse, FArrayBox& pi_hse,
                                   FArrayBox& r_hse,
                                   const Vector<FArrayBox>& NC_ALB_fab,
                                   const Vector<FArrayBox>& NC_PB_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
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
        const Real rdOcp = solverChoice.rdOcp;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            p_hse_arr(i,j,k)  = nc_pb_arr(i,j,k);
            pi_hse_arr(i,j,k) = getExnergivenP(p_hse_arr(i,j,k), rdOcp);
            r_hse_arr(i,j,k)  = 1.0 / alpha_arr(i,j,k);

        });
    } // idx
}

void
ERF::init_terrain_from_wrfinput(int lev, FArrayBox& z_phys,
                                const Vector<FArrayBox>& NC_PH_fab,
                                const Vector<FArrayBox>& NC_PHB_fab)
{
    for (int idx = 0; idx < num_boxes_at_level[lev]; idx++)
    {
        //
        // FArrayBox to FArrayBox copy does "copy on intersection"
        // This only works here because we have broadcast the FArrayBox of data from the netcdf file to all ranks
        //

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
        int khi = nodal_box.bigEnd()[2];

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
            } // k
        });
    } // idx
}
#endif // ERF_USE_NETCDF

void
ERF::init_from_input_sounding(int lev)
{
    // We only want to read the file once -- here we fill one FArrayBox (per variable) that spans the domain
    if (lev == 0) {
        if (input_sounding_file.empty())
            amrex::Error("input_sounding file name must be provided via input");
        Real ztop = geom[0].ProbHi(AMREX_SPACEDIM-1);
        input_sounding_data.read_from_file(input_sounding_file, ztop);

        if (init_sounding_ideal) input_sounding_data.calc_rho_p(ztop);
    }

    auto& lev_new = vars_new[lev];

    // update if init_sounding_ideal == true
    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.tilebox();
        const auto &cons_arr = lev_new[Vars::cons].array(mfi);
        const auto &xvel_arr = lev_new[Vars::xvel].array(mfi);
        const auto &yvel_arr = lev_new[Vars::yvel].array(mfi);
        const auto &zvel_arr = lev_new[Vars::zvel].array(mfi);
        Array4<Real> r_hse_arr = r_hse.array(mfi);
        Array4<Real> p_hse_arr = p_hse.array(mfi);
        Array4<Real> pi_hse_arr = pi_hse.array(mfi);

        if (init_sounding_ideal)
        {
            init_bx_scalars_from_input_sounding_hse(bx, cons_arr,
                                                    r_hse_arr, p_hse_arr, pi_hse_arr,
                                                    geom[lev].data(), input_sounding_data);
        }
        else
        {
            init_bx_scalars_from_input_sounding(bx, cons_arr,
                                                geom[lev].data(), input_sounding_data);
        }

        init_bx_velocities_from_input_sounding(bx, xvel_arr, yvel_arr, zvel_arr,
                                               geom[lev].data(), input_sounding_data);
    } //mfi
}

void
ERF::init_bx_scalars_from_input_sounding(
        const amrex::Box &bx,
        amrex::Array4<amrex::Real> const &state,
        amrex::GeometryData const &geomdata,
        InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d.dataPtr();
#ifdef ERF_USE_MOISTURE
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d.dataPtr();
#endif
    const int   inp_sound_size  = inputSoundingData.size();

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Geometry
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        amrex::Real rho_0 = 1.0;

        // Set the density
        state(i, j, k, Rho_comp) = rho_0;

        // Initial Rho0*Theta0
        state(i, j, k, RhoTheta_comp) = rho_0 * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state(i, j, k, RhoScalar_comp) = 0;

#ifdef ERF_USE_MOISTURE
        // total nonprecipitating water (Qt) == water vapor (Qv), i.e., there
        // is no cloud water or cloud ice
        state(i, j, k, RhoQt_comp) = rho_0 * interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);
#endif
    });
}

void
ERF::init_bx_scalars_from_input_sounding_hse(
        const amrex::Box &bx,
        amrex::Array4<amrex::Real> const &state,
        amrex::Array4<amrex::Real> const &r_hse,
        amrex::Array4<amrex::Real> const &p_hse,
        amrex::Array4<amrex::Real> const &pi_hse,
        amrex::GeometryData const &geomdata,
        InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* rho_inp_sound   = inputSoundingData.rho_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d.dataPtr();
#ifdef ERF_USE_MOISTURE
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d.dataPtr();
#endif
    const int   inp_sound_size  = inputSoundingData.size();

    amrex::Real l_gravity = solverChoice.gravity;

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    const Real rdOcp = solverChoice.rdOcp;

    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Geometry
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];
        int ktop = bx.bigEnd(2);

        Real rho_k, rhoTh_k;
        rho_k = interpolate_1d(z_inp_sound, rho_inp_sound, z, inp_sound_size);

        // Set the density
        state(i, j, k, Rho_comp) = rho_k;

        // Initial Rho0*Theta0
        rhoTh_k = rho_k * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);
        state(i, j, k, RhoTheta_comp) = rhoTh_k;

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state(i, j, k, RhoScalar_comp) = 0;

        // Update hse quantities with values calculated from InputSoundingData.calc_rho_p()
        r_hse(i, j, k) = rho_k;
        p_hse(i, j, k) = getPgivenRTh(rhoTh_k);
        pi_hse(i, j, k) = getExnergivenRTh(rhoTh_k, rdOcp);

        // Boundary treatment
        if (k==0)
        {
            // set the ghost cell with dz and rho at boundary
            amrex::Real rho_surf =
                interpolate_1d(z_inp_sound, rho_inp_sound, 0.0, inp_sound_size);
            amrex::Real rhoTh_surf =
                rho_surf * interpolate_1d(z_inp_sound, theta_inp_sound, 0.0, inp_sound_size);
             p_hse(i, j, k-1) = getPgivenRTh(rhoTh_surf) + dx[2]/2 * rho_surf * l_gravity;
            pi_hse(i, j, k-1) = getExnergivenP(p_hse(i, j, k-1), rdOcp);
        }
        else if (k==ktop)
        {
            // set the ghost cell with dz and rho at boundary
            amrex::Real rho_top =
                interpolate_1d(z_inp_sound, rho_inp_sound, z+dx[2]/2, inp_sound_size);
            amrex::Real rhoTh_top =
                rho_top * interpolate_1d(z_inp_sound, theta_inp_sound, z+dx[2]/2, inp_sound_size);
             p_hse(i, j, k+1) = getPgivenRTh(rhoTh_top) - dx[2]/2 * rho_top * l_gravity;
            pi_hse(i, j, k+1) = getExnergivenP(p_hse(i, j, k+1), rdOcp);
        }

#ifdef ERF_USE_MOISTURE
        // total nonprecipitating water (Qt) == water vapor (Qv), i.e., there
        // is no cloud water or cloud ice
        state(i, j, k, RhoQt_comp) = rho_k * interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);
#endif
    });
}

void
ERF::init_bx_velocities_from_input_sounding(
        const amrex::Box &bx,
        amrex::Array4<amrex::Real> const &x_vel,
        amrex::Array4<amrex::Real> const &y_vel,
        amrex::Array4<amrex::Real> const &z_vel,
        amrex::GeometryData const &geomdata,
        InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* U_inp_sound     = inputSoundingData.U_inp_sound_d.dataPtr();
    const Real* V_inp_sound     = inputSoundingData.V_inp_sound_d.dataPtr();
    const int   inp_sound_size  = inputSoundingData.size();

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    // Construct a box that is on x-faces
    const amrex::Box& xbx = amrex::surroundingNodes(gbx,0);
    // Construct a box that is on y-faces
    const amrex::Box& ybx = amrex::surroundingNodes(gbx,1);
    // Construct a box that is on z-faces
    const amrex::Box& zbx = amrex::surroundingNodes(gbx,2);

    // Set the x,y,z-velocities
    amrex::ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of x-faces
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel(i, j, k) = interpolate_1d(z_inp_sound, U_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of y-faces
        const amrex::Real* prob_lo = geomdata.ProbLo();
        const amrex::Real* dx = geomdata.CellSize();
        const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel(i, j, k) = interpolate_1d(z_inp_sound, V_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of z-faces
        // Set the z-velocity
        z_vel(i, j, k) = 0.0;
    });
}

void
ERF::init_custom(int lev)
{
    auto& lev_new = vars_new[lev];
#ifdef ERF_USE_MOISTURE
    auto& qv_new  = qv[lev];
    auto& qc_new  = qc[lev];
    auto& qi_new  = qi[lev];
#endif
    MultiFab r_hse(base_state[lev], make_alias, 0, 1); // r_0 is first  component
    MultiFab p_hse(base_state[lev], make_alias, 1, 1); // p_0 is second component

    MultiFab cons_pert(lev_new[Vars::cons].boxArray(), lev_new[Vars::cons].DistributionMap(),
                       lev_new[Vars::cons].nComp()   , lev_new[Vars::cons].nGrow());
    MultiFab xvel_pert(lev_new[Vars::xvel].boxArray(), lev_new[Vars::xvel].DistributionMap(), 1, lev_new[Vars::xvel].nGrowVect());
    MultiFab yvel_pert(lev_new[Vars::yvel].boxArray(), lev_new[Vars::yvel].DistributionMap(), 1, lev_new[Vars::yvel].nGrowVect());
    MultiFab zvel_pert(lev_new[Vars::zvel].boxArray(), lev_new[Vars::zvel].DistributionMap(), 1, lev_new[Vars::zvel].nGrowVect());

    // Default all perturbations to zero
    cons_pert.setVal(0.);
    xvel_pert.setVal(0.);
    yvel_pert.setVal(0.);
    zvel_pert.setVal(0.);

#ifdef ERF_USE_MOISTURE
    MultiFab qv_pert(qv[lev].boxArray(), qv[lev].DistributionMap(), 1, qv[lev].nGrow());
    MultiFab qc_pert(qc[lev].boxArray(), qc[lev].DistributionMap(), 1, qc[lev].nGrow());
    MultiFab qi_pert(qi[lev].boxArray(), qi[lev].DistributionMap(), 1, qi[lev].nGrow());
    qv_pert.setVal(0.);
    qc_pert.setVal(0.);
    qi_pert.setVal(0.);
#endif

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box &bx = mfi.tilebox();

        const auto &cons_pert_arr = cons_pert.array(mfi);
        const auto &xvel_pert_arr = xvel_pert.array(mfi);
        const auto &yvel_pert_arr = yvel_pert.array(mfi);
        const auto &zvel_pert_arr = zvel_pert.array(mfi);

        Array4<Real const> z_nd_arr = (solverChoice.use_terrain) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_cc_arr = (solverChoice.use_terrain) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};

        Array4<Real const> mf_m     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_u     = mapfac_m[lev]->array(mfi);
        Array4<Real const> mf_v     = mapfac_m[lev]->array(mfi);

        Array4<Real> r_hse_arr = r_hse.array(mfi);
        Array4<Real> p_hse_arr = p_hse.array(mfi);

#ifdef ERF_USE_MOISTURE
        const auto &qv_pert_arr = qv_pert.array(mfi);
        const auto &qc_pert_arr = qc_pert.array(mfi);
        const auto &qi_pert_arr = qi_pert.array(mfi);
#endif
        init_custom_prob(bx, cons_pert_arr, xvel_pert_arr, yvel_pert_arr, zvel_pert_arr,
                         r_hse_arr, p_hse_arr, z_nd_arr, z_cc_arr,
#ifdef ERF_USE_MOISTURE
                         qv_pert_arr, qc_pert_arr, qi_pert_arr,
#endif
                         geom[lev].data(), mf_m, mf_u, mf_v,
                         solverChoice);
    } //mfi

    // Add problem-specific perturbation to background flow
    MultiFab::Add(lev_new[Vars::cons], cons_pert, Rho_comp,      Rho_comp,      1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoTheta_comp, RhoTheta_comp, 1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoScalar_comp,RhoScalar_comp,1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQKE_comp,   RhoQKE_comp,   1, cons_pert.nGrow());
#ifdef ERF_USE_MOISTURE
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQt_comp,    RhoQt_comp,    1, cons_pert.nGrow());
    MultiFab::Add(lev_new[Vars::cons], cons_pert, RhoQp_comp,    RhoQp_comp,    1, cons_pert.nGrow());
    MultiFab::Add(             qv_new,   qv_pert, 0,             0,             1,   qv_pert.nGrow());
    MultiFab::Add(             qc_new,   qc_pert, 0,             0,             1,   qc_pert.nGrow());
    MultiFab::Add(             qi_new,   qi_pert, 0,             0,             1,   qi_pert.nGrow());
#endif
    MultiFab::Add(lev_new[Vars::xvel], xvel_pert, 0,             0,             1, xvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::yvel], yvel_pert, 0,             0,             1, yvel_pert.nGrowVect());
    MultiFab::Add(lev_new[Vars::zvel], zvel_pert, 0,             0,             1, zvel_pert.nGrowVect());
}
