/**
 * \file ERF_init_from_input_sounding.cpp
 */

#include <ERF.H>
#include <EOS.H>
#include <ERF_Constants.H>
#include <Utils.H>
#include <prob_common.H>

using namespace amrex;

void
init_bx_scalars_from_input_sounding (const amrex::Box &bx,
                                     amrex::Array4<amrex::Real> const &state,
                                     amrex::GeometryData const &geomdata,
                                     amrex::Array4<const amrex::Real> const &z_cc_arr,
                                     const bool& l_moist,
                                     InputSoundingData const &inputSoundingData);
void
init_bx_scalars_from_input_sounding_hse (const amrex::Box &bx,
                                         amrex::Array4<amrex::Real> const &state,
                                         amrex::Array4<amrex::Real> const &r_hse_arr,
                                         amrex::Array4<amrex::Real> const &p_hse_arr,
                                         amrex::Array4<amrex::Real> const &pi_hse_arr,
                                         amrex::GeometryData const &geomdata,
                                         amrex::Array4<const amrex::Real> const &z_cc_arr,
                                         const amrex::Real& l_gravity,
                                         const amrex::Real& l_rdOcp,
                                         const bool& l_moist,
                                         InputSoundingData const &inputSoundingData);

void
init_bx_velocities_from_input_sounding (const amrex::Box &bx,
                                        amrex::Array4<amrex::Real> const &x_vel,
                                        amrex::Array4<amrex::Real> const &y_vel,
                                        amrex::Array4<amrex::Real> const &z_vel,
                                        amrex::GeometryData const &geomdata,
                                        amrex::Array4<const amrex::Real> const &z_nd_arr,
                                        InputSoundingData const &inputSoundingData);

/**
 * High level wrapper for initializing scalar and velocity
 * level data from input sounding data.
 *
 * @param lev Integer specifying the current level
 */
void
ERF::init_from_input_sounding (int lev)
{
    // We only want to read the file once -- here we fill one FArrayBox (per variable) that spans the domain
    if (lev == 0) {
        if (input_sounding_file.empty())
            amrex::Error("input_sounding file name must be provided via input");

        // this will interpolate the input profiles to the nominal height levels
        // (ranging from 0 to the domain top)
        input_sounding_data.read_from_file(input_sounding_file, geom[lev], zlevels_stag);

        // this will calculate the hydrostatically balanced density and pressure
        // profiles following WRF ideal.exe
        if (init_sounding_ideal) input_sounding_data.calc_rho_p();
    }

    auto& lev_new = vars_new[lev];

    // update if init_sounding_ideal == true
    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    const Real l_gravity = solverChoice.gravity;
    const Real l_rdOcp   = solverChoice.rdOcp;
    const bool l_moist   = (solverChoice.moisture_type != MoistureType::None);

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

        Array4<Real const> z_cc_arr = (solverChoice.use_terrain) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_nd_arr = (solverChoice.use_terrain) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};

        if (init_sounding_ideal)
        {
            // HSE will be initialized here, interpolated from values previously
            // calculated by calc_rho_p()
            init_bx_scalars_from_input_sounding_hse(
                bx, cons_arr,
                r_hse_arr, p_hse_arr, pi_hse_arr,
                geom[lev].data(), z_cc_arr,
                l_gravity, l_rdOcp, l_moist, input_sounding_data);
        }
        else
        {
            // HSE will be calculated later with call to initHSE
            init_bx_scalars_from_input_sounding(
                bx, cons_arr,
                geom[lev].data(), z_cc_arr,
                l_moist, input_sounding_data);
        }

        init_bx_velocities_from_input_sounding(
            bx, xvel_arr, yvel_arr, zvel_arr,
            geom[lev].data(), z_nd_arr,
            input_sounding_data);

    } //mfi
}

/**
 * Box level wrapper for initializing scalar
 * data from input sounding data.
 *
 * @param bx Box specifying the indices we are initializing
 * @param state Array4 specifying the state data we are to initialize
 * @param geomdata GeometryData object specifying the domain geometry
 * @param inputSoundingData InputSoundingData object we are to initialize from
 */
void
init_bx_scalars_from_input_sounding (const amrex::Box &bx,
                                     amrex::Array4<amrex::Real> const &state,
                                     amrex::GeometryData const &geomdata,
                                     amrex::Array4<const amrex::Real> const &z_cc_arr,
                                     const bool& l_moist,
                                     InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d.dataPtr();
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d.dataPtr();
    const int   inp_sound_size  = inputSoundingData.size();

    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const amrex::Real z = (z_cc_arr) ? z_cc_arr(i,j,k)
                                         : prob_lo[2] + (k + 0.5) * dx[2];

        amrex::Real rho_0 = 1.0;

        // Set the density
        state(i, j, k, Rho_comp) = rho_0;

        // Initial Rho0*Theta0
        state(i, j, k, RhoTheta_comp) = rho_0 * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state(i, j, k, RhoScalar_comp) = 0;

        // total nonprecipitating water (Q1) == water vapor (Qv), i.e., there is no cloud water or cloud ice
        if (l_moist)
            state(i, j, k, RhoQ1_comp) = rho_0 * interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);

    });
}

/**
 * Box level wrapper for initializing scalar and hydrostatic
 * base state data from input sounding data.
 *
 * @param bx Box specifying the indices we are initializing
 * @param state Array4 specifying the state data we are to initialize
 * @param r_hse_arr Array4 specifying the density HSE base state data we are to initialize
 * @param p_hse_arr Array4 specifying the pressure HSE base state data we are to initialize
 * @param pi_hse_arr Array4 specifying the Exner pressure HSE base state data we are to initialize
 * @param geomdata GeometryData object specifying the domain geometry
 * @param l_gravity Real number specifying the gravitational acceleration constant
 * @param l_rdOcp Real number specifying the Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param inputSoundingData InputSoundingData object we are to initialize from
 */
void
init_bx_scalars_from_input_sounding_hse (const amrex::Box &bx,
                                         amrex::Array4<amrex::Real> const &state,
                                         amrex::Array4<amrex::Real> const &r_hse_arr,
                                         amrex::Array4<amrex::Real> const &p_hse_arr,
                                         amrex::Array4<amrex::Real> const &pi_hse_arr,
                                         amrex::GeometryData const &geomdata,
                                         amrex::Array4<const amrex::Real> const &z_cc_arr,
                                         const amrex::Real& l_gravity,
                                         const amrex::Real& l_rdOcp,
                                         const bool& l_moist,
                                         InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* rho_inp_sound   = inputSoundingData.rho_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d.dataPtr();
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d.dataPtr();
    const int   inp_sound_size  = inputSoundingData.size();

    // Geometry
    int ktop = bx.bigEnd(2);
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const amrex::Real z = (z_cc_arr) ? z_cc_arr(i,j,k)
                                         : prob_lo[2] + (k + 0.5) * dx[2];

        Real rho_k, rhoTh_k;

        // Set the density
        rho_k = interpolate_1d(z_inp_sound, rho_inp_sound, z, inp_sound_size);
        state(i, j, k, Rho_comp) = rho_k;

        // Initial Rho0*Theta0
        rhoTh_k = rho_k * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);
        state(i, j, k, RhoTheta_comp) = rhoTh_k;

        // Set scalar = A_0*exp(-10r^2), where r is distance from center of domain
        state(i, j, k, RhoScalar_comp) = 0;

        // Update hse quantities with values calculated from InputSoundingData.calc_rho_p()
        r_hse_arr (i, j, k) = rho_k;
        p_hse_arr (i, j, k) = getPgivenRTh(rhoTh_k); // NOTE WE ONLY USE THE DRY AIR PRESSURE
        pi_hse_arr(i, j, k) = getExnergivenRTh(rhoTh_k, l_rdOcp);

        // Boundary treatment
        if (k==0)
        {
            // set the ghost cell with dz and rho at boundary
            amrex::Real rho_surf =
                interpolate_1d(z_inp_sound, rho_inp_sound, 0.0, inp_sound_size);
            p_hse_arr (i, j, k-1) = p_hse_arr(i,j,k) + dx[2] * rho_surf * l_gravity;
            pi_hse_arr(i, j, k-1) = getExnergivenP(p_hse_arr(i, j, k-1), l_rdOcp);
        }
        else if (k==ktop)
        {
            // set the ghost cell with dz and rho at boundary
            amrex::Real rho_top =
                interpolate_1d(z_inp_sound, rho_inp_sound, z+dx[2]/2, inp_sound_size);
            p_hse_arr (i, j, k+1) = p_hse_arr(i,j,k) - dx[2] * rho_top * l_gravity;
            pi_hse_arr(i, j, k+1) = getExnergivenP(p_hse_arr(i, j, k+1), l_rdOcp);
        }

        // total nonprecipitating water (Q1) == water vapor (Qv), i.e., there
        // is no cloud water or cloud ice
        if (l_moist)
            state(i, j, k, RhoQ1_comp) = rho_k * interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);
    });
}

/**
 * Box level wrapper for initializing velocities from input sounding data.
 *
 * @param bx Box specifying the indices we are initializing
 * @param x_vel Array4 specifying the x-velocity data we are to initialize
 * @param y_vel Array4 specifying the y-velocity data we are to initialize
 * @param z_vel Array4 specifying the z-velocity data we are to initialize
 * @param geomdata GeometryData object specifying the domain geometry
 * @param inputSoundingData InputSoundingData object we are to initialize from
 */
void
init_bx_velocities_from_input_sounding (const amrex::Box &bx,
                                        amrex::Array4<amrex::Real> const &x_vel,
                                        amrex::Array4<amrex::Real> const &y_vel,
                                        amrex::Array4<amrex::Real> const &z_vel,
                                        amrex::GeometryData const &geomdata,
                                        amrex::Array4<const amrex::Real> const &z_nd_arr,
                                        InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d.dataPtr();
    const Real* U_inp_sound     = inputSoundingData.U_inp_sound_d.dataPtr();
    const Real* V_inp_sound     = inputSoundingData.V_inp_sound_d.dataPtr();
    const int   inp_sound_size  = inputSoundingData.size();

    // Geometry
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();

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
        const amrex::Real z = (z_nd_arr) ? 0.25*( z_nd_arr(i,j  ,k  )
                                                + z_nd_arr(i,j+1,k  )
                                                + z_nd_arr(i,j  ,k+1)
                                                + z_nd_arr(i,j+1,k+1))
                                         : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the x-velocity
        x_vel(i, j, k) = interpolate_1d(z_inp_sound, U_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of y-faces
        const amrex::Real z = (z_nd_arr) ? 0.25*( z_nd_arr(i  ,j,k  )
                                                + z_nd_arr(i+1,j,k  )
                                                + z_nd_arr(i  ,j,k+1)
                                                + z_nd_arr(i+1,j,k+1))
                                         : prob_lo[2] + (k + 0.5) * dx[2];

        // Set the y-velocity
        y_vel(i, j, k) = interpolate_1d(z_inp_sound, V_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of z-faces
        // Set the z-velocity
        z_vel(i, j, k) = 0.0;
    });
}
