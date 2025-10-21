/**
 * \file ERF_InitFromInputSounding.cpp
 */

#include <ERF.H>
#include <ERF_EOS.H>
#include <ERF_Constants.H>
#include <ERF_Utils.H>
#include <ERF_ProbCommon.H>

using namespace amrex;

void
init_bx_scalars_from_input_sounding (const Box &bx,
                                     Array4<Real> const &state,
                                     GeometryData const &geomdata,
                                     Array4<const Real> const &z_cc_arr,
                                     const bool& l_moist,
                                     InputSoundingData const &inputSoundingData);
void
init_bx_scalars_from_input_sounding_hse (const Box &bx,
                                         Array4<Real> const &state,
                                         Array4<Real> const &r_hse_arr,
                                         Array4<Real> const &p_hse_arr,
                                         Array4<Real> const &pi_hse_arr,
                                         Array4<Real> const &th_hse_arr,
                                         Array4<Real> const &qv_hse_arr,
                                         GeometryData const &geomdata,
                                         Array4<const Real> const &z_cc_arr,
                                         const Real& l_gravity,
                                         const Real& l_rdOcp,
                                         const bool& l_moist,
                                         InputSoundingData const &inputSoundingData,
                                         const bool& l_isentropic);

void
init_bx_velocities_from_input_sounding (const Box &bx,
                                        Array4<Real> const &x_vel,
                                        Array4<Real> const &y_vel,
                                        Array4<Real> const &z_vel,
                                        GeometryData const &geomdata,
                                        Array4<const Real> const &z_nd_arr,
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
        if (input_sounding_data.input_sounding_file.empty()) {
            Error("input_sounding file name must be provided via input");
        }

        input_sounding_data.resize_arrays();

        // this will interpolate the input profiles to the nominal height levels
        // (ranging from 0 to the domain top)
        for (int n = 0; n < input_sounding_data.n_sounding_files; n++) {
            input_sounding_data.read_from_file(geom[lev], zlevels_stag[lev], n);
        }

        // this will calculate the hydrostatically balanced density and pressure
        // profiles following WRF ideal.exe
        if (solverChoice.sounding_type == SoundingType::Ideal) {
            input_sounding_data.calc_rho_p(0);
        } else if (solverChoice.sounding_type == SoundingType::Isentropic ||
                   solverChoice.sounding_type == SoundingType::DryIsentropic) {
            input_sounding_data.assume_dry = (solverChoice.sounding_type == SoundingType::DryIsentropic);
            input_sounding_data.calc_rho_p_isentropic(0);
        }

    } else {
        //
        // We need to do this interp from coarse level in order to set the values of
        // the base state inside the domain but outside of the fine region
        //
        base_state[lev-1].FillBoundary(geom[lev-1].periodicity());
        //
        // NOTE: this interpolater assumes that ALL ghost cells of the coarse MultiFab
        //       have been pre-filled - this includes ghost cells both inside and outside
        //       the domain
        //
        InterpFromCoarseLevel(base_state[lev], base_state[lev].nGrowVect(),
                              IntVect(0,0,0), // do not fill ghost cells outside the domain
                              base_state[lev-1], 0, 0, base_state[lev].nComp(),
                              geom[lev-1], geom[lev],
                              refRatio(lev-1), &cell_cons_interp,
                              domain_bcs_type, BaseBCVars::rho0_bc_comp);

         // We need to do this here because the interpolation above may leave corners unfilled
         //    when the corners need to be filled by, for example, reflection of the fine ghost
         //    cell outside the fine region but inide the domain.
         (*physbcs_base[lev])(base_state[lev],0,base_state[lev].nComp(),base_state[lev].nGrowVect());
    }

    auto& lev_new = vars_new[lev];

    // updated if sounding is ideal (following WRF) or isentropic
    const bool l_isentropic = (solverChoice.sounding_type == SoundingType::Isentropic ||
                               solverChoice.sounding_type == SoundingType::DryIsentropic);
    const bool sounding_ideal_or_isentropic = (solverChoice.sounding_type == SoundingType::Ideal ||
                                               l_isentropic);
    MultiFab r_hse (base_state[lev], make_alias, BaseState::r0_comp, 1);
    MultiFab p_hse (base_state[lev], make_alias, BaseState::p0_comp, 1);
    MultiFab pi_hse(base_state[lev], make_alias, BaseState::pi0_comp, 1);
    MultiFab th_hse(base_state[lev], make_alias, BaseState::th0_comp, 1);
    MultiFab qv_hse(base_state[lev], make_alias, BaseState::qv0_comp, 1);

    const Real l_gravity = solverChoice.gravity;
    const Real l_rdOcp   = solverChoice.rdOcp;
    const bool l_moist   = (solverChoice.moisture_type != MoistureType::None);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(lev_new[Vars::cons], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box &bx = mfi.tilebox();
        const auto &cons_arr = lev_new[Vars::cons].array(mfi);
        const auto &xvel_arr = lev_new[Vars::xvel].array(mfi);
        const auto &yvel_arr = lev_new[Vars::yvel].array(mfi);
        const auto &zvel_arr = lev_new[Vars::zvel].array(mfi);
        Array4<Real>  r_hse_arr =  r_hse.array(mfi);
        Array4<Real>  p_hse_arr =  p_hse.array(mfi);
        Array4<Real> pi_hse_arr = pi_hse.array(mfi);
        Array4<Real> th_hse_arr = th_hse.array(mfi);
        Array4<Real> qv_hse_arr = qv_hse.array(mfi);

        Array4<Real const> z_cc_arr = (z_phys_cc[lev]) ? z_phys_cc[lev]->const_array(mfi) : Array4<Real const>{};
        Array4<Real const> z_nd_arr = (z_phys_nd[lev]) ? z_phys_nd[lev]->const_array(mfi) : Array4<Real const>{};

        if (sounding_ideal_or_isentropic)
        {
            // HSE will be initialized here, interpolated from values previously
            // calculated by calc_rho_p or calc_rho_p_isentropic
            init_bx_scalars_from_input_sounding_hse(
                bx, cons_arr,
                r_hse_arr, p_hse_arr, pi_hse_arr, th_hse_arr, qv_hse_arr,
                geom[lev].data(), z_cc_arr,
                l_gravity, l_rdOcp, l_moist, input_sounding_data,
                l_isentropic);
        }
        else
        {
            // This assumes rho_0 = 1.0
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
init_bx_scalars_from_input_sounding (const Box &bx,
                                     Array4<Real> const &state,
                                     GeometryData const &geomdata,
                                     Array4<const Real> const &z_cc_arr,
                                     const bool& l_moist,
                                     InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d[0].dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d[0].dataPtr();
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d[0].dataPtr();
    const int   inp_sound_size  = inputSoundingData.size(0);

    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real  z_lo = prob_lo[2];
    const Real  dz   = dx[2];

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const Real z = (z_cc_arr) ? z_cc_arr(i,j,k)
                                  : z_lo + (k + 0.5) * dz;

        Real rho_0 = 1.0;

        // Set the density
        state(i, j, k, Rho_comp) = rho_0;

        // Initial Rho0*Theta0
        state(i, j, k, RhoTheta_comp) = rho_0 * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);

        // Initialize all scalars to 0.
        for (int n = 0; n < NSCALARS; n++) {
            state(i, j, k, RhoScalar_comp+n) = 0;
        }

        // total nonprecipitating water (Q1) == water vapor (Qv), i.e., there is no cloud water or cloud ice
        if (l_moist) {
            state(i, j, k, RhoQ1_comp) = rho_0 * interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);
        }
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
 * @param th_hse_arr Array4 specifying the base state potential temperature we are to initialize
 * @param geomdata GeometryData object specifying the domain geometry
 * @param l_gravity Real number specifying the gravitational acceleration constant
 * @param l_rdOcp Real number specifying the Rhydberg constant ($R_d$) divided by specific heat at constant pressure ($c_p$)
 * @param inputSoundingData InputSoundingData object we are to initialize from
 */
void
init_bx_scalars_from_input_sounding_hse (const Box &bx,
                                         Array4<Real> const &state,
                                         Array4<Real> const &r_hse_arr,
                                         Array4<Real> const &p_hse_arr,
                                         Array4<Real> const &pi_hse_arr,
                                         Array4<Real> const &th_hse_arr,
                                         Array4<Real> const &qv_hse_arr,
                                         GeometryData const &geomdata,
                                         Array4<const Real> const &z_cc_arr,
                                         const Real& /*l_gravity*/,
                                         const Real& l_rdOcp,
                                         const bool& l_moist,
                                         InputSoundingData const &inputSoundingData,
                                         const bool& l_isentropic)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d[0].dataPtr();
    const Real* rho_inp_sound   = inputSoundingData.rho_inp_sound_d.dataPtr();
    const Real* theta_inp_sound = inputSoundingData.theta_inp_sound_d[0].dataPtr();
    const Real* qv_inp_sound    = inputSoundingData.qv_inp_sound_d[0].dataPtr();
    const int   inp_sound_size  = inputSoundingData.size(0);
    const bool anel_assume_dry  = inputSoundingData.assume_dry;

    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real  z_lo = prob_lo[2];
    const Real  dz   = dx[2];

    int kbot = geomdata.Domain().smallEnd(2);
    int ktop = geomdata.Domain().bigEnd(2);

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        const Real z = (z_cc_arr) ? z_cc_arr(i,j,k)
                                  : z_lo + (k + 0.5) * dz;

        Real rho_k = interpolate_1d(z_inp_sound, rho_inp_sound, z, inp_sound_size);
        Real rhoTh_k = rho_k * interpolate_1d(z_inp_sound, theta_inp_sound, z, inp_sound_size);

        Real rho_k_base = rho_k;
        if (l_isentropic) {
            // `rho_inp_sound` previously calculated in calc_rho_p_isentropic()
            // is in HSE and, when multiplied by the specified
            // `theta_input_sound`, give p, T, and theta that are consistent
            // with each other.
            //
            // Here, we do not require thermodynamic consistency between the
            // initial base state rho and the prognostic variables. Instead,
            // we calculate a `rho_hse` that is consistent with the isentropic
            // (constant theta) assumption.
            rho_k_base = rhoTh_k / theta_inp_sound[0];
        }

        // Set the density
        state(i, j, k, Rho_comp) = rho_k;

        // Initial Rho0*Theta0
        state(i, j, k, RhoTheta_comp) = rhoTh_k;

        // Initialize all scalars to 0.
        for (int n = 0; n < NSCALARS; n++) {
            state(i, j, k, RhoScalar_comp+n) = 0;
        }

        // total nonprecipitating water (Q1) == water vapor (Qv), i.e., there
        // is no cloud water or cloud ice
        Real qv_k = 0.0;
        if (l_moist) {
            qv_k = interpolate_1d(z_inp_sound, qv_inp_sound, z, inp_sound_size);
            state(i, j, k, RhoQ1_comp) = rho_k * qv_k;
        }

        // Update hse quantities with values calculated from InputSoundingData.calc_rho_p()
        if (anel_assume_dry) qv_k = 0;
        r_hse_arr (i,j,k) = rho_k_base * (1.0 + qv_k);
        p_hse_arr (i,j,k) = getPgivenRTh(rhoTh_k, qv_k);
        pi_hse_arr(i,j,k) = getExnergivenRTh(rhoTh_k, l_rdOcp, qv_k);
        th_hse_arr(i,j,k) = getRhoThetagivenP(p_hse_arr(i,j,k), qv_k) / rho_k_base;
        qv_hse_arr(i,j,k) = qv_k;
        if (l_isentropic) {
#if 0
            if (i==0 && j==0) {
                Print() << "HSE rho,p,T=pi*th,th at " << IntVect(i,j,k) << " : "
                    << r_hse_arr(i,j,k) << " "
                    << p_hse_arr(i,j,k) << " "
                    << pi_hse_arr(i,j,k)*th_hse_arr(i,j,k) << " "
                    << th_hse_arr(i,j,k)
                    << " with rho,rhotheta=" << rho_k << " " << rhoTh_k
                    << std::endl;
            }
#endif
            // If everything above is thermodynamically consistent, this should be constant
            AMREX_ALWAYS_ASSERT(std::abs(th_hse_arr(i,j,k) - theta_inp_sound[0]) < 1e-12);
        }

        // TODO: we should be setting this to the number of ghost cells of base_state[lev]
        //       instead of hard-wiring it here!
        int ng = 3;

        // FOEXTRAP hse arrays
        if (k==kbot)
        {
            for (int kk = 1; kk <= ng; kk++) {
                 r_hse_arr(i, j, k-kk) =  r_hse_arr(i,j,k);
                 p_hse_arr(i, j, k-kk) =  p_hse_arr(i,j,k);
                pi_hse_arr(i, j, k-kk) = pi_hse_arr(i,j,k);
                th_hse_arr(i, j, k-kk) = th_hse_arr(i,j,k);
                qv_hse_arr(i, j, k-kk) = qv_hse_arr(i,j,k);
            }
        }
        else if (k==ktop)
        {
            for (int kk = 1; kk <= ng; kk++) {
                 r_hse_arr(i, j, k+kk) =  r_hse_arr(i,j,k);
                 p_hse_arr(i, j, k+kk) =  p_hse_arr(i,j,k);
                pi_hse_arr(i, j, k+kk) = pi_hse_arr(i,j,k);
                qv_hse_arr(i, j, k+kk) = qv_hse_arr(i,j,k);
            }
        }
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
init_bx_velocities_from_input_sounding (const Box &bx,
                                        Array4<Real> const &x_vel,
                                        Array4<Real> const &y_vel,
                                        Array4<Real> const &z_vel,
                                        GeometryData const &geomdata,
                                        Array4<const Real> const &z_nd_arr,
                                        InputSoundingData const &inputSoundingData)
{
    const Real* z_inp_sound     = inputSoundingData.z_inp_sound_d[0].dataPtr();
    const Real* U_inp_sound     = inputSoundingData.U_inp_sound_d[0].dataPtr();
    const Real* V_inp_sound     = inputSoundingData.V_inp_sound_d[0].dataPtr();
    const int   inp_sound_size  = inputSoundingData.size(0);

    // Geometry
    const Real* prob_lo = geomdata.ProbLo();
    const Real* dx = geomdata.CellSize();
    const Real  z_lo = prob_lo[2];
    const Real  dz   = dx[2];

    // We want to set the lateral BC values, too
    Box gbx = bx; // Copy constructor
    gbx.grow(0,1); gbx.grow(1,1); // Grow by one in the lateral directions

    // Construct a box that is on x-faces
    const Box& xbx = surroundingNodes(gbx,0);
    // Construct a box that is on y-faces
    const Box& ybx = surroundingNodes(gbx,1);
    // Construct a box that is on z-faces
    const Box& zbx = surroundingNodes(gbx,2);

    // Set the x,y,z-velocities
    ParallelFor(xbx, ybx, zbx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of x-faces
        const Real z = (z_nd_arr) ? 0.25*( z_nd_arr(i,j  ,k  )
                                         + z_nd_arr(i,j+1,k  )
                                         + z_nd_arr(i,j  ,k+1)
                                         + z_nd_arr(i,j+1,k+1))
                                  : z_lo + (k + 0.5) * dz;

        // Set the x-velocity
        x_vel(i, j, k) = interpolate_1d(z_inp_sound, U_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of y-faces
        const Real z = (z_nd_arr) ? 0.25*( z_nd_arr(i  ,j,k  )
                                         + z_nd_arr(i+1,j,k  )
                                         + z_nd_arr(i  ,j,k+1)
                                         + z_nd_arr(i+1,j,k+1))
                                  : z_lo + (k + 0.5) * dz;

        // Set the y-velocity
        y_vel(i, j, k) = interpolate_1d(z_inp_sound, V_inp_sound, z, inp_sound_size);
    },
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
        // Note that this is called on a box of z-faces
        // Set the z-velocity
        z_vel(i, j, k) = 0.0;
    });
}
