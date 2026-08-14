#include "ERF_ImmersedForcing.H"
#include "ERF_TI_slow_headers.H"
#include "ERF_SrcHeaders.H"

using namespace amrex;

// helper function for immersed forcing wall model
/**
 * Compute the target velocity using Monin-Obukhov Similarity Theory for immersed forcing.
 *
 * @param[in] u1_2r First tangential velocity component.
 * @param[in] u2_2r Second tangential velocity component.
 * @param[in] delta Distance from the surface.
 * @param[in] z0 Roughness length.
 * @param[in] t_blank Volume fraction.
 * @param[in] theta_face Potential temperature at the face.
 * @param[in] theta_surf Potential temperature at the surface.
 * @param[in] tflux_in Surface heat flux.
 * @param[in] Olen_in Obukhov length.
 * @param[in] stability_correction Whether to apply stability corrections.
 * @return Target velocity component.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_if_most_target_vel(
    const amrex::Real u1_2r,
    const amrex::Real u2_2r,
    const amrex::Real delta,
    const amrex::Real z0,
    const amrex::Real t_blank,
    const amrex::Real theta_face,
    const amrex::Real theta_surf,
    const amrex::Real tflux_in,
    const amrex::Real Olen_in,
    const bool        stability_correction
)
{
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    Real psi_m                  = zero;
    Real psi_h                  = zero;
    Real tang_windspeed2r       = std::sqrt(u1_2r * u1_2r + u2_2r * u2_2r);

    Real ustar = tang_windspeed2r * KAPPA / (std::log(1.5 * delta / z0) - psi_m);
    Real tflux = (tflux_in != Real(1.e-8)) ? tflux_in : -(theta_face - theta_surf) * ustar * KAPPA / (std::log(1.5 * delta / z0) - psi_h);
    Real Olen  = (Olen_in != Real(1.e-8))  ? Olen_in  : -ustar * ustar * ustar * theta_face / (KAPPA * CONST_GRAV * tflux + tiny);
    Real zeta  = 1.5 * delta / Olen;

    // similarity functions
    similarity_funs sfuns;
    if (stability_correction){
        psi_m          = sfuns.calc_psi_m(zeta);
        psi_h          = sfuns.calc_psi_h(zeta);
    }
    ustar = tang_windspeed2r * KAPPA / (std::log(1.5 * delta / z0) - psi_m);

    // prevent some unphysical math
    if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
    if (!(ustar < 2.0  && !std::isnan(ustar))) { ustar = 2.0; }
    if (psi_m > std::log(myhalf * delta / z0)) { psi_m = std::log(myhalf * delta / z0); }

    Real uTarget      = (1 - t_blank) * ustar / KAPPA * (std::log(myhalf * delta / z0) - psi_m);
    Real u1Target     = uTarget * u1_2r / (tiny + tang_windspeed2r);

    return u1Target;
}

/**
 * Apply terrain immersed forcing to X-momentum
 */
void ImmersedForcingTerrain_Xmom (const Box& tbx,
                                  const Array4<const Real>& u,
                                  const Array4<const Real>& v,
                                  const Array4<const Real>& w,
                                  const Array4<const Real>& cell_data,
                                  const Array4<const Real>& t_blank_arr,
                                  const Array4<const Real>& t_blank_xface_arr,
                                  const Array4<const Real>& z_cc_arr,
                                  const Array4<      Real>& xmom_src_arr,
                                  const Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s = one; // unit velocity scale
    const bool l_implicit_drag = solverChoice.if_implicit_drag;

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg        = CONST_GRAV;
    const Real kappa      = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux_in           = solverChoice.if_surf_temp_flux;
    const Real Olen_in            = solverChoice.if_Olen_in;
    const bool l_use_most         = solverChoice.if_use_most;

    const Real small_volfrac = 0.005;

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux = u(i, j, k);
        const Real uy = fourth * ( v(i, j  , k  ) + v(i-1, j  , k  )
                               + v(i, j  , k+1) + v(i-1, j+1, k  ) );
        const Real uz = fourth * ( w(i, j  , k  ) + w(i-1, j  , k  )
                               + w(i, j  , k+1) + w(i-1, j  , k+1) );
        const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
        // Use face-centered terrain_blanking if available, otherwise average from cell centers
        Real t_blank_raw = (t_blank_xface_arr) ? t_blank_xface_arr(i, j, k) :
                           myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i-1, j, k));
        // Threshold: if averaged value is below small_volfrac, set to zero
        const Real t_blank = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        Real t_blank_above_raw = (t_blank_xface_arr) ? t_blank_xface_arr(i, j, k+1) :
                                 myhalf * (t_blank_arr(i, j, k+1) + t_blank_arr(i-1, j, k+1));
        const Real t_blank_above = (t_blank_above_raw < small_volfrac) ? zero : t_blank_above_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real rho_xface = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );

        if ((t_blank > 0 && (t_blank_above == zero)) && l_use_most) { // force to MOST value
            // calculate tangential velocity one cell above
            const Real ux2r = u(i, j, k+1) ;
            const Real uy2r = fourth * ( v(i, j  , k+1) + v(i-1, j  , k+1)
                               + v(i, j+1, k+1) + v(i-1, j+1, k+1) ) ;
            const Real h_windspeed2r = std::sqrt(ux2r * ux2r + uy2r * uy2r);

            // MOST
            const Real theta_xface = (myhalf * (cell_data(i,j,k  ,RhoTheta_comp) + cell_data(i-1,j,k, RhoTheta_comp))) / rho_xface;
            const Real rho_xface_below    = myhalf * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp) );
            const Real theta_xface_below  = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i-1,j,k-1, RhoTheta_comp))) / rho_xface_below;
            const Real theta_surf         = theta_xface_below;

            Real psi_m = zero;
            Real psi_h = zero;
            Real ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
            Real tflux = (tflux_in != Real(1e-8)) ? tflux_in : -(theta_xface - theta_surf) * ustar * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_h);
            Real Olen  = (Olen_in  != Real(1e-8)) ? Olen_in  : -ustar * ustar * ustar * theta_xface / (kappa * ggg * tflux + tiny);
            Real zeta  = Real(1.5) * dx_z / Olen;

            // similarity functions
            psi_m          = sfuns.calc_psi_m(zeta);
            psi_h          = sfuns.calc_psi_h(zeta);
            ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m);

            // prevent some unphysical math
            if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
            if (!(ustar < two && !std::isnan(ustar))) { ustar = two; }
            if (psi_m > std::log(myhalf * dx_z / z0)) { psi_m = std::log(myhalf * dx_z / z0); }

            // determine target velocity
            const Real uTarget  = ustar / kappa * (std::log(myhalf * dx_z / z0) - psi_m);
            Real uxTarget = uTarget * ux2r / (tiny + h_windspeed2r);
            const Real bc_forcing_x = -(uxTarget - ux); // BC forcing pushes nonrelative velocity toward target velocity
            const Real lambda = (1-t_blank) * CdM * U_s; // affine relaxation rate toward MOST target [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            xmom_src_arr(i, j, k) -= fac_local * rho_xface * bc_forcing_x; // if Vf low, force more strongly to MOST. If high, less forcing.
        } else {
            const Real lambda = t_blank * CdM * windspeed; // linear drag rate [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            xmom_src_arr(i, j, k) -= fac_local * rho_xface * ux;
        }
    });
}

/**
 * Apply terrain immersed forcing to Y-momentum
 */
void ImmersedForcingTerrain_Ymom (const Box& tby,
                                  const Array4<const Real>& u,
                                  const Array4<const Real>& v,
                                  const Array4<const Real>& w,
                                  const Array4<const Real>& cell_data,
                                  const Array4<const Real>& t_blank_arr,
                                  const Array4<const Real>& t_blank_yface_arr,
                                  const Array4<const Real>& z_cc_arr,
                                  const Array4<      Real>& ymom_src_arr,
                                  const Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s = one; // unit velocity scale
    const bool l_implicit_drag = solverChoice.if_implicit_drag;

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg        = CONST_GRAV;
    const Real kappa      = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux_in           = solverChoice.if_surf_temp_flux;
    const Real Olen_in            = solverChoice.if_Olen_in;
    const bool l_use_most         = solverChoice.if_use_most;

    const Real small_volfrac = 0.005;

    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux = fourth * ( u(i  , j  , k  ) + u(i  , j-1, k  )
                               + u(i+1, j  , k  ) + u(i+1, j-1, k  ) );
        const Real uy = v(i, j, k);
        const Real uz = fourth * ( w(i  , j  , k  ) + w(i  , j-1, k  )
                               + w(i  , j  , k+1) + w(i  , j-1, k+1) );
        const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
        // Use face-centered terrain_blanking if available, otherwise average from cell centers
        Real t_blank_raw = (t_blank_yface_arr) ? t_blank_yface_arr(i, j, k) :
                           myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i, j-1, k));
        const Real t_blank = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        Real t_blank_above_raw = (t_blank_yface_arr) ? t_blank_yface_arr(i, j, k+1) :
                                 myhalf * (t_blank_arr(i, j, k+1) + t_blank_arr(i, j-1, k+1));
        const Real t_blank_above = (t_blank_above_raw < small_volfrac) ? zero : t_blank_above_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real rho_yface =  myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );

        if ((t_blank > 0 && (t_blank_above == zero)) && l_use_most) { // force to MOST value
            // calculate tangential velocity one cell above
            const Real ux2r = fourth * ( u(i  , j  , k+1) + u(i  , j-1, k+1)
                               + u(i+1, j  , k+1) + u(i+1, j-1, k+1) );
            const Real uy2r = v(i, j, k+1) ;
            const Real h_windspeed2r = std::sqrt(ux2r * ux2r + uy2r * uy2r);

            // MOST
            const Real theta_yface = (myhalf * (cell_data(i,j,k  ,RhoTheta_comp) + cell_data(i,j-1,k, RhoTheta_comp))) / rho_yface;
            const Real rho_yface_below    =  myhalf * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
            const Real theta_yface_below  = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i,j-1,k-1, RhoTheta_comp))) / rho_yface_below;
            const Real theta_surf         = theta_yface_below;

            Real psi_m = zero;
            Real psi_h = zero;
            Real ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
            Real tflux = (tflux_in != Real(1e-8)) ? tflux_in : -(theta_yface - theta_surf) * ustar * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_h);
            Real Olen  = (Olen_in  != Real(1e-8)) ? Olen_in  : -ustar * ustar * ustar * theta_yface / (kappa * ggg * tflux + tiny);
            Real zeta  = Real(1.5) * dx_z / Olen;

            // similarity functions
            psi_m          = sfuns.calc_psi_m(zeta);
            psi_h          = sfuns.calc_psi_h(zeta);
            ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m);

            // prevent some unphysical math
            if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
            if (!(ustar < two && !std::isnan(ustar))) { ustar = two; }
            if (psi_m > std::log(myhalf * dx_z / z0)) { psi_m = std::log(myhalf * dx_z / z0); }

            // determine target velocity
            const Real uTarget  = ustar / kappa * (std::log(myhalf * dx_z / z0) - psi_m);
            Real uyTarget = uTarget * uy2r / (tiny + h_windspeed2r);
            const Real bc_forcing_y = -(uyTarget - uy);  // BC forcing pushes nonrelative velocity toward target velocity
            const Real lambda = (1 - t_blank) * CdM * U_s; // affine relaxation rate toward MOST target [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            ymom_src_arr(i, j, k) -= fac_local * rho_yface * bc_forcing_y; // if Vf low, force more strongly to MOST. If high, less forcing.
        } else {
            const Real lambda = t_blank * CdM * windspeed; // linear drag rate [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            ymom_src_arr(i, j, k) -= fac_local * rho_yface * uy;
        }
    });
}

/**
 * Apply terrain immersed forcing to Z-momentum
 */
void ImmersedForcingTerrain_Zmom (const Box& tbz,
                                  const Array4<const Real>& u,
                                  const Array4<const Real>& v,
                                  const Array4<const Real>& w,
                                  const Array4<const Real>& cell_data,
                                  const Array4<const Real>& t_blank_arr,
                                  const Array4<const Real>& t_blank_zface_arr,
                                  const Array4<const Real>& z_cc_arr,
                                  const Array4<      Real>& zmom_src_arr,
                                  const Geometry& geom,
                                  const SolverChoice& solverChoice,
                                  const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const bool l_implicit_drag = solverChoice.if_implicit_drag;

    const Real small_volfrac = 0.005;

    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux = fourth * ( u(i  , j  , k  ) + u(i+1, j  , k  )
                                 + u(i  , j  , k-1) + u(i+1, j  , k-1) );
        const Real uy = fourth * ( v(i  , j  , k  ) + v(i  , j+1, k  )
                                 + v(i  , j  , k-1) + v(i  , j+1, k-1) );
        const Real uz = w(i, j, k);
        const Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);
        // Use face-centered terrain_blanking if available, otherwise average from cell centers
        Real t_blank_raw = (t_blank_zface_arr) ? t_blank_zface_arr(i, j, k) :
                           myhalf * (t_blank_arr(i, j, k) + t_blank_arr(i, j, k-1));
        const Real t_blank = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real rho_zface =  myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
        const Real lambda = t_blank * CdM * windspeed; // linear drag rate [1/s]
        const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
        zmom_src_arr(i, j, k) -= fac_local * rho_zface * uz;
    });
}

/**
 * Apply buildings immersed forcing to X-momentum
 */
void ImmersedForcingBuildings_Xmom (const Box& tbx,
                                    const Array4<const Real>& u,
                                    const Array4<const Real>& v,
                                    const Array4<const Real>& w,
                                    const Array4<const Real>& cell_data,
                                    const Array4<const Real>& t_blank_arr,
                                    const Array4<const Real>& t_blank_xface_arr,
                                    const Array4<const Real>& z_cc_arr,
                                    const Array4<      Real>& xmom_src_arr,
                                    const Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;

    const Real alpha_m          = solverChoice.if_Cd_momentum;
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s              = one; // unit velocity scale

    // MOST parameters
    const Real z0                      = solverChoice.if_z0;
    const Real tflux_in                = solverChoice.if_surf_temp_flux;
    const Real Olen_in                 = solverChoice.if_Olen_in;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = 0.005;

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux   = u(i, j, k  );
        const Real uy   = fourth * ( v(i, j  , k  ) + v(i-1, j  , k  )
                                   + v(i, j+1, k  ) + v(i-1, j+1, k  ) );
        const Real uz   = fourth * ( w(i, j  , k  ) + w(i-1, j  , k  )
                                   + w(i, j  , k+1) + w(i-1, j  , k+1) );
        const amrex::Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

        const Real rho_xface   = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i-1,j,k,Rho_comp) );
        const Real theta_xface = (myhalf * (cell_data(i,j,k,RhoTheta_comp) + cell_data(i-1,j,k, RhoTheta_comp))) / rho_xface;

        // Use face-centered terrain_blanking if available, otherwise average from cell centers with threshold
        Real t_blank_raw       = (t_blank_xface_arr) ? t_blank_xface_arr(i, j  , k  ) :
                                 myhalf * (t_blank_arr(i, j  , k  ) + t_blank_arr(i-1, j  , k  ));
        const Real t_blank     = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        Real t_blank_below_raw = (k == 0) ? zero : (t_blank_xface_arr) ? t_blank_xface_arr(i, j  , k-1) :
                                 myhalf * (t_blank_arr(i, j  , k-1) + t_blank_arr(i-1, j  , k-1));
        const Real t_blank_below = (t_blank_below_raw < small_volfrac) ? zero : t_blank_below_raw;

        Real t_blank_above_raw = (t_blank_xface_arr) ? t_blank_xface_arr(i, j  , k+1) :
                                 myhalf * (t_blank_arr(i, j  , k+1) + t_blank_arr(i-1, j  , k+1));
        const Real t_blank_above = (t_blank_above_raw < small_volfrac) ? zero : t_blank_above_raw;

        Real t_blank_north_raw = (t_blank_xface_arr) ? t_blank_xface_arr(i, j+1, k  ) :
                                 myhalf * (t_blank_arr(i, j+1, k  ) + t_blank_arr(i-1, j+1, k  ));
        const Real t_blank_north = (t_blank_north_raw < small_volfrac) ? zero : t_blank_north_raw;

        Real t_blank_south_raw = (t_blank_xface_arr) ? t_blank_xface_arr(i, j-1, k  ) :
                                 myhalf * (t_blank_arr(i, j-1, k  ) + t_blank_arr(i-1, j-1, k  ));
        const Real t_blank_south = (t_blank_south_raw < small_volfrac) ? zero : t_blank_south_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real roof_mask     = (t_blank > zero && t_blank <  t_blank_below && t_blank_above == zero && l_use_most) ? one : zero; // roof cell
        const Real south_mask    = (t_blank > zero && t_blank <= t_blank_north && t_blank_south == zero && l_use_most) ? one : zero; // south wall cell
        const Real north_mask    = (t_blank > zero && t_blank <= t_blank_south && t_blank_north == zero && l_use_most) ? one : zero; // north wall cell
        const Real wall_mask     = (t_blank > zero && t_blank < one && !l_use_most) ? one : zero; // all walls when NOT using MOST
        const Real most_mask     = roof_mask + south_mask + north_mask; // cells getting MOST treatment
        const Real east_west_mask = (t_blank > zero && t_blank < one && l_use_most && most_mask == zero) ? one : zero; // partial cells not covered by MOST (east/west walls)
        const Real interior_mask = (t_blank == 1.0) ? one : zero; // interior cell

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_xface_inside = rho_xface;
        Real theta_surf       = theta_xface;
        Real bc_forcing_x     = zero;
        Real u_target         = zero;

        // roof forcing
        u1_cellaway         = u(i, j, k+1) ;
        u2_cellaway         = fourth * ( v(i, j  , k+1) + v(i-1, j  , k+1)
                                       + v(i, j+1, k+1) + v(i-1, j+1, k+1) ) ;
        rho_xface_inside    =  myhalf * (cell_data(i,j,k-1,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp));
        theta_surf          = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i-1,j,k-1, RhoTheta_comp))) / rho_xface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_z, z0, t_blank, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_x * roof_mask * rho_xface * CdM * U_s;

        // south wall forcing
        u1_cellaway         = u(i, j-1, k  );
        u2_cellaway         = fourth * ( w(i, j-1, k  ) + w(i-1, j-1, k  )
                                       + w(i, j-1, k+1) + w(i-1, j-1, k+1) ) ;
        rho_xface_inside    = myhalf * ( cell_data(i,j+1,k,Rho_comp) + cell_data(i-1,j+1,k,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i,j+1,k,RhoTheta_comp) + cell_data(i-1,j+1,k, RhoTheta_comp))) / rho_xface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_x * south_mask * rho_xface * CdM * U_s;

        // north wall forcing
        u1_cellaway         = u(i, j+1, k  ) ;
        u2_cellaway         = fourth * ( w(i, j+1, k  ) + w(i-1, j+1, k  )
                                       + w(i, j+1, k+1) + w(i-1, j+1, k+1) ) ;
        rho_xface_inside    = myhalf * ( cell_data(i,j-1,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i,j-1,k,RhoTheta_comp) + cell_data(i-1,j-1,k, RhoTheta_comp))) / rho_xface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_x * north_mask * rho_xface * CdM * U_s;

        // wall forcing (if not using most) or east/west walls when using MOST
        drag               += (wall_mask + east_west_mask) * t_blank * rho_xface * CdM * ux * windspeed;

        // interior cell forcing
        drag               += interior_mask * rho_xface * CdM * ux * windspeed;

        if (l_implicit_drag) {
            // point-implicit rescale of the aggregated drag
            const Real lambda = CdM * ( (roof_mask + south_mask + north_mask) * U_s
                                       + (wall_mask + east_west_mask) * t_blank * windspeed
                                       + interior_mask * windspeed );
            xmom_src_arr(i,j,k) -= drag / (one + lambda*dt);
        } else if (is_slow_step && !use_ImmersedForcing_fast) {
            // limit drag term for anelastic for numerical stability
            Real d_drag = dt * -drag; // time step * acceleration like tendency
            Real wsmax_change = damp_alpha * amrex::max(amrex::Math::abs(ux), ws_floor); // aims to prevent oscillations around 0.
            if (amrex::Math::abs(ux) < 0.1){ // no damping for smaller velocities
                wsmax_change =one * amrex::max(amrex::Math::abs(ux), ws_floor);
            }
            d_drag = amrex::min(amrex::max(d_drag, -wsmax_change), wsmax_change);
            xmom_src_arr(i,j,k) += d_drag / dt; // put back as limited tendency
        } else {
            xmom_src_arr(i, j, k) -= drag;
        }
    });
}

/**
 * Apply buildings immersed forcing to Y-momentum
 */
void ImmersedForcingBuildings_Ymom (const Box& tby,
                                    const Array4<const Real>& u,
                                    const Array4<const Real>& v,
                                    const Array4<const Real>& w,
                                    const Array4<const Real>& cell_data,
                                    const Array4<const Real>& t_blank_arr,
                                    const Array4<const Real>& t_blank_yface_arr,
                                    const Array4<const Real>& z_cc_arr,
                                    const Array4<      Real>& ymom_src_arr,
                                    const Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;

    const Real alpha_m          = solverChoice.if_Cd_momentum;
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s              = one; // unit velocity scale

    // MOST parameters
    const Real z0                      = solverChoice.if_z0;
    const Real tflux_in                = solverChoice.if_surf_temp_flux;
    const Real Olen_in                 = solverChoice.if_Olen_in;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = 0.005;

    ParallelFor(tby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux   = fourth * ( u(i  , j  , k  ) + u(i  , j-1, k  )
                                   + u(i+1, j  , k  ) + u(i+1, j-1, k  ) );
        const Real uy   = v(i, j, k);
        const Real uz   = fourth * ( w(i  , j  , k  ) + w(i  , j-1, k  )
                                   + w(i  , j  , k+1) + w(i  , j-1, k+1) );
        const amrex::Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

        const Real rho_yface   = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j-1,k,Rho_comp) );
        const Real theta_yface = (myhalf * (cell_data(i,j,k  ,RhoTheta_comp) + cell_data(i,j-1,k,RhoTheta_comp))) / rho_yface;

        // Use face-centered terrain_blanking if available, otherwise average from cell centers with threshold
        Real t_blank_raw       = (t_blank_yface_arr) ? t_blank_yface_arr(i  , j  , k  ) :
                                 myhalf * (t_blank_arr(i  , j  , k  ) + t_blank_arr(i  , j-1, k  ));
        const Real t_blank     = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        Real t_blank_below_raw = (k == 0) ? zero : (t_blank_yface_arr) ? t_blank_yface_arr(i  , j  , k-1) :
                                 myhalf * (t_blank_arr(i  , j  , k-1) + t_blank_arr(i  , j-1, k-1));
        const Real t_blank_below = (t_blank_below_raw < small_volfrac) ? zero : t_blank_below_raw;

        Real t_blank_above_raw = (t_blank_yface_arr) ? t_blank_yface_arr(i  , j  , k+1) :
                                 myhalf * (t_blank_arr(i  , j  , k+1) + t_blank_arr(i  , j-1, k+1));
        const Real t_blank_above = (t_blank_above_raw < small_volfrac) ? zero : t_blank_above_raw;

        Real t_blank_east_raw  = (t_blank_yface_arr) ? t_blank_yface_arr(i+1, j  , k  ) :
                                 myhalf * (t_blank_arr(i+1, j  , k  ) + t_blank_arr(i+1, j-1, k  ));
        const Real t_blank_east = (t_blank_east_raw < small_volfrac) ? zero : t_blank_east_raw;

        Real t_blank_west_raw  = (t_blank_yface_arr) ? t_blank_yface_arr(i-1, j  , k  ) :
                                 myhalf * (t_blank_arr(i-1, j  , k  ) + t_blank_arr(i-1, j-1, k  ));
        const Real t_blank_west = (t_blank_west_raw < small_volfrac) ? zero : t_blank_west_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real roof_mask     = (t_blank > zero && t_blank <  t_blank_below && t_blank_above == zero && l_use_most) ? one : zero; // roof cell
        const Real west_mask     = (t_blank > zero && t_blank <= t_blank_east  && t_blank_west  == zero && l_use_most) ? one : zero; // west wall cell
        const Real east_mask     = (t_blank > zero && t_blank <= t_blank_west  && t_blank_east  == zero && l_use_most) ? one : zero; // east wall cell
        const Real wall_mask     = (t_blank > zero && t_blank < one && !l_use_most) ? one : zero; // all walls when NOT using MOST
        const Real most_mask     = roof_mask + west_mask + east_mask; // cells getting MOST treatment
        const Real north_south_mask = (t_blank > zero && t_blank < one && l_use_most && most_mask == zero) ? one : zero; // partial cells not covered by MOST (north/south walls)
        const Real interior_mask = (t_blank == 1.0) ? one : zero; // interior cell

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_yface_inside = rho_yface;
        Real theta_surf       = theta_yface;
        Real bc_forcing_y     = zero;
        Real u_target         = zero;

        // roof forcing
        u1_cellaway         = fourth * ( u(i  , j  , k+1) + u(i  , j-1, k+1)
                                       + u(i+1, j  , k+1) + u(i+1, j-1, k+1) );
        u2_cellaway         = v(i, j, k+1);
        rho_yface_inside    = myhalf * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i,j-1,k-1,RhoTheta_comp))) / rho_yface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_z, z0, t_blank, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_y * roof_mask * rho_yface * CdM * U_s;

        // west wall forcing
        u1_cellaway         = v(i-1, j , k  );
        u2_cellaway         = fourth * ( w(i-1, j  , k  ) + w(i-1, j-1, k  )
                                       + w(i-1, j  , k+1) + w(i-1, j-1, k+1) );
        rho_yface_inside    = myhalf * ( cell_data(i+1,j,k,Rho_comp) + cell_data(i+1,j-1,k,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i+1,j,k,RhoTheta_comp) + cell_data(i+1,j-1,k,RhoTheta_comp))) / rho_yface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_y * west_mask * rho_yface * CdM * U_s;

        // east wall forcing
        u1_cellaway         = v(i+1, j , k  );
        u2_cellaway         = fourth * ( w(i+1, j  , k  ) + w(i+1, j-1, k  )
                                       + w(i+1, j  , k+1) + w(i+1, j-1, k+1) );
        rho_yface_inside    = myhalf * ( cell_data(i-1,j,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i-1,j,k,RhoTheta_comp) + cell_data(i-1,j-1,k,RhoTheta_comp))) / rho_yface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_y * east_mask * rho_yface * CdM * U_s;

        // wall forcing (if not using most) or north/south walls when using MOST
        drag               += (wall_mask + north_south_mask) * t_blank * rho_yface * CdM * uy * windspeed;

        // interior cell forcing
        drag               += interior_mask * rho_yface * CdM * uy * windspeed;

        if (l_implicit_drag) {
            // point-implicit rescale of the aggregated drag
            const Real lambda = CdM * ( (roof_mask + west_mask + east_mask) * U_s
                                       + (wall_mask + north_south_mask) * t_blank * windspeed
                                       + interior_mask * windspeed );
            ymom_src_arr(i,j,k) -= drag / (one + lambda*dt);
        } else if (is_slow_step && !use_ImmersedForcing_fast) {
            // limit drag term for anelastic for numerical stability
            Real d_drag = dt * -drag; // time step * acceleration like tendency
            Real wsmax_change = damp_alpha * amrex::max(amrex::Math::abs(uy), ws_floor); // aims to prevent oscillations around 0.
            if (amrex::Math::abs(uy) < 0.1){ // no damping for smaller velocities
                wsmax_change =one * amrex::max(amrex::Math::abs(uy), ws_floor);
            }
            d_drag = amrex::min(amrex::max(d_drag, -wsmax_change), wsmax_change);
            ymom_src_arr(i,j,k) += d_drag / dt; // put back as limited tendency
        } else {
            ymom_src_arr(i, j, k) -= drag;
        }
    });
}

/**
 * Apply buildings immersed forcing to Z-momentum
 */
void ImmersedForcingBuildings_Zmom (const Box& tbz,
                                    const Array4<const Real>& u,
                                    const Array4<const Real>& v,
                                    const Array4<const Real>& w,
                                    const Array4<const Real>& cell_data,
                                    const Array4<const Real>& t_blank_arr,
                                    const Array4<const Real>& t_blank_zface_arr,
                                    const Array4<const Real>& z_cc_arr,
                                    const Array4<      Real>& zmom_src_arr,
                                    const Geometry& geom,
                                    const SolverChoice& solverChoice,
                                    const Real fac)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;

    const Real alpha_m          = solverChoice.if_Cd_momentum;
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s              = one; // unit velocity scale

    // MOST parameters
    const Real z0                      = solverChoice.if_z0;
    const Real tflux_in                = solverChoice.if_surf_temp_flux;
    const Real Olen_in                 = solverChoice.if_Olen_in;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = 0.005;

    ParallelFor(tbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux   = fourth * ( u(i  , j  , k  ) + u(i+1, j  , k  )
                                   + u(i  , j  , k-1) + u(i+1, j  , k-1) );
        const Real uy   = fourth * ( v(i, j  , k  ) + v(i, j+1, k  )
                                   + v(i, j  , k-1) + v(i, j+1, k-1) );
        const Real uz   = w(i, j, k);
        const amrex::Real windspeed = std::sqrt(ux * ux + uy * uy + uz * uz);

        const Real rho_zface   = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i,j,k-1,Rho_comp) );
        const Real theta_zface = (myhalf * (cell_data(i,j,k,RhoTheta_comp) + cell_data(i,j,k-1,RhoTheta_comp))) / rho_zface;

        // Use face-centered terrain_blanking if available, otherwise average from cell centers with threshold
        Real t_blank_raw       = (t_blank_zface_arr) ? t_blank_zface_arr(i  ,j  , k  ) :
                                 myhalf * (t_blank_arr(i  ,j  , k)   + t_blank_arr(i  , j  , k-1));
        const Real t_blank     = (t_blank_raw < small_volfrac) ? zero : t_blank_raw;

        Real t_blank_below_raw = (k == 0) ? zero : (t_blank_zface_arr) ? t_blank_zface_arr(i  ,j  , k-1) :
                                 myhalf * (t_blank_arr(i  ,j  , k-1) + t_blank_arr(i  , j  , k-2));
        const Real t_blank_below = (t_blank_below_raw < small_volfrac) ? zero : t_blank_below_raw;

        Real t_blank_above_raw = (t_blank_zface_arr) ? t_blank_zface_arr(i  ,j  , k+1) :
                                 myhalf * (t_blank_arr(i  ,j  , k)   + t_blank_arr(i  , j  , k+1));
        const Real t_blank_above = (t_blank_above_raw < small_volfrac) ? zero : t_blank_above_raw;

        Real t_blank_north_raw = (t_blank_zface_arr) ? t_blank_zface_arr(i  , j+1, k  ) :
                                 myhalf * (t_blank_arr(i  ,j+1, k)   + t_blank_arr(i  , j+1, k-1));
        const Real t_blank_north = (t_blank_north_raw < small_volfrac) ? zero : t_blank_north_raw;

        Real t_blank_south_raw = (t_blank_zface_arr) ? t_blank_zface_arr(i  , j-1, k  ) :
                                 myhalf * (t_blank_arr(i  ,j-1, k)   + t_blank_arr(i  , j-1, k-1));
        const Real t_blank_south = (t_blank_south_raw < small_volfrac) ? zero : t_blank_south_raw;

        Real t_blank_east_raw  = (t_blank_zface_arr) ? t_blank_zface_arr(i+1, j  , k  ) :
                                 myhalf * (t_blank_arr(i+1,j  , k)   + t_blank_arr(i+1, j  , k-1));
        const Real t_blank_east = (t_blank_east_raw < small_volfrac) ? zero : t_blank_east_raw;

        Real t_blank_west_raw  = (t_blank_zface_arr) ? t_blank_zface_arr(i-1, j  , k  ) :
                                 myhalf * (t_blank_arr(i-1,j  , k)   + t_blank_arr(i-1, j  , k-1));
        const Real t_blank_west = (t_blank_west_raw < small_volfrac) ? zero : t_blank_west_raw;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);

        const Real south_mask    = (t_blank > zero && t_blank <= t_blank_north && t_blank_south == zero && l_use_most && k >= 1) ? one : zero; // south wall cell
        const Real north_mask    = (t_blank > zero && t_blank <= t_blank_south && t_blank_north == zero && l_use_most && k >= 1) ? one : zero; // north wall cell
        const Real west_mask     = (t_blank > zero && t_blank <= t_blank_east  && t_blank_west  == zero && l_use_most && k >= 1) ? one : zero; // west wall cell
        const Real east_mask     = (t_blank > zero && t_blank <= t_blank_west  && t_blank_east  == zero && l_use_most && k >= 1) ? one : zero; // east wall cell
        const Real wall_mask     = (t_blank > zero && t_blank < one && !l_use_most) ? one : zero; // all walls when NOT using MOST
        const Real most_mask     = south_mask + north_mask + west_mask + east_mask; // cells getting MOST treatment
        const Real roof_mask     = (t_blank > zero && t_blank_above == zero && l_use_most) ? one : zero; // roof cell (horizontal surface) - uses simple drag
        const Real interior_mask = (t_blank == 1.0) ? one : zero; // interior cell

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_zface_inside = rho_zface;
        Real theta_surf       = theta_zface;
        Real bc_forcing_z     = zero;
        Real u_target         = zero;

        // south wall forcing
        u1_cellaway         = fourth * ( u(i  , j-1, k  ) + u(i+1, j-1, k  )
                                       + u(i  , j-1, k-1) + u(i+1, j-1, k-1) );
        u2_cellaway         = w(i, j-1, k);
        rho_zface_inside    = myhalf * ( cell_data(i,j+1,k,Rho_comp) + cell_data(i,j+1,k-1,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i,j+1,k,RhoTheta_comp) + cell_data(i,j+1,k-1,RhoTheta_comp))) / rho_zface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_z * south_mask * rho_zface * CdM * U_s;

        // north wall forcing
        u1_cellaway         = fourth * ( u(i  , j+1, k  ) + u(i+1, j+1, k  )
                                       + u(i  , j+1, k-1) + u(i+1, j+1, k-1) );
        u2_cellaway         = w(i, j+1, k);
        rho_zface_inside    = myhalf * ( cell_data(i,j-1,k,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i,j-1,k,RhoTheta_comp) + cell_data(i,j-1,k-1,RhoTheta_comp))) / rho_zface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_z * north_mask * rho_zface * CdM * U_s;

        // west wall forcing
        u1_cellaway         = fourth * ( v(i-1, j  , k  ) + v(i-1, j+1, k  )
                                       + v(i-1, j  , k-1) + v(i-1, j+1, k-1) );
        u2_cellaway         = w(i-1, j, k);
        rho_zface_inside    = myhalf * ( cell_data(i+1,j,k,Rho_comp) + cell_data(i+1,j,k-1,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i+1,j,k,RhoTheta_comp) + cell_data(i+1,j,k-1,RhoTheta_comp))) / rho_zface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_z * west_mask * rho_zface * CdM * U_s;

        // east wall forcing
        u1_cellaway         = fourth * ( v(i+1, j  , k  ) + v(i+1, j+1, k  )
                                       + v(i+1, j  , k-1) + v(i+1, j+1, k-1) );
        u2_cellaway         = w(i+1, j, k);
        rho_zface_inside    = myhalf * ( cell_data(i-1,j,k,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp) );
        theta_surf          = (myhalf * (cell_data(i-1,j,k,RhoTheta_comp) + cell_data(i-1,j,k-1,RhoTheta_comp))) / rho_zface_inside;
        u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction);
        bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
        drag               += bc_forcing_z * east_mask * rho_zface * CdM * U_s;

        // wall forcing (if not using most) or roof when using MOST
        drag               += (wall_mask + roof_mask) * t_blank * rho_zface * CdM * uz * windspeed;

        // interior cell forcing
        drag               += interior_mask * rho_zface * CdM * uz * windspeed;

        if (l_implicit_drag) {
            // point-implicit rescale of the aggregated drag
            const Real lambda = CdM * ( (south_mask + north_mask + west_mask + east_mask) * U_s
                                       + (wall_mask + roof_mask) * t_blank * windspeed
                                       + interior_mask * windspeed );
            zmom_src_arr(i,j,k) -= drag / (one + lambda*dt);
        } else if (is_slow_step && !use_ImmersedForcing_fast) {
            // limit drag term for anelastic for numerical stability
            Real d_drag = dt * -drag; // time step * acceleration like tendency
            Real wsmax_change = damp_alpha * amrex::max(amrex::Math::abs(uz), ws_floor); // aims to prevent oscillations around 0.
            if (amrex::Math::abs(uz) < 0.1){ // no damping for smaller velocities
                wsmax_change = one * amrex::max(amrex::Math::abs(uz), ws_floor);
            }
            d_drag = amrex::min(amrex::max(d_drag, -wsmax_change), wsmax_change);
            zmom_src_arr(i,j,k) += d_drag / dt; // put back as limited tendency
        } else {
            zmom_src_arr(i, j, k) -= drag;
        }
    });
}

/**
 * Apply terrain immersed forcing to scalars (Rho, RhoTheta)
 */
void ImmersedForcingTerrain_Scalar (const Box& bx,
                                   const Array4<const Real>& u,
                                   const Array4<const Real>& v,
                                   const Array4<const Real>& cell_data,
                                   const Array4<const Real>& t_blank_arr,
                                   const Array4<const Real>& z_cc_arr,
                                   const Array4<      Real>& cell_src,
                                   const Geometry& geom,
                                   const SolverChoice& solverChoice,
                                   const Table1D<Real>& r_avg,
                                   const Table1D<Real>& t_avg,
                                   const Real time)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];

    const Real alpha_h          = solverChoice.if_Cd_scalar;
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s              = one; // unit velocity scale

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg                = CONST_GRAV;
    const Real kappa              = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux              = solverChoice.if_surf_temp_flux;
    const Real init_surf_temp     = solverChoice.if_init_surf_temp;

    // Note this has been converted to K / s when it was read in;
    const Real surf_heating_rate  = solverChoice.if_surf_heating_rate;

    const Real Olen_in            = solverChoice.if_Olen_in;

    ParallelFor(bx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, one/three);

        const Real t_blank       = t_blank_arr(i, j, k);
        const Real t_blank_above = t_blank_arr(i, j, k+1);
        const Real ux_cc_2r = myhalf * (u(i  ,j  ,k+1) + u(i+1,j  ,k+1));
        const Real uy_cc_2r = myhalf * (v(i  ,j  ,k+1) + v(i  ,j+1,k+1));
        const Real h_windspeed2r  = std::sqrt(ux_cc_2r * ux_cc_2r + uy_cc_2r * uy_cc_2r);

        const Real theta          = cell_data(i,j,k  ,RhoTheta_comp) / cell_data(i,j,k  ,Rho_comp);
        const Real theta_neighbor = cell_data(i,j,k+1,RhoTheta_comp) / cell_data(i,j,k+1,Rho_comp);

        // SURFACE TEMP AND HEATING/COOLING RATE
        if (init_surf_temp > zero) {
            if (t_blank > 0 && (t_blank_above == zero)) { // force to MOST value
                const Real surf_temp    = init_surf_temp + surf_heating_rate*time;
                const Real bc_forcing_rt_srf = -(cell_data(i,j,k-1,Rho_comp) * surf_temp - cell_data(i,j,k-1,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
            }
        }

        // SURFACE HEAT FLUX
        if (tflux != Real(1e-8)){
            if (t_blank > 0 && (t_blank_above == zero)) { // force to MOST value
                Real psi_m           = zero;
                Real psi_h           = zero;
                Real psi_h_neighbor  = zero;
                Real ustar = h_windspeed2r * kappa / (std::log((Real(1.5)) * dx_z / z0) - psi_m);
                const Real Olen  = -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny);
                const Real zeta          = (myhalf) * dx_z / Olen;
                const Real zeta_neighbor = (Real(1.5)) * dx_z / Olen;

                // similarity functions
                psi_m          = sfuns.calc_psi_m(zeta);
                psi_h          = sfuns.calc_psi_h(zeta);
                psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                ustar = h_windspeed2r * kappa / (std::log((Real(1.5)) * dx_z / z0) - psi_m);

                // prevent some unphysical math
                if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
                if (!(ustar < two && !std::isnan(ustar))) { ustar = two; }
                if (psi_h_neighbor > std::log(Real(1.5) * dx_z / z0)) { psi_h_neighbor = std::log(Real(1.5) * dx_z / z0); }
                if (psi_h > std::log(myhalf * dx_z / z0)) { psi_h = std::log(myhalf * dx_z / z0); }

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((Real(1.5)) * dx_z / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * dx_z / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
            }
        }

        // OBUKHOV LENGTH
        if (Olen_in != Real(1e-8)){
            if (t_blank > 0 && (t_blank_above == zero)) { // force to MOST value
                const Real Olen  = Olen_in;
                const Real zeta          = (myhalf) * dx_z / Olen;
                const Real zeta_neighbor = (Real(1.5)) * dx_z / Olen;

                // similarity functions
                const Real psi_m          = sfuns.calc_psi_m(zeta);
                const Real psi_h          = sfuns.calc_psi_h(zeta);
                const Real psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                const Real ustar = h_windspeed2r * kappa / (std::log((Real(1.5)) * dx_z / z0) - psi_m);

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((Real(1.5)) * dx_z / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * dx_z / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
            }
        }

        // Force fully immersed cells to planar average rho and theta
        if (t_blank == one && r_avg && t_avg) {
            const Real rho_avg = r_avg(k);
            const Real theta_avg = t_avg(k) / rho_avg;  // Convert from RhoTheta to Theta
            const Real rho_cell = cell_data(i,j,k,Rho_comp);
            const Real theta_cell = cell_data(i,j,k,RhoTheta_comp) / rho_cell;
            const Real bc_forcing_r = -(rho_avg - rho_cell);
            const Real bc_forcing_rt = -(rho_avg * theta_avg - cell_data(i,j,k,RhoTheta_comp));

            cell_src(i, j, k, Rho_comp) -= drag_coefficient * U_s * bc_forcing_r;
            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
        }
    });
}

/**
 * Apply buildings immersed forcing to scalars (Rho, RhoTheta)
 */
void ImmersedForcingBuildings_Scalar (const Box& bx,
                                     const Array4<const Real>& u,
                                     const Array4<const Real>& v,
                                     const Array4<const Real>& w,
                                     const Array4<const Real>& cell_data,
                                     const Array4<const Real>& t_blank_arr,
                                     const Array4<const Real>& z_cc_arr,
                                     const Array4<      Real>& cell_src,
                                     const Geometry& geom,
                                     const SolverChoice& solverChoice,
                                     const Table1D<Real>& r_avg,
                                     const Table1D<Real>& t_avg,
                                     const Real time)
{
    // geometric properties
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];

    const Real alpha_h          = solverChoice.if_Cd_scalar;
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s              = one; // unit velocity scale

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg                = CONST_GRAV;
    const Real kappa              = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux              = solverChoice.if_surf_temp_flux;
    const Real init_surf_temp     = solverChoice.if_init_surf_temp;
    const Real surf_heating_rate  = solverChoice.if_surf_heating_rate;
    const Real Olen_in            = solverChoice.if_Olen_in;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real t_blank       = t_blank_arr(i, j, k);
        const Real t_blank_below = t_blank_arr(i, j, k-1);
        const Real t_blank_above = t_blank_arr(i, j, k+1);
        const Real t_blank_north = t_blank_arr(i  , j+1, k);
        const Real t_blank_south = t_blank_arr(i  , j-1, k);
        const Real t_blank_east  = t_blank_arr(i+1, j  , k);
        const Real t_blank_west  = t_blank_arr(i-1, j  , k);

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, one/three);

        // SURFACE TEMP AND HEATING/COOLING RATE
        if (init_surf_temp > zero) {
            const Real surf_temp    = init_surf_temp + surf_heating_rate*time;
            if (t_blank > 0 && (t_blank_above == zero) && (t_blank_below == one)) { // building roof
                const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;

            } else if (((t_blank > zero && t_blank < t_blank_west && t_blank_east == zero) ||
                        (t_blank > zero && t_blank < t_blank_east && t_blank_west == zero) ||
                        (t_blank > zero && t_blank < t_blank_north && t_blank_south == zero) ||
                        (t_blank > zero && t_blank < t_blank_south && t_blank_north == zero))) {
                // this should enter for just building walls
                // walls are currently separated to allow for flexibility in the future to heat walls differently

                // south face
                if ((t_blank < t_blank_north) && (t_blank_north == one)) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // north face
                if ((t_blank < t_blank_south) && (t_blank_south == one)) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // west face
                if ((t_blank < t_blank_east) && (t_blank_east == one)) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // east face
                if ((t_blank < t_blank_west) && (t_blank_west == one)) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

            }
        }

        // SURFACE HEAT FLUX
        if (tflux != Real(1.e-8)){
            const Real ux_cc_2r = myhalf * (u(i  ,j  ,k+1) + u(i+1,j  ,k+1));
            const Real uy_cc_2r = myhalf * (v(i  ,j  ,k+1) + v(i  ,j+1,k+1));
            const Real h_windspeed2r  = std::sqrt(ux_cc_2r * ux_cc_2r + uy_cc_2r * uy_cc_2r);

            const Real theta          = cell_data(i,j,k  ,RhoTheta_comp) / cell_data(i,j,k  ,Rho_comp);
            Real theta_neighbor       = cell_data(i,j,k+1,RhoTheta_comp) / cell_data(i,j,k+1,Rho_comp);

            if (t_blank > zero && (t_blank_above == zero)) { // building roof
                Real psi_m           = zero;
                Real psi_h           = zero;
                Real psi_h_neighbor  = zero;
                Real ustar           = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);
                Real Olen            = (Olen_in  != Real(1e-8)) ? Olen_in  : -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny);

                for (int iter = 0; iter < 2; ++iter) {
                    if (iter > 0) { Olen  = -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny); }
                    Real zeta          = (myhalf) * dx_z / Olen;
                    Real zeta_neighbor = (1.5)    * dx_z / Olen;

                    // similarity functions
                    psi_m          = sfuns.calc_psi_m(zeta);
                    psi_h          = sfuns.calc_psi_h(zeta);
                    psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                    ustar = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);
                }

                // prevent some unphysical math
                if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
                if (!(ustar < 2.0  && !std::isnan(ustar))) { ustar = 2.0; }
                if (psi_h_neighbor > std::log(1.5 * dx_z / z0)) { psi_h_neighbor = std::log(1.5 * dx_z / z0); }
                if (psi_h > std::log(myhalf * dx_z / z0)) { psi_h = std::log(myhalf * dx_z / z0); }

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * dx_z / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * dx_z / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;

            } else if (((t_blank > zero && t_blank < t_blank_west && t_blank_east == zero) ||
                        (t_blank > zero && t_blank < t_blank_east && t_blank_west == zero) ||
                        (t_blank > zero && t_blank < t_blank_north && t_blank_south == zero) ||
                        (t_blank > zero && t_blank < t_blank_south && t_blank_north == zero))) { // this should enter for just building walls

                Real ux_cellaway = zero;
                Real uy_cellaway = zero;
                Real uz_cellaway = zero;
                Real u1          = zero;
                Real u2          = zero;
                Real delta       = zero;

                // south face
                if (t_blank > zero && t_blank < t_blank_north && t_blank_south == zero) {
                    ux_cellaway = myhalf * (u(i  ,j-1,k) + u(i+1,j-1,k  ));
                    uz_cellaway = myhalf * (w(i  ,j-1,k) + w(i  ,j-1,k+1));
                    u1 = ux_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_y;

                    // MOST
                    theta_neighbor = cell_data(i,j-1,k,RhoTheta_comp) / cell_data(i,j-1,k,Rho_comp);
                }

                // north face
                if (t_blank > zero && t_blank < t_blank_south && t_blank_north == zero) {
                    ux_cellaway = myhalf * (u(i  ,j+1,k) + u(i+1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i  ,j+1,k) + w(i  ,j+1,k+1));
                    u1 = ux_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_y;

                    // MOST
                    theta_neighbor = cell_data(i,j+1,k,RhoTheta_comp) / cell_data(i,j+1,k,Rho_comp);
                }

                // west face
                if (t_blank > zero && t_blank < t_blank_east && t_blank_west == zero) {
                    uy_cellaway = myhalf * (v(i-1,j  ,k) + v(i-1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i-1,j  ,k) + w(i-1,j  ,k+1));
                    u1 = uy_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_x;

                    // MOST
                    theta_neighbor = cell_data(i-1,j,k,RhoTheta_comp) / cell_data(i-1,j,k,Rho_comp);
                }

                // east face
                if (t_blank > zero && t_blank < t_blank_west && t_blank_east == zero) {
                    uy_cellaway = myhalf * (v(i+1,j  ,k) + v(i+1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i+1,j  ,k) + w(i+1,j  ,k+1));
                    u1 = uy_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_x;

                    // MOST
                    theta_neighbor = cell_data(i+1,j,k,RhoTheta_comp) / cell_data(i+1,j,k,Rho_comp);
                }

                Real tan_wspd = std::sqrt(u1 * u1 + u2 * u2);

                Real psi_m           = zero;
                Real psi_h           = zero;
                Real psi_h_neighbor  = zero;
                Real ustar           = tan_wspd * kappa / (std::log(1.5 * delta / z0) - psi_m);
                Real Olen            = (Olen_in  != Real(1e-8)) ? Olen_in  : -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny);

                for (int iter = 0; iter < 2; ++iter) {
                    if (iter > 0) { Olen  = -ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny); }
                    Real zeta          = (myhalf) * delta / Olen;
                    Real zeta_neighbor = (1.5)    * delta / Olen;

                    // similarity functions
                    psi_m          = sfuns.calc_psi_m(zeta);
                    psi_h          = sfuns.calc_psi_h(zeta);
                    psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                    ustar = tan_wspd * kappa / (std::log((1.5) * delta / z0) - psi_m);
                }

                // prevent some unphysical math
                if (!(ustar > zero && !std::isnan(ustar))) { ustar = zero; }
                if (!(ustar < 2.0  && !std::isnan(ustar))) { ustar = 2.0; }
                if (psi_h_neighbor > std::log(1.5 * delta / z0)) { psi_h_neighbor = std::log(1.5 * delta / z0); }
                if (psi_h > std::log(myhalf * delta / z0)) { psi_h = std::log(myhalf * delta / z0); }

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * delta / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * delta / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
            }
        }

        // Force fully immersed cells to planar average rho and theta
        if (t_blank == 1.0) {
            const Real rho_avg = r_avg(k);
            const Real theta_avg = t_avg(k) / rho_avg;  // Convert from RhoTheta to Theta
            const Real rho_cell = cell_data(i,j,k,Rho_comp);
            const Real theta_cell = cell_data(i,j,k,RhoTheta_comp) / rho_cell;
            const Real bc_forcing_r = -(rho_avg - rho_cell);
            const Real bc_forcing_rt = -(rho_avg * theta_avg - cell_data(i,j,k,RhoTheta_comp));

            cell_src(i, j, k, Rho_comp) -= drag_coefficient * U_s * bc_forcing_r;
            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
        }
    });
}
