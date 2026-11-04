#include "ERF_ImmersedForcing.H"
#include "ERF_Constants.H"
#include "ERF_TI_slow_headers.H"
#include "ERF_SrcHeaders.H"
#include "ERF_ImmersedWallStability.H"
#include "ERF_ImmersedWallCell.H"
#include "ERF_ImmersedBuildingMasks.H"

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
 * @param[in] wind_floor Wind-speed floor [m/s] behind the friction velocity of the stability estimate (erf.if_stability_wind_floor).
 * @param[in] psi_cap_factor Cap on psi_m as a fraction of ln(z/z0) (erf.if_psi_cap_factor).
 * @return Target velocity component.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
compute_if_most_target_vel (
    const amrex::Real u1_2r,
    const amrex::Real u2_2r,
    const amrex::Real delta,
    const amrex::Real z0,
    const amrex::Real t_blank,
    const amrex::Real theta_face,
    const amrex::Real theta_surf,
    const amrex::Real tflux_in,
    const amrex::Real Olen_in,
    const bool        stability_correction,
    const amrex::Real wind_floor,
    const amrex::Real psi_cap_factor
)
{
    const Real tiny             = std::numeric_limits<amrex::Real>::epsilon();
    Real psi_m                  = zero;
    Real psi_h                  = zero;
    Real tang_windspeed2r       = std::sqrt(u1_2r * u1_2r + u2_2r * u2_2r);

    // The stability estimate uses the wind floored at wind_floor, a derived Obukhov length bounded
    // so that |zeta| <= ib_stability::zeta_max() as on flat ground, and a prescribed one as given
    // with its zeta held there; the target below uses the actual wind
    Real ustar = ib_stability::floored_wind(tang_windspeed2r, wind_floor) * KAPPA / (std::log(1.5 * delta / z0) - psi_m);
    Real tflux = (tflux_in != Real(1.e-8)) ? tflux_in : -(theta_face - theta_surf) * ustar * KAPPA / (std::log(1.5 * delta / z0) - psi_h);
    Real Olen  = (Olen_in != Real(1.e-8)) ? Olen_in
               : ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta_face / (KAPPA * CONST_GRAV * tflux + tiny), Real(1.5) * delta);
    Real zeta  = ib_stability::bounded_zeta(Real(1.5) * delta, Olen);

    // similarity functions
    similarity_funs sfuns;
    if (stability_correction){
        psi_m          = sfuns.calc_psi_m(zeta);
        psi_h          = sfuns.calc_psi_h(zeta);
    }
    ustar = tang_windspeed2r * KAPPA / (std::log(1.5 * delta / z0) - ib_stability::psi_m_for_ustar(psi_m, Real(1.5) * delta, z0, psi_cap_factor));

    // prevent some unphysical math
    ustar = ib_stability::clamped_ustar(ustar);
    psi_m = ib_stability::capped_psi(psi_m, myhalf * delta, z0, psi_cap_factor);

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
    const auto dx_arr = geom.CellSizeArray();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s = one; // unit velocity scale
    const bool l_implicit_drag = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the drag of the solid is applied inside the anelastic projection
    const bool l_drag = !solverChoice.if_implicit_projection;

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg        = CONST_GRAV;
    const Real kappa      = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux_in           = solverChoice.if_surf_temp_flux;
    const Real Olen_in            = solverChoice.if_Olen_in;
    const Real stab_wind_floor    = solverChoice.if_stability_wind_floor;
    const Real psi_cap_factor     = solverChoice.if_psi_cap_factor;
    const bool l_use_most         = solverChoice.if_use_most;
    const bool l_stability_correction = solverChoice.if_stability_correction;

    const Real small_volfrac = solverChoice.if_small_volfrac;
    // The terrain kernels keep the raw fractions: erf.if_snap_partial_cells
    // applies to the buildings kernels only (the wall law below is weighted
    // by the fluid fraction of the face, which a snapped face has not).

    ParallelFor(tbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real ux = u(i, j, k);
        const Real uy = fourth * ( v(i, j  , k  ) + v(i-1, j  , k  )
                                 + v(i, j+1, k  ) + v(i-1, j+1, k  ) );
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

        // With erf.if_implicit_projection a fully solid face is held by the drag in the projection
        // (the explicit path zeroes it after the projection instead), so the wall law is not applied there
        if ((t_blank > 0 && (t_blank_above == zero)) && l_use_most && (l_drag || t_blank < one)) { // force to MOST value
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
            // The stability estimate uses the wind floored at erf.if_stability_wind_floor, a derived
            // Obukhov length bounded so that |zeta| <= ib_stability::zeta_max() as on flat ground,
            // and a prescribed one as given with its zeta held there; the target uses the actual wind
            Real ustar = ib_stability::floored_wind(h_windspeed2r, stab_wind_floor) * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
            Real tflux = (tflux_in != Real(1e-8)) ? tflux_in : -(theta_xface - theta_surf) * ustar * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_h);
            Real Olen  = (Olen_in  != Real(1e-8)) ? Olen_in
                       : ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta_xface / (kappa * ggg * tflux + tiny), Real(1.5) * dx_z);
            Real zeta  = ib_stability::bounded_zeta(Real(1.5) * dx_z, Olen);

            // similarity functions (off by default: the plain log law)
            if (l_stability_correction) {
                psi_m          = sfuns.calc_psi_m(zeta);
                psi_h          = sfuns.calc_psi_h(zeta);
            }
            ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - ib_stability::psi_m_for_ustar(psi_m, Real(1.5) * dx_z, z0, psi_cap_factor));

            // prevent some unphysical math
            ustar = ib_stability::clamped_ustar(ustar);
            psi_m = ib_stability::capped_psi(psi_m, myhalf * dx_z, z0, psi_cap_factor);

            // determine target velocity
            const Real uTarget  = ustar / kappa * (std::log(myhalf * dx_z / z0) - psi_m);
            Real uxTarget = uTarget * ux2r / (tiny + h_windspeed2r);
            const Real bc_forcing_x = -(uxTarget - ux); // BC forcing pushes nonrelative velocity toward target velocity
            const Real lambda = (1-t_blank) * CdM * U_s; // affine relaxation rate toward MOST target [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            xmom_src_arr(i, j, k) -= fac_local * rho_xface * bc_forcing_x; // if Vf low, force more strongly to MOST. If high, less forcing.
        } else if (l_drag) {
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
    const auto dx_arr = geom.CellSizeArray();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const Real U_s = one; // unit velocity scale
    const bool l_implicit_drag = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the drag of the solid is applied inside the anelastic projection
    const bool l_drag = !solverChoice.if_implicit_projection;

    // MOST parameters
    similarity_funs sfuns;
    const Real ggg        = CONST_GRAV;
    const Real kappa      = KAPPA;
    const Real z0                 = solverChoice.if_z0;
    const Real tflux_in           = solverChoice.if_surf_temp_flux;
    const Real Olen_in            = solverChoice.if_Olen_in;
    const Real stab_wind_floor    = solverChoice.if_stability_wind_floor;
    const Real psi_cap_factor     = solverChoice.if_psi_cap_factor;
    const bool l_use_most         = solverChoice.if_use_most;
    const bool l_stability_correction = solverChoice.if_stability_correction;

    const Real small_volfrac = solverChoice.if_small_volfrac;
    // The terrain kernels keep the raw fractions: erf.if_snap_partial_cells
    // applies to the buildings kernels only (the wall law below is weighted
    // by the fluid fraction of the face, which a snapped face has not).

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

        // With erf.if_implicit_projection a fully solid face is held by the drag in the projection
        // (the explicit path zeroes it after the projection instead), so the wall law is not applied there
        if ((t_blank > 0 && (t_blank_above == zero)) && l_use_most && (l_drag || t_blank < one)) { // force to MOST value
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
            // The stability estimate uses the wind floored at erf.if_stability_wind_floor, a derived
            // Obukhov length bounded so that |zeta| <= ib_stability::zeta_max() as on flat ground,
            // and a prescribed one as given with its zeta held there; the target uses the actual wind
            Real ustar = ib_stability::floored_wind(h_windspeed2r, stab_wind_floor) * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_m); // calculated from bottom of cell. Maintains flexibility for different Vf values
            Real tflux = (tflux_in != Real(1e-8)) ? tflux_in : -(theta_yface - theta_surf) * ustar * kappa / (std::log(Real(1.5) * dx_z / z0) - psi_h);
            Real Olen  = (Olen_in  != Real(1e-8)) ? Olen_in
                       : ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta_yface / (kappa * ggg * tflux + tiny), Real(1.5) * dx_z);
            Real zeta  = ib_stability::bounded_zeta(Real(1.5) * dx_z, Olen);

            // similarity functions (off by default: the plain log law)
            if (l_stability_correction) {
                psi_m          = sfuns.calc_psi_m(zeta);
                psi_h          = sfuns.calc_psi_h(zeta);
            }
            ustar = h_windspeed2r * kappa / (std::log(Real(1.5) * dx_z / z0) - ib_stability::psi_m_for_ustar(psi_m, Real(1.5) * dx_z, z0, psi_cap_factor));

            // prevent some unphysical math
            ustar = ib_stability::clamped_ustar(ustar);
            psi_m = ib_stability::capped_psi(psi_m, myhalf * dx_z, z0, psi_cap_factor);

            // determine target velocity
            const Real uTarget  = ustar / kappa * (std::log(myhalf * dx_z / z0) - psi_m);
            Real uyTarget = uTarget * uy2r / (tiny + h_windspeed2r);
            const Real bc_forcing_y = -(uyTarget - uy);  // BC forcing pushes nonrelative velocity toward target velocity
            const Real lambda = (1 - t_blank) * CdM * U_s; // affine relaxation rate toward MOST target [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            ymom_src_arr(i, j, k) -= fac_local * rho_yface * bc_forcing_y; // if Vf low, force more strongly to MOST. If high, less forcing.
        } else if (l_drag) {
            const Real lambda = t_blank * CdM * windspeed; // linear drag rate [1/s]
            const Real fac_local    = l_implicit_drag ? lambda / (one + lambda*dt) : lambda; // point-implicit rescale (else explicit)
            ymom_src_arr(i, j, k) -= fac_local * rho_yface * uy;
        }
    });
}

/**
 * Surface condition of the fraction-stress wall law from the erf.if_* inputs at this time: a given
 * Obukhov length, a given heat flux, a given surface temperature (erf.if_init_surf_temp, changing
 * at erf.if_surf_heating_rate, which is stored in K/s), or none (neutral). The inputs guarantee
 * that at most one is set.
 */
static ib_wall::WallCond
make_wall_cond (const SolverChoice& sc, const Real time)
{
    ib_wall::WallCond c;
    if (sc.if_Olen_in != Real(1e-8)) {
        c.type = ib_wall::WallCond::obukhov;
        c.Linv = one / sc.if_Olen_in;
    } else if (sc.if_surf_temp_flux != Real(1e-8)) {
        c.type = ib_wall::WallCond::heat_flux;
        c.q = sc.if_surf_temp_flux;
    } else if (sc.if_init_surf_temp > zero) {
        c.type = ib_wall::WallCond::surface_temp;
        c.theta_s = sc.if_init_surf_temp + sc.if_surf_heating_rate * time;
    }
    return c;
}

namespace {

/// Solid fraction on a face column at level kk, clamped to the domain in k, fluid below small
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real if_face_beta (const Array4<const Real>& tb, const Array4<const Real>& tb_face,
                   int i, int j, int kk, int io, int jo, int klo, int khi, Real small) noexcept
{
    const int kc = amrex::min(amrex::max(kk, klo), khi);
    const Real b = (tb_face) ? tb_face(i, j, kc) : myhalf * (tb(i, j, kc) + tb(i-io, j-jo, kc));
    return (b < small) ? zero : b;
}

/// Solid fraction of a cell column at level kk, clamped to the domain in k, fluid below small
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real if_cell_beta (const Array4<const Real>& tb, int i, int j, int kk, int klo, int khi, Real small) noexcept
{
    const Real b = tb(i, j, amrex::min(amrex::max(kk, klo), khi));
    return (b < small) ? zero : b;
}

/// Horizontal speed on a face of direction dir at level kk (the cross component averaged to it)
template <int dir>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real if_face_speed (const Array4<const Real>& u, const Array4<const Real>& v, int i, int j, int kk) noexcept
{
    const Real un = (dir == 0) ? u(i, j, kk) : v(i, j, kk);
    const Real uc = (dir == 0) ? fourth * ( v(i, j, kk) + v(i-1, j, kk) + v(i, j+1, kk) + v(i-1, j+1, kk) )
                               : fourth * ( u(i, j, kk) + u(i+1, j, kk) + u(i, j-1, kk) + u(i+1, j-1, kk) );
    return std::sqrt(un * un + uc * uc);
}

/// Potential temperature on a face of direction dir at level kk
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real if_face_theta (const Array4<const Real>& cell_data, int i, int j, int kk, int io, int jo) noexcept
{
    return myhalf * ( cell_data(i,j,kk,RhoTheta_comp) / cell_data(i,j,kk,Rho_comp)
                    + cell_data(i-io,j-jo,kk,RhoTheta_comp) / cell_data(i-io,j-jo,kk,Rho_comp) );
}

} // namespace

/**
 * Fraction-stress wall law (erf.if_wall_form = fraction_stress) for a horizontal momentum
 * component, see ERF_ImmersedWallCell.H: the wall stress in the wall cell of each face column,
 * the drag of the terrain kernels elsewhere (solid cells, and partial cells under the wall cell).
 * dir = 0 forces x-momentum, dir = 1 y-momentum.
 */
template <int dir>
void ImmersedForcingTerrain_HorizMom_FractionStress (const Box& tb,
                                                    const Array4<const Real>& u,
                                                    const Array4<const Real>& v,
                                                    const Array4<const Real>& w,
                                                    const Array4<const Real>& cell_data,
                                                    const Array4<const Real>& t_blank_arr,
                                                    const Array4<const Real>& t_blank_face_arr,
                                                    const Array4<      Real>& mom_src_arr,
                                                    const Geometry& geom,
                                                    const SolverChoice& solverChoice,
                                                    const Real dt,
                                                    const Real time)
{
    const Real* dx_arr = geom.CellSize();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dz   = dx_arr[2];   // fraction_stress runs on a constant-dz mesh
    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real z0      = solverChoice.if_z0;
    const bool l_implicit_drag = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the drag of the solid is applied inside the anelastic projection
    const bool l_drag = !solverChoice.if_implicit_projection;
    const Real tiny  = std::numeric_limits<Real>::epsilon();
    const Real small = solverChoice.if_small_volfrac;
    const ib_wall::WallCond cond = make_wall_cond(solverChoice, time);
    constexpr int io = (dir == 0) ? 1 : 0;
    constexpr int jo = (dir == 1) ? 1 : 0;
    const int klo = geom.Domain().smallEnd(2);
    const int khi = geom.Domain().bigEnd(2);

    ParallelFor(tb, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // solid fraction on this face column, as the terrain kernels read it, with the column
        // ends clamped to the domain as make_ib_wall_face_masks does
        const Real b_km1 = if_face_beta(t_blank_arr, t_blank_face_arr, i, j, k-1, io, jo, klo, khi, small);
        const Real b_k   = if_face_beta(t_blank_arr, t_blank_face_arr, i, j, k,   io, jo, klo, khi, small);
        const Real b_kp1 = if_face_beta(t_blank_arr, t_blank_face_arr, i, j, k+1, io, jo, klo, khi, small);
        const Real rho_face = myhalf * ( cell_data(i,j,k,Rho_comp) + cell_data(i-io,j-jo,k,Rho_comp) );
        const Real un = (dir == 0) ? u(i, j, k) : v(i, j, k);

        if (ib_wall::is_wall_cell(b_km1, b_k, b_kp1, small)) {
            const Real ut     = if_face_speed<dir>(u, v, i, j, k);
            const Real ut_ref = if_face_speed<dir>(u, v, i, j, k+1);
            const ib_wall::WallState ws = ib_wall::wall_state(ut, ut_ref,
                                                              if_face_theta(cell_data, i, j, k,   io, jo),
                                                              if_face_theta(cell_data, i, j, k+1, io, jo),
                                                              ib_wall::wall_beta(b_k, small), dz, z0, cond);
            mom_src_arr(i, j, k) -= ib_wall::stress_rate(ws.ustar, ut, dz, dt) * rho_face * un;
        } else if (l_drag) {
            const Real wz = (dir == 0) ? fourth * ( w(i, j, k) + w(i-1, j, k) + w(i, j, k+1) + w(i-1, j, k+1) )
                                       : fourth * ( w(i, j, k) + w(i, j-1, k) + w(i, j, k+1) + w(i, j-1, k+1) );
            const Real uh = if_face_speed<dir>(u, v, i, j, k);
            const Real windspeed = std::sqrt(uh * uh + wz * wz);
            const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dz, one/three);
            const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);
            const Real lambda = b_k * CdM * windspeed;   // linear drag rate [1/s]; zero in fluid cells
            const Real fac_local = l_implicit_drag ? lambda / (one + lambda*dt) : lambda;
            mom_src_arr(i, j, k) -= fac_local * rho_face * un;
        }
    });
}

void ImmersedForcingTerrain_Xmom_FractionStress (const Box& tbx,
                                                 const Array4<const Real>& u,
                                                 const Array4<const Real>& v,
                                                 const Array4<const Real>& w,
                                                 const Array4<const Real>& cell_data,
                                                 const Array4<const Real>& t_blank_arr,
                                                 const Array4<const Real>& t_blank_xface_arr,
                                                 const Array4<      Real>& xmom_src_arr,
                                                 const Geometry& geom,
                                                 const SolverChoice& solverChoice,
                                                 const Real fac,
                                                 const Real time)
{
    ImmersedForcingTerrain_HorizMom_FractionStress<0>(tbx, u, v, w, cell_data, t_blank_arr, t_blank_xface_arr,
                                                      xmom_src_arr, geom, solverChoice, fac, time);
}

void ImmersedForcingTerrain_Ymom_FractionStress (const Box& tby,
                                                 const Array4<const Real>& u,
                                                 const Array4<const Real>& v,
                                                 const Array4<const Real>& w,
                                                 const Array4<const Real>& cell_data,
                                                 const Array4<const Real>& t_blank_arr,
                                                 const Array4<const Real>& t_blank_yface_arr,
                                                 const Array4<      Real>& ymom_src_arr,
                                                 const Geometry& geom,
                                                 const SolverChoice& solverChoice,
                                                 const Real fac,
                                                 const Real time)
{
    ImmersedForcingTerrain_HorizMom_FractionStress<1>(tby, u, v, w, cell_data, t_blank_arr, t_blank_yface_arr,
                                                      ymom_src_arr, geom, solverChoice, fac, time);
}

/**
 * Fraction-stress wall law for the scalars (erf.if_wall_form = fraction_stress), see
 * ERF_ImmersedWallCell.H: the wall heat flux -u* theta* in the wall cell of each column, with the
 * u* and theta* of the momentum law, and the relaxation of fully immersed cells to the planar
 * average of the legacy kernel. The kinematic wall flux -u* theta* is stored in wall_hfx at the
 * wall cell (zero elsewhere), for the TKE buoyancy on the wall face; with wall_tke the wall value
 * of the TKE is stored in wall_ktarget at the wall cell (-1 elsewhere), and a time factor of 1
 * leaves the pin to erf_slow_rhs_post and the implicit TKE solve.
 */
void ImmersedForcingTerrain_Scalar_FractionStress (const Box& bx,
                                                   const Array4<const Real>& u,
                                                   const Array4<const Real>& v,
                                                   const Array4<const Real>& cell_data,
                                                   const Array4<const Real>& t_blank_arr,
                                                   const Array4<      Real>& cell_src,
                                                   const Array4<      Real>& wall_hfx,
                                                   const Array4<      Real>& wall_ktarget,
                                                   const Geometry& geom,
                                                   const SolverChoice& solverChoice,
                                                   const Table1D<Real>& r_avg,
                                                   const Table1D<Real>& t_avg,
                                                   const Real dt,
                                                   const Real time,
                                                   const bool wall_tke,
                                                   const Real Cmu0)
{
    const Real* dx_arr = geom.CellSize();
    const Real dz = dx_arr[2];   // fraction_stress runs on a constant-dz mesh
    const Real tke_time_factor = solverChoice.if_wall_tke_time_factor;
    // a time factor of 1 pins the wall cell at the wall value after every stage (erf_slow_rhs_post
    // and the implicit TKE solve read wall_ktarget), as erf.dirichlet_k does on flat ground;
    // larger factors relax toward it here as a source
    const bool tke_pin = (tke_time_factor <= one);
    const Real drag_coefficient = solverChoice.if_Cd_scalar / std::pow(dx_arr[0]*dx_arr[1]*dz, one/three);
    const Real z0    = solverChoice.if_z0;
    const Real small = solverChoice.if_small_volfrac;
    const ib_wall::WallCond cond = make_wall_cond(solverChoice, time);
    const int klo = geom.Domain().smallEnd(2);
    const int khi = geom.Domain().bigEnd(2);
    const bool have_avg = (r_avg && t_avg);

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real b_k = if_cell_beta(t_blank_arr, i, j, k, klo, khi, small);
        if (wall_hfx) { wall_hfx(i, j, k) = zero; }
        if (wall_ktarget) { wall_ktarget(i, j, k) = -one; }   // no wall value: not a wall cell

        if (ib_wall::is_wall_cell(if_cell_beta(t_blank_arr, i, j, k-1, klo, khi, small), b_k,
                                  if_cell_beta(t_blank_arr, i, j, k+1, klo, khi, small), small)) {
            // horizontal speed at the cell centres of the wall cell and the cell above
            const Real uc0 = myhalf * (u(i, j, k)   + u(i+1, j, k));
            const Real vc0 = myhalf * (v(i, j, k)   + v(i, j+1, k));
            const Real uc1 = myhalf * (u(i, j, k+1) + u(i+1, j, k+1));
            const Real vc1 = myhalf * (v(i, j, k+1) + v(i, j+1, k+1));
            const Real rho       = cell_data(i, j, k,   Rho_comp);
            const Real theta     = cell_data(i, j, k,   RhoTheta_comp) / rho;
            const Real theta_ref = cell_data(i, j, k+1, RhoTheta_comp) / cell_data(i, j, k+1, Rho_comp);
            const ib_wall::WallState ws = ib_wall::wall_state(std::sqrt(uc0*uc0 + vc0*vc0), std::sqrt(uc1*uc1 + vc1*vc1),
                                                              theta, theta_ref,
                                                              ib_wall::wall_beta(b_k, small), dz, z0, cond);
            cell_src(i, j, k, RhoTheta_comp) += rho * ib_wall::heat_source(ws, theta, cond, dz, dt);
            if (wall_hfx) { wall_hfx(i, j, k) = -ws.ustar * ws.thetastar; }
            if (wall_tke) {
                const Real k_wall = ib_wall::wall_tke(ws, ib_wall::wall_beta(b_k, small), dz, z0, theta, Cmu0);
                if (wall_ktarget) { wall_ktarget(i, j, k) = k_wall; }
                if (!tke_pin) {
                    // relax the wall-cell TKE toward its wall value over tke_time_factor steps
                    const Real k_now  = cell_data(i, j, k, RhoKE_comp) / rho;
                    cell_src(i, j, k, RhoKE_comp) += rho * (k_wall - k_now) / (tke_time_factor * dt);
                }
            }
        } else if (b_k == one && have_avg) {
            // fully immersed cells: relax to the planar average, as the legacy kernel does
            const Real rho_avg   = r_avg(k);
            const Real theta_avg = t_avg(k) / rho_avg;
            cell_src(i, j, k, Rho_comp)      -= drag_coefficient * (cell_data(i,j,k,Rho_comp) - rho_avg);
            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * (cell_data(i,j,k,RhoTheta_comp) - rho_avg * theta_avg);
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
    const auto dx_arr = geom.CellSizeArray();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];
    const Real dt = fac;  // fac is actually dt in the calling code

    const Real alpha_m = solverChoice.if_Cd_momentum;
    const Real tiny = std::numeric_limits<amrex::Real>::epsilon();
    const bool l_implicit_drag = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the drag of the solid is applied inside the anelastic projection
    const bool l_drag = !solverChoice.if_implicit_projection;

    const Real small_volfrac = solverChoice.if_small_volfrac;
    // The terrain kernels keep the raw fractions: erf.if_snap_partial_cells
    // applies to the buildings kernels only (the wall law below is weighted
    // by the fluid fraction of the face, which a snapped face has not).

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
        if (l_drag) { zmom_src_arr(i, j, k) -= fac_local * rho_zface * uz; }
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
    const auto dx_arr = geom.CellSizeArray();
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
    const Real stab_wind_floor         = solverChoice.if_stability_wind_floor;
    const Real psi_cap_factor          = solverChoice.if_psi_cap_factor;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the linear drag of the solid is applied inside the anelastic
    // projection (ERF::make_if_projection_sigma), not here
    const bool l_ip               = solverChoice.if_implicit_projection;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = solverChoice.if_small_volfrac;
    // erf.if_snap_partial_cells: read the cell blanking snapped to solid (1)
    // or fluid (0) at half, so a height-map building becomes the same
    // staircase of whole cells an exact box is. A face is solid when either
    // cell it joins is; a face between a solid and a fluid cell is
    // wall-normal and gets the interior drag toward zero (no penetration); a
    // face between two solid cells carries the roof or wall law of its row
    // (the full log-law target, not the partial-cell weighted one) or the
    // interior drag, never both; the partial-cell branches (wall_mask,
    // east_west_mask and the like, faces with 0 < t_blank < 1) do not arise
    // under the snap. Off (the default), the raw fractions
    // are used and nothing below changes: a boundary face of an exact box
    // carries the wall law and the interior drag together as before. See
    // SolverChoice::if_snap_partial_cells.
    const bool l_snap = solverChoice.if_snap_partial_cells;
    // The face blanking and its neighbours are read through ERF_ImmersedBuildingMasks.H, as the
    // implicit projection reads them
    const int dom_klo = geom.Domain().smallEnd(2);


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

        // the face, its neighbours and the raw blanking the explicit path zeroes on (ERF_ImmersedBuildingMasks.H)
        const if_bld::HorizFace fc = if_bld::horiz_face<0>(t_blank_arr, t_blank_xface_arr, i, j, k, dom_klo, l_snap, small_volfrac);
        const Real t_blank       = fc.t;
        // With the snap the boundary solid face stands for the wall layer and
        // takes the full log-law target; the partial-cell weight (1 - t_blank)
        // of compute_if_most_target_vel() would make it zero (a no-slip
        // staircase), so the snap path passes a zero blanking to it.
        const Real t_blank_law = l_snap ? zero : t_blank;
        const Real t_blank_below = fc.below;
        const Real t_blank_above = fc.above;
        const Real t_blank_south = fc.lo;
        const Real t_blank_north = fc.hi;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);
        // With the snap a face joining a solid and a fluid cell is wall-normal:
        // it gets the interior drag (no penetration), not a wall law, whatever
        // row it lies in.
        const bool normal_face = fc.normal;

        // With the snap every solid face has t_blank = 1: a roof face lies in
        // the top solid row (t_blank <= t_blank_below, the face above fluid),
        // and a face carrying a wall law is not also an interior face.
        // roof, south/north wall law, partial-wall and interior drag (ERF_ImmersedBuildingMasks.H)
        if_bld::HorizMasks mk = if_bld::horiz_masks(t_blank, t_blank_below, t_blank_above, t_blank_south, t_blank_north,
                                                    normal_face, l_use_most, l_snap);
        // erf.if_implicit_projection: the linear drag is in the projection, and a face the explicit
        // path would zero after the projection (solid on both sides) is held there, with no wall law
        if (l_ip) {
            mk.wall = zero; mk.side = zero; mk.interior = zero;
            if (fc.raw == one) { mk.roof = zero; mk.lo = zero; mk.hi = zero; }
        }
        const Real roof_mask      = mk.roof;
        const Real south_mask     = mk.lo;
        const Real north_mask     = mk.hi;
        const Real wall_mask      = mk.wall;
        const Real east_west_mask = mk.side;
        const Real interior_mask  = mk.interior;

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_xface_inside = rho_xface;
        Real theta_surf       = theta_xface;
        Real bc_forcing_x     = zero;
        Real u_target         = zero;

        // roof forcing
        if (roof_mask == one) {
            u1_cellaway         = u(i, j, k+1) ;
            u2_cellaway         = fourth * ( v(i, j  , k+1) + v(i-1, j  , k+1)
                                           + v(i, j+1, k+1) + v(i-1, j+1, k+1) ) ;
            rho_xface_inside    =  myhalf * (cell_data(i,j,k-1,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp));
            theta_surf          = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i-1,j,k-1, RhoTheta_comp))) / rho_xface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_z, z0, t_blank_law, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_x * roof_mask * rho_xface * CdM * U_s;
        }

        // south wall forcing
        if (south_mask == one) {
            u1_cellaway         = u(i, j-1, k  );
            u2_cellaway         = fourth * ( w(i, j-1, k  ) + w(i-1, j-1, k  )
                                           + w(i, j-1, k+1) + w(i-1, j-1, k+1) ) ;
            rho_xface_inside    = myhalf * ( cell_data(i,j+1,k,Rho_comp) + cell_data(i-1,j+1,k,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i,j+1,k,RhoTheta_comp) + cell_data(i-1,j+1,k, RhoTheta_comp))) / rho_xface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank_law, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_x * south_mask * rho_xface * CdM * U_s;
        }

        // north wall forcing
        if (north_mask == one) {
            u1_cellaway         = u(i, j+1, k  ) ;
            u2_cellaway         = fourth * ( w(i, j+1, k  ) + w(i-1, j+1, k  )
                                           + w(i, j+1, k+1) + w(i-1, j+1, k+1) ) ;
            rho_xface_inside    = myhalf * ( cell_data(i,j-1,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i,j-1,k,RhoTheta_comp) + cell_data(i-1,j-1,k, RhoTheta_comp))) / rho_xface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank_law, theta_xface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_x        = -(u_target - ux); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_x * north_mask * rho_xface * CdM * U_s;
        }

        // wall forcing (if not using most) or east/west walls when using MOST
        if (wall_mask == one || east_west_mask == one) {
            drag               += (wall_mask + east_west_mask) * t_blank * rho_xface * CdM * ux * windspeed;
        }

        // interior cell forcing
        if (interior_mask == one) {
            drag               += interior_mask * rho_xface * CdM * ux * windspeed;
        }

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
    const auto dx_arr = geom.CellSizeArray();
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
    const Real stab_wind_floor         = solverChoice.if_stability_wind_floor;
    const Real psi_cap_factor          = solverChoice.if_psi_cap_factor;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the linear drag of the solid is applied inside the anelastic
    // projection (ERF::make_if_projection_sigma), not here
    const bool l_ip               = solverChoice.if_implicit_projection;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = solverChoice.if_small_volfrac;
    // erf.if_snap_partial_cells: read the cell blanking snapped to solid (1)
    // or fluid (0) at half, so a height-map building becomes the same
    // staircase of whole cells an exact box is. A face is solid when either
    // cell it joins is; a face between a solid and a fluid cell is
    // wall-normal and gets the interior drag toward zero (no penetration); a
    // face between two solid cells carries the roof or wall law of its row
    // (the full log-law target, not the partial-cell weighted one) or the
    // interior drag, never both; the partial-cell branches (wall_mask,
    // east_west_mask and the like, faces with 0 < t_blank < 1) do not arise
    // under the snap. Off (the default), the raw fractions
    // are used and nothing below changes: a boundary face of an exact box
    // carries the wall law and the interior drag together as before. See
    // SolverChoice::if_snap_partial_cells.
    const bool l_snap = solverChoice.if_snap_partial_cells;
    // The face blanking and its neighbours are read through ERF_ImmersedBuildingMasks.H, as the
    // implicit projection reads them
    const int dom_klo = geom.Domain().smallEnd(2);


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

        // the face, its neighbours and the raw blanking the explicit path zeroes on (ERF_ImmersedBuildingMasks.H)
        const if_bld::HorizFace fc = if_bld::horiz_face<1>(t_blank_arr, t_blank_yface_arr, i, j, k, dom_klo, l_snap, small_volfrac);
        const Real t_blank       = fc.t;
        // With the snap the boundary solid face stands for the wall layer and
        // takes the full log-law target; the partial-cell weight (1 - t_blank)
        // of compute_if_most_target_vel() would make it zero (a no-slip
        // staircase), so the snap path passes a zero blanking to it.
        const Real t_blank_law = l_snap ? zero : t_blank;
        const Real t_blank_below = fc.below;
        const Real t_blank_above = fc.above;
        const Real t_blank_west  = fc.lo;
        const Real t_blank_east  = fc.hi;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);
        // With the snap a face joining a solid and a fluid cell is wall-normal:
        // it gets the interior drag (no penetration), not a wall law, whatever
        // row it lies in.
        const bool normal_face = fc.normal;

        // As in the x-momentum: with the snap a roof face lies in the top
        // solid row and a wall-law face is not also an interior face.
        // roof, west/east wall law, partial-wall and interior drag (ERF_ImmersedBuildingMasks.H)
        if_bld::HorizMasks mk = if_bld::horiz_masks(t_blank, t_blank_below, t_blank_above, t_blank_west, t_blank_east,
                                                    normal_face, l_use_most, l_snap);
        // erf.if_implicit_projection: as in the x-momentum
        if (l_ip) {
            mk.wall = zero; mk.side = zero; mk.interior = zero;
            if (fc.raw == one) { mk.roof = zero; mk.lo = zero; mk.hi = zero; }
        }
        const Real roof_mask        = mk.roof;
        const Real west_mask        = mk.lo;
        const Real east_mask        = mk.hi;
        const Real wall_mask        = mk.wall;
        const Real north_south_mask = mk.side;
        const Real interior_mask    = mk.interior;

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_yface_inside = rho_yface;
        Real theta_surf       = theta_yface;
        Real bc_forcing_y     = zero;
        Real u_target         = zero;

        // roof forcing
        if (roof_mask == one) {
            u1_cellaway         = fourth * ( u(i  , j  , k+1) + u(i  , j-1, k+1)
                                           + u(i+1, j  , k+1) + u(i+1, j-1, k+1) );
            u2_cellaway         = v(i, j, k+1);
            rho_yface_inside    = myhalf * ( cell_data(i,j,k-1,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i,j,k-1,RhoTheta_comp) + cell_data(i,j-1,k-1,RhoTheta_comp))) / rho_yface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_z, z0, t_blank_law, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_y * roof_mask * rho_yface * CdM * U_s;
        }

        // west wall forcing
        if (west_mask == one) {
            u1_cellaway         = v(i-1, j , k  );
            u2_cellaway         = fourth * ( w(i-1, j  , k  ) + w(i-1, j-1, k  )
                                           + w(i-1, j  , k+1) + w(i-1, j-1, k+1) );
            rho_yface_inside    = myhalf * ( cell_data(i+1,j,k,Rho_comp) + cell_data(i+1,j-1,k,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i+1,j,k,RhoTheta_comp) + cell_data(i+1,j-1,k,RhoTheta_comp))) / rho_yface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank_law, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_y * west_mask * rho_yface * CdM * U_s;
        }

        // east wall forcing
        if (east_mask == one) {
            u1_cellaway         = v(i+1, j , k  );
            u2_cellaway         = fourth * ( w(i+1, j  , k  ) + w(i+1, j-1, k  )
                                           + w(i+1, j  , k+1) + w(i+1, j-1, k+1) );
            rho_yface_inside    = myhalf * ( cell_data(i-1,j,k,Rho_comp) + cell_data(i-1,j-1,k,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i-1,j,k,RhoTheta_comp) + cell_data(i-1,j-1,k,RhoTheta_comp))) / rho_yface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank_law, theta_yface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_y        = -(u_target - uy); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_y * east_mask * rho_yface * CdM * U_s;
        }

        // wall forcing (if not using most) or north/south walls when using MOST
        if (wall_mask == one || north_south_mask == one) {
            drag               += (wall_mask + north_south_mask) * t_blank * rho_yface * CdM * uy * windspeed;
        }

        // interior cell forcing
        if (interior_mask == one) {
            drag               += interior_mask * rho_yface * CdM * uy * windspeed;
        }

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
    const auto dx_arr = geom.CellSizeArray();
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
    const Real stab_wind_floor         = solverChoice.if_stability_wind_floor;
    const Real psi_cap_factor          = solverChoice.if_psi_cap_factor;
    const bool l_use_most              = solverChoice.if_use_most;
    const bool l_stability_correction  = solverChoice.if_stability_correction;

    // To limit stiffness of drag when using anelastic
    const Real ws_floor           = solverChoice.if_ws_floor;
    const Real damp_alpha         = solverChoice.if_damp_alpha;
    // Point-implicit alternative to the clamp above; stabilizes both compressible and anelastic
    const bool l_implicit_drag    = solverChoice.if_implicit_drag;
    // erf.if_implicit_projection: the linear drag of the solid is applied inside the anelastic
    // projection (ERF::make_if_projection_sigma), not here
    const bool l_ip               = solverChoice.if_implicit_projection;

    const bool is_slow_step = true;  // This is determined by calling context
    const bool use_ImmersedForcing_fast = solverChoice.immersed_forcing_substep;
    const Real small_volfrac = solverChoice.if_small_volfrac;
    // erf.if_snap_partial_cells: read the cell blanking snapped to solid (1)
    // or fluid (0) at half, so a height-map building becomes the same
    // staircase of whole cells an exact box is. A face is solid when either
    // cell it joins is; a face between a solid and a fluid cell is
    // wall-normal and gets the interior drag toward zero (no penetration); a
    // face between two solid cells carries the roof or wall law of its row
    // (the full log-law target, not the partial-cell weighted one) or the
    // interior drag, never both; the partial-cell branches (wall_mask,
    // east_west_mask and the like, faces with 0 < t_blank < 1) do not arise
    // under the snap. Off (the default), the raw fractions
    // are used and nothing below changes: a boundary face of an exact box
    // carries the wall law and the interior drag together as before. See
    // SolverChoice::if_snap_partial_cells.
    const bool l_snap = solverChoice.if_snap_partial_cells;
    // The face blanking and its neighbours are read through ERF_ImmersedBuildingMasks.H, as the
    // implicit projection reads them


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

        // the face, its neighbours and the raw blanking the explicit path zeroes on (ERF_ImmersedBuildingMasks.H)
        const if_bld::VertFace fc = if_bld::vert_face(t_blank_arr, t_blank_zface_arr, i, j, k, l_snap, small_volfrac);
        const Real t_blank       = fc.t;
        // With the snap the boundary solid face stands for the wall layer and
        // takes the full log-law target; the partial-cell weight (1 - t_blank)
        // of compute_if_most_target_vel() would make it zero (a no-slip
        // staircase), so the snap path passes a zero blanking to it.
        const Real t_blank_law = l_snap ? zero : t_blank;
        const Real t_blank_above = fc.above;
        const Real t_blank_south = fc.south;
        const Real t_blank_north = fc.north;
        const Real t_blank_west  = fc.west;
        const Real t_blank_east  = fc.east;

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_m / std::pow(dx_x*dx_y*dx_z, one/three);
        const Real CdM = std::min(drag_coefficient / (windspeed + tiny), drag_coefficient);
        // With the snap a face joining a solid and a fluid cell is wall-normal:
        // it gets the interior drag (no penetration), not a wall law, whatever
        // row it lies in.
        const bool normal_face = fc.normal;

        // wall law on the four sides, partial-wall, roof and interior drag (ERF_ImmersedBuildingMasks.H);
        // with the snap a face carrying a wall law or the roof drag is not also an interior face
        if_bld::VertMasks mk = if_bld::vert_masks(t_blank, t_blank_above, t_blank_south, t_blank_north,
                                                  t_blank_west, t_blank_east, k >= 1, normal_face, l_use_most, l_snap);
        // erf.if_implicit_projection: as in the x-momentum
        if (l_ip) {
            mk.wall = zero; mk.roof = zero; mk.interior = zero;
            if (fc.raw == one) { mk.south = zero; mk.north = zero; mk.west = zero; mk.east = zero; }
        }
        const Real south_mask    = mk.south;
        const Real north_mask    = mk.north;
        const Real west_mask     = mk.west;
        const Real east_mask     = mk.east;
        const Real wall_mask     = mk.wall;
        const Real roof_mask     = mk.roof;
        const Real interior_mask = mk.interior;

        Real drag             = zero;
        Real u1_cellaway      = zero;
        Real u2_cellaway      = zero;
        Real rho_zface_inside = rho_zface;
        Real theta_surf       = theta_zface;
        Real bc_forcing_z     = zero;
        Real u_target         = zero;

        // south wall forcing
        if (south_mask == one) {
            u1_cellaway         = fourth * ( u(i  , j-1, k  ) + u(i+1, j-1, k  )
                                           + u(i  , j-1, k-1) + u(i+1, j-1, k-1) );
            u2_cellaway         = w(i, j-1, k);
            rho_zface_inside    = myhalf * ( cell_data(i,j+1,k,Rho_comp) + cell_data(i,j+1,k-1,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i,j+1,k,RhoTheta_comp) + cell_data(i,j+1,k-1,RhoTheta_comp))) / rho_zface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank_law, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_z * south_mask * rho_zface * CdM * U_s;
        }

        // north wall forcing
        if (north_mask == one) {
            u1_cellaway         = fourth * ( u(i  , j+1, k  ) + u(i+1, j+1, k  )
                                           + u(i  , j+1, k-1) + u(i+1, j+1, k-1) );
            u2_cellaway         = w(i, j+1, k);
            rho_zface_inside    = myhalf * ( cell_data(i,j-1,k,Rho_comp) + cell_data(i,j-1,k-1,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i,j-1,k,RhoTheta_comp) + cell_data(i,j-1,k-1,RhoTheta_comp))) / rho_zface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_y, z0, t_blank_law, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_z * north_mask * rho_zface * CdM * U_s;
        }

        // west wall forcing
        if (west_mask == one) {
            u1_cellaway         = fourth * ( v(i-1, j  , k  ) + v(i-1, j+1, k  )
                                           + v(i-1, j  , k-1) + v(i-1, j+1, k-1) );
            u2_cellaway         = w(i-1, j, k);
            rho_zface_inside    = myhalf * ( cell_data(i+1,j,k,Rho_comp) + cell_data(i+1,j,k-1,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i+1,j,k,RhoTheta_comp) + cell_data(i+1,j,k-1,RhoTheta_comp))) / rho_zface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank_law, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_z * west_mask * rho_zface * CdM * U_s;
        }

        // east wall forcing
        if (east_mask == one) {
            u1_cellaway         = fourth * ( v(i+1, j  , k  ) + v(i+1, j+1, k  )
                                           + v(i+1, j  , k-1) + v(i+1, j+1, k-1) );
            u2_cellaway         = w(i+1, j, k);
            rho_zface_inside    = myhalf * ( cell_data(i-1,j,k,Rho_comp) + cell_data(i-1,j,k-1,Rho_comp) );
            theta_surf          = (myhalf * (cell_data(i-1,j,k,RhoTheta_comp) + cell_data(i-1,j,k-1,RhoTheta_comp))) / rho_zface_inside;
            u_target            = compute_if_most_target_vel(u1_cellaway, u2_cellaway, dx_x, z0, t_blank_law, theta_zface, theta_surf, tflux_in, Olen_in, l_stability_correction, stab_wind_floor, psi_cap_factor);
            bc_forcing_z        = -(u_target - uz); // BC forcing pushes nonrelative velocity toward target velocity
            drag               += bc_forcing_z * east_mask * rho_zface * CdM * U_s;
        }

        // wall forcing (if not using most) or roof when using MOST
        if (wall_mask == one || roof_mask == one) {
            drag               += (wall_mask + roof_mask) * t_blank * rho_zface * CdM * uz * windspeed;
        }

        // interior cell forcing
        if (interior_mask == one) {
            drag               += interior_mask * rho_zface * CdM * uz * windspeed;
        }

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
    const auto dx_arr = geom.CellSizeArray();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];

    const Real alpha_h          = solverChoice.if_Cd_scalar;
    // The terrain kernels keep the raw fractions: erf.if_snap_partial_cells
    // applies to the buildings kernels only (the wall law below is weighted
    // by the fluid fraction of the face, which a snapped face has not).
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
    const Real stab_wind_floor    = solverChoice.if_stability_wind_floor;
    // erf.if_stability_correction (off by default): the similarity functions of the temperature
    // kernels, as of the momentum kernels; off, the log law (erf-model/ERF#4931)
    const bool l_stability_correction = solverChoice.if_stability_correction;

    ParallelFor(bx, [=]
                AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        const Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, one/three);

        const Real t_blank       = t_blank_arr(i, j, k);
        const Real t_blank_above = t_blank_arr(i, j, k+1);
        const Real ux_cc_2r = myhalf * (u(i  ,j  ,k+1) + u(i+1,j  ,k+1));
        const Real uy_cc_2r = myhalf * (v(i  ,j  ,k+1) + v(i  ,j+1,k+1));
        // The heat transfer uses the wind floored at erf.if_stability_wind_floor and an Obukhov
        // length bounded so that |zeta| <= ib_stability::zeta_max(), as on flat ground
        const Real h_windspeed2r  = ib_stability::floored_wind(std::sqrt(ux_cc_2r * ux_cc_2r + uy_cc_2r * uy_cc_2r), stab_wind_floor);

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
                const Real Olen  = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny),
                                                                        Real(1.5) * dx_z);
                const Real zeta          = ib_stability::bounded_zeta(myhalf * dx_z, Olen);
                const Real zeta_neighbor = ib_stability::bounded_zeta(Real(1.5) * dx_z, Olen);

                // similarity functions (off by default: the plain log law, as in the momentum kernels)
                if (l_stability_correction) {
                    psi_m          = sfuns.calc_psi_m(zeta);
                    psi_h          = sfuns.calc_psi_h(zeta);
                    psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                }
                ustar = h_windspeed2r * kappa / (std::log((Real(1.5)) * dx_z / z0) - psi_m);

                // prevent some unphysical math
                ustar = ib_stability::clamped_ustar(ustar);
                // caps of ln(z / z0): erf.if_psi_cap_factor stays off the temperature forcing, whose
                // theta* comes from an L not made consistent with u* (the legacy law, kept as it is;
                // the fraction-stress law of erf-model/ERF#4134 forms them consistently)
                // (with the flag off psi_h is zero and stays so; a cap of ln(z / z0) <= 0, for z0 >= z, is not applied)
                if (l_stability_correction) {
                    psi_h_neighbor = ib_stability::capped_psi(psi_h_neighbor, Real(1.5) * dx_z, z0, one);
                    psi_h = ib_stability::capped_psi(psi_h, myhalf * dx_z, z0, one);
                }

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
                // the prescribed length as given, with its zeta held within +-zeta_max()
                const Real Olen  = Olen_in;
                const Real zeta          = ib_stability::bounded_zeta(myhalf * dx_z, Olen);
                const Real zeta_neighbor = ib_stability::bounded_zeta(Real(1.5) * dx_z, Olen);

                // similarity functions (off by default: the log law, theta* from the given L); psi_h
                // capped at ln(z / z0) as on the other temperature branches
                const Real psi_m          = l_stability_correction ? sfuns.calc_psi_m(zeta) : zero;
                const Real psi_h          = l_stability_correction ? ib_stability::capped_psi(sfuns.calc_psi_h(zeta), myhalf * dx_z, z0, one) : zero;
                const Real psi_h_neighbor = l_stability_correction ? ib_stability::capped_psi(sfuns.calc_psi_h(zeta_neighbor), Real(1.5) * dx_z, z0, one) : zero;
                // psi_m at 0.5 dz against ln(1.5 dz / z0) can cancel; u* is held in [0, 2] m/s as on the
                // other branches. In unstable air psi_h >= psi_m and psi_h grows with |zeta|, so wherever
                // psi_m reaches ln(1.5 dz / z0) both psi_h caps above bind, the two log brackets close
                // and the target is the cell above whatever u* is: the clamp only guards the arithmetic.
                const Real ustar = ib_stability::clamped_ustar(h_windspeed2r * kappa / (std::log((Real(1.5)) * dx_z / z0) - psi_m));

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
    const auto dx_arr = geom.CellSizeArray();
    const Real dx_x = dx_arr[0];
    const Real dx_y = dx_arr[1];

    const Real alpha_h          = solverChoice.if_Cd_scalar;
    // erf.if_snap_partial_cells: read the cell blanking snapped to solid (1)
    // or fluid (0) at half, so a height-map building becomes the same
    // staircase of whole cells an exact box is: the thermal conditions sit
    // on the boundary solid cells, roofs included, found from the neighbour
    // blanking. Off (the default), the raw fractions are used and nothing
    // below changes. See SolverChoice::if_snap_partial_cells.
    const bool l_snap = solverChoice.if_snap_partial_cells;
    auto snapb = [=] AMREX_GPU_DEVICE (amrex::Real b) noexcept -> amrex::Real {
        return l_snap ? ((b >= myhalf) ? one : zero) : b;
    };
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
    const Real stab_wind_floor    = solverChoice.if_stability_wind_floor;
    // erf.if_stability_correction (off by default): the similarity functions of the temperature
    // kernels, as of the momentum kernels; off, the log law (erf-model/ERF#4931)
    const bool l_stability_correction = solverChoice.if_stability_correction;

    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        const Real t_blank       = snapb(t_blank_arr(i, j, k));
        const Real t_blank_below = snapb(t_blank_arr(i, j, k-1));
        const Real t_blank_above = snapb(t_blank_arr(i, j, k+1));
        const Real t_blank_north = snapb(t_blank_arr(i  , j+1, k));
        const Real t_blank_south = snapb(t_blank_arr(i  , j-1, k));
        const Real t_blank_east  = snapb(t_blank_arr(i+1, j  , k));
        const Real t_blank_west  = snapb(t_blank_arr(i-1, j  , k));

        const Real dx_z = (z_cc_arr) ? (z_cc_arr(i,j,k) - z_cc_arr(i,j,k-1)) : dx_arr[2];
        Real drag_coefficient = alpha_h / std::pow(dx_x*dx_y*dx_z, one/three);

        // Wall cells of the building, named by the face they carry: a
        // partially blanked cell with a more solid neighbour behind it and a
        // fluid one in front, or, with the snap (every solid cell at 1), a
        // solid cell with a solid neighbour behind and a fluid one in front.
        const bool south_face = l_snap ? (t_blank == one && t_blank_north == one && t_blank_south == zero)
                                       : (t_blank > zero && t_blank < t_blank_north && t_blank_south == zero);
        const bool north_face = l_snap ? (t_blank == one && t_blank_south == one && t_blank_north == zero)
                                       : (t_blank > zero && t_blank < t_blank_south && t_blank_north == zero);
        const bool west_face  = l_snap ? (t_blank == one && t_blank_east  == one && t_blank_west  == zero)
                                       : (t_blank > zero && t_blank < t_blank_east  && t_blank_west  == zero);
        const bool east_face  = l_snap ? (t_blank == one && t_blank_west  == one && t_blank_east  == zero)
                                       : (t_blank > zero && t_blank < t_blank_west  && t_blank_east  == zero);

        // SURFACE TEMP AND HEATING/COOLING RATE
        if (init_surf_temp > zero) {
            const Real surf_temp    = init_surf_temp + surf_heating_rate*time;
            if (t_blank > 0 && (t_blank_above == zero) && (t_blank_below == one)) { // building roof
                const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;

            } else if (east_face || west_face || south_face || north_face) {
                // this should enter for just building walls
                // walls are currently separated to allow for flexibility in the future to heat walls differently

                // south face
                if (l_snap ? south_face : ((t_blank < t_blank_north) && (t_blank_north == one))) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // north face
                if (l_snap ? north_face : ((t_blank < t_blank_south) && (t_blank_south == one))) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // west face
                if (l_snap ? west_face : ((t_blank < t_blank_east) && (t_blank_east == one))) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

                // east face
                if (l_snap ? east_face : ((t_blank < t_blank_west) && (t_blank_west == one))) {
                    const Real bc_forcing_rt_srf = -(cell_data(i,j,k,Rho_comp) * surf_temp - cell_data(i,j,k,RhoTheta_comp));
                    cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt_srf;
                }

            }
        }

        // SURFACE HEAT FLUX
        if (tflux != Real(1.e-8)){
            const Real ux_cc_2r = myhalf * (u(i  ,j  ,k+1) + u(i+1,j  ,k+1));
            const Real uy_cc_2r = myhalf * (v(i  ,j  ,k+1) + v(i  ,j+1,k+1));
            // The heat transfer uses the wind floored at erf.if_stability_wind_floor and an Obukhov
            // length bounded so that |zeta| <= ib_stability::zeta_max(), as on flat ground
            const Real h_windspeed2r  = ib_stability::floored_wind(std::sqrt(ux_cc_2r * ux_cc_2r + uy_cc_2r * uy_cc_2r), stab_wind_floor);

            const Real theta          = cell_data(i,j,k  ,RhoTheta_comp) / cell_data(i,j,k  ,Rho_comp);
            Real theta_neighbor       = cell_data(i,j,k+1,RhoTheta_comp) / cell_data(i,j,k+1,Rho_comp);

            if (t_blank > zero && (t_blank_above == zero)) { // building roof
                Real psi_m           = zero;
                Real psi_h           = zero;
                Real psi_h_neighbor  = zero;
                Real ustar           = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);
                // (start-up refuses erf.if_Olen with a surface flux, so L always comes from the flux here)
                Real Olen            = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny), Real(1.5) * dx_z);

                for (int iter = 0; iter < 2; ++iter) {
                    if (iter > 0) { Olen  = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny), Real(1.5) * dx_z); }
                    Real zeta          = ib_stability::bounded_zeta(myhalf * dx_z, Olen);
                    Real zeta_neighbor = ib_stability::bounded_zeta(Real(1.5) * dx_z, Olen);

                    // similarity functions (off by default: the plain log law, as in the momentum kernels)
                    if (l_stability_correction) {
                        psi_m          = sfuns.calc_psi_m(zeta);
                        psi_h          = sfuns.calc_psi_h(zeta);
                        psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                    }
                    ustar = h_windspeed2r * kappa / (std::log((1.5) * dx_z / z0) - psi_m);
                }

                // prevent some unphysical math
                ustar = ib_stability::clamped_ustar(ustar);
                // caps of ln(z / z0): erf.if_psi_cap_factor stays off the temperature forcing, whose
                // theta* comes from an L not made consistent with u* (the legacy law, kept as it is;
                // the fraction-stress law of erf-model/ERF#4134 forms them consistently)
                // (with the flag off psi_h is zero and stays so; a cap of ln(z / z0) <= 0, for z0 >= z, is not applied)
                if (l_stability_correction) {
                    psi_h_neighbor = ib_stability::capped_psi(psi_h_neighbor, Real(1.5) * dx_z, z0, one);
                    psi_h = ib_stability::capped_psi(psi_h, myhalf * dx_z, z0, one);
                }

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * dx_z / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * dx_z / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;

            } else if (east_face || west_face || south_face || north_face) { // this should enter for just building walls

                Real ux_cellaway = zero;
                Real uy_cellaway = zero;
                Real uz_cellaway = zero;
                Real u1          = zero;
                Real u2          = zero;
                Real delta       = zero;

                // south face
                if (south_face) {
                    ux_cellaway = myhalf * (u(i  ,j-1,k) + u(i+1,j-1,k  ));
                    uz_cellaway = myhalf * (w(i  ,j-1,k) + w(i  ,j-1,k+1));
                    u1 = ux_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_y;

                    // MOST
                    theta_neighbor = cell_data(i,j-1,k,RhoTheta_comp) / cell_data(i,j-1,k,Rho_comp);
                }

                // north face
                if (north_face) {
                    ux_cellaway = myhalf * (u(i  ,j+1,k) + u(i+1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i  ,j+1,k) + w(i  ,j+1,k+1));
                    u1 = ux_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_y;

                    // MOST
                    theta_neighbor = cell_data(i,j+1,k,RhoTheta_comp) / cell_data(i,j+1,k,Rho_comp);
                }

                // west face
                if (west_face) {
                    uy_cellaway = myhalf * (v(i-1,j  ,k) + v(i-1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i-1,j  ,k) + w(i-1,j  ,k+1));
                    u1 = uy_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_x;

                    // MOST
                    theta_neighbor = cell_data(i-1,j,k,RhoTheta_comp) / cell_data(i-1,j,k,Rho_comp);
                }

                // east face
                if (east_face) {
                    uy_cellaway = myhalf * (v(i+1,j  ,k) + v(i+1,j+1,k  ));
                    uz_cellaway = myhalf * (w(i+1,j  ,k) + w(i+1,j  ,k+1));
                    u1 = uy_cellaway;
                    u2 = uz_cellaway;
                    delta = dx_x;

                    // MOST
                    theta_neighbor = cell_data(i+1,j,k,RhoTheta_comp) / cell_data(i+1,j,k,Rho_comp);
                }

                Real tan_wspd = ib_stability::floored_wind(std::sqrt(u1 * u1 + u2 * u2), stab_wind_floor);

                Real psi_m           = zero;
                Real psi_h           = zero;
                Real psi_h_neighbor  = zero;
                Real ustar           = tan_wspd * kappa / (std::log(1.5 * delta / z0) - psi_m);
                Real Olen            = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny), Real(1.5) * delta);

                for (int iter = 0; iter < 2; ++iter) {
                    if (iter > 0) { Olen  = ib_stability::bounded_obukhov_length(-ustar * ustar * ustar * theta / (kappa * ggg * tflux + tiny), Real(1.5) * delta); }
                    Real zeta          = ib_stability::bounded_zeta(myhalf * delta, Olen);
                    Real zeta_neighbor = ib_stability::bounded_zeta(Real(1.5) * delta, Olen);

                    // similarity functions (off by default: the plain log law, as in the momentum kernels)
                    if (l_stability_correction) {
                        psi_m          = sfuns.calc_psi_m(zeta);
                        psi_h          = sfuns.calc_psi_h(zeta);
                        psi_h_neighbor = sfuns.calc_psi_h(zeta_neighbor);
                    }
                    ustar = tan_wspd * kappa / (std::log((1.5) * delta / z0) - psi_m);
                }

                // prevent some unphysical math
                ustar = ib_stability::clamped_ustar(ustar);
                // caps of ln(z / z0), as above (erf-model/ERF#4134)
                // (with the flag off psi_h is zero and stays so; a cap of ln(z / z0) <= 0, for z0 >= z, is not applied)
                if (l_stability_correction) {
                    psi_h_neighbor = ib_stability::capped_psi(psi_h_neighbor, Real(1.5) * delta, z0, one);
                    psi_h = ib_stability::capped_psi(psi_h, myhalf * delta, z0, one);
                }

                // We do not know the actual temperature so use cell above
                const Real thetastar    = theta * ustar * ustar / (kappa * ggg * Olen);
                const Real surf_temp    = theta_neighbor - thetastar / kappa * (std::log((1.5) * delta / z0) - psi_h_neighbor);
                const Real tTarget      = surf_temp + thetastar / kappa * (std::log((myhalf) * delta / z0) - psi_h);

                const Real bc_forcing_rt = -(cell_data(i,j,k,Rho_comp) * tTarget - cell_data(i,j,k,RhoTheta_comp));
                cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
            }
        }

        // Force fully immersed cells to planar average rho and theta
        if (t_blank == 1.0 && r_avg && t_avg) {
            const Real rho_avg = r_avg(k);
            const Real theta_avg = t_avg(k) / rho_avg;  // Convert from RhoTheta to Theta
            const Real rho_cell = cell_data(i,j,k,Rho_comp);
            const Real bc_forcing_r = -(rho_avg - rho_cell);
            const Real bc_forcing_rt = -(rho_avg * theta_avg - cell_data(i,j,k,RhoTheta_comp));

            cell_src(i, j, k, Rho_comp) -= drag_coefficient * U_s * bc_forcing_r;
            cell_src(i, j, k, RhoTheta_comp) -= drag_coefficient * U_s * bc_forcing_rt;
        }
    });
}
