#ifndef ERF_EBMOSTStress_H
#define ERF_EBMOSTStress_H

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>


/**
 * Adiabatic with constant roughness
 */
struct adiabatic_eb
{
    adiabatic_eb (amrex::Real Tflux,
               amrex::Real Qvflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& /*max_iters*/,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& /*q_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        olen_arr(i,j,k)   = bogus_large_value;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,0) / std::log(zref_arr(i,j,0) / z0_arr(i,j,k));
        t_star_arr(i,j,k) = zero;
        q_star_arr(i,j,k) = zero;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
};


/**
 * Surface temperature with constant roughness
 */
struct surface_temp_eb
{
    surface_temp_eb (amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = zero;
        amrex::Real zeta     = zero;
        amrex::Real zeta_old = zero;
        amrex::Real psi_m    = zero;
        amrex::Real psi_h    = zero;
        amrex::Real num      = zero;
        amrex::Real den      = zero;
        amrex::Real zref     = zref_arr(i,j,0);
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,0), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == bogus_large_value) {
            olen_arr(i,j,k) = amrex::Real(1000);
        }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > zero) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,0);
        }
        qv_a = qvm_arr(i,j,0);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
            amrex::Real ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,0) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(1 + amrex::Real(0.61)*qvm_arr(i,j,0)) + amrex::Real(0.61)*tm_arr(i,j,0)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,0));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,0)*umm_arr(i,j,0) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (one + amrex::Real(0.61)*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,0) * (one + amrex::Real(0.61)*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,0) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-amrex::Real(4.0)),amrex::Real(4.0));

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Limiting
            num = std::max(C - psi_m, amrex::Real(1.0));
            den = std::max(C - psi_h, amrex::Real(1.0));

            // Update with under relaxation
            zeta = (one - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,0) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,0);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,0) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = amrex::Real(1.0e-3);
    const amrex::Real alpha = myhalf;
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};

/**
 * Surface flux with constant roughness
 */
struct surface_flux_eb
{
    surface_flux_eb (amrex::Real Tflux,
                  amrex::Real Qvflux,
                  bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        spec_qflux = cons_qflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = zero;
        amrex::Real wstar = zero;
        amrex::Real tflux = zero;
        amrex::Real qflux = zero;
        amrex::Real zeta  = zero;
        amrex::Real psi_m = zero;
        amrex::Real psi_h = zero;
        amrex::Real Olen  = zero;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == bogus_large_value) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,0) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + amrex::Real(0.61)*qvm_arr(i,j,0)) + qflux*amrex::Real(0.61)*tm_arr(i,j,0);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,0));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,0)*umm_arr(i,j,0) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,0) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,0);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,0);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,0) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = amrex::Real(1.0e-5);
    const amrex::Real WSMIN = amrex::Real(0.1); // minimum wind speed
};

/**
 * Moeng flux formulation
 */
struct moeng_flux_eb
{
    moeng_flux_eb () {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& qvm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& q_surf_arr,
                    const amrex::Array4<const amrex::Real>& u_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& v_vfrac_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real qv    = cons_arr(i,j,k,RhoQ1_comp) / rho;

        // Volume-weighted average of x-face velocities to cell center
        amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k);
        amrex::Real velx = (u_vfrac_sum > eps) ?
            (velx_arr(i,j,k) * u_vfrac_arr(i,j,k) + velx_arr(i+1,j,k) * u_vfrac_arr(i+1,j,k))
            / u_vfrac_sum : zero;

        // Volume-weighted average of y-face velocities to cell center
        amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k);
        amrex::Real vely = (v_vfrac_sum > eps) ?
            (vely_arr(i,j,k) * v_vfrac_arr(i,j,k) + vely_arr(i,j+1,k) * v_vfrac_arr(i,j+1,k))
            / v_vfrac_sum : zero;

        amrex::Real qv_mean  = qvm_arr(i,j,0);
        amrex::Real ustar    = u_star_arr(i,j,k);
        amrex::Real qstar    = q_star_arr(i,j,k);
        amrex::Real qv_surf  = q_surf_arr(i,j,k);
        amrex::Real wsp_mean = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = std::sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * (qv_mean-qv_surf);
        amrex::Real num2    = wsp_mean * (qv-qv_mean);

        // NOTE: this is rho*<Qv'w'> = -K dQvdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ?
                              -rho*qstar*ustar*(num1+num2)/((qv_mean-qv_surf)*wsp_mean) : zero;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& velz_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr,
                    const amrex::Array4<const amrex::Real>& u_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& v_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& w_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& bnorm_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real theta = cons_arr(i,j,k,RhoTheta_comp) / rho;

        // Volume-weighted average of x-face velocities to cell center
        amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k);
        amrex::Real velx = (u_vfrac_sum > eps) ?
            (velx_arr(i,j,k) * u_vfrac_arr(i,j,k) + velx_arr(i+1,j,k) * u_vfrac_arr(i+1,j,k))
            / u_vfrac_sum : zero;

        // Volume-weighted average of y-face velocities to cell center
        amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k);
        amrex::Real vely = (v_vfrac_sum > eps) ?
            (vely_arr(i,j,k) * v_vfrac_arr(i,j,k) + vely_arr(i,j+1,k) * v_vfrac_arr(i,j+1,k))
            / v_vfrac_sum : zero;

        // Volume-weighted average of z-face velocities to cell center
        amrex::Real w_vfrac_sum = w_vfrac_arr(i,j,k) + w_vfrac_arr(i,j,k+1);
        amrex::Real velz = (w_vfrac_sum > eps) ?
            (velz_arr(i,j,k) * w_vfrac_arr(i,j,k) + velz_arr(i,j,k+1) * w_vfrac_arr(i,j,k+1))
            / w_vfrac_sum : zero;

        // Get boundary normal components
        amrex::Real nx = bnorm_arr(i,j,k,0);
        amrex::Real ny = bnorm_arr(i,j,k,1);
        amrex::Real nz = bnorm_arr(i,j,k,2);

        // Project velocity onto tangent plane (remove normal component)
        amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
        amrex::Real velx_tangent = velx - v_dot_n * nx;
        amrex::Real vely_tangent = vely - v_dot_n * ny;

        amrex::Real theta_mean  = tm_arr(i,j,0);
        amrex::Real ustar       = u_star_arr(i,j,k);
        amrex::Real tstar       = t_star_arr(i,j,k);
        amrex::Real theta_surf  = t_surf_arr(i,j,k);
        amrex::Real wsp_mean    = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Use tangential velocity magnitude instead of Cartesian
        amrex::Real wsp     = std::sqrt(velx_tangent*velx_tangent+vely_tangent*vely_tangent);
        amrex::Real num1    = wsp * (theta_mean-theta_surf);
        amrex::Real num2    = wsp_mean * (theta-theta_mean);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ?
                              -rho*tstar*ustar*(num1+num2)/((theta_mean-theta_surf)*wsp_mean) : zero;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (int i,
                    int j,
                    int k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& velz_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& um_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& u_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& v_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& w_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& cc_vfrac_arr,
                    const amrex::Array4<const amrex::EBCellFlag>& cc_flag_arr,
                    const amrex::Array4<const amrex::Real>& bnorm_arr,
                    int idir = 0) const
    {
        amrex::Real velx, vely, rho, ustar, wsp_mean;
        amrex::Real velx_tangent, vely_tangent;

        if (idir == 0) {
            // x-face: average to x-face
            velx = velx_arr(i,j,k);

            // Volume-weighted average of y-face velocities to x-face
            amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k) +
                                      v_vfrac_arr(i-1,j,k) + v_vfrac_arr(i-1,j+1,k);
            vely = (v_vfrac_sum > eps) ?
                (vely_arr(i,j,k)   * v_vfrac_arr(i,j,k)   + vely_arr(i,j+1,k) * v_vfrac_arr(i,j+1,k) +
                 vely_arr(i-1,j,k) * v_vfrac_arr(i-1,j,k) + vely_arr(i-1,j+1,k) * v_vfrac_arr(i-1,j+1,k))
                 / v_vfrac_sum : zero;

            // Volume-weighted average of z-face velocities to x-face
            amrex::Real w_vfrac_sum = w_vfrac_arr(i,j,k) + w_vfrac_arr(i,j,k+1) +
                                      w_vfrac_arr(i-1,j,k) + w_vfrac_arr(i-1,j,k+1);
            amrex::Real velz = (w_vfrac_sum > eps) ?
                (velz_arr(i,j,k)   * w_vfrac_arr(i,j,k)   + velz_arr(i,j,k+1)   * w_vfrac_arr(i,j,k+1) +
                 velz_arr(i-1,j,k) * w_vfrac_arr(i-1,j,k) + velz_arr(i-1,j,k+1) * w_vfrac_arr(i-1,j,k+1))
                 / w_vfrac_sum : zero;

            // Get boundary normal at x-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to x-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i-1,j,k) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i-1,j,k,Rho_comp) * cc_vfrac_arr(i-1,j,k) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to x-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i-1,j,k).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i-1,j,k) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i-1,j,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i-1,j,k);
                wsp_mean = umm_arr(i-1,j,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }

        } else if (idir == 1) {
            // y-face: average to y-face
            vely = vely_arr(i,j,k);

            // Volume-weighted average of x-face velocities to y-face
            amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k) +
                                      u_vfrac_arr(i,j-1,k) + u_vfrac_arr(i+1,j-1,k);
            velx = (u_vfrac_sum > eps) ?
                (velx_arr(i,j,k)   * u_vfrac_arr(i,j,k)   + velx_arr(i+1,j,k)   * u_vfrac_arr(i+1,j,k) +
                 velx_arr(i,j-1,k) * u_vfrac_arr(i,j-1,k) + velx_arr(i+1,j-1,k) * u_vfrac_arr(i+1,j-1,k))
                 / u_vfrac_sum : zero;

            // Volume-weighted average of z-face velocities to y-face
            amrex::Real w_vfrac_sum = w_vfrac_arr(i,j,k) + w_vfrac_arr(i,j,k+1) +
                                      w_vfrac_arr(i,j-1,k) + w_vfrac_arr(i,j-1,k+1);
            amrex::Real velz = (w_vfrac_sum > eps) ?
                (velz_arr(i,j,k)   * w_vfrac_arr(i,j,k)   + velz_arr(i,j,k+1)   * w_vfrac_arr(i,j,k+1) +
                 velz_arr(i,j-1,k) * w_vfrac_arr(i,j-1,k) + velz_arr(i,j-1,k+1) * w_vfrac_arr(i,j-1,k+1))
                 / w_vfrac_sum : zero;

            // Get boundary normal at y-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to y-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i,j-1,k) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i,j-1,k,Rho_comp) * cc_vfrac_arr(i,j-1,k) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to y-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i,j-1,k).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i,j-1,k) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i,j-1,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i,j-1,k);
                wsp_mean = umm_arr(i,j-1,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }

        } else {
            // z-face: average to z-face
            // Volume-weighted average of x-face velocities to z-face
            amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k-1) + u_vfrac_arr(i+1,j,k-1) +
                                      u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k);
            velx = (u_vfrac_sum > eps) ?
                (velx_arr(i,j,k-1)   * u_vfrac_arr(i,j,k-1)   + velx_arr(i+1,j,k-1)   * u_vfrac_arr(i+1,j,k-1) +
                 velx_arr(i,j,k) * u_vfrac_arr(i,j,k) + velx_arr(i+1,j,k) * u_vfrac_arr(i+1,j,k))
                 / u_vfrac_sum : zero;

            // Volume-weighted average of y-face velocities to z-face
            amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k-1) + v_vfrac_arr(i,j+1,k-1) +
                                      v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k);
            vely = (v_vfrac_sum > eps) ?
                (vely_arr(i,j,k-1)   * v_vfrac_arr(i,j,k-1)   + vely_arr(i,j+1,k-1)   * v_vfrac_arr(i,j+1,k-1) +
                 vely_arr(i,j,k) * v_vfrac_arr(i,j,k) + vely_arr(i,j+1,k) * v_vfrac_arr(i,j+1,k))
                 / v_vfrac_sum : zero;

            // z-velocity is already at z-face
            amrex::Real velz = velz_arr(i,j,k);

            // Get boundary normal at z-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to z-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i,j,k-1) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i,j,k-1,Rho_comp) * cc_vfrac_arr(i,j,k-1) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to z-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i,j,k-1).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i,j,k-1) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i,j,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i,j,k-1);
                wsp_mean = umm_arr(i,j,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }
        }

        wsp_mean = std::max(wsp_mean, WSMIN);
        amrex::Real umean = um_arr(i,j,0);

        // Note: The surface mean shear stress is decomposed into tau_xz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of umean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
        amrex::Real wsp     = std::sqrt(velx_tangent*velx_tangent+vely_tangent*vely_tangent);
        amrex::Real num1    = wsp * umean;
        amrex::Real num2    = wsp_mean * (velx_tangent-umean);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressx;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (int i,
                    int j,
                    int k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& velz_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& vm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& u_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& v_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& w_vfrac_arr,
                    const amrex::Array4<const amrex::Real>& cc_vfrac_arr,
                    const amrex::Array4<const amrex::EBCellFlag>& cc_flag_arr,
                    const amrex::Array4<const amrex::Real>& bnorm_arr,
                    int idir = 0) const
    {
        amrex::Real velx, vely, rho, ustar, wsp_mean;
        amrex::Real velx_tangent, vely_tangent;

        if (idir == 0) {
            // x-face: average from cells (i-1) and (i)
            // Volume-weighted average of y-face velocities to x-face
            amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k) +
                                      v_vfrac_arr(i-1,j,k) + v_vfrac_arr(i-1,j+1,k);
            vely = (v_vfrac_sum > eps) ?
                (vely_arr(i,j,k)   * v_vfrac_arr(i,j,k)   + vely_arr(i,j+1,k)   * v_vfrac_arr(i,j+1,k) +
                 vely_arr(i-1,j,k) * v_vfrac_arr(i-1,j,k) + vely_arr(i-1,j+1,k) * v_vfrac_arr(i-1,j+1,k))
                 / v_vfrac_sum : zero;

            velx = velx_arr(i,j,k);

            // Volume-weighted average of z-face velocities to x-face
            amrex::Real w_vfrac_sum = w_vfrac_arr(i,j,k) + w_vfrac_arr(i,j,k+1) +
                                      w_vfrac_arr(i-1,j,k) + w_vfrac_arr(i-1,j,k+1);
            amrex::Real velz = (w_vfrac_sum > eps) ?
                (velz_arr(i,j,k)   * w_vfrac_arr(i,j,k)   + velz_arr(i,j,k+1)   * w_vfrac_arr(i,j,k+1) +
                 velz_arr(i-1,j,k) * w_vfrac_arr(i-1,j,k) + velz_arr(i-1,j,k+1) * w_vfrac_arr(i-1,j,k+1))
                 / w_vfrac_sum : zero;

            // Get boundary normal at x-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to x-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i-1,j,k) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i-1,j,k,Rho_comp) * cc_vfrac_arr(i-1,j,k) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to x-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i-1,j,k).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i-1,j,k) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i-1,j,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i-1,j,k);
                wsp_mean = umm_arr(i-1,j,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }

        } else if (idir == 1) {
            // y-face: average from cells (j-1) and (j)
            // Volume-weighted average of x-face velocities to y-face
            amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k) +
                                      u_vfrac_arr(i,j-1,k) + u_vfrac_arr(i+1,j-1,k);
            velx = (u_vfrac_sum > eps) ?
                (velx_arr(i,j,k)   * u_vfrac_arr(i,j,k)   + velx_arr(i+1,j,k)   * u_vfrac_arr(i+1,j,k) +
                 velx_arr(i,j-1,k) * u_vfrac_arr(i,j-1,k) + velx_arr(i+1,j-1,k) * u_vfrac_arr(i+1,j-1,k))
                 / u_vfrac_sum : zero;

            vely = vely_arr(i,j,k);

            // Volume-weighted average of z-face velocities to y-face
            amrex::Real w_vfrac_sum = w_vfrac_arr(i,j,k) + w_vfrac_arr(i,j,k+1) +
                                      w_vfrac_arr(i,j-1,k) + w_vfrac_arr(i,j-1,k+1);
            amrex::Real velz = (w_vfrac_sum > eps) ?
                (velz_arr(i,j,k)   * w_vfrac_arr(i,j,k)   + velz_arr(i,j,k+1)   * w_vfrac_arr(i,j,k+1) +
                 velz_arr(i,j-1,k) * w_vfrac_arr(i,j-1,k) + velz_arr(i,j-1,k+1) * w_vfrac_arr(i,j-1,k+1))
                 / w_vfrac_sum : zero;

            // Get boundary normal at y-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to y-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i,j-1,k) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i,j-1,k,Rho_comp) * cc_vfrac_arr(i,j-1,k) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to y-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i,j-1,k).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i,j-1,k) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i,j-1,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i,j-1,k);
                wsp_mean = umm_arr(i,j-1,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }

        } else {
            // z-face: average from cells (k-1) and (k)
            // Volume-weighted average of x-face velocities to z-face
            amrex::Real u_vfrac_sum = u_vfrac_arr(i,j,k-1) + u_vfrac_arr(i+1,j,k-1) +
                                      u_vfrac_arr(i,j,k) + u_vfrac_arr(i+1,j,k);
            velx = (u_vfrac_sum > eps) ?
                (velx_arr(i,j,k-1)   * u_vfrac_arr(i,j,k-1)   + velx_arr(i+1,j,k-1)   * u_vfrac_arr(i+1,j,k-1) +
                 velx_arr(i,j,k)     * u_vfrac_arr(i,j,k)     + velx_arr(i+1,j,k)     * u_vfrac_arr(i+1,j,k))
                 / u_vfrac_sum : zero;

            // Volume-weighted average of y-face velocities to z-face
            amrex::Real v_vfrac_sum = v_vfrac_arr(i,j,k-1) + v_vfrac_arr(i,j+1,k-1) +
                                      v_vfrac_arr(i,j,k) + v_vfrac_arr(i,j+1,k);
            vely = (v_vfrac_sum > eps) ?
                (vely_arr(i,j,k-1)   * v_vfrac_arr(i,j,k-1)   + vely_arr(i,j+1,k-1)   * v_vfrac_arr(i,j+1,k-1) +
                 vely_arr(i,j,k)     * v_vfrac_arr(i,j,k)     + vely_arr(i,j+1,k)     * v_vfrac_arr(i,j+1,k))
                 / v_vfrac_sum : zero;

            // z-velocity is already at z-face
            amrex::Real velz = velz_arr(i,j,k);

            // Get boundary normal at z-face (already at correct staggered location)
            amrex::Real nx = bnorm_arr(i,j,k,0);
            amrex::Real ny = bnorm_arr(i,j,k,1);
            amrex::Real nz = bnorm_arr(i,j,k,2);

            // Project velocity onto tangent plane
            amrex::Real v_dot_n = velx*nx + vely*ny + velz*nz;
            velx_tangent = velx - v_dot_n * nx;
            vely_tangent = vely - v_dot_n * ny;

            // Volume-weighted average of cell-centered density to z-face
            amrex::Real cc_vfrac_sum = cc_vfrac_arr(i,j,k-1) + cc_vfrac_arr(i,j,k);
            rho = (cc_vfrac_sum > eps) ?
                (cons_arr(i,j,k-1,Rho_comp) * cc_vfrac_arr(i,j,k-1) + cons_arr(i,j,k,Rho_comp) * cc_vfrac_arr(i,j,k))
                / cc_vfrac_sum : zero;

            // Average cell-centered u_star and wsp_mean to z-face, using only valid (SingleValued) cells
            bool low_valid  = cc_flag_arr(i,j,k-1).isSingleValued();
            bool high_valid = cc_flag_arr(i,j,k).isSingleValued();

            if (low_valid && high_valid) {
                ustar    = myhalf * (u_star_arr(i,j,k-1) + u_star_arr(i,j,k));
                wsp_mean = myhalf * (umm_arr(i,j,0) + umm_arr(i,j,0));
            } else if (low_valid) {
                ustar    = u_star_arr(i,j,k-1);
                wsp_mean = umm_arr(i,j,0);
            } else if (high_valid) {
                ustar    = u_star_arr(i,j,k);
                wsp_mean = umm_arr(i,j,0);
            } else {
                ustar    = zero;
                wsp_mean = WSMIN;
            }
        }

        wsp_mean = std::max(wsp_mean, WSMIN);
        amrex::Real vmean = vm_arr(i,j,0);

        // Note: The surface mean shear stress is decomposed into tau_yz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of vmean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
        amrex::Real wsp     = std::sqrt(velx_tangent*velx_tangent+vely_tangent*vely_tangent);
        amrex::Real num1    = wsp * vmean;
        amrex::Real num2    = wsp_mean * (vely_tangent-vmean);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressy;
    }

private:
#ifdef AMREX_USE_FLOAT
    const amrex::Real     eps = amrex::Real(1e-6);
#else
    const amrex::Real     eps = amrex::Real(1e-12);
#endif
    const amrex::Real   WSMIN = amrex::Real(0.1); // minimum wind speed
};
#endif
