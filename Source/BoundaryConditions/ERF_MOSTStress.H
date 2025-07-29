#ifndef ERF_MOSTStress_H
#define ERF_MOSTStress_H

#include <ERF_MOST.H>
#include <ERF_IndexDefines.H>
#include <ERF_Wstar.H>

/**
 * Moeng flux formulation
 */
struct moeng_flux
{
    moeng_flux () {}

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
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real qv    = cons_arr(i,j,k,RhoQ1_comp) / rho;
        amrex::Real velx  = 0.5 * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = 0.5 * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real qv_mean  = qvm_arr(i,j,k);
        amrex::Real ustar    = u_star_arr(i,j,k);
        amrex::Real qstar    = q_star_arr(i,j,k);
        amrex::Real qv_surf  = q_surf_arr(i,j,k);
        amrex::Real wsp_mean = umm_arr(i,j,k);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * (qv_mean-qv_surf);
        amrex::Real num2    = wsp_mean * (qv-qv_mean);

        // NOTE: this is rho*<Qv'w'> = -K dQvdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ?
                              -rho*qstar*ustar*(num1+num2)/((qv_mean-qv_surf)*wsp_mean) : 0.0;

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
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real theta = cons_arr(i,j,k,RhoTheta_comp) / rho;
        amrex::Real velx  = 0.5 * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = 0.5 * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real theta_mean  = tm_arr(i,j,k);
        amrex::Real ustar       = u_star_arr(i,j,k);
        amrex::Real tstar       = t_star_arr(i,j,k);
        amrex::Real theta_surf  = t_surf_arr(i,j,k);
        amrex::Real wsp_mean    = umm_arr(i,j,k);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * (theta_mean-theta_surf);
        amrex::Real num2    = wsp_mean * (theta-theta_mean);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ?
                              -rho*tstar*ustar*(num1+num2)/((theta_mean-theta_surf)*wsp_mean) : 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& um_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = 0.25 * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                   + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = 0.5  * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real umean    = um_arr(i,j,k);
        amrex::Real ustar    = 0.5 * ( u_star_arr(i-1,j,k) + u_star_arr(i,j,k) );
        amrex::Real wsp_mean = 0.5 * (    umm_arr(i-1,j,k) +    umm_arr(i,j,k) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_xz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of umean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng 1984.
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * umean;
        amrex::Real num2    = wsp_mean * (velx-umean);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressx;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& vm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = 0.25 * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                   + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = 0.5  * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real vmean     = vm_arr(i,j,k);
        amrex::Real ustar     = 0.5 * ( u_star_arr(i,j-1,k) + u_star_arr(i,j,k) );
        amrex::Real wsp_mean  = 0.5 * (    umm_arr(i,j-1,k) +    umm_arr(i,j,k) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_yz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of vmean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng 1984.
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * vmean;
        amrex::Real num2    = wsp_mean * (vely-vmean);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressy;
    }

private:
    const amrex::Real   eps = 1e-15;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Donelan flux formulation
 */
struct donelan_flux
{
    donelan_flux () {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& /*i*/,
                    const int& /*j*/,
                    const int& /*k*/,
                    const amrex::Array4<const amrex::Real>& /*cons_arr*/,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*q_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*q_surf_arr*/) const
    {
        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*t_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real Ch          = 0.0012;
        amrex::Real rho         = cons_arr(i,j,k,Rho_comp);
        amrex::Real theta_surf  = t_surf_arr(i,j,k);
        amrex::Real theta_mean  =  tm_arr(i,j,k);
        amrex::Real wsp_mean    = umm_arr(i,j,k);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux = -rho * Ch * wsp_mean * (theta_mean - theta_surf);

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = 0.25 * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                   + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = 0.5  * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real Cd       = 0.001;
        const amrex::Real c  = 7e-5;
        amrex::Real wsp      = sqrt(velx*velx+vely*vely);
        amrex::Real wsp_mean = 0.5 * ( umm_arr(i-1,j,k) + umm_arr(i,j,k) );
        if (wsp_mean <= 5.0) {
            Cd = 0.001;
        } else if (wsp_mean < 25.0 && wsp_mean > 5.0) {
            Cd = 0.001 + c * (wsp_mean - 5.0);
        } else {
            Cd = 0.0024;
        }

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho * Cd * velx * wsp;

        return stressx;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real velx  = 0.25 * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                   + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = 0.5  * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real Cd       = 0.001;
        const amrex::Real c  = 7e-5;
        amrex::Real wsp      = sqrt(velx*velx+vely*vely);
        amrex::Real wsp_mean = 0.5 * ( umm_arr(i,j-1,k) + umm_arr(i,j,k) );
        if (wsp_mean <= 5.0) {
            Cd = 0.001;
        } else if (wsp_mean < 25.0 && wsp_mean > 5.0) {
            Cd = 0.001 + c * (wsp_mean - 5.0);
        } else {
            Cd = 0.0024;
        }

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho * Cd * vely * wsp;

        return stressy;
    }

private:

};


/**
 * Custom flux formulation
 */
struct custom_flux
{
    custom_flux (bool specified_rho_surf)
    : fluxes_include_rho(specified_rho_surf)
    {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& /*q_surf_arr*/) const
    {
        amrex::Real rho   = (fluxes_include_rho) ? 1.0 : cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar = q_star_arr(i,j,k);

        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = (std::abs(qstar) > eps) ? rho * qstar : 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        amrex::Real rho   = (fluxes_include_rho) ? 1.0 : cons_arr(i,j,k,Rho_comp);
        amrex::Real tstar = t_star_arr(i,j,k);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = (std::abs(tstar) > eps) ? rho * tstar : 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = 0.25 * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                   + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = (fluxes_include_rho) ? 1.0 : 0.5 * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = 0.5 * ( u_star_arr(i-1,j,k) + u_star_arr(i,j,k) );
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<u'w'> = -K dudz
        amrex::Real stressx = -rho * ustar * ustar * velx / wsp;

        return stressx;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        amrex::Real velx  = 0.25 * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                   + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = (fluxes_include_rho) ? 1.0 : 0.5 * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = 0.5 * ( u_star_arr(i,j-1,k) + u_star_arr(i,j,k) );
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho * ustar * ustar * vely / wsp;

        return stressy;
    }

private:
    const amrex::Real eps = 1e-15;
    const bool fluxes_include_rho{false};

};

/**
 * Rotate flux formulation
 */
struct rotate_flux
{
    rotate_flux () {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& qvm_arr,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_star_arr,
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real qv_mean = qvm_arr(i,j,k);
        amrex::Real qv_surf = q_surf_arr(i,j,k);
        amrex::Real rho     = cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar   = q_star_arr(i,j,k);
        amrex::Real moflux  = (std::abs(qstar) > eps) ? -rho * qstar * (qv_mean-qv_surf): 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real theta_mean = tm_arr(i,j,k);
        amrex::Real theta_surf = t_surf_arr(i,j,k);
        amrex::Real rho        = cons_arr(i,j,k,Rho_comp);
        amrex::Real tstar      = t_star_arr(i,j,k);
        amrex::Real moflux = (std::abs(tstar) > eps) ? -rho * tstar * (theta_mean-theta_surf) : 0.0;

        return moflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*um_arr*/,
                    const amrex::Array4<const amrex::Real>& u_star_arr) const
    {
        // NOTE: this is the total stress
        amrex::Real rho     = 0.5 *( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i  ,j,k,Rho_comp) );
        amrex::Real ustar   = 0.5 * ( u_star_arr(i-1,j,k) + u_star_arr(i,j,k) );
        amrex::Real stressx = -rho * ustar * ustar;

        return stressx;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_v_flux (const int& /*i*/,
                    const int& /*j*/,
                    const int& /*k*/,
                    const amrex::Array4<const amrex::Real>& /*cons_arr*/,
                    const amrex::Array4<const amrex::Real>& /*velx_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vely_arr*/,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*vm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        // NOTE: this is the total stress
        amrex::Real stressy = 0.0;

        return stressy;
    }

private:
    const amrex::Real eps = 1e-15;
};
#endif
