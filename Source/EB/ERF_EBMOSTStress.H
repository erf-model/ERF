#ifndef ERF_EBMOSTStress_H
#define ERF_EBMOSTStress_H

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>

/**
 * Moeng flux formulation
 */
struct eb_moeng_flux
{
    eb_moeng_flux () {}

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
        amrex::Real velx  = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real qv_mean  = qvm_arr(i,j,0);
        amrex::Real ustar    = u_star_arr(i,j,0);
        amrex::Real qstar    = q_star_arr(i,j,0);
        amrex::Real qv_surf  = q_surf_arr(i,j,0);
        amrex::Real wsp_mean = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
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
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real rho   = cons_arr(i,j,k,Rho_comp);
        amrex::Real theta = cons_arr(i,j,k,RhoTheta_comp) / rho;
        amrex::Real velx  = myhalf * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely  = myhalf * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );

        amrex::Real theta_mean  = tm_arr(i,j,0);
        amrex::Real ustar       = u_star_arr(i,j,0);
        amrex::Real tstar       = t_star_arr(i,j,0);
        amrex::Real theta_surf  = t_surf_arr(i,j,0);
        amrex::Real wsp_mean    = umm_arr(i,j,0);
        wsp_mean = std::max(wsp_mean, WSMIN);

        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
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
        amrex::Real vely  = fourth * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                   + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real rho   = myhalf  * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real umean    = um_arr(i,j,0);
        amrex::Real ustar    = myhalf * ( u_star_arr(i-1,j,0) + u_star_arr(i,j,0) );
        amrex::Real wsp_mean = myhalf * (    umm_arr(i-1,j,0) +    umm_arr(i,j,0) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_xz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of umean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
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
        amrex::Real velx  = fourth * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                   + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real rho   = myhalf  * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real vmean     = vm_arr(i,j,0);
        amrex::Real ustar     = myhalf * ( u_star_arr(i,j-1,0) + u_star_arr(i,j,0) );
        amrex::Real wsp_mean  = myhalf * (    umm_arr(i,j-1,0) +    umm_arr(i,j,0) );
        wsp_mean = std::max(wsp_mean, WSMIN);

        // Note: The surface mean shear stress is decomposed into tau_yz by
        //       multiplying the modeled shear stress (rho*ustar^2) with
        //       a factor of vmean/wsp_mean for directionality; this factor
        //       modifies the denominator from what is in Moeng amrex::Real(1984.)
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = wsp * vmean;
        amrex::Real num2    = wsp_mean * (vely-vmean);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho*ustar*ustar * (num1+num2)/(wsp_mean*wsp_mean);

        return stressy;
    }

private:
    const amrex::Real     eps = 1e-15;
    const amrex::Real   WSMIN = amrex::Real(0.1); // minimum wind speed
};
#endif