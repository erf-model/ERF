#ifndef MOSTStress_H
#define MOSTStress_H

#include <ERF_Constants.H>
#include <IndexDefines.H>

struct most_data
{
public:
    amrex::Real zref{10.0};          ///< Reference height (m)
    amrex::Real z0_const{0.1};       ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};        ///< von Karman constant
    amrex::Real gravity{CONST_GRAV}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0}; ///< Heat flux
    amrex::Real surf_temp;           ///< Surface temperature
};


struct similarity_funs
{
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(1.0 - gamma_m * zeta));
            return 2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1.0 + x * x)) -
                   2.0 * std::atan(x) + PIoTwo;
        }
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_h * zeta;
        } else {
            amrex::Real x = std::sqrt(1.0 - gamma_h * zeta);
            return 2.0 * std::log(0.5 * (1.0 + x));
        }
    }

private:
    amrex::Real beta_m{5.0};         ///< Constants from Dyer, BLM, 1974
    amrex::Real beta_h{5.0};         ///< https://doi.org/10.1007/BF00240838
    amrex::Real gamma_m{16.0};
    amrex::Real gamma_h{16.0};
};


struct adiabatic
{
    adiabatic (amrex::Real temp,
               amrex::Real flux)
    {mdata.surf_temp = temp; mdata.surf_temp_flux = flux;}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& /*max_iters*/,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr) const
    {
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(mdata.zref / z0_arr(i,j,k));
        t_star_arr(i,j,k) = 0.0;
        olen_arr(i,j,k)   = 1.0e16;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
};


struct surface_flux
{
    surface_flux (amrex::Real temp,
                  amrex::Real flux)
    {mdata.surf_temp = temp; mdata.surf_temp_flux = flux;}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            Olen = -ustar * ustar * ustar * tm_arr(i,j,k) /
                   (mdata.kappa * mdata.gravity * mdata.surf_temp_flux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        olen_arr(i,j,k)   = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    amrex::Real tol = 1.0e-5;
};


struct surface_temp
{
    surface_temp (amrex::Real temp,
                  amrex::Real flux)
    {mdata.surf_temp = temp; mdata.surf_temp_flux = flux;}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
            Olen = -ustar * ustar * ustar * tm_arr(i,j,k) /
                    (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
        olen_arr(i,j,k)   = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    amrex::Real tol = 1.0e-5;
};


struct moeng_flux
{
    moeng_flux(int l_zlo,
               amrex::Real l_dz)
      :  zlo(l_zlo), dz(l_dz) {}


    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    compute_t_flux (const int& i,
                    const int& j,
                    const int& k,
                    const int& n,
                    const int& icomp,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& eta_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& tm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<const amrex::Real>& t_star_arr,
                    const amrex::Array4<const amrex::Real>& t_surf_arr,
                    const amrex::Array4<amrex::Real>& dest_arr) const
    {
        amrex::Real velx, vely, rho, theta, eta;
        int ix, jx, iy, jy, ie, je, ic, jc;

        ix = i < lbound(velx_arr).x    ? lbound(velx_arr).x   : i;
        jx = j < lbound(velx_arr).y    ? lbound(velx_arr).y   : j;
        ix = ix > ubound(velx_arr).x-1 ? ubound(velx_arr).x-1 : ix;
        jx = jx > ubound(velx_arr).y   ? ubound(velx_arr).y   : jx;

        iy = i  < lbound(vely_arr).x   ? lbound(vely_arr).x   : i;
        jy = j  < lbound(vely_arr).y   ? lbound(vely_arr).y   : j;
        iy = iy > ubound(vely_arr).x   ? ubound(vely_arr).x   : iy;
        jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;

        ic = i  < lbound(cons_arr).x ? lbound(cons_arr).x : i;
        jc = j  < lbound(cons_arr).y ? lbound(cons_arr).y : j;
        ic = ic > ubound(cons_arr).x ? ubound(cons_arr).x : ic;
        jc = jc > ubound(cons_arr).y ? ubound(cons_arr).y : jc;

        velx  = 0.5 *( velx_arr(ix,jx,zlo) + velx_arr(ix+1,jx  ,zlo) );
        vely  = 0.5 *( vely_arr(iy,jy,zlo) + vely_arr(iy  ,jy+1,zlo) );
        rho   = cons_arr(ic,jc,zlo,Rho_comp);
        theta = cons_arr(ic,jc,zlo,RhoTheta_comp) / rho;
        eta   = eta_arr(ie,je,zlo,EddyDiff::Theta_v); // == rho * alpha [kg/m^3 * m^2/s]

        amrex::Real d_thM  =  tm_arr(ic,jc,zlo);
        amrex::Real d_vmM  = umm_arr(ic,jc,zlo);
        amrex::Real d_utau = u_star_arr(ic,jc,zlo);
        amrex::Real d_ttau = t_star_arr(ic,jc,zlo);
        amrex::Real d_sfcT = t_surf_arr(ic,jc,zlo);

        amrex::Real vmag    = sqrt(velx*velx+vely*vely);
        amrex::Real num1    = (theta-d_thM)*d_vmM;
        amrex::Real num2    = (d_thM-d_sfcT)*vmag;
        amrex::Real moflux  = d_ttau*d_utau*(num1+num2)/((d_thM-d_sfcT)*d_vmM);
        amrex::Real deltaz  = dz * (zlo - k);

        dest_arr(i,j,k,icomp+n) = rho*(theta - moflux*rho/eta*deltaz);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    compute_u_flux (const int& i,
                    const int& j,
                    const int& k,
                    const int& icomp,
                    const int& var_idx,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& eta_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& um_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<amrex::Real>& dest_arr) const
    {
        amrex::Real velx, vely, rho, eta;
        int jy, ie, je, ic, jc;

        int iylo = i <= lbound(vely_arr).x ? lbound(vely_arr).x : i-1;
        int iyhi = i >  ubound(vely_arr).x ? ubound(vely_arr).x : i;

        jy = j  < lbound(vely_arr).y   ? lbound(vely_arr).y   : j;
        jy = jy > ubound(vely_arr).y-1 ? ubound(vely_arr).y-1 : jy;

        ie = i  < lbound(eta_arr).x+1 ? lbound(eta_arr).x+1 : i;
        je = j  < lbound(eta_arr).y   ? lbound(eta_arr).y   : j;
        ie = ie > ubound(eta_arr).x   ? ubound(eta_arr).x   : ie;
        je = je > ubound(eta_arr).y   ? ubound(eta_arr).y   : je;

        ic = i  < lbound(cons_arr).x+1 ? lbound(cons_arr).x+1 : i;
        jc = j  < lbound(cons_arr).y   ? lbound(cons_arr).y   : j;
        ic = ic > ubound(cons_arr).x   ? ubound(cons_arr).x   : ic;
        jc = jc > ubound(cons_arr).y   ? ubound(cons_arr).y   : jc;

        velx  = velx_arr(i,j,zlo);
        vely  = 0.25*( vely_arr(iyhi,jy,zlo)+vely_arr(iyhi,jy+1,zlo)
                       +vely_arr(iylo,jy,zlo)+vely_arr(iylo,jy+1,zlo));
        rho   = 0.5 *( cons_arr(ic-1,jc,zlo,Rho_comp)+
                       cons_arr(ic  ,jc,zlo,Rho_comp));
        eta   = 0.5 *( eta_arr(ie-1,je,zlo,EddyDiff::Mom_v)+
                       eta_arr(ie  ,je,zlo,EddyDiff::Mom_v));

        amrex::Real d_vxM  = um_arr(i,j,zlo);
        amrex::Real d_vmM  = 0.5 * ( umm_arr(ic-1,jc,zlo) + umm_arr(ic,jc,zlo) );
        amrex::Real d_utau = 0.5 * ( u_star_arr(ic-1,jc,zlo) + u_star_arr(ic,jc,zlo) );

        amrex::Real vmag    = sqrt(velx*velx+vely*vely);
        amrex::Real stressx = ( (velx-d_vxM)*d_vmM + vmag*d_vxM ) / (d_vmM*d_vmM) * d_utau*d_utau;
        amrex::Real deltaz  = dz * (zlo - k);

        if (var_idx == Vars::xmom) {
            dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho*rho/eta*deltaz;
        } else {
            AMREX_ALWAYS_ASSERT(var_idx == Vars::xvel);
            dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressx*rho/eta*deltaz;
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    compute_v_flux (const int& i,
                    const int& j,
                    const int& k,
                    const int& icomp,
                    const int& var_idx,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& eta_arr,
                    const amrex::Array4<const amrex::Real>& umm_arr,
                    const amrex::Array4<const amrex::Real>& vm_arr,
                    const amrex::Array4<const amrex::Real>& u_star_arr,
                    const amrex::Array4<amrex::Real>& dest_arr) const
    {
        amrex::Real velx, vely, rho, eta;
        int ix, ie, je, ic, jc;

        ix = i  < lbound(velx_arr).x ? lbound(velx_arr).x : i;
        ix = ix > ubound(velx_arr).x ? ubound(velx_arr).x : ix;

        int jxlo = j <= lbound(velx_arr).y ? lbound(velx_arr).y : j-1;
        int jxhi = j >  ubound(velx_arr).y ? ubound(velx_arr).y : j;

        ie = i  < lbound(eta_arr).x   ? lbound(eta_arr).x   : i;
        je = j  < lbound(eta_arr).y+1 ? lbound(eta_arr).y+1 : j;
        ie = ie > ubound(eta_arr).x   ? ubound(eta_arr).x   : ie;
        je = je > ubound(eta_arr).y   ? ubound(eta_arr).y   : je;

        ic = i  < lbound(cons_arr).x   ? lbound(cons_arr).x   : i;
        jc = j  < lbound(cons_arr).y+1 ? lbound(cons_arr).y+1 : j;
        ic = ic > ubound(cons_arr).x   ? ubound(cons_arr).x   : ic;
        jc = jc > ubound(cons_arr).y   ? ubound(cons_arr).y   : jc;

        velx  = 0.25*( velx_arr(ix,jxhi,zlo)+velx_arr(ix+1,jxhi,zlo)
                      +velx_arr(ix,jxlo,zlo)+velx_arr(ix+1,jxlo,zlo));
        vely  = vely_arr(i,j,zlo);
        rho   = 0.5*(cons_arr(ic,jc-1,zlo,Rho_comp)+
                     cons_arr(ic,jc  ,zlo,Rho_comp));
        eta   = 0.5*(eta_arr(ie,je-1,zlo,EddyDiff::Mom_v)+
                     eta_arr(ie,je  ,zlo,EddyDiff::Mom_v));

        amrex::Real d_vyM  = vm_arr(i,j,zlo);
        amrex::Real d_vmM  = 0.5 * ( umm_arr(ic,jc-1,zlo) + umm_arr(ic,jc,zlo) );
        amrex::Real d_utau = 0.5 * ( u_star_arr(ic,jc-1,zlo) + u_star_arr(ic,jc,zlo) );

        amrex::Real vmag    = sqrt(velx*velx+vely*vely);
        amrex::Real stressy = ( (vely-d_vyM)*d_vmM + vmag*d_vyM ) / (d_vmM*d_vmM)*d_utau*d_utau;
        amrex::Real deltaz  = dz * (zlo - k);

        if (var_idx == Vars::ymom) {
            dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho*rho/eta*deltaz;
        } else {
            AMREX_ALWAYS_ASSERT(var_idx == Vars::yvel);
            dest_arr(i,j,k,icomp) = dest_arr(i,j,zlo,icomp) - stressy*rho/eta*deltaz;
        }
    }

private:
    int zlo;
    amrex::Real dz;
};
#endif
