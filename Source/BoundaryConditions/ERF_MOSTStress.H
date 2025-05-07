#ifndef ERF_MOSTStress_H
#define ERF_MOSTStress_H

#include <ERF_Constants.H>
#include <ERF_IndexDefines.H>
#include <ERF_MOSTRoughness.H>
#include <ERF_Wstar.H>

/**
 * Structure of plain old data relevant to MOST BCs
 */
struct most_data
{
public:
    amrex::Real zref{10.0};          ///< Reference height (m)
    amrex::Real z0_const{0.1};       ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};        ///< von Karman constant
    amrex::Real gravity{CONST_GRAV}; ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0}; ///< Heat flux

    amrex::Real Cnk_a{0.0185};       ///< Standard Charnock constant https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b1{1.0/30.0};    ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b2{1260.0};      ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_d{30.0};         ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b;
    bool visc{false};                ///< Use viscous Charnock formulation

    const amrex::Real Bjr_beta = 1.2; // Empirical parameter from Beljaars 1995 QJRMS
};


/**
 * Structure of similarity functions for Moeng formulation
 */
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


/**
 * Empirical kinematic viscosity [m2/s] formula from Andreas (1989) CRREL Rep.
 * 89-11, valid between -173 and 277 deg C.
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
air_viscosity (amrex::Real T_degK)
{
    amrex::Real TC = T_degK - 273.15;
    return 1.326e-5*(1. + 6.542e-3*TC + 8.301e-6*TC*TC - 4.84e-9*TC*TC*TC);
}


/**
 * Adiabatic with constant roughness
 */
struct adiabatic
{
    adiabatic (amrex::Real zref,
               amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

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
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*q_star_arr*/,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(mdata.zref / z0_arr(i,j,k));
        t_star_arr(i,j,k) = 0.0;
        olen_arr(i,j,k)   = 1.0e16;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
};


/**
 * Adiabatic with charnock roughness
 */
struct adiabatic_charnock
{
    adiabatic_charnock (amrex::Real zref,
                        amrex::Real flux,
                        amrex::Real cnk_a,
                        bool cnk_visc)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_a = cnk_a;
        mdata.visc = cnk_visc;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*q_star_arr*/,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real z0    = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            if (mdata.Cnk_a > 0) {
                z0 = (mdata.Cnk_a / mdata.gravity) * ustar * ustar;
                if (mdata.visc) {
                    z0 += air_viscosity(tm_arr(i,j,k)) / std::max(ustar, 0.05);
                }
            } else {
                z0 = COARE3_roughness(mdata.zref, umm, ustar);
            }
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_star_arr(i,j,k) = 0.0;
          olen_arr(i,j,k) = 1.0e16;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Adiabatic with modified charnock roughness
 */
struct adiabatic_mod_charnock
{
    adiabatic_mod_charnock (amrex::Real zref,
                            amrex::Real flux,
                            amrex::Real depth)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*q_star_arr*/,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real z0    = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_star_arr(i,j,k) = 0.0;
          olen_arr(i,j,k) = 1.0e16;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Adiabatic with Donelan roughness
 */
struct adiabatic_donelan
{
    adiabatic_donelan (amrex::Real zref,
                       amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*q_star_arr*/,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real z0    = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_star_arr(i,j,k) = 0.0;
          olen_arr(i,j,k) = 1.0e16;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Adiabatic with wave-coupled roughness
 */
struct adiabatic_wave_coupled
{
    adiabatic_wave_coupled (amrex::Real zref,
                            amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*tvm_arr*/,
                  const amrex::Array4<const amrex::Real>& /*qvm_arr*/,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& /*w_star_arr*/,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& /*q_star_arr*/,
                  const amrex::Array4<amrex::Real>& /*t_surf_arr*/,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& /*pblh_arr*/,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr) const
    {
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real z0    = 0.0;
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);

        t_star_arr(i,j,k) = 0.0;
          olen_arr(i,j,k) = 1.0e16;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with constant roughness
 */
struct surface_flux
{
    surface_flux (amrex::Real zref,
                  amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

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
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) - 0.61*tm_arr(i,j,k)*ustar*q_star_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        olen_arr(i,j,k)   = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with charnock roughness
 */
struct surface_flux_charnock
{
    surface_flux_charnock (amrex::Real zref,
                           amrex::Real flux,
                           amrex::Real cnk_a,
                           bool cnk_visc)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_a = cnk_a;
        mdata.visc = cnk_visc;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            if (mdata.Cnk_a > 0) {
                z0 = (mdata.Cnk_a / mdata.gravity) * ustar * ustar;
                if (mdata.visc) {
                    z0 += air_viscosity(tm_arr(i,j,k)) / std::max(ustar, 0.05);
                }
            } else {
                z0 = COARE3_roughness(mdata.zref, umm, ustar);
            }
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) - 0.61*tm_arr(i,j,k)*ustar*q_star_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
          olen_arr(i,j,k) = Olen;
           z0_arr(i,j,k)  = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with modified charnock roughness
 */
struct surface_flux_mod_charnock
{
    surface_flux_mod_charnock (amrex::Real zref,
                               amrex::Real flux,
                               amrex::Real depth)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) - 0.61*tm_arr(i,j,k)*ustar*q_star_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
          olen_arr(i,j,k) = Olen;
           z0_arr(i,j,k)  = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with donelan roughness
 */
struct surface_flux_donelan
{
    surface_flux_donelan (amrex::Real zref,
                          amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) - 0.61*tm_arr(i,j,k)*ustar*q_star_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
          olen_arr(i,j,k) = Olen;
           z0_arr(i,j,k)  = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with wave-coupled roughness
 */
struct surface_flux_wave_coupled
{
    surface_flux_wave_coupled (amrex::Real zref,
                               amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) - 0.61*tm_arr(i,j,k)*ustar*q_star_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(mdata.zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
          olen_arr(i,j,k) = Olen;
           z0_arr(i,j,k)  = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with constant roughness
 */
struct surface_temp
{
    surface_temp (amrex::Real zref,
                  amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

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
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            tflux += 0.61*tm_arr(i,j,k) * -ustar*q_star_arr(i,j,k); // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( ((Olen >= 0.0) && (Oleno <= 0.0)) ||
                 ((Olen <= 0.0) && (Oleno >= 0.0)) ) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0_arr(i,j,k)) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0_arr(i,j,k)) - psi_h);
        olen_arr(i,j,k)   = Olen;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with charnock roughness
 */
struct surface_temp_charnock
{
    surface_temp_charnock (amrex::Real zref,
                           amrex::Real flux,
                           amrex::Real cnk_a,
                           bool cnk_visc)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_a = cnk_a;
        mdata.visc = cnk_visc;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            if (mdata.Cnk_a > 0) {
                z0 = (mdata.Cnk_a / mdata.gravity) * ustar * ustar;
                if (mdata.visc) {
                    z0 += air_viscosity(tm_arr(i,j,k)) / std::max(ustar, 0.05);
                }
            } else {
                z0 = COARE3_roughness(mdata.zref, umm, ustar);
            }
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            tflux += 0.61*tm_arr(i,j,k) * -ustar*q_star_arr(i,j,k); // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( ((Olen >= 0.0) && (Oleno <= 0.0)) ||
                 ((Olen <= 0.0) && (Oleno >= 0.0)) ) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0) - psi_h);
          olen_arr(i,j,k) = Olen;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with modified charnock roughness
 */
struct surface_temp_mod_charnock
{
    surface_temp_mod_charnock (amrex::Real zref,
                               amrex::Real flux,
                               amrex::Real depth)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            tflux += 0.61*tm_arr(i,j,k) * -ustar*q_star_arr(i,j,k); // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( ((Olen >= 0.0) && (Oleno <= 0.0)) ||
                 ((Olen <= 0.0) && (Oleno >= 0.0)) ) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0) - psi_h);
          olen_arr(i,j,k) = Olen;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with donelan roughness
 */
struct surface_temp_donelan
{
    surface_temp_donelan (amrex::Real zref,
                          amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real umm = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            tflux += 0.61*tm_arr(i,j,k) * -ustar*q_star_arr(i,j,k); // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( ((Olen >= 0.0) && (Oleno <= 0.0)) ||
                 ((Olen <= 0.0) && (Oleno >= 0.0)) ) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0) - psi_h);
          olen_arr(i,j,k) = Olen;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with wave-coupled roughness
 */
struct surface_temp_wave_coupled
{
    surface_temp_wave_coupled (amrex::Real zref,
                               amrex::Real flux)
    {
        mdata.zref = zref;
        mdata.surf_temp_flux = flux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
                  const amrex::Array4<const amrex::Real>& tvm_arr,
                  const amrex::Array4<const amrex::Real>& qvm_arr,
                  const amrex::Array4<amrex::Real>& u_star_arr,
                  const amrex::Array4<amrex::Real>& w_star_arr,
                  const amrex::Array4<amrex::Real>& t_star_arr,
                  const amrex::Array4<amrex::Real>& q_star_arr,
                  const amrex::Array4<amrex::Real>& t_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr) const
    {
        int iter = 0;
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(mdata.zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(mdata.zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            tflux += 0.61*tm_arr(i,j,k) * -ustar*q_star_arr(i,j,k); // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( ((Olen >= 0.0) && (Oleno <= 0.0)) ||
                 ((Olen <= 0.0) && (Oleno >= 0.0)) ) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = mdata.zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(mdata.zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ASSERT_WITH_MESSAGE(iter < max_iters, "Maximum number of MOST iterations reached.");

        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(mdata.zref / z0) - psi_h);
          olen_arr(i,j,k) = Olen;
            z0_arr(i,j,k) = z0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Moeng flux formulation
 */
struct moeng_flux
{
    moeng_flux () {}

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
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
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
    const amrex::Real     eps = 1e-15;
    const amrex::Real   WSMIN = 0.1; // minimum wind speed
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
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
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
    custom_flux () {}

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
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        amrex::Real rho   =   cons_arr(i,j,k,Rho_comp);
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
        amrex::Real rho   =   cons_arr(i,j,k,Rho_comp);
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
        amrex::Real rho   = 0.5  * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

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
        amrex::Real rho   = 0.5  * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );

        amrex::Real ustar   = 0.5 * ( u_star_arr(i,j-1,k) + u_star_arr(i,j,k) );
        amrex::Real wsp     = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        amrex::Real stressy = -rho * ustar * ustar * vely / wsp;

        return stressy;
    }

private:
    const amrex::Real     eps = 1e-15;

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
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*q_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*t_surf_arr*/) const
    {
        // NOTE: this is the total stress
        amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        amrex::Real qstar  = 0.0; //q_star_arr(i,j,k);
        amrex::Real moflux = (std::abs(qstar) > eps) ? -rho * qstar : 0.0;

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
