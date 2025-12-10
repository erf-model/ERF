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
    amrex::Real z0_const{0.1};        ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};         ///< von Karman constant
    amrex::Real gravity{CONST_GRAV};  ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{0.0};  ///< Heat flux TODO: decide whether this is <θ'w'> or <θv'w'> under moist conditions
    amrex::Real surf_moist_flux{0.0}; ///< Moisture flux

    amrex::Real Cnk_a{0.0185};        ///< Standard Charnock constant https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b1{1.0/30.0};     ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b2{1260.0};       ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_d{30.0};          ///< Modified Charnock Eq (4) https://doi.org/10.1175/JAMC-D-17-0137.1
    amrex::Real Cnk_b;
    bool visc{false};                 ///< Use viscous Charnock formulation

    amrex::Real Cd{0.001};
    amrex::Real Ch{0.001};
    amrex::Real Cq{0.001};

    const amrex::Real Bjr_beta = 1.2; // Empirical parameter from Beljaars 1995 QJRMS
};


/**
 * Structure of similarity functions for Moeng formulation
 */
struct similarity_funs
{
    //
    // Similarity functions from Jimenez (2012)
    //
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(1.0 + std::pow(zeta, 2.5), 1.0/2.5);
            return ( -6.1*std::log(zeta + x) );
        } else {
            amrex::Real x = std::pow(1.0 - 16.0*zeta, 0.25);
            amrex::Real psi_k_m =  2.0 * std::log(0.5 * (1.0 + x)) + log(0.5 * (1.0 + x * x)) -
                                   2.0 * std::atan(x) + PIoTwo;
            amrex::Real y = std::pow(1.0 - 10.0*zeta, 1.0/3.0);
            amrex::Real psi_c_m = (3.0/2.0)*std::log((y*y + y + 1.0)/3.0)
                                - std::sqrt(3.0)*std::atan((2.0*y + 1.0)/std::sqrt(3.0))
                                + PI/std::sqrt(3.0);
            return ( (psi_k_m + zeta*zeta*psi_c_m) / (1. + zeta*zeta) );
        }
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(1.0 + std::pow(zeta, 1.1), 1.0/1.1);
            return ( -6.1*std::log(zeta + x) );
        } else {
            amrex::Real x = std::sqrt(1.0 - 16.0*zeta);
            amrex::Real psi_k_h = 2.0 * std::log(0.5 * (1.0 + x));
            amrex::Real y = std::pow(1.0 - 34.0*zeta, 1.0/3.0);
            amrex::Real psi_c_h = (3.0/2.0)*std::log((y*y + y + 1.0)/3.0)
                                - std::sqrt(3.0)*std::atan((2.0*y + 1.0)/std::sqrt(3.0))
                                + PI/std::sqrt(3.0);
            return ( (psi_k_h + zeta*zeta*psi_c_h) / (1. + zeta*zeta) );
        }
    }


    //
    // Similarity functions from Businger & Dyer (1966)
    //
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
    adiabatic (amrex::Real Tflux,
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
        olen_arr(i,j,k)   = 1.0e16;
        u_star_arr(i,j,k) = mdata.kappa * umm_arr(i,j,k) / std::log(zref_arr(i,j,k) / z0_arr(i,j,k));
        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
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
    adiabatic_charnock (amrex::Real Tflux,
                        amrex::Real Qvflux,
                        amrex::Real cnk_a,
                        bool cnk_visc)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
                  const amrex::Array4<const amrex::Real>& umm_arr,
                  const amrex::Array4<const amrex::Real>& tm_arr,
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
        amrex::Real ustar  = 0.0;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = 0.0;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            if (mdata.Cnk_a > 0) {
                z0 = (mdata.Cnk_a / mdata.gravity) * ustar * ustar;
                if (mdata.visc) {
                    z0 += air_viscosity(tm_arr(i,j,k)) / std::max(ustar, 0.05);
                }
            } else {
                z0 = COARE3_roughness(zref, umm, ustar);
            }

            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = 1.0e16;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol_z = 1.0e-10;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Adiabatic with modified charnock roughness
 */
struct adiabatic_mod_charnock
{
    adiabatic_mod_charnock (amrex::Real Tflux,
                            amrex::Real Qvflux,
                            amrex::Real depth)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
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
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
        amrex::Real ustar  = 0.0;
        amrex::Real zref   = zref_arr(i,j,k);
        amrex::Real z0     = z0_arr(i,j,k);
        amrex::Real z0_old = z0;
        amrex::Real psi_m  = 0.0;
        amrex::Real umm    = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C      = std::log(zref / z0);

        // Fixed point iteration on roughness
        int iter_z = 0;
        do {
            // Transfer curr to old
            z0_old = z0;

            // Update
            C = std::log(zref / z0_old);
            ustar = mdata.kappa * umm / (C - psi_m);
            z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );

            ++iter_z;
        } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                         "Maximum number of MOST roughness iterations reached.");
        C = std::log(zref / z0);

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = 1.0e16;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol_z = 1.0e-10;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Adiabatic with Donelan roughness
 */
struct adiabatic_donelan
{
    adiabatic_donelan (amrex::Real Tflux,
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
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real z0    = 0.0;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = 1.0e16;
        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
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
    adiabatic_wave_coupled (amrex::Real Tflux,
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
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
                  const amrex::Array4<amrex::Real>& z0_arr,
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
                  const amrex::Array4<amrex::Real>& Hwave_arr,
                  const amrex::Array4<amrex::Real>& Lwave_arr,
                  const amrex::Array4<amrex::Real>& eta_arr) const
    {
        int iter = 0;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real ustar = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real z0    = 0.0;
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = 1.0e16;
        t_star_arr(i,j,k) = 0.0;
        q_star_arr(i,j,k) = 0.0;
    }

private:
    most_data mdata;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN  = 0.1; // minimum wind speed
};


/**
 * Surface flux with constant roughness
 */
struct surface_flux
{
    surface_flux (amrex::Real Tflux,
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
        amrex::Real ustar = 0.0;
        amrex::Real wstar = 0.0;
        amrex::Real tflux = 0.0;
        amrex::Real qflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
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
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
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
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with charnock roughness
 */
struct surface_flux_charnock
{
    surface_flux_charnock (amrex::Real Tflux,
                           amrex::Real Qvflux,
                           amrex::Real cnk_a,
                           bool cnk_visc,
                           bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
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
                z0 = COARE3_roughness(zref, umm, ustar);
            }
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with modified charnock roughness
 */
struct surface_flux_mod_charnock
{
    surface_flux_mod_charnock (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               amrex::Real depth,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
        spec_qflux  = cons_qflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with donelan roughness
 */
struct surface_flux_donelan
{
    surface_flux_donelan (amrex::Real Tflux,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface flux with wave-coupled roughness
 */
struct surface_flux_wave_coupled
{
    surface_flux_wave_coupled (amrex::Real Tflux,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real z0    = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux = mdata.surf_temp_flux*(1 + 0.61*qvm_arr(i,j,k)) + qflux*0.61*tm_arr(i,j,k);
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen  = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_surf_arr(i,j,k) = mdata.surf_temp_flux * (std::log(zref / z0) - psi_h) /
                            (u_star_arr(i,j,k) * mdata.kappa) + tm_arr(i,j,k);
        t_star_arr(i,j,k) = -mdata.surf_temp_flux / u_star_arr(i,j,k);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN  = 0.1; // minimum wind speed
};


/**
 * Surface temperature with constant roughness
 */
struct surface_temp
{
    surface_temp (amrex::Real Tflux,
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
        amrex::Real Rib      = 0.0;
        amrex::Real zeta     = 0.0;
        amrex::Real zeta_old = 0.0;
        amrex::Real psi_m    = 0.0;
        amrex::Real psi_h    = 0.0;
        amrex::Real num      = 0.0;
        amrex::Real den      = 0.0;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == 1.E34) { olen_arr(i,j,k) = 1.0e3; }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > 0.) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
            amrex::Real ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(1 + 0.61*qvm_arr(i,j,k)) + 0.61*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (1.0 + 0.61*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (1.0 + 0.61*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-4.0),4.0);

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Limiting
            num = std::max(C - psi_m, 1.0);
            den = std::max(C - psi_h, 1.0);

            // Update with under relaxation
            zeta = (1.0 - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = 1.0e-3;
    const amrex::Real alpha = 0.5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with charnock roughness
 */
struct surface_temp_charnock
{
    surface_temp_charnock (amrex::Real Tflux,
                           amrex::Real Qvflux,
                           amrex::Real cnk_a,
                           bool cnk_visc,
                           bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_a = cnk_a;
        mdata.visc  = cnk_visc;
        spec_qflux  = cons_qflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = 0.0;
        amrex::Real zeta     = 0.0;
        amrex::Real zeta_old = 0.0;
        amrex::Real psi_m    = 0.0;
        amrex::Real psi_h    = 0.0;
        amrex::Real num      = 0.0;
        amrex::Real den      = 0.0;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = 0.0;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == 1.E34) { olen_arr(i,j,k) = 1.0e3; }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > 0.) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
                        ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(1 + 0.61*qvm_arr(i,j,k)) + 0.61*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (1.0 + 0.61*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (1.0 + 0.61*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-4.0),4.0);

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to old
                z0_old = z0;

                // Update
                C = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                if (mdata.Cnk_a > 0) {
                    z0 = (mdata.Cnk_a / mdata.gravity) * ustar * ustar;
                    if (mdata.visc) {
                        z0 += air_viscosity(tm_arr(i,j,k)) / std::max(ustar, 0.05);
                    }
                } else {
                    z0 = COARE3_roughness(zref, umm, ustar);
                }

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, 1.0);
            den = std::max(C - psi_h, 1.0);

            // Update with under relaxation
            zeta = (1.0 - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = 1.0e-3;
    const amrex::Real tol_z = 1.0e-10;
    const amrex::Real alpha = 0.5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with modified charnock roughness
 */
struct surface_temp_mod_charnock
{
    surface_temp_mod_charnock (amrex::Real Tflux,
                               amrex::Real Qvflux,
                               amrex::Real depth,
                               bool cons_qflux)
    {
        mdata.surf_temp_flux  = Tflux;
        mdata.surf_moist_flux = Qvflux;
        mdata.Cnk_d = depth;
        mdata.Cnk_b = mdata.Cnk_b1 * std::log(mdata.Cnk_b2 / mdata.Cnk_d);
        spec_qflux  = cons_qflux;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    iterate_flux (const int& i,
                  const int& j,
                  const int& k,
                  const int& max_iters,
                  const amrex::Array4<const amrex::Real>& zref_arr,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
                  const amrex::Array4<amrex::Real>& olen_arr,
                  const amrex::Array4<amrex::Real>& pblh_arr,
                  const amrex::Array4<amrex::Real>& /*Hwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*Lwave_arr*/,
                  const amrex::Array4<amrex::Real>& /*eta_arr*/) const
    {
        amrex::Real Rib      = 0.0;
        amrex::Real zeta     = 0.0;
        amrex::Real zeta_old = 0.0;
        amrex::Real psi_m    = 0.0;
        amrex::Real psi_h    = 0.0;
        amrex::Real num      = 0.0;
        amrex::Real den      = 0.0;
        amrex::Real z0       = z0_arr(i,j,k);
        amrex::Real z0_old   = z0;
        amrex::Real ustar    = 0.0;
        amrex::Real zref     = zref_arr(i,j,k);
        amrex::Real umm      = std::max(umm_arr(i,j,k), WSMIN);
        amrex::Real C        = std::log(zref / z0);

        // First iteration we assume neutral (L -> inf)
        if (u_star_arr(i,j,k) == 1.E34) { olen_arr(i,j,k) = 1.0e3; }
        zeta  = zref / olen_arr(i,j,k);

        // Water vapor in atmos and surface
        amrex::Real qv_s, qv_a;
        if (q_surf_arr(i,j,k) > 0.) {
            qv_s = q_surf_arr(i,j,k);
        } else {
            // First iteration and no qv_surf was specified
            // Use mean since there will be no flux
            qv_s = qvm_arr(i,j,k);
        }
        qv_a = qvm_arr(i,j,k);

        // update w* and Umagmean from Beljaars (1995)
        if (w_star_arr) {
            // NOTE: Thv flux is lagged, similar to WRF
            psi_m = sfuns.calc_psi_m2(zeta);
            psi_h = sfuns.calc_psi_h2(zeta);
                        ustar = mdata.kappa * umm / (C - psi_m);
            amrex::Real tstar = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
            amrex::Real qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -ustar * mdata.kappa * (qv_a - qv_s) / (C - psi_h);
            amrex::Real tflux = -ustar*tstar*(1 + 0.61*qvm_arr(i,j,k)) + 0.61*tm_arr(i,j,k)*qflux;
            w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
            amrex::Real wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
            umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
            umm = std::max(umm, WSMIN);
        }

        // Bulk Richardson number w/ moisture
        amrex::Real thv_s = t_surf_arr(i,j,k) * (1.0 + 0.61*qv_s);
        amrex::Real thv_a =     tm_arr(i,j,k) * (1.0 + 0.61*qv_a);
        Rib = ( (mdata.gravity * zref) / tm_arr(i,j,k) ) *
              ( (thv_a - thv_s) / (umm * umm) );
        Rib = std::min(std::max(Rib,-4.0),4.0);

        // Fixed point iteration on zeta
        int iter = 0;
        do {
            // Transfer curr to old
            zeta_old = zeta;

            // Stability functions
            psi_m = sfuns.calc_psi_m2(zeta_old);
            psi_h = sfuns.calc_psi_h2(zeta_old);

            // Fixed point iteration on roughness
            int iter_z = 0;
            do {
                // Transfer curr to old
                z0_old = z0;

                // Update
                C = std::log(zref / z0_old);
                ustar = mdata.kappa * umm / (C - psi_m);
                z0    = std::exp( (2.7*ustar - 1.8/mdata.Cnk_b) / (ustar + 0.17/mdata.Cnk_b) );

                ++iter_z;
            } while ( (std::abs(z0 - z0_old) > tol_z) && (iter_z <= max_iters) );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter_z < max_iters,
                                             "Maximum number of MOST roughness iterations reached.");
            C = std::log(zref / z0);

            // Limiting
            num = std::max(C - psi_m, 1.0);
            den = std::max(C - psi_h, 1.0);

            // Update with under relaxation
            zeta = (1.0 - alpha) * zeta_old + alpha * Rib * num * num / den;

            ++iter;
        } while ( (std::abs(zeta - zeta_old) > tol) && (iter <= max_iters) );
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = zref / zeta;
        u_star_arr(i,j,k) = mdata.kappa * umm / (C - psi_m);
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) / (C - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (C - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) / (C - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol   = 1.0e-3;
    const amrex::Real tol_z = 1.0e-10;
    const amrex::Real alpha = 0.5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with donelan roughness
 */
struct surface_temp_donelan
{
    surface_temp_donelan (amrex::Real Tflux,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        amrex::Real umm = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = Donelan_roughness(ustar);
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux += 0.61*tm_arr(i,j,k) * qflux; // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( (((Olen >= 0.0) && (Oleno <= 0.0)) ||
                  ((Olen <= 0.0) && (Oleno >= 0.0))) &&
                 std::fabs(Olen) + std::fabs(Oleno) < 1.) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(zref / z0) - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real WSMIN = 0.1; // minimum wind speed
};


/**
 * Surface temperature with wave-coupled roughness
 */
struct surface_temp_wave_coupled
{
    surface_temp_wave_coupled (amrex::Real Tflux,
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
                  const amrex::Array4<amrex::Real>& q_surf_arr,
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
        amrex::Real qflux = 0.0;
        amrex::Real zeta  = 0.0;
        amrex::Real psi_m = 0.0;
        amrex::Real psi_h = 0.0;
        amrex::Real Olen  = 0.0;
        amrex::Real Oleno = 0.0;
        amrex::Real zref  = zref_arr(i,j,k);
        int ie, je;
        ie = i  < lbound(eta_arr).x ? lbound(eta_arr).x : i;
        je = j  < lbound(eta_arr).y ? lbound(eta_arr).y : j;
        ie = ie > ubound(eta_arr).x ? ubound(eta_arr).x : ie;
        je = je > ubound(eta_arr).y ? ubound(eta_arr).y : je;
        amrex::Real umm   = std::max(umm_arr(i,j,k), WSMIN);
        if (u_star_arr(i,j,k) == 1.E34) {
            u_star_arr(i,j,k) = mdata.kappa * umm / std::log(zref / z0_arr(i,j,k));
        } else {
            Olen  = olen_arr(i,j,k);
            Oleno = Olen;
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
        }
        do {
            ustar = u_star_arr(i,j,k);
            z0    = std::min( std::max(1200.0 * Hwave_arr(i,j,k) * std::pow( Hwave_arr(i,j,k)/(Lwave_arr(i,j,k)+eps), 4.5 )
                                      + 0.11 * eta_arr(ie,je,k,EddyDiff::Mom_v) / ustar, z0_eps), z0_max );
            tflux = -(tm_arr(i,j,k) - t_surf_arr(i,j,k)) * ustar * mdata.kappa /
                     (std::log(zref / z0) - psi_h); // <w'T'>
            tflux *= (1 + 0.61*qvm_arr(i,j,k));
            qflux = (spec_qflux) ? mdata.surf_moist_flux :
                                 -(qvm_arr(i,j,k) - q_surf_arr(i,j,k)) * ustar * mdata.kappa /
                                  (std::log(zref / z0_arr(i,j,k)) - psi_h); // <w'Qv'>
            tflux += 0.61*tm_arr(i,j,k) * qflux; // ~= <w'Tv'>
            if (w_star_arr) {
                // update w* and Umagmean
                w_star_arr(i,j,k) = calc_wstar(tflux, pblh_arr(i,j,k), tvm_arr(i,j,k));
                wstar = mdata.Bjr_beta * w_star_arr(i,j,k);
                umm = std::sqrt(umm_arr(i,j,k)*umm_arr(i,j,k) + wstar*wstar);
                umm = std::max(umm, WSMIN);
            }
            Olen = -ustar * ustar * ustar * tvm_arr(i,j,k) / (mdata.kappa * mdata.gravity * tflux);
            if ( (((Olen >= 0.0) && (Oleno <= 0.0)) ||
                  ((Olen <= 0.0) && (Oleno >= 0.0))) &&
                 std::fabs(Olen) + std::fabs(Oleno) < 1.) {
                Olen = 0.5 * (Olen + Oleno);
            }
            Oleno = Olen;
            zeta  = zref / Olen;
            psi_m = sfuns.calc_psi_m(zeta);
            psi_h = sfuns.calc_psi_h(zeta);
            u_star_arr(i,j,k) = mdata.kappa * umm / (std::log(zref / z0) - psi_m);
            ++iter;
        } while ((std::abs(u_star_arr(i,j,k) - ustar) > tol) && iter <= max_iters);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(iter < max_iters,
                                         "Maximum number of MOST iterations reached.");

        // Populate the stored MOST arrays
        z0_arr(i,j,k)     = z0;
        olen_arr(i,j,k)   = Olen;
        t_star_arr(i,j,k) = mdata.kappa * (tm_arr(i,j,k) - t_surf_arr(i,j,k)) /
                            (std::log(zref / z0) - psi_h);
        if (spec_qflux) {
            q_surf_arr(i,j,k) = mdata.surf_moist_flux * (std::log(zref / z0_arr(i,j,k)) - psi_h) /
                                (u_star_arr(i,j,k) * mdata.kappa) + qvm_arr(i,j,k);
            q_star_arr(i,j,k) = -mdata.surf_moist_flux / u_star_arr(i,j,k);
        } else {
            q_star_arr(i,j,k) = mdata.kappa * (qvm_arr(i,j,k) - q_surf_arr(i,j,k)) /
                                (std::log(zref / z0_arr(i,j,k)) - psi_h);
        }
    }

private:
    most_data mdata;
    bool spec_qflux;
    similarity_funs sfuns;
    const amrex::Real tol = 1.0e-5;
    const amrex::Real eps = 1e-15;
    const amrex::Real z0_eps = 1.0e-6;
    const amrex::Real z0_max = 0.1;
    const amrex::Real WSMIN  = 0.1; // minimum wind speed
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
 * Bulk coefficient flux formulation
 */
struct bulk_coeff_flux
{
    bulk_coeff_flux (amrex::Real m_Cd,
                     amrex::Real m_Ch,
                     amrex::Real m_Cq)
    {
        mdata.Cd = m_Cd;
        mdata.Ch = m_Ch;
        mdata.Cq = m_Cq;
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    compute_q_flux (const int& i,
                    const int& j,
                    const int& k,
                    const amrex::Array4<const amrex::Real>& cons_arr,
                    const amrex::Array4<const amrex::Real>& velx_arr,
                    const amrex::Array4<const amrex::Real>& vely_arr,
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*qm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*q_star_arr*/,
                    const amrex::Array4<const amrex::Real>& q_surf_arr) const
    {
        amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        amrex::Real qv     = cons_arr(i,j,k,RhoQ1_comp)/cons_arr(i,j,k,Rho_comp);
        amrex::Real qvsurf = q_surf_arr(i,j,k);
        amrex::Real velx   = 0.5 * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely   = 0.5 * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );
        amrex::Real wsp    = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<q'w'> = -K dqdz
        amrex::Real moflux  = rho * mdata.Cq * wsp * (qvsurf - qv);

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
                    const amrex::Array4<const amrex::Real>& /*umm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*tm_arr*/,
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/,
                    const amrex::Array4<const amrex::Real>& /*t_star_arr*/,
                    const amrex::Array4<const amrex::Real>& t_surf_arr) const
    {
        amrex::Real rho    = cons_arr(i,j,k,Rho_comp);
        amrex::Real th     = cons_arr(i,j,k,RhoTheta_comp)/cons_arr(i,j,k,Rho_comp);
        amrex::Real thsurf = t_surf_arr(i,j,k);
        amrex::Real velx   = 0.5 * ( velx_arr(i,j,k) + velx_arr(i+1,j  ,k) );
        amrex::Real vely   = 0.5 * ( vely_arr(i,j,k) + vely_arr(i  ,j+1,k) );
        amrex::Real wsp    = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<T'w'> = -K dTdz
        amrex::Real moflux  = rho * mdata.Ch * wsp * (thsurf - th);

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
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real rho   = 0.5 * ( cons_arr(i-1,j,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        amrex::Real velx  = velx_arr(i,j,k);
        amrex::Real vely  = 0.25 * ( vely_arr(i  ,j,k) + vely_arr(i  ,j+1,k)
                                   + vely_arr(i-1,j,k) + vely_arr(i-1,j+1,k) );
        amrex::Real wsp   = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<u'w'> = -K dudz
        // NOTE: tau_tot = rho * Cd * wsp^2, multiply by u/wsp to get tau_13
        amrex::Real stressx = -rho * mdata.Cd * wsp * velx;

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
                    const amrex::Array4<const amrex::Real>& /*u_star_arr*/) const
    {
        amrex::Real rho   = 0.5 * ( cons_arr(i,j-1,k,Rho_comp) + cons_arr(i,j,k,Rho_comp) );
        amrex::Real velx  = 0.25 * ( velx_arr(i,j  ,k) + velx_arr(i+1,j  ,k)
                                   + velx_arr(i,j-1,k) + velx_arr(i+1,j-1,k) );
        amrex::Real vely  = vely_arr(i,j,k);
        amrex::Real wsp   = sqrt(velx*velx+vely*vely);

        // NOTE: this is rho*<v'w'> = -K dvdz
        // NOTE: tau_tot = rho * Cd * wsp^2, multiply by v/wsp to get tau_13
        amrex::Real stressy = -rho * mdata.Cd * wsp * vely;

        return stressy;
    }

private:
    most_data mdata;
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
