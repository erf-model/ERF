#ifndef ERF_MOSTUtils_H
#define ERF_MOSTUtils_H

#include "ERF_Constants.H"

/**
 * Structure of plain old data relevant to MOST BCs
 */
struct most_data
{
public:
    amrex::Real z0_const{amrex::Real(0.1)};    ///< Roughness height -- default constant value(m)
    amrex::Real kappa{KAPPA};                  ///< von Karman constant
    amrex::Real gravity{CONST_GRAV};           ///< Acceleration due to gravity (m/s^2)
    amrex::Real surf_temp_flux{zero};          ///< Heat flux TODO: decide whether this is <θ'w'> or <θv'w'> under moist conditions
    amrex::Real surf_moist_flux{zero};         ///< Moisture flux

    amrex::Real Cnk_a{amrex::Real(0.0185)};    ///< Standard Charnock constant https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b1{one/amrex::Real(30.0)}; ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b2{amrex::Real(1260.0)};   ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_d{amrex::Real(30.0)};      ///< Modified Charnock Eq (4) https://doi.org/amrex::Real(10.1175)/JAMC-D-17-amrex::Real(0137.1)
    amrex::Real Cnk_b;                         ///< Modified Charnock derived quantity
    bool visc{true};                           ///< Use viscous Charnock formulation

    amrex::Real Cd{amrex::Real(0.001)};
    amrex::Real Ch{amrex::Real(0.001)};
    amrex::Real Cq{amrex::Real(0.001)};

    const amrex::Real Bjr_beta = amrex::Real(1.2); // Empirical parameter from Beljaars 1995 QJRMS
};


/**
 * Structure of similarity functions for Moeng formulation
 */
struct similarity_funs
{
    //
    // Similarity functions from Jimenez (2012)
    //
    /**
     * Compute the Jimenez momentum stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(one + std::pow(zeta, amrex::Real(2.5)), one/amrex::Real(2.5));
            return ( -amrex::Real(6.1)*std::log(zeta + x) );
        } else {
            amrex::Real x = std::pow(one - amrex::Real(16.0)*zeta, fourth);
            amrex::Real psi_k_m =  two * std::log(myhalf * (one + x)) + std::log(myhalf * (one + x * x)) -
                                   two * std::atan(x) + PIoTwo;
            amrex::Real y = std::pow(one - amrex::Real(10.0)*zeta, one/three);
            amrex::Real psi_c_m = (three/two)*std::log((y*y + y + one)/three)
                                - std::sqrt(three)*std::atan((two*y + one)/std::sqrt(three))
                                + PI/std::sqrt(three);
            return ( (psi_k_m + zeta*zeta*psi_c_m) / (one + zeta*zeta) );
        }
    }

    /**
     * Compute the Jimenez heat stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h2 (amrex::Real zeta) const
    {
        if (zeta > 0) {
            amrex::Real x = std::pow(one + std::pow(zeta, amrex::Real(1.1)), one/amrex::Real(1.1));
            return ( -amrex::Real(5.3)*std::log(zeta + x) );
        } else {
            amrex::Real x = std::sqrt(one - amrex::Real(16.0)*zeta);
            amrex::Real psi_k_h = two * std::log(myhalf * (one + x));
            amrex::Real y = std::pow(one - amrex::Real(34.0)*zeta, one/three);
            amrex::Real psi_c_h = (three/two)*std::log((y*y + y + one)/three)
                                - std::sqrt(three)*std::atan((two*y + one)/std::sqrt(three))
                                + PI/std::sqrt(three);
            return ( (psi_k_h + zeta*zeta*psi_c_h) / (one + zeta*zeta) );
        }
    }


    //
    // Similarity functions from Businger & Dyer (1966)
    //
    /**
     * Compute the Businger-Dyer momentum stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_m (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_m * zeta;
        } else {
            amrex::Real x = std::sqrt(std::sqrt(one - gamma_m * zeta));
            return two * std::log(myhalf * (one + x)) + std::log(myhalf * (one + x * x)) -
                   two * std::atan(x) + PIoTwo;
        }
    }

    /**
     * Compute the Businger-Dyer heat stability correction.
     *
     * @param[in] zeta stability coordinate z/L
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real
    calc_psi_h (amrex::Real zeta) const
    {
        if (zeta > 0) {
            return -beta_h * zeta;
        } else {
            amrex::Real x = std::sqrt(one - gamma_h * zeta);
            return two * std::log(myhalf * (one + x));
        }
    }

private:
    amrex::Real beta_m{amrex::Real(5.0)};         ///< Constants from Dyer, BLM, 1974
    amrex::Real beta_h{amrex::Real(5.0)};         ///< https://doi.org/amrex::Real(10.1007)/BF00240838
    amrex::Real gamma_m{amrex::Real(16.0)};
    amrex::Real gamma_h{amrex::Real(16.0)};
};


/**
 * Empirical kinematic viscosity [m2/s] formula from Andreas (1989) CRREL Rep.
 * 89-11, valid between -173 and 277 deg C.
 *
 * @param[in] T_degK air temperature in Kelvin
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
air_viscosity (amrex::Real T_degK)
{
    amrex::Real TC = T_degK - amrex::Real(273.15);
    return amrex::Real(1.326e-5)*(one + amrex::Real(6.542e-3)*TC + amrex::Real(8.301e-6)*TC*TC - amrex::Real(4.84e-9)*TC*TC*TC);
}

/**
 * Charnock roughness with optional aerodynamically-smooth contribution
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] ustar friction velocity
 * @param[in] cnk_a charnock constant
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Charnock_roughness (amrex::Real ustar,
                    amrex::Real Cnk_a,
                    amrex::Real nu,
                    bool visc)
{
    amrex::Real z0 = (Cnk_a / CONST_GRAV) * ustar * ustar;
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}

/**
 * Modified Charnock roughness
 *
 * Smooth-flow limit is intrinsic.
 *
 * @param[in] ustar friction velocity
 * @param[in] cnk_b modified charnock constant
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mod_Charnock_roughness (amrex::Real ustar,
                        amrex::Real Cnk_b)
{
    amrex::Real z0 = std::exp( (amrex::Real(2.7)*ustar - amrex::Real(1.8)/Cnk_b) / (ustar + amrex::Real(0.17)/Cnk_b) );
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(1.0e-2));
}

/**
 * COARE3.0, following WRF phys/module_sf_mynn.F
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] zref  reference height for the wind speed
 * @param[in] umm   horizontal wind speed magnitude at zref
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
COARE3_roughness (amrex::Real zref,
                  amrex::Real umm,
                  amrex::Real ustar,
                  amrex::Real nu,
                  bool visc)
{
    amrex::Real ws10m = umm*std::log(amrex::Real(10.)/amrex::Real(1e-4)) / std::log(zref/amrex::Real(1e-4));
    amrex::Real Czc = amrex::Real(0.011) + amrex::Real(0.007)*std::min(std::max((ws10m-amrex::Real(10.))/amrex::Real(8.), amrex::Real(0.0)), amrex::Real(1.0));
    amrex::Real z0 = (Czc / CONST_GRAV) * ustar * ustar;
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}

/**
 * Updated Donelan 2004 as found in AHW
 *
 * Matches implementations in WRF phys/physics_mmm/sf_sfclayrev.F90 and
 * phys/module_sf_mynn.F
 *
 * Smooth flow limit is intrinsic in zn2.
 *
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Donelan_roughness (amrex::Real ustar,
                   amrex::Real nu,
                   bool visc)
{
    constexpr amrex::Real ozo{amrex::Real(1.59e-5)};
    amrex::Real zw  = std::min(std::pow(ustar/amrex::Real(1.06), amrex::Real(0.3)),amrex::Real(1.0));
    amrex::Real zn1 = amrex::Real(0.011)*ustar*ustar/CONST_GRAV + ozo;
    amrex::Real zn2 = amrex::Real(10.)*std::exp(-amrex::Real(9.5)*std::pow(ustar, -amrex::Real(.3333)));
    if (visc) {
        zn2 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    amrex::Real z0 = (amrex::Real(1.0)-zw) * zn1 + zw * zn2;
    // From Davis et al 2008 MWR
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}


/**
 * Wave coupled roughness
 *
 * Smooth-flow limit follows z0 = 0.11 nu / u* (Smith 1988; ECMWF; WRF)
 *
 * @param[in] Hwave wave height
 * @param[in] Lwave wave length
 * @param[in] eps   avoid divide by zero
 * @param[in] ustar friction velocity
 * @param[in] nu    intrinsic fluid viscosity
 * @param[in] visc  flag to apply smooth flow limit
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
WaveCoupled_roughness (amrex::Real Hwave,
                       amrex::Real Lwave,
                       amrex::Real eps,
                       amrex::Real ustar,
                       amrex::Real nu,
                       bool visc)
{

    amrex::Real z0 = amrex::Real(1200.0) * Hwave * std::pow( Hwave/(Lwave+eps), amrex::Real(4.5) );
    if (visc) {
        z0 += amrex::Real(0.11) * nu / std::max(ustar, amrex::Real(0.01));
    }
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(1.0e-2));
}


/**
 * Invert log-law to get roughness
 *
 * @param[in] zref  reference height for the wind speed
 * @param[in] Olen  Obukhov length
 * @param[in] umm   horizontal wind speed magnitude at zref
 * @param[in] ustar friction velocity
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Compute_roughness (amrex::Real zref,
                   amrex::Real Olen,
                   amrex::Real umm,
                   amrex::Real ustar)
{
    similarity_funs sfuns;
    amrex::Real zeta  = zref / Olen;
    amrex::Real psi_m = sfuns.calc_psi_m(zeta);
    amrex::Real z0    = zref / std::exp(KAPPA*umm/ustar + psi_m);
    return std::max(z0, amrex::Real(1.0e-7));
}
#endif
