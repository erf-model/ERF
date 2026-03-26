#ifndef ERF_STORM_DIAGNOSTICS_H_
#define ERF_STORM_DIAGNOSTICS_H_

#include "ERF_DataStruct.H"

/**
 * Routines to compute strom diagnostics
 */


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real compute_max_reflectivity_dbz (amrex::Real rho_air, amrex::Real tmk,
              amrex::Real qra, amrex::Real qsn, amrex::Real qgr,
              int in0r, int in0s, int in0g, int iliqskin)
{
    using namespace amrex;


    // ----------------------------------------------------------------
    // Constants
    // ----------------------------------------------------------------
    constexpr amrex::Real pi = Real(3.14159265358979323846);
    constexpr amrex::Real rhowat = Real(1000.0);   // liquid water density [kg/m^3]
    constexpr amrex::Real rho_r  = Real(1000.0);   // rain
    constexpr amrex::Real rho_s  = Real(100.0);    // snow
    constexpr amrex::Real rho_g  = Real(400.0);    // graupel
    constexpr amrex::Real celkel = Real(273.15);   // 0°C in Kelvin
    constexpr amrex::Real alpha  = Real(0.224);
    constexpr amrex::Real gamma_seven = Real(720.0);

    // Constant intercepts
    constexpr amrex::Real rn0_r = Real(8.e6);
    constexpr amrex::Real rn0_s = Real(2.e7);
    constexpr amrex::Real rn0_g = Real(4.e6);

    // Variable intercept constants
    constexpr amrex::Real r1 = Real(1.e-15);
    constexpr amrex::Real ron2 = Real(1.e10);
    constexpr amrex::Real gon  = Real(5.e7);
    constexpr amrex::Real ron_min    = Real(8.e6);
    constexpr amrex::Real ron_qr0    = Real(1.0e-4);
    constexpr amrex::Real ron_delqr0 = fourth * ron_qr0;
    constexpr amrex::Real ron_const1r = (ron2 - ron_min) * Real(0.5);
    constexpr amrex::Real ron_const2r = (ron2 + ron_min) * Real(0.5);

    // Precompute factors
    const amrex::Real factor_r = gamma_seven * Real(1.e18) * std::pow(one / (pi * rho_r), Real(1.75));
    const amrex::Real factor_s = gamma_seven * Real(1.e18) * std::pow(one / (pi * rho_s), Real(1.75))
                          * std::pow(rho_s / rhowat, two) * alpha;
    const amrex::Real factor_g = gamma_seven * Real(1.e18) * std::pow(one / (pi * rho_g), Real(1.75))
                          * std::pow(rho_g / rhowat, two) * alpha;

    // ----------------------------------------------------------------
    // Adjust for bright band (liquid skin on frozen particles)
    // ----------------------------------------------------------------
    amrex::Real factorb_s = factor_s;
    amrex::Real factorb_g = factor_g;
    if (iliqskin == 1 && tmk > celkel) {
        factorb_s = factor_s / alpha;
        factorb_g = factor_g / alpha;
    }

    // ----------------------------------------------------------------
    // Variable intercepts
    // ----------------------------------------------------------------
    amrex::Real sonv = rn0_s;
    if (in0s == 1) { // Thompson for snow
        amrex::Real temp_c = std::min(-Real(0.001), tmk - celkel);
        sonv = std::min(Real(2.0e8), Real(2.0e6) * std::exp(-Real(0.12) * temp_c));
    }

    amrex::Real gonv = rn0_g;
    if (in0g == 1) { // Thompson for graupel
        if (qgr > r1) {
            gonv = Real(2.38) * std::pow(pi * rho_g / (rho_air * qgr), Real(0.92));
            gonv = std::max(Real(1.e4), std::min(gonv, gon));
        }
    }

    amrex::Real ronv = rn0_r;
    if (in0r == 1) { // Thompson for rain
        ronv = ron2;
        if (qra > r1) {
            ronv = ron_const1r * std::tanh((ron_qr0 - qra) / ron_delqr0) + ron_const2r;
        }
    }

    // ----------------------------------------------------------------
    // Equivalent reflectivity factor Z_e [mm^6 m^-3]
    // ----------------------------------------------------------------
    amrex::Real z_e = factor_r  * std::pow(rho_air * qra, Real(1.75)) / std::pow(ronv, Real(0.75))
             + factorb_s * std::pow(rho_air * qsn, Real(1.75)) / std::pow(sonv, Real(0.75))
             + factorb_g * std::pow(rho_air * qgr, Real(1.75)) / std::pow(gonv, Real(0.75));

    // Clamp small values
    z_e = std::max(z_e, Real(0.01));

    // ----------------------------------------------------------------
    // Convert to dBZ
    // ----------------------------------------------------------------
    return Real(10.0) * std::log10(z_e);
}
#endif
