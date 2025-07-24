#ifndef ERF_MOSTRoughness_H
#define ERF_MOSTRoughness_H

#include <ERF_Constants.H>


/**
 * Modified to include water depth effects
 * See Jimenez and Dudhia (2018) https://doi.org/10.1175/JAMC-D-17-0137.1
 */
struct modified_Charnock_roughness {
    modified_Charnock_roughness(amrex::Real /*zref*/, amrex::Real depth)
    : Cnk_b(Cnk_b1 * std::log(Cnk_b2 / depth))
    {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real /*umm*/,
                amrex::Real ustar)
    {
        z0 = std::exp( (2.7*ustar - 1.8/Cnk_b) / (ustar + 0.17/Cnk_b) );
    }

    const amrex::Real Cnk_b;

    static constexpr amrex::Real Cnk_b1{1.0/30.0};  ///< Modified Charnock constants, Eq (4)
    static constexpr amrex::Real Cnk_b2{1260.0};    ///< Modified Charnock constants, Eq (4)
};

/**
 * COARE3.0, following WRF phys/module_sf_mynn.F
 */
struct COARE3_roughness {
    COARE3_roughness(amrex::Real _zref, amrex::Real /*depth*/) : zref(_zref) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real umm,
                amrex::Real ustar)
    {
        amrex::Real ws10m = umm * std::log(10./z0ref) / std::log(zref/z0ref);
        amrex::Real Czc = 0.011 + 0.007 * std::clamp((ws10m-10.)/8., 0., 1.0);
        z0 = (Czc / CONST_GRAV) * ustar * ustar;
        // From Davis et al 2008 MWR, but their lower limit is 0.125e-6 -- WRF typo?
        z0 = std::clamp(z0, 1.27e-7, 2.85e-3);
    }

    const amrex::Real zref;

    static constexpr amrex::Real z0ref{1e-4};
};

/**
 * Updated Donelan 2004 as found in AHW
 *
 * Matches implementations in WRF phys/physics_mmm/sf_sfclayrev.F90 and
 * phys/module_sf_mynn.F
 */
struct Donelan_roughness {
    Donelan_roughness(amrex::Real /*zref*/, amrex::Real /*depth*/) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real /*umm*/,
                amrex::Real ustar)
    {
        amrex::Real zw  = std::min(std::pow(ustar/1.06, 0.3), 1.0);
        amrex::Real zn1 = 0.011 * ustar * ustar / CONST_GRAV + ozo;
        amrex::Real zn2 = 10. * std::exp(-9.5 * std::pow(ustar, -.3333))
                        + 0.11 * 1.5e-5 / std::max(ustar, 0.01);
        z0 = (1.0-zw) * zn1 + zw * zn2;
        // From Davis et al 2008 MWR, but their lower limit is 0.125e-6 -- WRF typo?
        z0 = std::clamp(z0, 1.27e-7, 2.85e-3);
    }

    static constexpr amrex::Real ozo{1.59e-5};
};
#endif
