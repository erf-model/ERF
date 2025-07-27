#ifndef ERF_MOSTRoughness_H
#define ERF_MOSTRoughness_H

#include <ERF_Constants.H>
#include <ERF_EOS.H>

using RoughnessModel = std::function<void(amrex::Real& /*z0*/,
                                          amrex::Real /*ustar*/,
                                          amrex::Real /*Umagmean*/,
                                          amrex::Real /*Tmean*/,
                                          amrex::Real /*depth*/,
                                          amrex::Real /*Kmv*/,
                                          amrex::Real /*Hwave*/,
                                          amrex::Real /*Lwave*/)>;

struct ConstantRoughness {

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    virtual void
    operator() (amrex::Real& /*z0*/,
                amrex::Real /*ustar*/,
                amrex::Real /*Umagmean*/,
                amrex::Real /*Tmean*/,
                amrex::Real /*depth*/,
                amrex::Real /*Kmv*/,
                amrex::Real /*Hwave*/,
                amrex::Real /*Lwave*/) const
    {
        // no change to z0
    }
};

/**
 * See Fairall et al. 2003
 * https://doi.org/10.1175/1520-0442(2003)016<0571:BPOASF>2.0.CO;2
 */
struct CharnockRoughness {

    CharnockRoughness(amrex::Real _Cnk_a, bool _visc)
    : Cnk_a(_Cnk_a)
    , visc(_visc)
    {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real ustar,
                amrex::Real /*Umagmean*/,
                amrex::Real Tmean,
                amrex::Real /*depth*/,
                amrex::Real /*Kmv*/,
                amrex::Real /*Hwave*/,
                amrex::Real /*Lwave*/) const
    {
        z0 = (Cnk_a / CONST_GRAV) * ustar * ustar;
        if (visc) {
            z0 += compute_air_viscosity(Tmean) / std::max(ustar, 0.05);
        }
    }

    const amrex::Real Cnk_a{0.0185}; // standard Charnock constant
    const bool visc{false};          // use viscous Charnock formulation
};

/**
 * Modified to include water depth effects
 * See Jimenez and Dudhia (2018) https://doi.org/10.1175/JAMC-D-17-0137.1
 */
struct ModifiedCharnockRoughness {

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real ustar,
                amrex::Real /*Umagmean*/,
                amrex::Real /*Tmean*/,
                amrex::Real depth,
                amrex::Real /*Kmv*/,
                amrex::Real /*Hwave*/,
                amrex::Real /*Lwave*/) const
    {
        const amrex::Real Cnk_b = Cnk_b1 * std::log(Cnk_b2 / depth);
        z0 = std::exp( (2.7*ustar - 1.8/Cnk_b) / (ustar + 0.17/Cnk_b) );
    }

    // Modified Charnock constant (Eqn. 4)
    static constexpr amrex::Real Cnk_b1{1.0/30.0};
    static constexpr amrex::Real Cnk_b2{1260.0};
};

/**
 * COARE3.0, following WRF phys/module_sf_mynn.F
 */
struct COARE3Roughness {

    COARE3Roughness(amrex::Real _zref) : zref(_zref) {}

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real ustar,
                amrex::Real Umagmean,
                amrex::Real /*Tmean*/,
                amrex::Real /*depth*/,
                amrex::Real /*Kmv*/,
                amrex::Real /*Hwave*/,
                amrex::Real /*Lwave*/) const
    {
        amrex::Real ws10m = Umagmean * std::log(10./z0ref) / std::log(zref/z0ref);
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
struct DonelanRoughness {

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real ustar,
                amrex::Real /*Umagmean*/,
                amrex::Real /*Tmean*/,
                amrex::Real /*depth*/,
                amrex::Real /*Kmv*/,
                amrex::Real /*Hwave*/,
                amrex::Real /*Lwave*/) const
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

/**
 * Wave Coupled
 */
struct WaveCoupledRoughness {

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void
    operator() (amrex::Real& z0,
                amrex::Real ustar,
                amrex::Real /*Umagmean*/,
                amrex::Real /*Tmean*/,
                amrex::Real /*depth*/,
                amrex::Real Kmv,
                amrex::Real Hwave,
                amrex::Real Lwave) const
    {
        z0 = 1200.0 * Hwave * std::pow( Hwave / (Lwave + eps), 4.5 )
           + 0.11 * Kmv / ustar;
        z0 = std::clamp(z0, z0_eps, z0_max);
    }

    static constexpr amrex::Real eps{1e-15};
    static constexpr amrex::Real z0_eps{1.0e-6};
    static constexpr amrex::Real z0_max{0.1};
};
#endif
