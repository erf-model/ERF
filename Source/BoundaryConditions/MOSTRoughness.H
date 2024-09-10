#include <ERF_Constants.H>

/**
 * COARE3.0, following WRF phys/module_sf_mynn.F
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
COARE3_roughness (amrex::Real zref,
                  amrex::Real umm,
                  amrex::Real ustar)
{
    amrex::Real ws10m = umm*std::log(10./1e-4) / std::log(zref/1e-4);
    amrex::Real Czc = 0.011 + 0.007*std::min(std::max((ws10m-10.)/8., 0.), 1.0);
    amrex::Real z0 = (Czc / CONST_GRAV) * ustar * ustar;
    // From Davis et al 2008 MWR, but their lower limit is 0.125e-6 -- WRF typo?
    return std::min(std::max(z0, 1.27e-7), 2.85e-3);
}

/**
 * Updated Donelan 2004 as found in AHW
 *
 * Matches implementations in WRF phys/physics_mmm/sf_sfclayrev.F90 and
 * phys/module_sf_mynn.F
 */
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Donelan_roughness (amrex::Real ustar)
{
    constexpr amrex::Real ozo{1.59e-5};
    amrex::Real zw  = std::min(std::pow(ustar/1.06, 0.3),1.0);
    amrex::Real zn1 = 0.011*ustar*ustar/CONST_GRAV + ozo;
    amrex::Real zn2 = 10.*std::exp(-9.5*std::pow(ustar, -.3333))
                    + 0.11*1.5e-5/std::max(ustar,0.01);
    amrex::Real z0 = (1.0-zw) * zn1 + zw * zn2;
    // From Davis et al 2008 MWR, but their lower limit is 0.125e-6 -- WRF typo?
    return std::min(std::max(z0, 1.27e-7), 2.85e-3);
}
