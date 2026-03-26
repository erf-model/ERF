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
    amrex::Real ws10m = umm*std::log(amrex::Real(10.)/1e-4) / std::log(zref/1e-4);
    amrex::Real Czc = amrex::Real(0.011) + amrex::Real(0.007)*std::min(std::max((ws10m-amrex::Real(10.))/amrex::Real(8.), amrex::Real(0.0)), amrex::Real(1.0));
    amrex::Real z0 = (Czc / CONST_GRAV) * ustar * ustar;
    // From Davis et al 2008 MWR, but their lower limit is amrex::Real(0.125e-6) -- WRF typo?
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
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
    constexpr amrex::Real ozo{amrex::Real(1.59e-5)};
    amrex::Real zw  = std::min(std::pow(ustar/amrex::Real(1.06), amrex::Real(0.3)),amrex::Real(1.0));
    amrex::Real zn1 = amrex::Real(0.011)*ustar*ustar/CONST_GRAV + ozo;
    amrex::Real zn2 = amrex::Real(10.)*std::exp(-amrex::Real(9.5)*std::pow(ustar, -amrex::Real(.3333)))
                    + amrex::Real(0.11)*amrex::Real(1.5e-5)/std::max(ustar,amrex::Real(0.01));
    amrex::Real z0 = (amrex::Real(1.0)-zw) * zn1 + zw * zn2;
    // From Davis et al 2008 MWR, but their lower limit is amrex::Real(0.125e-6) -- WRF typo?
    return std::min(std::max(z0, amrex::Real(1.27e-7)), amrex::Real(2.85e-3));
}
