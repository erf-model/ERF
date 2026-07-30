#ifndef ERF_RADIATION_CLOUD_FRACTION_H
#define ERF_RADIATION_CLOUD_FRACTION_H

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

struct RadiationCloudFractions
{
    amrex::Real liquid;
    amrex::Real ice;
    amrex::Real total;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RadiationCloudFractions
radiation_cloud_fractions (const amrex::Real qc,
                           const amrex::Real qi,
                           const amrex::Real supplied_liquid_fraction,
                           const bool use_supplied_liquid_fraction) noexcept
{
    constexpr amrex::Real min_cloud_fraction = amrex::Real(1.0e-4);

    amrex::Real liquid_cf = amrex::Real(0.0);
    if (qc > amrex::Real(0.0)) {
        if (use_supplied_liquid_fraction) {
            liquid_cf = supplied_liquid_fraction;
            liquid_cf = (liquid_cf < amrex::Real(0.0)) ? amrex::Real(0.0) : liquid_cf;
            liquid_cf = (liquid_cf > amrex::Real(1.0)) ? amrex::Real(1.0) : liquid_cf;
            liquid_cf = (liquid_cf < min_cloud_fraction) ? min_cloud_fraction : liquid_cf;
        } else {
            liquid_cf = amrex::Real(1.0);
        }
    }

    const amrex::Real ice_cf = (qi > amrex::Real(0.0)) ? amrex::Real(1.0) : amrex::Real(0.0);
    return {liquid_cf, ice_cf,
            (liquid_cf > ice_cf) ? liquid_cf : ice_cf};
}

#endif
