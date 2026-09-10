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

enum class RadiationCloudSamplingMode
{
    MCICA,
    DeterministicBinary
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
RadiationCloudSamplingMode
radiation_cloud_sampling_mode (const bool do_subcol_sampling) noexcept
{
    return do_subcol_sampling ? RadiationCloudSamplingMode::MCICA
                              : RadiationCloudSamplingMode::DeterministicBinary;
}

inline bool
radiation_sampling_config_unsupported (const bool do_subcol_sampling,
                                       const bool native_shoc,
                                       const bool use_shoc_cldfrac) noexcept
{
    return !do_subcol_sampling && native_shoc && use_shoc_cldfrac;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real
radiation_cloud_mass (const amrex::Real mixing_ratio,
                      const amrex::Real cloud_fraction,
                      const amrex::Real rho,
                      const amrex::Real dz) noexcept
{
    constexpr amrex::Real min_cloud_fraction = amrex::Real(1.0e-4);
    constexpr amrex::Real max_incloud_mixing_ratio = amrex::Real(0.005);

    if (cloud_fraction <= amrex::Real(0.0)) {
        return amrex::Real(0.0);
    }

    const amrex::Real effective_fraction =
        (cloud_fraction < min_cloud_fraction) ? min_cloud_fraction : cloud_fraction;
    const amrex::Real in_cloud_mixing_ratio = mixing_ratio / effective_fraction;
    const amrex::Real capped_mixing_ratio =
        (in_cloud_mixing_ratio < max_incloud_mixing_ratio)
            ? in_cloud_mixing_ratio : max_incloud_mixing_ratio;
    return capped_mixing_ratio * rho * dz;
}

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
