#ifndef ERF_NEAR_SURFACE_DIAGNOSTICS_H_
#define ERF_NEAR_SURFACE_DIAGNOSTICS_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <cmath>

#include "ERF_Constants.H"
#include "ERF_MOSTStress.H"
#include "ERF_SurfaceDiagnosticSource.H"

namespace near_surface_diagnostics
{

struct NativeBundle
{
    amrex::Real temperature_vegetated;
    amrex::Real temperature_bare;
    amrex::Real mixing_ratio_vegetated;
    amrex::Real mixing_ratio_bare;
    amrex::Real vegetation_fraction;
};

struct NativeAggregate
{
    amrex::Real temperature;
    amrex::Real mixing_ratio;
};

struct MostState
{
    amrex::Real theta_surface;
    amrex::Real theta_star;
    amrex::Real mixing_ratio_surface;
    amrex::Real mixing_ratio_star;
    amrex::Real roughness_height;
    amrex::Real obukhov_length;
};

struct MostProfile
{
    amrex::Real potential_temperature;
    amrex::Real mixing_ratio;
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_real (amrex::Real value) noexcept
{
    return std::isfinite(value);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool valid_noahmp_component (amrex::Real value) noexcept
{
    return valid_real(value) && value > amrex::Real(-9990.0) &&
           value < amrex::Real(0.5) * lsm_undefined;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_temperature_is_valid (const NativeBundle& bundle) noexcept;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_mixing_ratio_is_valid (const NativeBundle& bundle) noexcept;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_vegetation_fraction_is_valid (const NativeBundle& bundle) noexcept;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_is_valid (const NativeBundle& bundle) noexcept
{
    return native_bundle_temperature_is_valid(bundle) &&
           native_bundle_mixing_ratio_is_valid(bundle) &&
           native_bundle_vegetation_fraction_is_valid(bundle);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_temperature_is_valid (const NativeBundle& bundle) noexcept
{
    return valid_noahmp_component(bundle.temperature_vegetated) &&
           valid_noahmp_component(bundle.temperature_bare) &&
           bundle.temperature_vegetated > amrex::Real(0.0) &&
           bundle.temperature_bare > amrex::Real(0.0) &&
           valid_noahmp_component(bundle.vegetation_fraction) &&
           bundle.vegetation_fraction >= amrex::Real(0.0) &&
           bundle.vegetation_fraction <= amrex::Real(1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_mixing_ratio_is_valid (const NativeBundle& bundle) noexcept
{
    return valid_noahmp_component(bundle.mixing_ratio_vegetated) &&
           valid_noahmp_component(bundle.mixing_ratio_bare) &&
           valid_noahmp_component(bundle.vegetation_fraction) &&
           bundle.mixing_ratio_vegetated >= amrex::Real(0.0) &&
           bundle.mixing_ratio_bare >= amrex::Real(0.0) &&
           bundle.vegetation_fraction >= amrex::Real(0.0) &&
           bundle.vegetation_fraction <= amrex::Real(1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool native_bundle_vegetation_fraction_is_valid (const NativeBundle& bundle) noexcept
{
    return valid_noahmp_component(bundle.vegetation_fraction) &&
           bundle.vegetation_fraction >= amrex::Real(0.0) &&
           bundle.vegetation_fraction <= amrex::Real(1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
NativeAggregate aggregate_native_bundle (const NativeBundle& bundle) noexcept
{
    const amrex::Real f = bundle.vegetation_fraction;
    return {
        f * bundle.temperature_vegetated + (amrex::Real(1.0)-f) * bundle.temperature_bare,
        f * bundle.mixing_ratio_vegetated + (amrex::Real(1.0)-f) * bundle.mixing_ratio_bare
    };
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool evaluate_most_factor (const MostState& state,
                           amrex::Real height,
                           amrex::Real& factor) noexcept
{
    if (!valid_real(state.roughness_height) ||
        !valid_real(state.obukhov_length) ||
        state.roughness_height <= amrex::Real(0.0) ||
        height <= state.roughness_height ||
        state.obukhov_length == amrex::Real(0.0)) {
        return false;
    }
    const similarity_funs sfuns{};
    factor = std::log(height / state.roughness_height) -
             sfuns.calc_psi_h(height / state.obukhov_length);
    return valid_real(factor);
}

// Evaluate ERF's integrated scalar profile without changing the production
// surface-layer state. The inputs use theta and mixing-ratio scales.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool evaluate_most_profile (const MostState& state,
                            amrex::Real height,
                            MostProfile& profile) noexcept
{
    if (!valid_real(state.theta_surface) || !valid_real(state.theta_star) ||
        !valid_real(state.mixing_ratio_surface) || !valid_real(state.mixing_ratio_star)) {
        return false;
    }

    amrex::Real profile_factor = amrex::Real(0.0);
    if (!evaluate_most_factor(state, height, profile_factor)) {
        return false;
    }

    profile.potential_temperature = state.theta_surface +
        state.theta_star / KAPPA * profile_factor;
    profile.mixing_ratio = state.mixing_ratio_surface +
        state.mixing_ratio_star / KAPPA * profile_factor;
    return valid_real(profile.potential_temperature) &&
           valid_real(profile.mixing_ratio) &&
           profile.potential_temperature > amrex::Real(0.0) &&
           profile.mixing_ratio >= amrex::Real(0.0);
}

struct Sources
{
    const amrex::MultiFab* native_temperature_vegetated = nullptr;
    const amrex::MultiFab* native_temperature_bare = nullptr;
    const amrex::MultiFab* native_mixing_ratio_vegetated = nullptr;
    const amrex::MultiFab* native_mixing_ratio_bare = nullptr;
    const amrex::MultiFab* native_vegetation_fraction = nullptr;

    const amrex::MultiFab* theta_surface = nullptr;
    const amrex::MultiFab* theta_star = nullptr;
    const amrex::MultiFab* mixing_ratio_surface = nullptr;
    const amrex::MultiFab* mixing_ratio_star = nullptr;
    const amrex::MultiFab* roughness_height = nullptr;
    const amrex::MultiFab* obukhov_length = nullptr;
    const amrex::MultiFab* source_mask = nullptr;
    const amrex::iMultiFab* land_mask = nullptr;
    const amrex::MultiFab* cons = nullptr;
    const amrex::MultiFab* z_phys_nd = nullptr;

    amrex::Real dz = amrex::Real(0.0);
    int klo = 0;
    bool moist = false;
    bool has_lsm = false;
};

// Fill any requested unified 2-m components. A component index below zero
// means that the caller did not select that field. A source-only request uses
// the temperature path; provenance paired with humidity describes humidity
// without adding a temperature prerequisite. Dry runs retain a meaningful
// temperature source while humidity remains unavailable.
void fill (amrex::MultiFab& dst,
           int temperature_comp,
           int mixing_ratio_comp,
           int source_comp,
           const Sources& sources,
           amrex::Real missing_value = amrex::Real(-999.0));

} // namespace near_surface_diagnostics

#endif
