/**
 * \file ERF_Plotfile2DCloudOptics.H
 */
#ifndef ERF_PLOTFILE2DCLOUDOPTICS_H_
#define ERF_PLOTFILE2DCLOUDOPTICS_H_

#include <string>

#include <AMReX_Geometry.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include "ERF_DataStruct.H"
#include "ERF_NumericalConstants.H"

namespace plotfile2d {

// Illustrative mass extinction [m^2/kg], in qc, qi, qr, qs, qg order.
// These fixed bulk proxies are not scheme-calibrated optical properties.
struct CloudOpticsParameters
{
    amrex::GpuArray<amrex::Real, 5> kext{{amrex::Real(150.), amrex::Real(65.4),
                                          amrex::Real(3.), amrex::Real(30.),
                                          amrex::Real(7.5)}};
    // Two-stream-inspired default for conservative scattering (g = 0.85).
    // User-tunable; not calibrated to the active microphysics or radiation scheme.
    amrex::Real reflectance_scale = amrex::Real(7.7);
};

//! Value written where the optical diagnostics are missing or invalid
constexpr amrex::Real cloud_optics_missing = amrex::Real(-999.);

std::string
cloud_optics_parameter_error (const CloudOpticsParameters& parameters);
// Reads and validates erf.plot2d.cloud_optics.*, aborting on an invalid value.
CloudOpticsParameters parse_cloud_optics_parameters ();
// Prints the resolved settings as one "Plot2D cloud optics:" line.
void print_cloud_optics_parameters (const CloudOpticsParameters& parameters);
bool cloud_optics_available (const MoistureComponentIndices& indices) noexcept;

// Reference-illumination hydrometeor-only proxy, not a radiative-transfer solution.
// The two branches avoid overflow in tau + scale for large finite values.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cloud_reflectance (amrex::Real tau, amrex::Real scale) noexcept
{
    if (!amrex::Math::isfinite(tau) || tau < zero ||
        !amrex::Math::isfinite(scale) || scale <= zero) {
        return cloud_optics_missing;
    }
    if (tau == zero) {
        return zero;
    }
    const amrex::Real value = tau > scale ? one / (one + scale / tau)
                                          : (tau / scale) / (one + tau / scale);
    return amrex::min(amrex::Real(one), amrex::max(amrex::Real(zero), value));
}

// dst is the collapsed slab of cons with the same box ordering and mapping.
// Explicit destination components may be -1 for unrequested outputs. Computes
// one tau slab, sums z-stacked boxes, then transforms the completed level sum.
// This is the sum over this level's covered cells, even for an incomplete
// bottom-reaching column; no coarse-level mass is added. The writer retains
// its existing missing-value policy for columns with no ground source.
// Nonfinite present mass (even with kext=0), bad metrics, or overflow -> -999.
void fill_cloud_optics (amrex::MultiFab& dst,
                        int tau_component,
                        int albedo_component,
                        const amrex::MultiFab& cons,
                        const MoistureComponentIndices& indices,
                        const CloudOpticsParameters& parameters,
                        const amrex::Geometry& geom,
                        const amrex::MultiFab& detJ);

} // namespace plotfile2d
#endif
