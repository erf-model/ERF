#ifndef ERF_AUXILIARY_DIFFUSION_H_
#define ERF_AUXILIARY_DIFFUSION_H_

#include "ERF_AuxiliaryDiffusionFace.H"
#include "ERF_AuxiliaryMappedTransport.H"
#include "ERF_ScalarDiffusion.H"

#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>

#include <array>
#include <string>
#include <vector>

namespace erf_auxiliary {

struct MappedFaceTransmissibilityTag {};

/**
 * Nonnegative face exchange coefficients T_f in ERF's canonical mapped face
 * convention: the low-order transfer of an intensive field z through face f
 * is -T_f (z_R - z_L), and its divergence is applied with dx_inv / omega.
 */
using MappedFaceTransmissibility = MappedFaceField<MappedFaceTransmissibilityTag>;

/**
 * Host meshes. ConstantDz and StretchedDz are orthogonal, with a two-point
 * host diffusion stencil. VariableDz is static terrain-fitted: its low-order
 * exchange is the diagonal of the mapped diffusion tensor, which the host's
 * centred cross-diffusion stencil is not.
 */
enum class DiffusionMesh { ConstantDz, StretchedDz, VariableDz };

/**
 * Host metric views. Map factors are ERF's k=0 fields: u_x/u_y on x faces,
 * v_x/v_y on y faces, m_x/m_y at cell centers. ax/ay are read for StretchedDz
 * and VariableDz, the device pointer to the stretched vertical spacing only
 * for StretchedDz, and z_nd, ERF's nodal physical heights on the cell layout
 * with at least one ghost node in z, only for VariableDz.
 */
struct MappedDiffusionMetrics
{
    DiffusionMesh mesh{DiffusionMesh::ConstantDz};
    const amrex::MultiFab* mf_ux{nullptr};
    const amrex::MultiFab* mf_uy{nullptr};
    const amrex::MultiFab* mf_vx{nullptr};
    const amrex::MultiFab* mf_vy{nullptr};
    const amrex::MultiFab* mf_mx{nullptr};
    const amrex::MultiFab* mf_my{nullptr};
    const amrex::MultiFab* ax{nullptr};
    const amrex::MultiFab* ay{nullptr};
    const amrex::Real* stretched_dz{nullptr};
    const amrex::MultiFab* z_nd{nullptr};
};

/**
 * Face coefficients under the host's rule: an eddy coefficient is already
 * density weighted, and a kinematic molecular coefficient is multiplied by
 * the face-averaged density. Horizontal faces read `horizontal.eddy_h_comp`
 * and vertical faces `vertical.eddy_v_comp`.
 */
struct MappedDiffusionCoefficients
{
    ScalarDiffusionCoefficients horizontal;
    ScalarDiffusionCoefficients vertical;
    ScalarDiffusionCoefficientMode mode;
};

/**
 * Build T_f equal to the host N/S scalar-diffusion stencil at every interior
 * and periodic face, and zero at every nonperiodic domain face. Fails, with a
 * diagnostic, on any nonfinite or negative coefficient, nonpositive density
 * read by a kinematic coefficient, or nonfinite or negative result.
 *
 * On VariableDz, T_f is the diagonal of the mapped diffusion tensor with the
 * host's face averages: side faces take a (m_num / m_den) mu_h / dxi, and a
 * k face (mu_v + mu_h [(m_x h_xi)^2 + (m_y h_eta)^2]) / (m_x m_y J dzeta),
 * where J, h_xi and h_eta are ERF's k-face metrics and mu_h on a k face
 * averages the horizontal channel of its two cells. The slope term keeps a
 * slope's horizontal diffusion in the vertical exchange, which is nonzero
 * with mu_v = 0. Fails, with a diagnostic, on a periodic or open z face, or
 * on a nonfinite or nonpositive side area, k-face Jacobian or map factor, or
 * a nonfinite slope.
 */
bool BuildMappedDiffusionTransmissibility (
    MappedFaceTransmissibility& transmissibility,
    const amrex::Geometry& geom,
    const MappedDiffusionMetrics& metrics,
    const MappedDiffusionCoefficients& coefficients,
    const amrex::MultiFab& density,
    int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    std::string& diagnostic);

/**
 * Diffusive kind of each domain face, indexed by amrex::Orientation (Closed
 * by default; periodic directions ignore it), and the dry-air carrier of each
 * direction on the transmissibility's face layout, which a
 * DirichletWhereInward face reads for its sign.
 */
struct MappedDiffusionBoundary
{
    std::array<DiffusionFaceKind, 2 * AMREX_SPACEDIM> kind{};
    std::array<const amrex::MultiFab*, AMREX_SPACEDIM> carrier{};
};

/**
 * As above, with open domain faces. Where a face is open, T_f is the host
 * face coefficient, read from the ghost and the edge cell, over half the edge
 * cell: the distance from its center to the face, which holds the ghost
 * value. Where it is closed, T_f is zero. The validation also reads the ghost
 * cells beyond open faces.
 */
bool BuildMappedDiffusionTransmissibility (
    MappedFaceTransmissibility& transmissibility,
    const amrex::Geometry& geom,
    const MappedDiffusionMetrics& metrics,
    const MappedDiffusionCoefficients& coefficients,
    const amrex::MultiFab& density,
    int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    const MappedDiffusionBoundary& boundary,
    std::string& diagnostic);

/**
 * Low-order diffusive face rates -T_f (z_R - z_L) of the selected intensive
 * components; rate slot n holds component components[n]. A domain face with
 * T_f = 0 carries nothing and reads no ghost; an open one exchanges with the
 * ghost value. Fails, with a diagnostic, if any rate is nonfinite.
 */
bool BuildLowOrderDiffusionRate (
    MappedFaceFluxRate& rate,
    const MappedFaceTransmissibility& transmissibility,
    const amrex::Geometry& geom,
    const amrex::MultiFab& intensive,
    const std::vector<int>& components,
    std::string& diagnostic);

/**
 * Outgoing donor rate of one cell for low-order advection and diffusion from
 * the same anchor: the advective outflow plus every face exchange, over the
 * cell's carrier inventory. Separate advective and diffusive bounds are not
 * sufficient; only this sum bounds the shared positive update.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
MappedCombinedOutgoingDemandRate (const amrex::Real fx_hi,
                                  const amrex::Real fx_lo,
                                  const amrex::Real fy_hi,
                                  const amrex::Real fy_lo,
                                  const amrex::Real fz_hi,
                                  const amrex::Real fz_lo,
                                  const amrex::Real tx_hi,
                                  const amrex::Real tx_lo,
                                  const amrex::Real ty_hi,
                                  const amrex::Real ty_lo,
                                  const amrex::Real tz_hi,
                                  const amrex::Real tz_lo,
                                  const amrex::Real measure,
                                  const amrex::Real density,
                                  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv) noexcept
{
    const amrex::Real exchange = (tx_hi + tx_lo) * dx_inv[0] +
                                 (ty_hi + ty_lo) * dx_inv[1] +
                                 (tz_hi + tz_lo) * dx_inv[2];
    return MappedOutgoingDemandRate(fx_hi, fx_lo, fy_hi, fy_lo, fz_hi, fz_lo,
                                    measure, density, dx_inv) +
           exchange / (measure * density);
}

/**
 * Reduce the largest combined outgoing rate over a cell-centered domain. On
 * VariableDz, T's slope term carries the stiffness of the host stencil's
 * cross diffusion, which a slope-blind rate misses; pass the full T, never one
 * reduced by a limiter.
 */
bool ComputeMaxMappedCombinedOutgoingRate (
    const MappedFaceFluxRate& carrier,
    const MappedFaceTransmissibility& transmissibility,
    const amrex::MultiFab& measure,
    const amrex::MultiFab& conserved,
    int density_component,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
    amrex::Real& max_rate,
    std::string& diagnostic);

} // namespace erf_auxiliary

#endif
