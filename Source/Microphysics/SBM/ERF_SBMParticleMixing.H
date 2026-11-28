#ifndef ERF_SBM_PARTICLE_MIXING_H_
#define ERF_SBM_PARTICLE_MIXING_H_

#include "ERF_AuxiliaryDiffusion.H"
#include "ERF_DataStruct.H"
#include "ERF_SBMPhysicalBoundary.H"

#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace erf_sbm {

/**
 * Turbulent gradient mixing of spectral particles. One coefficient acts on
 * every member of every bin group, and particle molecular (Brownian)
 * diffusion is zero under every policy.
 *
 * HostPassiveScalar consumes ERF's Scalar_h/Scalar_v eddy channels, which are
 * already density weighted and Schmidt scaled. Under ERF's LES closures they
 * equal the Q channels that diffuse moisture, but the PBL schemes other than
 * SHOC fill Q_v and not Scalar_v, so this policy is for LES closures only and
 * enabling it must reject a PBL scheme. ConstantReference multiplies
 * kinematic reference diffusivities by the face-averaged dry-air density.
 */
enum class ParticleMixingPolicy { None, ConstantReference, HostPassiveScalar };

struct ParticleMixingConfig
{
    ParticleMixingPolicy policy{ParticleMixingPolicy::None};
    // Kinematic horizontal and vertical diffusivities [m^2/s], ConstantReference only.
    amrex::Real reference_kh{0};
    amrex::Real reference_kv{0};

    /** Exact policy and coefficient-source line for a restart identity. */
    [[nodiscard]] std::string identity () const;
};

[[nodiscard]] const char* particle_mixing_name (ParticleMixingPolicy policy) noexcept;

[[nodiscard]] bool parse_particle_mixing_policy (const std::string& name,
                                                 ParticleMixingPolicy& policy);

[[nodiscard]] bool validate_particle_mixing_config (const ParticleMixingConfig& config,
                                                    std::string& reason);

/** The policy named by erf.sbm_particle_mixing and its reference diffusivities. */
[[nodiscard]] bool make_particle_mixing_config (const std::string& name,
                                                amrex::Real reference_kh,
                                                amrex::Real reference_kv,
                                                ParticleMixingConfig& config,
                                                std::string& reason);

/**
 * The ERF configuration a policy is qualified with: no molecular, PBL or RANS
 * diffusion; ERF turbulence only as the same Smagorinsky or Deardorff LES on
 * every level, and only with HostPassiveScalar, which requires it; and a
 * constant- or stretched-dz mesh. SBM's own AMR envelope (two levels, a
 * constant-dz mesh) bounds the levels.
 */
[[nodiscard]] bool validate_particle_mixing_host (const ParticleMixingConfig& config,
                                                  const std::vector<TurbChoice>& turbulence,
                                                  MolecDiffType molecular,
                                                  MeshType mesh,
                                                  std::string& reason);

/**
 * As above. With `terrain_qualification`, for qualification only, a policy
 * other than None also accepts a terrain-fitted (VariableDz) mesh on one
 * level; the stage still checks the geometry admission below.
 */
[[nodiscard]] bool validate_particle_mixing_host (const ParticleMixingConfig& config,
                                                  const std::vector<TurbChoice>& turbulence,
                                                  MolecDiffType molecular,
                                                  MeshType mesh,
                                                  bool terrain_qualification,
                                                  std::string& reason);

/**
 * Mixing accepts periodic, wall, outflow and inflow faces. An inflow_outflow
 * z face is rejected: ERF's own vertical scalar diffusion has no inward-only
 * rule there to follow, and ERF accepts inflow_outflow on x and y faces only.
 */
[[nodiscard]] bool validate_particle_mixing_boundaries (
    const ParticleMixingConfig& config,
    const SBMPhysicalBoundaryConfig& boundaries,
    std::string& reason);

/** Restart-schema lines of a mixing policy; empty for None. */
[[nodiscard]] std::string particle_mixing_schema (const ParticleMixingConfig& config);

/**
 * As above, followed with a line of the diffusive semantics of every face
 * that is neither periodic nor a wall, when there is one.
 */
[[nodiscard]] std::string particle_mixing_schema (const ParticleMixingConfig& config,
                                                  const SBMPhysicalBoundaryConfig& boundaries);

/** The host coefficient semantics of a mixing policy other than None. */
[[nodiscard]] bool particle_mixing_coefficients (
    const ParticleMixingConfig& config,
    erf_auxiliary::MappedDiffusionCoefficients& coefficients,
    std::string& reason);

/**
 * The low-order operator mesh; terrain-fitted (VariableDz) meshes are rejected
 * unless `terrain_qualification` is set, for qualification only.
 */
[[nodiscard]] bool particle_mixing_mesh (MeshType mesh_type,
                                         erf_auxiliary::DiffusionMesh& mesh,
                                         std::string& reason,
                                         bool terrain_qualification = false);

/**
 * The terrain geometry terrain mixing is qualified on: in every valid cell,
 * the drop sqrt(dx_drop^2 + dy_drop^2) of ComputeTerrainCellDrops, unclamped,
 * at most `max_drop` cell thicknesses, and a finite, positive thickness.
 * Fails with the worst cell and its drop otherwise. An engineering envelope,
 * not a stability bound.
 */
[[nodiscard]] bool admit_terrain_mixing_geometry (const amrex::MultiFab& z_nd,
                                                  const amrex::BoxArray& cells,
                                                  amrex::Real max_drop,
                                                  std::string& reason);

/**
 * The face transmissibility of a policy; None gives zero on every face.
 * `density` must be the dry-air density (ERF's Rho_comp), which
 * ConstantReference multiplies by.
 */
[[nodiscard]] bool build_particle_mixing_transmissibility (
    const ParticleMixingConfig& config,
    const amrex::Geometry& geom,
    const erf_auxiliary::MappedDiffusionMetrics& metrics,
    const amrex::MultiFab& density,
    int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    erf_auxiliary::MappedFaceTransmissibility& transmissibility,
    std::string& diagnostic);

/** As above, with the diffusive faces of `boundary` open. */
[[nodiscard]] bool build_particle_mixing_transmissibility (
    const ParticleMixingConfig& config,
    const amrex::Geometry& geom,
    const erf_auxiliary::MappedDiffusionMetrics& metrics,
    const amrex::MultiFab& density,
    int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    const erf_auxiliary::MappedDiffusionBoundary& boundary,
    erf_auxiliary::MappedFaceTransmissibility& transmissibility,
    std::string& diagnostic);

/**
 * Terrain-fitted mixing only: the fraction theta in [0, 1] of one matched
 * coordinate's nominal low-order diffusive rate `low` that a face carries,
 * given that coordinate's native rate `high`: 1 for low = 0, 0 where high is
 * zero or opposes low, |high| / |low| where it is smaller, and 1 otherwise,
 * so theta low lies between zero and high. A nonfinite rate gives NaN, which
 * the caller rejects. A face's group takes the least fraction of its
 * coordinates.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
matched_diffusion_fraction (const amrex::Real high, const amrex::Real low) noexcept
{
    if (!amrex::Math::isfinite(high) || !amrex::Math::isfinite(low)) {
        return std::numeric_limits<amrex::Real>::quiet_NaN();
    }
    if (low == amrex::Real(0.0)) { return amrex::Real(1.0); }
    if (high == amrex::Real(0.0) || (high > amrex::Real(0.0)) != (low > amrex::Real(0.0))) {
        return amrex::Real(0.0);
    }
    const amrex::Real h = amrex::Math::abs(high);
    const amrex::Real l = amrex::Math::abs(low);
    return h < l ? h / l : amrex::Real(1.0);
}

/**
 * As particle_mixing_schema with boundaries, for the mixing mesh: a
 * terrain-fitted (VariableDz) mesh adds a line naming its operator and
 * boundary policy; every other mesh leaves the schema unchanged.
 */
[[nodiscard]] std::string particle_mixing_schema (const ParticleMixingConfig& config,
                                                  const SBMPhysicalBoundaryConfig& boundaries,
                                                  erf_auxiliary::DiffusionMesh mesh);

/** The diffusive faces of a boundary configuration, with the dry-air carrier of each direction. */
[[nodiscard]] erf_auxiliary::MappedDiffusionBoundary particle_mixing_boundary (
    const SBMPhysicalBoundaryConfig& boundaries,
    const std::array<const amrex::MultiFab*, AMREX_SPACEDIM>& carrier);

} // namespace erf_sbm

#endif
