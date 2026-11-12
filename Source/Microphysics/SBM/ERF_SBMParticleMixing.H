#ifndef ERF_SBM_PARTICLE_MIXING_H_
#define ERF_SBM_PARTICLE_MIXING_H_

#include "ERF_AuxiliaryDiffusion.H"
#include "ERF_DataStruct.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <string>

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

/** The host coefficient semantics of a mixing policy other than None. */
[[nodiscard]] bool particle_mixing_coefficients (
    const ParticleMixingConfig& config,
    erf_auxiliary::MappedDiffusionCoefficients& coefficients,
    std::string& reason);

/** The low-order operator mesh; terrain-fitted (VariableDz) meshes are rejected. */
[[nodiscard]] bool particle_mixing_mesh (MeshType mesh_type,
                                         erf_auxiliary::DiffusionMesh& mesh,
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

} // namespace erf_sbm

#endif
