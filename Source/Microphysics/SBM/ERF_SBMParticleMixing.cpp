#include "ERF_SBMParticleMixing.H"
#include "ERF_SBMCanonicalIdentity.H"
#include "ERF_IndexDefines.H"

#include <cmath>

namespace erf_sbm {

const char*
particle_mixing_name (const ParticleMixingPolicy policy) noexcept
{
    switch (policy) {
    case ParticleMixingPolicy::None: return "None";
    case ParticleMixingPolicy::ConstantReference: return "ConstantReference";
    case ParticleMixingPolicy::HostPassiveScalar: return "HostPassiveScalar";
    }
    return "Unknown";
}

bool
parse_particle_mixing_policy (const std::string& name, ParticleMixingPolicy& policy)
{
    for (const auto candidate : {ParticleMixingPolicy::None,
                                 ParticleMixingPolicy::ConstantReference,
                                 ParticleMixingPolicy::HostPassiveScalar}) {
        if (name == particle_mixing_name(candidate)) {
            policy = candidate;
            return true;
        }
    }
    return false;
}

std::string
ParticleMixingConfig::identity () const
{
    switch (policy) {
    case ParticleMixingPolicy::None:
        return "particle-mixing=none-v1";
    case ParticleMixingPolicy::ConstantReference:
        return "particle-mixing=constant-reference-v1[kh=" + canonical_real(reference_kh) +
               ",kv=" + canonical_real(reference_kv) + "]";
    case ParticleMixingPolicy::HostPassiveScalar:
        return "particle-mixing=host-passive-scalar-v1[Scalar_h,Scalar_v]";
    }
    return "particle-mixing=unknown";
}

bool
validate_particle_mixing_config (const ParticleMixingConfig& config, std::string& reason)
{
    reason.clear();
    if (config.policy == ParticleMixingPolicy::ConstantReference) {
        for (const amrex::Real value : {config.reference_kh, config.reference_kv}) {
            if (!(std::isfinite(value) && value >= amrex::Real(0.0))) {
                reason = "ConstantReference diffusivities must be finite and nonnegative";
                return false;
            }
        }
        return true;
    }
    if (config.reference_kh != amrex::Real(0.0) || config.reference_kv != amrex::Real(0.0)) {
        reason = std::string("reference diffusivities apply only to ConstantReference, not ") +
                 particle_mixing_name(config.policy);
        return false;
    }
    return true;
}

bool
particle_mixing_coefficients (const ParticleMixingConfig& config,
                              erf_auxiliary::MappedDiffusionCoefficients& coefficients,
                              std::string& reason)
{
    if (!validate_particle_mixing_config(config, reason)) { return false; }
    coefficients = {};
    switch (config.policy) {
    case ParticleMixingPolicy::None:
        reason = "particle mixing policy None has no diffusion coefficients";
        return false;
    case ParticleMixingPolicy::ConstantReference:
        coefficients.horizontal.molecular_coeff = config.reference_kh;
        coefficients.vertical.molecular_coeff = config.reference_kv;
        coefficients.mode.molecular_is_kinematic = true;
        coefficients.mode.add_turbulence = false;
        return true;
    case ParticleMixingPolicy::HostPassiveScalar:
        for (auto* channel : {&coefficients.horizontal, &coefficients.vertical}) {
            channel->molecular_coeff = amrex::Real(0.0);
            channel->eddy_h_comp = EddyDiff::Scalar_h;
            channel->eddy_v_comp = EddyDiff::Scalar_v;
        }
        coefficients.mode.molecular_is_kinematic = false;
        coefficients.mode.add_turbulence = true;
        return true;
    }
    reason = "unknown particle mixing policy";
    return false;
}

bool
particle_mixing_mesh (const MeshType mesh_type,
                      erf_auxiliary::DiffusionMesh& mesh,
                      std::string& reason)
{
    reason.clear();
    if (mesh_type == MeshType::ConstantDz) {
        mesh = erf_auxiliary::DiffusionMesh::ConstantDz;
        return true;
    }
    if (mesh_type == MeshType::StretchedDz) {
        mesh = erf_auxiliary::DiffusionMesh::StretchedDz;
        return true;
    }
    reason = "SBM particle mixing has no qualified monotone low-order operator on a "
             "terrain-fitted (VariableDz) mesh: its centred cross-diffusion stencil "
             "is not positive, and the terrain design gate is open";
    return false;
}

bool
build_particle_mixing_transmissibility (
    const ParticleMixingConfig& config,
    const amrex::Geometry& geom,
    const erf_auxiliary::MappedDiffusionMetrics& metrics,
    const amrex::MultiFab& density,
    const int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    erf_auxiliary::MappedFaceTransmissibility& transmissibility,
    std::string& diagnostic)
{
    if (!validate_particle_mixing_config(config, diagnostic)) { return false; }
    if (config.policy == ParticleMixingPolicy::None) {
        if (!transmissibility.is_defined() || transmissibility.nComp() != 1) {
            diagnostic = "particle mixing transmissibility requires one face component";
            return false;
        }
        transmissibility.setVal(amrex::Real(0.0));
        return true;
    }
    erf_auxiliary::MappedDiffusionCoefficients coefficients;
    if (!particle_mixing_coefficients(config, coefficients, diagnostic)) { return false; }
    if (!erf_auxiliary::BuildMappedDiffusionTransmissibility(
            transmissibility, geom, metrics, coefficients, density, density_component,
            eddy_diffusivity, diagnostic)) {
        diagnostic = std::string("SBM particle mixing ") +
                     particle_mixing_name(config.policy) + ": " + diagnostic;
        return false;
    }
    return true;
}

} // namespace erf_sbm
