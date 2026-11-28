#include "ERF_SBMParticleMixing.H"
#include "ERF_SBMCanonicalIdentity.H"
#include "ERF_IndexDefines.H"
#include "ERF_TerrainDiffusionLimits.H"

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>
#include <sstream>

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
make_particle_mixing_config (const std::string& name,
                             const amrex::Real reference_kh,
                             const amrex::Real reference_kv,
                             ParticleMixingConfig& config,
                             std::string& reason)
{
    reason.clear();
    config = {};
    if (!parse_particle_mixing_policy(name, config.policy)) {
        reason = "erf.sbm_particle_mixing must be None, ConstantReference or "
                 "HostPassiveScalar, not '" + name + "'";
        return false;
    }
    config.reference_kh = reference_kh;
    config.reference_kv = reference_kv;
    if (!validate_particle_mixing_config(config, reason)) {
        reason = "erf.sbm_particle_mixing: " + reason;
        return false;
    }
    return true;
}

bool
validate_particle_mixing_host (const ParticleMixingConfig& config,
                               const std::vector<TurbChoice>& turbulence,
                               const MolecDiffType molecular,
                               const MeshType mesh,
                               std::string& reason)
{
    return validate_particle_mixing_host(config, turbulence, molecular, mesh, false, reason);
}

bool
validate_particle_mixing_host (const ParticleMixingConfig& config,
                               const std::vector<TurbChoice>& turbulence,
                               const MolecDiffType molecular,
                               const MeshType mesh,
                               const bool terrain_qualification,
                               std::string& reason)
{
    if (!validate_particle_mixing_config(config, reason)) { return false; }
    const std::string policy = particle_mixing_name(config.policy);
    if (molecular != MolecDiffType::None) {
        reason = "SBM does not support molecular diffusion (erf.molec_diff_type must be "
                 "None): particle molecular diffusion is zero, and molecular vapor "
                 "diffusion has not been qualified with SBM";
        return false;
    }
    for (std::size_t lev = 0; lev < turbulence.size(); ++lev) {
        const TurbChoice& turb = turbulence[lev];
        const std::string at_level = " (level " + std::to_string(lev) + ")";
        if (turb.pbl_type != PBLType::None) {
            reason = "SBM does not support a PBL closure (erf.pbl_type = " +
                     amrex::getEnumNameString(turb.pbl_type) + at_level +
                     "): ERF's PBL schemes other than SHOC set the moisture eddy "
                     "diffusivity but not the passive-scalar one particles use, and "
                     "several add nonlocal transport a local coefficient does not represent";
            return false;
        }
        if (turb.rans_type != RANSType::None) {
            reason = "SBM does not support a RANS closure (erf.rans_type)" + at_level;
            return false;
        }
        if (config.policy == ParticleMixingPolicy::HostPassiveScalar) {
            if (turb.les_type != LESType::Smagorinsky && turb.les_type != LESType::Deardorff) {
                reason = "SBM particle mixing HostPassiveScalar requires erf.les_type = "
                         "Smagorinsky or Deardorff on every level" + at_level;
                return false;
            }
            // ERF reads Smagorinsky2D as Smagorinsky with smag2d set
            if (turb.smag2d) {
                reason = "SBM particle mixing HostPassiveScalar does not support erf.les_type = "
                         "Smagorinsky2D" + at_level + ": its horizontal-only closure "
                         "leaves vertical mixing to a PBL scheme, which SBM does not support";
                return false;
            }
            if (turb.les_type != turbulence.front().les_type) {
                reason = "SBM particle mixing HostPassiveScalar requires the same "
                         "erf.les_type on every level";
                return false;
            }
        } else if (turb.les_type != LESType::None) {
            reason = "SBM particle mixing " + policy + " does not mix the spectrum with "
                     "ERF's LES closure (erf.les_type)" + at_level +
                     "; set erf.sbm_particle_mixing = HostPassiveScalar";
            return false;
        }
    }
    if (config.policy == ParticleMixingPolicy::None) { return true; }
    erf_auxiliary::DiffusionMesh diffusion_mesh = erf_auxiliary::DiffusionMesh::ConstantDz;
    if (!particle_mixing_mesh(mesh, diffusion_mesh, reason, terrain_qualification)) {
        return false;
    }
    if (diffusion_mesh == erf_auxiliary::DiffusionMesh::VariableDz && turbulence.size() != 1) {
        reason = "SBM particle mixing on a terrain-fitted (VariableDz) mesh is qualified on "
                 "one level only";
        return false;
    }
    return true;
}

bool
validate_particle_mixing_boundaries (const ParticleMixingConfig& config,
                                     const SBMPhysicalBoundaryConfig& boundaries,
                                     std::string& reason)
{
    reason.clear();
    for (int face = 0; face < physical_face_count; ++face) {
        const auto kind = boundaries.kind[static_cast<std::size_t>(face)];
        const auto diffusion = boundaries.diffusion[static_cast<std::size_t>(face)];
        if (diffusion != erf_auxiliary::DiffusionFaceKind::Closed &&
            kind != AdvectionFaceKind::PrescribedInflow) {
            reason = std::string("SBM particle mixing can exchange only through an inflow "
                                 "face, not the ") + face_kind_name(kind) + " face " +
                     face_name(face);
            return false;
        }
    }
    if (config.policy == ParticleMixingPolicy::None) { return true; }
    for (const int face : {2, 2 + AMREX_SPACEDIM}) {
        if (boundaries.diffusion[static_cast<std::size_t>(face)] ==
            erf_auxiliary::DiffusionFaceKind::DirichletWhereInward) {
            reason = std::string("SBM particle mixing ") + particle_mixing_name(config.policy) +
                     " does not support an inflow_outflow " + face_name(face) +
                     " face: ERF's vertical scalar diffusion has no inward-only rule to follow";
            return false;
        }
    }
    return true;
}

std::string
particle_mixing_schema (const ParticleMixingConfig& config)
{
    if (config.policy == ParticleMixingPolicy::None) { return {}; }
    return config.identity() + "\ndiffusion=two-point-low-host-high-group-fct-v1";
}

std::string
particle_mixing_schema (const ParticleMixingConfig& config,
                        const SBMPhysicalBoundaryConfig& boundaries)
{
    std::string schema = particle_mixing_schema(config);
    if (schema.empty() || !boundaries.has_inflow_or_outflow_face()) { return schema; }
    schema += "\nboundary-diffusion=open-faces-v1";
    for (const int face : {0, 3, 1, 4, 2, 5}) {
        const auto kind = boundaries.kind[static_cast<std::size_t>(face)];
        if (kind == AdvectionFaceKind::Periodic || kind == AdvectionFaceKind::ImpermeableWall) {
            continue;
        }
        schema += std::string(";") + face_name(face) + '=' +
                  diffusion_face_kind_name(boundaries.diffusion[static_cast<std::size_t>(face)]);
    }
    return schema;
}

std::string
particle_mixing_schema (const ParticleMixingConfig& config,
                        const SBMPhysicalBoundaryConfig& boundaries,
                        const erf_auxiliary::DiffusionMesh mesh)
{
    std::string schema = particle_mixing_schema(config, boundaries);
    if (schema.empty() || mesh != erf_auxiliary::DiffusionMesh::VariableDz) { return schema; }
    return schema + "\nterrain-diffusion=diagonal-low-matched-native-high-v1;"
                    "walls=complete-mapped-closure,quadratic-height-ghosts;"
                    "sides=face-trace";
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
                      std::string& reason,
                      const bool terrain_qualification)
{
    reason.clear();
    if (mesh_type == MeshType::VariableDz && terrain_qualification) {
        mesh = erf_auxiliary::DiffusionMesh::VariableDz;
        return true;
    }
    if (mesh_type == MeshType::ConstantDz) {
        mesh = erf_auxiliary::DiffusionMesh::ConstantDz;
        return true;
    }
    if (mesh_type == MeshType::StretchedDz) {
        mesh = erf_auxiliary::DiffusionMesh::StretchedDz;
        return true;
    }
    reason = "SBM particle mixing on a terrain-fitted (VariableDz) mesh has not passed its "
             "release gates; erf.sbm_particle_mixing_terrain_qualification enables it for "
             "qualification only";
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
    return build_particle_mixing_transmissibility(
        config, geom, metrics, density, density_component, eddy_diffusivity,
        erf_auxiliary::MappedDiffusionBoundary{}, transmissibility, diagnostic);
}

bool
build_particle_mixing_transmissibility (
    const ParticleMixingConfig& config,
    const amrex::Geometry& geom,
    const erf_auxiliary::MappedDiffusionMetrics& metrics,
    const amrex::MultiFab& density,
    const int density_component,
    const amrex::MultiFab* eddy_diffusivity,
    const erf_auxiliary::MappedDiffusionBoundary& boundary,
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
            eddy_diffusivity, boundary, diagnostic)) {
        diagnostic = std::string("SBM particle mixing ") +
                     particle_mixing_name(config.policy) + ": " + diagnostic;
        return false;
    }
    return true;
}

namespace {

// The unclamped drop of cell (i,j,k) in cell thicknesses, or infinity where
// the thickness is not positive or the drop not finite.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
terrain_cell_drop (const amrex::Array4<const amrex::Real>& z,
                   const int i, const int j, const int k) noexcept
{
    const TerrainCellDrops d = ComputeTerrainCellDrops(i, j, k, z);
    const amrex::Real drop = std::sqrt(d.dx_drop * d.dx_drop + d.dy_drop * d.dy_drop) / d.dz;
    return (amrex::Math::isfinite(drop) && d.dz > amrex::Real(0.0))
               ? drop : std::numeric_limits<amrex::Real>::infinity();
}

} // namespace

bool
admit_terrain_mixing_geometry (const amrex::MultiFab& z_nd,
                               const amrex::BoxArray& cells,
                               const amrex::Real max_drop,
                               std::string& reason)
{
    reason.clear();
    if (z_nd.boxArray() != amrex::convert(cells, amrex::IntVect(1))) {
        reason = "the terrain admission requires nodal heights on the cell layout";
        return false;
    }
    // The worst drop, with -1 for a nonfinite or nonpositive thickness, then
    // the first cell (in a fixed linear order) attaining it.
    const amrex::Box domain = cells.minimalBox();
    const auto lo = amrex::lbound(domain);
    const auto len = amrex::length(domain);
    amrex::ReduceOps<amrex::ReduceOpMax> max_op;
    amrex::ReduceData<amrex::Real> max_data(max_op);
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        const auto z = z_nd.const_array(mfi);
        max_op.eval(amrex::enclosedCells(mfi.validbox()), max_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real> {
                return {terrain_cell_drop(z, i, j, k)};
            });
    }
    amrex::Real worst = amrex::get<0>(max_data.value(max_op));
    amrex::ParallelAllReduce::Max(worst, amrex::ParallelContext::CommunicatorSub());
    if (worst <= max_drop) { return true; }

    amrex::ReduceOps<amrex::ReduceOpMin> min_op;
    amrex::ReduceData<amrex::Long> min_data(min_op);
    for (amrex::MFIter mfi(z_nd, false); mfi.isValid(); ++mfi) {
        const auto z = z_nd.const_array(mfi);
        min_op.eval(amrex::enclosedCells(mfi.validbox()), min_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Long> {
                const amrex::Real drop = terrain_cell_drop(z, i, j, k);
                const amrex::Long index = (i - lo.x) + amrex::Long(len.x) *
                                          ((j - lo.y) + amrex::Long(len.y) * (k - lo.z));
                return {drop >= worst ? index : std::numeric_limits<amrex::Long>::max()};
            });
    }
    amrex::Long first = amrex::get<0>(min_data.value(min_op));
    amrex::ParallelAllReduce::Min(first, amrex::ParallelContext::CommunicatorSub());
    const amrex::Long i = first % len.x;
    const amrex::Long j = (first / len.x) % len.y;
    const amrex::Long k = first / (amrex::Long(len.x) * len.y);
    std::ostringstream message;
    message.precision(6);
    message << "SBM terrain particle mixing is qualified for terrain dropping at most "
            << max_drop << " cell thickness across a cell; cell (" << i + lo.x << ","
            << j + lo.y << "," << k + lo.z << ") drops ";
    if (std::isfinite(worst)) {
        message << worst;
    } else {
        message << "with a nonfinite or nonpositive thickness";
    }
    reason = message.str();
    return false;
}

erf_auxiliary::MappedDiffusionBoundary
particle_mixing_boundary (const SBMPhysicalBoundaryConfig& boundaries,
                          const std::array<const amrex::MultiFab*, AMREX_SPACEDIM>& carrier)
{
    erf_auxiliary::MappedDiffusionBoundary boundary;
    boundary.kind = boundaries.diffusion;
    boundary.carrier = carrier;
    return boundary;
}

} // namespace erf_sbm
