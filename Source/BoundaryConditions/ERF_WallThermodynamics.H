#ifndef ERF_WALL_THERMODYNAMICS_H_
#define ERF_WALL_THERMODYNAMICS_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include <type_traits>

namespace erf_wall_thermodynamics {

// These compact enums and the face record are copied into GPU kernels.  Keep
// parsing strings and all other host-only state outside this representation.
enum class ThermalMode : unsigned char {
    LegacyNumerical,
    FixedPhysicalTemperature
};

enum class MoistureMode : unsigned char {
    LegacyNumerical,
    DryImpermeable,
    WetEquilibrium
};

enum class MomentumModel : unsigned char {
    ResolvedNoSlip,
    BulkAero,
    NeutralRoughnessLog
};

enum class ScalarModel : unsigned char {
    ResolvedMolecular,
    BulkAero,
    NeutralRoughnessLog
};

enum class CoefficientProvider : unsigned char {
    None,
    Fixed,
    MOST
};

struct ThermalBoundary {
    ThermalMode mode = ThermalMode::LegacyNumerical;
    amrex::Real temperature_K = amrex::Real(0.0);
};

struct MomentumTransfer {
    MomentumModel model = MomentumModel::ResolvedNoSlip;
    CoefficientProvider provider = CoefficientProvider::None;
    amrex::Real C_D = amrex::Real(0.0);
    amrex::Real z0_m = amrex::Real(0.0);
};

struct ScalarTransfer {
    ScalarModel model = ScalarModel::ResolvedMolecular;
    CoefficientProvider provider = CoefficientProvider::None;
    amrex::Real coefficient = amrex::Real(0.0);
    amrex::Real z0 = amrex::Real(0.0);
};

struct FaceWall {
    ThermalBoundary thermal;
    MoistureMode moisture = MoistureMode::LegacyNumerical;
    MomentumTransfer momentum;
    ScalarTransfer heat;
    ScalarTransfer vapor;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wall_velocity{};
};

using Boundary = amrex::GpuArray<FaceWall, 2 * AMREX_SPACEDIM>;

static_assert(std::is_trivially_copyable_v<ThermalBoundary>);
static_assert(std::is_standard_layout_v<ThermalBoundary>);
static_assert(std::is_trivially_copyable_v<MomentumTransfer>);
static_assert(std::is_standard_layout_v<MomentumTransfer>);
static_assert(std::is_trivially_copyable_v<ScalarTransfer>);
static_assert(std::is_standard_layout_v<ScalarTransfer>);
static_assert(std::is_trivially_copyable_v<FaceWall>);
static_assert(std::is_standard_layout_v<FaceWall>);

} // namespace erf_wall_thermodynamics

#endif
