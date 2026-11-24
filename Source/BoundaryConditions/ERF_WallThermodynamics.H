#ifndef ERF_WALL_THERMODYNAMICS_H_
#define ERF_WALL_THERMODYNAMICS_H_

#include <AMReX_Array.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <type_traits>

namespace erf_wall_thermodynamics {

// These compact enums and the face record are copied into GPU kernels.  Keep
// parsing strings and all other host-only state outside this representation.
enum class ThermalMode : unsigned char {
    LegacyNumerical,
    FixedPhysicalTemperature,
    Adiabatic
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

// Non-owning device views for sidewall scalar fluxes. MOST is restricted to
// the horizontal plates; supported sidewall drag and rate laws need no panel
// temperature. Null offsets select the original uniform-temperature wall law.
struct TemperaturePanel {
    amrex::Real temperature_K;
    amrex::Real fraction;
};

struct TemperatureBandView {
    const int* offsets = nullptr;
    const TemperaturePanel* panels = nullptr;
    int k_lo = 0;
};

using SidewallBandViews = amrex::GpuArray<TemperatureBandView, 4>; // xlo, xhi, ylo, yhi
static_assert(std::is_trivially_copyable_v<TemperatureBandView>);

struct MomentumTransfer {
    MomentumModel model = MomentumModel::ResolvedNoSlip;
    CoefficientProvider provider = CoefficientProvider::None;
    amrex::Real C_D = amrex::Real(0.0);
    amrex::Real z0_m = amrex::Real(0.0);
    amrex::Real molecular_speed = amrex::Real(0.0); // nu / half-cell distance [m/s]
    amrex::Real resolved_row_rate = amrex::Real(0.0); // ordinary no-slip stencil [1/s]; guard only
};

struct ScalarTransfer {
    ScalarModel model = ScalarModel::ResolvedMolecular;
    CoefficientProvider provider = CoefficientProvider::None;
    amrex::Real coefficient = amrex::Real(0.0);
    amrex::Real z0 = amrex::Real(0.0);
    // Both speeds have units of m/s (alpha / half-cell distance); at most one is active.
    amrex::Real molecular_speed = amrex::Real(0.0); // additive part of a combined bulk law only
    // Resolved-only reference for the corner timestep. The resolved flux already uses alpha.
    amrex::Real resolved_speed = amrex::Real(0.0);
};

struct FaceWall {
    ThermalBoundary thermal;
    MoistureMode moisture = MoistureMode::LegacyNumerical;
    amrex::Real wetness = amrex::Real(1.0); // Prescribed RH over liquid on wet walls.
    MomentumTransfer momentum;
    ScalarTransfer heat;
    ScalarTransfer vapor;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> wall_velocity{};

    // Activity, not ownership: an adiabatic face still replaces the heat flux with zero.
    AMREX_GPU_HOST_DEVICE constexpr bool heat_is_conductive () const noexcept
    { return thermal.mode == ThermalMode::FixedPhysicalTemperature; }

    AMREX_GPU_HOST_DEVICE constexpr bool vapor_is_wet () const noexcept
    { return moisture == MoistureMode::WetEquilibrium; }

    AMREX_GPU_HOST_DEVICE constexpr bool momentum_uses_most () const noexcept
    {
        return momentum.model == MomentumModel::BulkAero &&
               momentum.provider == CoefficientProvider::MOST;
    }

    AMREX_GPU_HOST_DEVICE constexpr bool heat_uses_most () const noexcept
    {
        return heat_is_conductive() && heat.model == ScalarModel::BulkAero &&
               heat.provider == CoefficientProvider::MOST;
    }

    AMREX_GPU_HOST_DEVICE constexpr bool vapor_uses_most () const noexcept
    {
        return vapor_is_wet() && vapor.model == ScalarModel::BulkAero &&
               vapor.provider == CoefficientProvider::MOST;
    }

    AMREX_GPU_HOST_DEVICE constexpr bool momentum_requires_speed () const noexcept
    { return momentum.model != MomentumModel::ResolvedNoSlip; }

    AMREX_GPU_HOST_DEVICE constexpr bool heat_requires_speed () const noexcept
    { return heat_is_conductive() && heat.model != ScalarModel::ResolvedMolecular; }

    AMREX_GPU_HOST_DEVICE constexpr bool vapor_requires_speed () const noexcept
    { return vapor_is_wet() && vapor.model != ScalarModel::ResolvedMolecular; }

    AMREX_GPU_HOST_DEVICE constexpr bool has_wall_rate_channel () const noexcept
    { return momentum_requires_speed() || heat_requires_speed() || vapor_requires_speed(); }
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
