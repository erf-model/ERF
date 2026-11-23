#ifndef ERF_CLOUD_CHAMBER_WALL_FLUX_H_
#define ERF_CLOUD_CHAMBER_WALL_FLUX_H_

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_Math.H>
#include <AMReX_MultiFab.H>

#include "ERF_Constants.H"
#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"
#include "BoundaryConditions/ERF_MOSTUtils.H"
#include "BoundaryConditions/ERF_WallThermodynamics.H"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <type_traits>

namespace erf_cloud_chamber_wall_flux {

enum OwnedChannel : unsigned char {
    OwnNone = 0,
    OwnMomentum = 1 << 0,
    OwnHeat = 1 << 1,
    OwnVapor = 1 << 2,
    OwnCloudWater = 1 << 3
};

struct WallFlux {
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> traction_on_fluid{};
    amrex::Real rhoTheta_in = amrex::Real(0.0);
    amrex::Real rhoQv_in = amrex::Real(0.0);
    amrex::Real rhoQc_in = amrex::Real(0.0);
    unsigned char owned_channels = OwnNone;
};

struct FaceRuntimeCoefficients {
    amrex::Real C_D = amrex::Real(0.0);
    amrex::Real C_H = amrex::Real(0.0);
    amrex::Real C_E = amrex::Real(0.0);
    amrex::Real u_star = amrex::Real(0.0);
    amrex::Real zeta = amrex::Real(0.0);
    unsigned char valid = 0;
};

struct MostStabilityResult {
    amrex::Real zeta = amrex::Real(0.0);
    int iterations = 0;
};

enum class ScalarChannel : unsigned char {
    None,
    Heat,
    Vapor,
    CloudWater
};

struct ScalarWallSample {
    amrex::Real rho = amrex::Real(0.0);
    amrex::Real scalar_air = amrex::Real(0.0);
    amrex::Real p_hse = amrex::Real(0.0);
    amrex::Real U_t = amrex::Real(0.0);
    amrex::Real alpha_T = amrex::Real(0.0);
    amrex::Real alpha_C = amrex::Real(0.0);
    amrex::Real dx_inv = amrex::Real(0.0);
    amrex::Real rdOcp = amrex::Real(0.0);
    // Explicit distance of this scalar sample from the physical wall.
    amrex::Real wall_distance = amrex::Real(0.0);
    amrex::Real theta_air = amrex::Real(0.0);
    amrex::Real qv_air = amrex::Real(0.0);
    int gravity_sign = 0;
};

struct NeutralLogMomentumState {
    amrex::Real log_m = amrex::Real(0.0);
    amrex::Real C_D = amrex::Real(0.0);
    amrex::Real u_star = amrex::Real(0.0);
};

struct MomentumWallSample {
    amrex::Real rho = amrex::Real(0.0);
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u_t{};
    amrex::Real U_t = amrex::Real(0.0);
    amrex::Real wall_distance = amrex::Real(0.0);
    amrex::Real theta_air = amrex::Real(0.0);
    amrex::Real qv_air = amrex::Real(0.0);
    amrex::Real p_hse = amrex::Real(0.0);
    amrex::Real rdOcp = amrex::Real(0.0);
    int gravity_sign = 0;
};

struct ScalarFluxReplacement {
    amrex::Real coordinate_flux = amrex::Real(0.0);
    amrex::Real rhs_delta = amrex::Real(0.0);
};

static_assert(std::is_trivially_copyable_v<WallFlux>);
static_assert(std::is_standard_layout_v<WallFlux>);
static_assert(std::is_trivially_copyable_v<FaceRuntimeCoefficients>);
static_assert(std::is_standard_layout_v<FaceRuntimeCoefficients>);
static_assert(std::is_trivially_copyable_v<MostStabilityResult>);
static_assert(std::is_standard_layout_v<MostStabilityResult>);
static_assert(std::is_trivially_copyable_v<ScalarWallSample>);
static_assert(std::is_standard_layout_v<ScalarWallSample>);
static_assert(std::is_trivially_copyable_v<NeutralLogMomentumState>);
static_assert(std::is_standard_layout_v<NeutralLogMomentumState>);
static_assert(std::is_trivially_copyable_v<MomentumWallSample>);
static_assert(std::is_standard_layout_v<MomentumWallSample>);
static_assert(std::is_trivially_copyable_v<ScalarFluxReplacement>);
static_assert(std::is_standard_layout_v<ScalarFluxReplacement>);

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_uses_most_heat (const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return wall.heat_uses_most();
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_uses_most_vapor (const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return wall.vapor_uses_most();
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_uses_most (const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return wall.momentum_uses_most() || wall.heat_uses_most() || wall.vapor_uses_most();
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
ScalarChannel
scalar_channel_from_quantity (int quantity) noexcept
{
    if (quantity == RhoTheta_comp) { return ScalarChannel::Heat; }
    if (quantity == RhoQ1_comp) { return ScalarChannel::Vapor; }
    if (quantity == RhoQ2_comp) { return ScalarChannel::CloudWater; }
    return ScalarChannel::None;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
int
primitive_component (ScalarChannel channel) noexcept
{
    switch (channel) {
    case ScalarChannel::Heat:
        return PrimTheta_comp;
    case ScalarChannel::Vapor:
        return PrimQ1_comp;
    case ScalarChannel::CloudWater:
        return PrimQ2_comp;
    case ScalarChannel::None:
        break;
    }
    return -1;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_owns_scalar_channel (
    const erf_wall_thermodynamics::FaceWall& wall,
    ScalarChannel channel) noexcept
{
    using namespace erf_wall_thermodynamics;

    switch (channel) {
    case ScalarChannel::Heat:
        return wall.thermal.mode != ThermalMode::LegacyNumerical;
    case ScalarChannel::Vapor:
    case ScalarChannel::CloudWater:
        return wall.moisture != MoistureMode::LegacyNumerical;
    case ScalarChannel::None:
        return false;
    }

    return false;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
unsigned char
scalar_channel_mask (ScalarChannel channel) noexcept
{
    switch (channel) {
    case ScalarChannel::Heat:
        return OwnHeat;
    case ScalarChannel::Vapor:
        return OwnVapor;
    case ScalarChannel::CloudWater:
        return OwnCloudWater;
    case ScalarChannel::None:
        return OwnNone;
    }

    return OwnNone;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
owns (const WallFlux& flux, ScalarChannel channel) noexcept
{
    return (flux.owned_channels & scalar_channel_mask(channel)) != 0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
NeutralLogMomentumState
neutral_log_momentum_state (amrex::Real U_t,
                           amrex::Real wall_distance,
                           amrex::Real z0_m) noexcept
{
    const amrex::Real log_m = std::log(wall_distance / z0_m);
    const amrex::Real k_over_log = KAPPA / log_m;
    return {log_m, k_over_log * k_over_log, k_over_log * U_t};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
neutral_log_scalar_coefficient (amrex::Real wall_distance,
                                amrex::Real z0_m,
                                amrex::Real z0_scalar) noexcept
{
    const amrex::Real log_m = std::log(wall_distance / z0_m);
    const amrex::Real log_s = std::log(wall_distance / z0_scalar);
    return (KAPPA * KAPPA) / (log_m * log_s);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
most_stability_require (bool condition, const char* message) noexcept
{
    if (!condition) {
        // AMReX device assertions and Abort can be disabled by NDEBUG.
        AMREX_IF_ON_DEVICE((AMREX_DEVICE_PRINTF("%s\n", message);))
#if defined(AMREX_USE_CUDA)
        AMREX_IF_ON_DEVICE((__trap();))
#else
        AMREX_IF_ON_DEVICE((__builtin_trap();))
#endif
        AMREX_IF_ON_HOST((amrex::Abort(message); std::abort();))
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
MostStabilityResult
most_stability_coordinate (amrex::Real Ri,
                           amrex::Real D_m,
                           amrex::Real D_s) noexcept
{
    similarity_funs sfuns;
    constexpr amrex::Real tolerance = amrex::Real(1.0e-3);
    if (Ri > amrex::Real(0.0)) {
        int iterations = 0;
        const auto state =
            [=] AMREX_GPU_HOST_DEVICE(
                amrex::Real z) noexcept -> amrex::GpuArray<amrex::Real, 3> {
            const amrex::Real psi_m = sfuns.calc_psi_m2(z);
            const amrex::Real psi_s = sfuns.calc_psi_h2(z);
            const amrex::Real A_m = amrex::max(D_m - psi_m, amrex::Real(1.0));
            const amrex::Real A_s = amrex::max(D_s - psi_s, amrex::Real(1.0));
            const amrex::Real rhs = Ri * A_m * A_m / A_s;
            most_stability_require(
                amrex::Math::isfinite(static_cast<double>(psi_m)) &&
                    amrex::Math::isfinite(static_cast<double>(psi_s)) &&
                    amrex::Math::isfinite(static_cast<double>(A_m)) &&
                    amrex::Math::isfinite(static_cast<double>(A_s)) &&
                    amrex::Math::isfinite(static_cast<double>(rhs)) &&
                    amrex::Math::isfinite(static_cast<double>(z - rhs)),
                "Cloud Chamber MOST produced a nonfinite stable residual");
            return {A_m, A_s, z - rhs};
        };
        most_stability_require(
            amrex::Math::isfinite(static_cast<double>(Ri)) &&
                amrex::Math::isfinite(static_cast<double>(D_m)) &&
                amrex::Math::isfinite(static_cast<double>(D_s)) &&
                D_m > amrex::Real(0.0) && D_s > amrex::Real(0.0),
            "Cloud Chamber MOST requires finite positive stable resistances");

        const auto interpolate = [=] AMREX_GPU_HOST_DEVICE(
                                     amrex::Real a, amrex::Real b,
                                     amrex::Real fa, amrex::Real fb) noexcept {
            most_stability_require(
                b > a && b - a <= tolerance && fa <= amrex::Real(0.0) &&
                    fb >= amrex::Real(0.0) && fb > fa &&
                    amrex::Math::isfinite(static_cast<double>(fb - fa)),
                "Cloud Chamber MOST requires a certified stable root bracket");
            const amrex::Real z = a + (b - a) * (-fa / (fb - fa));
            most_stability_require(
                amrex::Math::isfinite(static_cast<double>(z)) && z >= a &&
                    z <= b,
                "Cloud Chamber MOST produced an invalid bracket interpolation");
            return z;
        };

        amrex::Real seed = amrex::Real(0.0);
        for (int iter = 0; iter < 100; ++iter) {
            const auto value = state(seed);
            ++iterations;
            if (value[2] == amrex::Real(0.0)) {
                return {seed, iterations};
            }
            const amrex::Real next = seed - amrex::Real(0.5) * value[2];
            most_stability_require(
                amrex::Math::isfinite(static_cast<double>(next)),
                "Cloud Chamber MOST produced a nonfinite stability coordinate");
            const bool small_update =
                std::abs(next - seed) <= amrex::Real(0.125) * tolerance;
            seed = next;
            if (small_update) {
                break;
            }
        }

        // Keep the iteration's branch only if a root really lies nearby.
        amrex::Real a = amrex::max(amrex::Real(0.0), seed - tolerance);
        amrex::Real b = seed + tolerance;
        auto va = state(a);
        auto vb = state(b);
        if (va[2] <= amrex::Real(0.0) && vb[2] >= amrex::Real(0.0)) {
            for (int iter = 0; iter < 64; ++iter) {
                ++iterations;
                if (va[2] == amrex::Real(0.0)) {
                    return {a, iterations};
                }
                if (vb[2] == amrex::Real(0.0)) {
                    return {b, iterations};
                }
                const amrex::Real mid = a + amrex::Real(0.5) * (b - a);
                if (b - a <= tolerance) {
                    return {interpolate(a, b, va[2], vb[2]), iterations};
                }
                most_stability_require(
                    mid > a && mid < b,
                    "Cloud Chamber MOST exhausted stable root precision");
                const auto vm = state(mid);
                if (vm[2] < amrex::Real(0.0)) {
                    a = mid;
                    va = vm;
                } else {
                    b = mid;
                    vb = vm;
                }
            }
        }

        // Search from neutral in order: endpoint signs alone can miss a
        // narrow first-root pair. Stable A_m and A_s increase with z, giving
        // F(z) <= b - Ri*A_m(a)^2/A_s(b) throughout [a,b].
        constexpr int max_depth = 64;
        amrex::GpuArray<amrex::Real, max_depth> pending{};
        int depth = 0;
        amrex::Real unresolved = amrex::Real(-1.0);
        a = amrex::Real(0.0);
        b = amrex::max(amrex::Real(1.0), seed);
        va = state(a);
        vb = state(b);
        for (int search = 0; search < 8192; ++search) {
            ++iterations;
            most_stability_require(
                amrex::Math::isfinite(static_cast<double>(a)) &&
                    amrex::Math::isfinite(static_cast<double>(b)) && b > a &&
                    va[2] <= amrex::Real(0.0),
                "Cloud Chamber MOST exhausted stable bracket range");
            if (unresolved >= amrex::Real(0.0)) {
                const amrex::Real limit = unresolved + tolerance;
                most_stability_require(
                    a < limit,
                    "Cloud Chamber MOST cannot resolve the stable root branch");
                if (b > limit) {
                    b = limit;
                    vb = state(b);
                }
            }
            if (va[2] == amrex::Real(0.0)) {
                return {a, iterations};
            }
            const amrex::Real mid = a + amrex::Real(0.5) * (b - a);
            if (vb[2] >= amrex::Real(0.0) && b - a <= tolerance) {
                return {interpolate(a, b, va[2], vb[2]), iterations};
            }
            const amrex::Real lower_rhs = Ri * va[0] * va[0] / vb[1];
            const amrex::Real roundoff =
                amrex::Real(16.0) *
                std::numeric_limits<amrex::Real>::epsilon() *
                amrex::max(b, lower_rhs);
            if (vb[2] < amrex::Real(0.0) && b - lower_rhs < -roundoff) {
                a = b;
                va = vb;
                b = depth > 0 ? pending[--depth] : amrex::Real(2.0) * b;
            } else if (depth < max_depth && mid > a && mid < b) {
                pending[depth++] = b;
                b = mid;
            } else {
                // Retain a rounding-ambiguous interval in the location error
                // bound; never step past it to a distant root.
                most_stability_require(
                    b - a <= tolerance && vb[2] < amrex::Real(0.0),
                    "Cloud Chamber MOST exhausted stable root precision");
                if (unresolved < amrex::Real(0.0)) {
                    unresolved = a;
                }
                a = b;
                va = vb;
                b = depth > 0 ? pending[--depth] : amrex::Real(2.0) * b;
            }
            most_stability_require(
                amrex::Math::isfinite(static_cast<double>(b)) && b > a,
                "Cloud Chamber MOST exhausted stable bracket range");
            vb = state(b);
        }
        most_stability_require(
            false, "Cloud Chamber MOST exhausted stable root search budget");
        return {std::numeric_limits<amrex::Real>::quiet_NaN(), -iterations};
    }

    amrex::Real zeta = amrex::Real(0.0);
    bool converged = false;
    int iterations = 0;
    constexpr int max_iters = 100;
    for (int iter = 0; iter < max_iters; ++iter) {
        const amrex::Real zeta_old = zeta;
        const amrex::Real psi_m = sfuns.calc_psi_m2(zeta_old);
        const amrex::Real psi_s = sfuns.calc_psi_h2(zeta_old);
        const amrex::Real A_m = amrex::max(D_m - psi_m, amrex::Real(1.0));
        const amrex::Real A_s = amrex::max(D_s - psi_s, amrex::Real(1.0));
        const amrex::Real zeta_new = amrex::Real(0.5) * zeta_old +
            amrex::Real(0.5) * Ri * A_m * A_m / A_s;
        most_stability_require(
            amrex::Math::isfinite(static_cast<double>(zeta_new)),
            "Cloud Chamber MOST produced a nonfinite stability coordinate");
        zeta = zeta_new;
        iterations = iter + 1;
        if (std::abs(zeta - zeta_old) <= amrex::Real(1.0e-3)) {
            converged = true;
            break;
        }
    }
    most_stability_require(
        converged, "Cloud Chamber MOST fixed-point iteration did not converge");
    return {zeta, iterations};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wet_wall_qv (const erf_wall_thermodynamics::FaceWall& wall, amrex::Real p_hse) noexcept
{
    // Apply RH to vapor pressure (hPa) before the capped mixing-ratio conversion.
    // At RH=1 this is the same calculation as erf_qsatw.
    const amrex::Real e = wall.wetness * erf_esatw(wall.thermal.temperature_K);
    return erf_qsat_from_esat(e, p_hse * amrex::Real(0.01));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
FaceRuntimeCoefficients
most_wall_coefficients (
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real theta_air,
    amrex::Real qv_air,
    amrex::Real p_hse,
    amrex::Real rdOcp,
    amrex::Real U_t,
    amrex::Real wall_distance,
    int gravity_sign) noexcept
{
    using namespace erf_wall_thermodynamics;
    FaceRuntimeCoefficients result{};
    const bool most_momentum = wall.momentum_uses_most();
    const bool most_heat = wall_uses_most_heat(wall);
    const bool most_vapor = wall_uses_most_vapor(wall);
    if (!(most_momentum || most_heat || most_vapor)) {
        result.valid = 1;
        return result;
    }

    amrex::Real qv_wall = qv_air;
    if (wall.vapor_is_wet()) {
        qv_wall = wet_wall_qv(wall, p_hse);
    }
    const bool adiabatic = wall.thermal.mode == ThermalMode::Adiabatic;
    amrex::Real delta_theta_v;
    if (adiabatic) {
        // theta_wall = theta_air: only wet-wall saturation can drive buoyancy.
        // Difference the mixing ratios first so equal humidities are exactly neutral.
        delta_theta_v = theta_air * epsv * (qv_air - qv_wall);
    } else {
        const amrex::Real theta_wall = wall.thermal.temperature_K *
            std::pow(p_0 / p_hse, rdOcp);
        const amrex::Real theta_v_air = theta_air * (amrex::Real(1.0) + epsv * qv_air);
        const amrex::Real theta_v_wall = theta_wall * (amrex::Real(1.0) + epsv * qv_wall);
        delta_theta_v = theta_v_air - theta_v_wall;
    }
    const amrex::Real U_most = amrex::max(U_t, amrex::Real(0.01));
    const amrex::Real Ri_b = gravity_sign == 0 ? amrex::Real(0.0) :
        amrex::Real(gravity_sign) * CONST_GRAV * wall_distance / theta_air *
        delta_theta_v / (U_most * U_most);
    const amrex::Real Ri = amrex::max(amrex::Real(-4.0),
        amrex::min(amrex::Real(4.0), Ri_b));

    const amrex::Real z0_m = wall.momentum.z0_m;
    const amrex::Real z0_q = wall.vapor.z0 > amrex::Real(0.0) ?
        wall.vapor.z0 : z0_m;
    const amrex::Real D_m = std::log(wall_distance / z0_m);
    // Keep returning all coefficients, but never sample an adiabatic heat roughness.
    const amrex::Real D_h = wall.heat_is_conductive() && wall.heat.z0 > amrex::Real(0.0) ?
        std::log(wall_distance / wall.heat.z0) : D_m;
    const amrex::Real D_q = std::log(wall_distance / z0_q);

    // A coupled MOST heat/vapor face is host-constrained to have equal scalar
    // roughness.  For a single active scalar channel, that channel supplies
    // the scalar resistance used by the bulk-Richardson fixed point.
    const amrex::Real D_stability = most_heat ? D_h : (most_vapor ? D_q : D_m);
    const auto stability = most_stability_coordinate(Ri, D_m, D_stability);
    const amrex::Real zeta = stability.zeta;
    most_stability_require(
        stability.iterations > 0 &&
            amrex::Math::isfinite(static_cast<double>(zeta)),
        "Cloud Chamber MOST returned an uncertified stability coordinate");

    similarity_funs sfuns;
    const amrex::Real psi_m = sfuns.calc_psi_m2(zeta);
    const amrex::Real psi_h = sfuns.calc_psi_h2(zeta);
    const amrex::Real A_m = amrex::max(D_m - psi_m, amrex::Real(1.0));
    const amrex::Real A_h = amrex::max(D_h - psi_h, amrex::Real(1.0));
    const amrex::Real A_q = amrex::max(D_q - psi_h, amrex::Real(1.0));
    result.C_D = (KAPPA / A_m) * (KAPPA / A_m);
    result.C_H = KAPPA * KAPPA / (A_m * A_h);
    result.C_E = KAPPA * KAPPA / (A_m * A_q);
    result.u_star = KAPPA * U_most / A_m;
    result.zeta = zeta;
    result.valid = 1;
    most_stability_require(
        amrex::Math::isfinite(static_cast<double>(result.C_D)) &&
            amrex::Math::isfinite(static_cast<double>(result.C_H)) &&
            amrex::Math::isfinite(static_cast<double>(result.C_E)),
        "Cloud Chamber MOST produced nonfinite transfer coefficients");
    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
neutral_log_traction_on_fluid (const MomentumWallSample& sample,
                               const NeutralLogMomentumState& neutral) noexcept
{
    const amrex::Real factor =
        -sample.rho * neutral.C_D * sample.U_t;
    return {factor * sample.u_t[0],
            factor * sample.u_t[1],
            factor * sample.u_t[2]};
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
neutral_scalar_rate (amrex::Real coefficient,
                     amrex::Real U_t,
                     amrex::Real dx_inv) noexcept
{
    return coefficient * U_t * dx_inv;
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
momentum_infinity_row_sum_factor () noexcept
{
    // For a fixed C_D, the maximum absolute row sum of the two-component
    // tangential Jacobian of C_D |u_t| u_t is (3 + sqrt(2))/2 * C_D |u_t|.
    // Cloud Chamber MOST applies this factor to its frozen current C_D.
    return (amrex::Real(3.0) + std::sqrt(amrex::Real(2.0))) / amrex::Real(2.0);
}

// Totals for incident faces with the same normal (both sides for a one-cell direction).
struct DirectionalMomentumRate {
    amrex::Real modeled = 0;
    amrex::Real resolved = 0;
    bool combined = false;
};

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
momentum_row_sum_rate (int component, const amrex::GpuArray<DirectionalMomentumRate,
                       AMREX_SPACEDIM>& rates) noexcept
{
    // Only walls tangential to this component enter its row. A combined momentum
    // wall activates ordinary no-slip shear on that same row, not on all three.
    // Otherwise keep the modeled-only sum; never add a second stress to the RHS.
    amrex::Real modeled = 0, resolved = 0;
    bool combined = false;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (dir != component) {
            modeled += rates[dir].modeled;
            resolved += rates[dir].resolved;
            combined |= rates[dir].combined;
        }
    }
    return combined ? modeled + resolved : modeled;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
momentum_rate_for_face (
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real U_t,
    amrex::Real dx_inv,
    const FaceRuntimeCoefficients& runtime = FaceRuntimeCoefficients{}) noexcept
{
    using namespace erf_wall_thermodynamics;
    amrex::Real C_D = amrex::Real(0.0);
    const amrex::Real wall_distance = amrex::Real(0.5) / dx_inv;
    if (wall.momentum.model !=
        MomentumModel::ResolvedNoSlip) {
        if (wall.momentum.model == MomentumModel::NeutralRoughnessLog) {
            C_D = neutral_log_momentum_state(
                U_t, wall_distance, wall.momentum.z0_m).C_D;
        } else if (wall.momentum.model == MomentumModel::BulkAero) {
            if (wall.momentum.provider == CoefficientProvider::MOST) {
                AMREX_ALWAYS_ASSERT(runtime.valid != 0);
                C_D = runtime.C_D;
            } else {
                C_D = wall.momentum.C_D;
            }
        }
    }
    return momentum_infinity_row_sum_factor() * C_D * U_t * dx_inv +
        wall.momentum.molecular_speed * dx_inv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
requires_tangential_speed (
    const erf_wall_thermodynamics::FaceWall& wall,
    ScalarChannel channel) noexcept
{
    switch (channel) {
    case ScalarChannel::Heat:
        return wall.heat_requires_speed();
    case ScalarChannel::Vapor:
        return wall.vapor_requires_speed();
    case ScalarChannel::CloudWater:
    case ScalarChannel::None:
        return false;
    }

    return false;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
requires_tangential_momentum (
    const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return wall.momentum_requires_speed();
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_rate_requires_tangential_speed (
    const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return wall.has_wall_rate_channel();
}

// One scalar row at one cell. A resolved wall contributes 2*alpha/dx^2 to
// the row, but is not itself a wall-law timestep constraint. Only an active
// combined bulk face selects the sum in place of the legacy face maximum.
struct CellScalarWallRate {
    amrex::Real max_wall_law = 0;
    amrex::Real sum = 0;
    bool combined = false;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (const erf_wall_thermodynamics::ScalarTransfer& scalar,
              amrex::Real U_t, amrex::Real dx_inv, amrex::Real z0_m,
              amrex::Real most_coefficient) noexcept
    {
        using namespace erf_wall_thermodynamics;
        amrex::Real rate = 0;
        switch (scalar.model) {
        case ScalarModel::ResolvedMolecular:
            sum += scalar.resolved_speed * dx_inv;
            return;
        case ScalarModel::BulkAero: {
            const auto C = scalar.provider == CoefficientProvider::MOST ?
                most_coefficient : scalar.coefficient;
            rate = neutral_scalar_rate(C, U_t, dx_inv) + scalar.molecular_speed * dx_inv;
            // Setup requires positive diffusivity for bulk_aero_plus_molecular;
            // plain bulk has zero molecular speed.
            combined |= scalar.molecular_speed > 0;
            break;
        }
        case ScalarModel::NeutralRoughnessLog:
            rate = neutral_scalar_rate(neutral_log_scalar_coefficient(
                amrex::Real(0.5) / dx_inv, z0_m, scalar.z0), U_t, dx_inv);
            break;
        }
        sum += rate;
        max_wall_law = amrex::max(max_wall_law, rate);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real max_rate () const noexcept
    {
        return combined ? sum : max_wall_law;
    }
};

template <int DIR, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
to_coordinate_flux (amrex::Real inward) noexcept
{
    // Wall evaluators are inward-positive. ERF stores coordinate-positive
    // face fluxes, so this is the one low/high sign adapter for every scalar.
    amrex::ignore_unused(DIR);
    return HIGH ? -inward : inward;
}

template <int DIR, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
traction_from_local_inward (amrex::Real normal_in,
                            amrex::Real tangent_a,
                            amrex::Real tangent_b) noexcept
{
    // Physical traction result/orientation support. The stress adapter maps
    // this pure result into ERF tau_ij arrays for active wall models.
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> result{};
    result[DIR] = to_coordinate_flux<DIR,HIGH>(normal_in);
    if constexpr (DIR == 0) {
        result[1] = tangent_a;
        result[2] = tangent_b;
    } else if constexpr (DIR == 1) {
        result[0] = tangent_a;
        result[2] = tangent_b;
    } else {
        result[0] = tangent_a;
        result[1] = tangent_b;
    }
    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
bulk_theta_flux_in (amrex::Real rho,
                    amrex::Real C_H,
                    amrex::Real U_t,
                    amrex::Real wall_temperature_K,
                    amrex::Real theta_air,
                    amrex::Real p_hse,
                    amrex::Real rdOcp,
                    amrex::Real molecular_speed = amrex::Real(0)) noexcept
{
    const amrex::Real theta_wall =
        wall_temperature_K * std::pow(p_0 / p_hse, rdOcp);
    return (rho * C_H * U_t + rho * molecular_speed) * (theta_wall - theta_air);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
bulk_vapor_flux_in (amrex::Real rho,
                    amrex::Real C_E,
                    amrex::Real U_t,
                    const erf_wall_thermodynamics::FaceWall& wall,
                    amrex::Real qv_air,
                    amrex::Real p_hse,
                    amrex::Real molecular_speed = amrex::Real(0)) noexcept
{
    const amrex::Real qv_wall = wet_wall_qv(wall, p_hse);
    return (rho * C_E * U_t + rho * molecular_speed) * (qv_wall - qv_air);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>
relative_tangential_velocity_cell_centered (
    int dir, int i, int j, int k,
    const amrex::Array4<const amrex::Real>& u,
    const amrex::Array4<const amrex::Real>& v,
    const amrex::Array4<const amrex::Real>& w,
    const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    const amrex::Real uc =
        amrex::Real(0.5) * (u(i,j,k) + u(i+1,j,k)) - wall.wall_velocity[0];
    const amrex::Real vc =
        amrex::Real(0.5) * (v(i,j,k) + v(i,j+1,k)) - wall.wall_velocity[1];
    const amrex::Real wc =
        amrex::Real(0.5) * (w(i,j,k) + w(i,j,k+1)) - wall.wall_velocity[2];
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> result{uc, vc, wc};
    result[dir] = amrex::Real(0.0);
    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
tangential_magnitude (
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& u_t) noexcept
{
    return std::sqrt(u_t[0]*u_t[0] + u_t[1]*u_t[1] + u_t[2]*u_t[2]);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
tangential_speed (int dir,
                  int i, int j, int k,
                  const amrex::Array4<const amrex::Real>& u,
                  const amrex::Array4<const amrex::Real>& v,
                  const amrex::Array4<const amrex::Real>& w,
                  const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return tangential_magnitude(
        relative_tangential_velocity_cell_centered(dir, i, j, k, u, v, w, wall));
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
tangential_speed_cell_centered (int dir,
                                int i, int j, int k,
                                const amrex::Array4<const amrex::Real>& velocity,
                                const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    const amrex::Real uc = velocity(i,j,k,0) - wall.wall_velocity[0];
    const amrex::Real vc = velocity(i,j,k,1) - wall.wall_velocity[1];
    const amrex::Real wc = velocity(i,j,k,2) - wall.wall_velocity[2];
    const amrex::Real normal = (dir == 0) ? uc : ((dir == 1) ? vc : wc);
    const amrex::Real t0 = uc - ((dir == 0) ? normal : amrex::Real(0.0));
    const amrex::Real t1 = vc - ((dir == 1) ? normal : amrex::Real(0.0));
    const amrex::Real t2 = wc - ((dir == 2) ? normal : amrex::Real(0.0));
    return std::sqrt(t0*t0 + t1*t1 + t2*t2);
}

// Use both incident faces of a one-cell direction. For a MOST face freeze
// its current coefficients once for all channels (not a nonlinear Jacobian).
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wall_rate_for_cell (int i, int j, int k, const amrex::Box& domain,
                    const erf_wall_thermodynamics::Boundary& walls,
                    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                    const amrex::Array4<const amrex::Real>& velocity) noexcept
{
    using namespace erf_wall_thermodynamics;
    const amrex::IntVect cell(i,j,k);
    amrex::GpuArray<DirectionalMomentumRate, AMREX_SPACEDIM> momentum{};
    CellScalarWallRate heat{}, vapor{};
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        for (int side = 0; side < 2; ++side) {
            if (cell[dir] != (side ? domain.bigEnd(dir) : domain.smallEnd(dir))) { continue; }
            const auto& wall = walls[2*dir+side];
            const auto inv = dx_inv[dir];
            amrex::Real U_t = 0;
            FaceRuntimeCoefficients runtime{};
            if (wall_rate_requires_tangential_speed(wall)) {
                U_t = tangential_speed_cell_centered(dir, i, j, k, velocity, wall);
                if (wall_uses_most(wall)) {
                    runtime = most_wall_coefficients(
                        wall, velocity(i,j,k,3), velocity(i,j,k,4), velocity(i,j,k,5),
                        velocity(i,j,k,6), U_t, amrex::Real(0.5) / inv, dir == 2 ? 1-2*side : 0);
                    AMREX_ALWAYS_ASSERT(runtime.valid != 0);
                }
            }
            if (requires_tangential_momentum(wall)) {
                momentum[dir].modeled += momentum_rate_for_face(wall, U_t, inv, runtime);
                momentum[dir].combined |= wall.momentum.molecular_speed > 0;
            } else {
                momentum[dir].resolved += wall.momentum.resolved_row_rate;
            }
            if (wall.heat_is_conductive()) {
                heat.add(wall.heat, U_t, inv, wall.momentum.z0_m, runtime.C_H);
            }
            if (wall.vapor_is_wet()) {
                vapor.add(wall.vapor, U_t, inv, wall.momentum.z0_m, runtime.C_E);
            }
        }
    }
    auto rate = amrex::max(heat.max_rate(), vapor.max_rate());
    for (int component = 0; component < AMREX_SPACEDIM; ++component) {
        rate = amrex::max(rate, momentum_row_sum_rate(component, momentum));
    }
    return rate;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
resolved_theta_flux_in (
    amrex::Real rho,
    amrex::Real alpha_T,
    amrex::Real wall_temperature_K,
    amrex::Real theta_air,
    amrex::Real p_hse,
    amrex::Real dx_inv,
    amrex::Real rdOcp) noexcept
{
    const amrex::Real theta_wall =
        wall_temperature_K * std::pow(p_0 / p_hse, rdOcp);
    return rho * alpha_T * (theta_wall - theta_air) *
        amrex::Real(2.0) * dx_inv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
resolved_vapor_flux_in (
    amrex::Real rho,
    amrex::Real alpha_C,
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real qv_air,
    amrex::Real p_hse,
    amrex::Real dx_inv) noexcept
{
    const amrex::Real qv_wall = wet_wall_qv(wall, p_hse);
    return rho * alpha_C * (qv_wall - qv_air) *
        amrex::Real(2.0) * dx_inv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
WallFlux
evaluate_momentum_traction (
    const erf_wall_thermodynamics::FaceWall& wall,
    const MomentumWallSample& sample) noexcept
{
    using namespace erf_wall_thermodynamics;

    WallFlux result{};
    switch (wall.momentum.model) {
    case MomentumModel::ResolvedNoSlip:
        return result;
    case MomentumModel::NeutralRoughnessLog: {
        result.owned_channels = OwnMomentum;
        const auto neutral = neutral_log_momentum_state(
            sample.U_t, sample.wall_distance, wall.momentum.z0_m);
        result.traction_on_fluid =
            neutral_log_traction_on_fluid(sample, neutral);
        return result;
    }
    case MomentumModel::BulkAero: {
        result.owned_channels = OwnMomentum;
        if (sample.U_t == amrex::Real(0.0)) { return result; }
        amrex::Real C_D = wall.momentum.C_D;
        if (wall.momentum.provider == CoefficientProvider::MOST) {
            const auto runtime = most_wall_coefficients(
                wall, sample.theta_air, sample.qv_air, sample.p_hse,
                sample.rdOcp, sample.U_t, sample.wall_distance,
                sample.gravity_sign);
            AMREX_ALWAYS_ASSERT(runtime.valid != 0);
            C_D = runtime.C_D;
        }
        const amrex::Real factor = -sample.rho * C_D * sample.U_t - sample.rho * wall.momentum.molecular_speed;
        result.traction_on_fluid = {
            factor * sample.u_t[0], factor * sample.u_t[1], factor * sample.u_t[2]};
        return result;
    }
    }
    AMREX_ALWAYS_ASSERT(false);
    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
WallFlux
evaluate_scalar_flux_in (
    const erf_wall_thermodynamics::FaceWall& wall,
    ScalarChannel channel,
    const ScalarWallSample& sample) noexcept
{
    using namespace erf_wall_thermodynamics;

    WallFlux result{};
    if (!wall_owns_scalar_channel(wall, channel)) {
        return result;
    }

    result.owned_channels = scalar_channel_mask(channel);

    switch (channel) {
    case ScalarChannel::Heat:
        // An owned zero replaces ordinary diffusion and removes its RHS contribution.
        if (!wall.heat_is_conductive()) { return result; }
        if (sample.U_t == amrex::Real(0.0) &&
            wall.heat.model != ScalarModel::ResolvedMolecular && wall.heat.molecular_speed == 0) {
            return result;
        }
        switch (wall.heat.model) {
        case ScalarModel::ResolvedMolecular:
            result.rhoTheta_in = resolved_theta_flux_in(
                sample.rho,
                sample.alpha_T,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.dx_inv,
                sample.rdOcp);
            return result;

        case ScalarModel::BulkAero: {
            amrex::Real C_H = wall.heat.coefficient;
            if (wall.heat.provider == CoefficientProvider::MOST) {
                const auto runtime = most_wall_coefficients(
                    wall, sample.theta_air, sample.qv_air, sample.p_hse,
                    sample.rdOcp, sample.U_t, sample.wall_distance,
                    sample.gravity_sign);
                AMREX_ALWAYS_ASSERT(runtime.valid != 0);
                C_H = runtime.C_H;
            }
            result.rhoTheta_in = bulk_theta_flux_in(
                sample.rho,
                C_H,
                sample.U_t,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.rdOcp, wall.heat.molecular_speed);
            return result;
        }
        case ScalarModel::NeutralRoughnessLog:
            result.rhoTheta_in = bulk_theta_flux_in(
                sample.rho,
                neutral_log_scalar_coefficient(
                    sample.wall_distance,
                    wall.momentum.z0_m,
                    wall.heat.z0),
                sample.U_t,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.rdOcp);
            return result;
        }
        AMREX_ALWAYS_ASSERT(false);
        return result;

    case ScalarChannel::Vapor:
        if (!wall.vapor_is_wet()) { return result; }
        if (sample.U_t == amrex::Real(0.0) &&
            wall.vapor.model != ScalarModel::ResolvedMolecular && wall.vapor.molecular_speed == 0) {
            return result;
        }
        switch (wall.vapor.model) {
        case ScalarModel::ResolvedMolecular:
            result.rhoQv_in = resolved_vapor_flux_in(
                sample.rho,
                sample.alpha_C,
                wall,
                sample.scalar_air,
                sample.p_hse,
                sample.dx_inv);
            return result;

        case ScalarModel::BulkAero: {
            amrex::Real C_E = wall.vapor.coefficient;
            if (wall.vapor.provider == CoefficientProvider::MOST) {
                const auto runtime = most_wall_coefficients(
                    wall, sample.theta_air, sample.qv_air, sample.p_hse,
                    sample.rdOcp, sample.U_t, sample.wall_distance,
                    sample.gravity_sign);
                AMREX_ALWAYS_ASSERT(runtime.valid != 0);
                C_E = runtime.C_E;
            }
            result.rhoQv_in = bulk_vapor_flux_in(
                sample.rho,
                C_E,
                sample.U_t,
                wall,
                sample.scalar_air,
                sample.p_hse, wall.vapor.molecular_speed);
            return result;
        }
        case ScalarModel::NeutralRoughnessLog:
            result.rhoQv_in = bulk_vapor_flux_in(
                sample.rho,
                neutral_log_scalar_coefficient(
                    sample.wall_distance,
                    wall.momentum.z0_m,
                    wall.vapor.z0),
                sample.U_t,
                wall,
                sample.scalar_air,
                sample.p_hse);
            return result;
        }
        AMREX_ALWAYS_ASSERT(false);
        return result;

    case ScalarChannel::CloudWater:
        return result;

    case ScalarChannel::None:
        return result;
    }

    return result;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
scalar_flux_in (const WallFlux& flux, ScalarChannel channel) noexcept
{
    switch (channel) {
    case ScalarChannel::Heat:
        return flux.rhoTheta_in;
    case ScalarChannel::Vapor:
        return flux.rhoQv_in;
    case ScalarChannel::CloudWater:
        return flux.rhoQc_in;
    case ScalarChannel::None:
        return amrex::Real(0.0);
    }

    return amrex::Real(0.0);
}

// Each panel sees the same adjacent air and pressure. In particular evaluate
// wet-wall saturation per panel, not at an area-averaged temperature.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real evaluate_face_scalar_flux_in (
    const erf_wall_thermodynamics::FaceWall& wall, ScalarChannel channel,
    const ScalarWallSample& sample, const erf_wall_thermodynamics::TemperatureBandView& bands,
    int k) noexcept
{
    const bool uses_temperature = (channel == ScalarChannel::Heat && wall.heat_is_conductive()) ||
                                  (channel == ScalarChannel::Vapor && wall.vapor_is_wet());
    if (!bands.offsets || !uses_temperature) {
        return scalar_flux_in(evaluate_scalar_flux_in(wall, channel, sample), channel);
    }
    const int row = k - bands.k_lo;
    auto panel_wall = wall;
    amrex::Real inward = 0;
    for (int n = bands.offsets[row]; n < bands.offsets[row+1]; ++n) {
        const auto panel = bands.panels[n];
        panel_wall.thermal.temperature_K = panel.temperature_K;
        inward += panel.fraction * scalar_flux_in(evaluate_scalar_flux_in(panel_wall, channel, sample), channel);
    }
    return inward;
}

template <int DIR, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
ScalarFluxReplacement
scalar_flux_replacement (
    amrex::Real old_flux,
    amrex::Real inward_flux,
    amrex::Real dx_inv) noexcept
{
    const amrex::Real new_flux =
        to_coordinate_flux<DIR,HIGH>(inward_flux);
    const amrex::Real correction =
        (new_flux - old_flux) * dx_inv;
    return {
        new_flux,
        HIGH ? -correction : correction
    };
}

template <int DIR, bool HIGH>
AMREX_FORCE_INLINE
void apply_face (const amrex::Box& bx,
                 const amrex::Box& domain,
                 const int quantity,
                 const int flux_comp,
                 const amrex::Array4<const amrex::Real>& state,
                 const amrex::Array4<const amrex::Real>& prim,
                 const amrex::Array4<const amrex::Real>& base_state,
                 const amrex::Array4<const amrex::Real>& u,
                 const amrex::Array4<const amrex::Real>& v,
                 const amrex::Array4<const amrex::Real>& w,
                 const amrex::Array4<amrex::Real>& rhs,
                 const amrex::Array4<amrex::Real>& xflux,
                 const amrex::Array4<amrex::Real>& yflux,
                 const amrex::Array4<amrex::Real>& zflux,
                 const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                 const erf_wall_thermodynamics::Boundary& walls,
                 const amrex::Real alpha_T,
                 const amrex::Real alpha_C,
                 const amrex::Real rdOcp,
                 const bool cloudy,
                 const erf_wall_thermodynamics::TemperatureBandView bands = {})
{
    // A tile may belong to an interior BoxArray box. Do not relocate that
    // tile to a global boundary: its Array4 views do not own data there.
    if constexpr (HIGH) {
        if (bx.bigEnd(DIR) != domain.bigEnd(DIR)) { return; }
    } else {
        if (bx.smallEnd(DIR) != domain.smallEnd(DIR)) { return; }
    }

    const ScalarChannel channel =
        scalar_channel_from_quantity(quantity);
    if (channel == ScalarChannel::None) { return; }

    const auto wall = walls[2*DIR + (HIGH ? 1 : 0)];
    if (!wall_owns_scalar_channel(wall, channel)) { return; }

    const int prim_comp = primitive_component(channel);
    AMREX_ALWAYS_ASSERT(prim_comp >= 0);

    const int cell_face = HIGH ? domain.bigEnd(DIR) : domain.smallEnd(DIR);
    amrex::Box face_cells = bx;
    face_cells.setSmall(DIR, cell_face);
    face_cells.setBig(DIR, cell_face);
    face_cells &= domain;
    if (face_cells.isEmpty()) { return; }
    AMREX_ALWAYS_ASSERT(bx.contains(face_cells));

    const bool needs_speed = requires_tangential_speed(wall, channel);
    const auto& flux =
        (DIR == 0) ? xflux : ((DIR == 1) ? yflux : zflux);

    amrex::ParallelFor(face_cells,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const int fi = HIGH && DIR == 0 ? i + 1 : i;
        const int fj = HIGH && DIR == 1 ? j + 1 : j;
        const int fk = HIGH && DIR == 2 ? k + 1 : k;

        const amrex::Real old_flux = flux(fi,fj,fk,flux_comp);
        const ScalarWallSample sample{
            state(i,j,k,Rho_comp),
            prim(i,j,k,prim_comp),
            base_state(i,j,k,BaseState::p0_comp),
            needs_speed ?
                tangential_speed(DIR, i, j, k, u, v, w, wall) :
                amrex::Real(0.0),
            alpha_T,
            alpha_C,
            dx_inv[DIR],
            rdOcp,
            amrex::Real(0.5) / dx_inv[DIR],
            prim(i,j,k,PrimTheta_comp),
            cloudy ? prim(i,j,k,PrimQ1_comp) : amrex::Real(0.0),
            DIR == 2 ? (HIGH ? -1 : 1) : 0
        };

        // Ownership was checked for this face/channel before the launch.
        const auto inward = evaluate_face_scalar_flux_in(wall, channel, sample, bands, k);
        const ScalarFluxReplacement update =
            scalar_flux_replacement<DIR,HIGH>(old_flux, inward, dx_inv[DIR]);

        flux(fi,fj,fk,flux_comp) = update.coordinate_flux;
        rhs(i,j,k,quantity) += update.rhs_delta;
    });
}

inline void apply (const amrex::Box& bx, const amrex::Box& domain,
                   const int quantity, const int flux_comp,
                   const amrex::Array4<const amrex::Real>& state,
                   const amrex::Array4<const amrex::Real>& prim,
                   const amrex::Array4<const amrex::Real>& base_state,
                   const amrex::Array4<const amrex::Real>& u,
                   const amrex::Array4<const amrex::Real>& v,
                   const amrex::Array4<const amrex::Real>& w,
                   const amrex::Array4<amrex::Real>& rhs,
                   const amrex::Array4<amrex::Real>& xflux,
                   const amrex::Array4<amrex::Real>& yflux,
                   const amrex::Array4<amrex::Real>& zflux,
                   const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx_inv,
                   const erf_wall_thermodynamics::Boundary& walls,
                   const amrex::Real alpha_T, const amrex::Real alpha_C,
                   const amrex::Real rdOcp,
                   const bool cloudy = false,
                   const erf_wall_thermodynamics::SidewallBandViews& bands = {})
{
    apply_face<0,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy,bands[0]);
    apply_face<0,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy,bands[1]);
    apply_face<1,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy,bands[2]);
    apply_face<1,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy,bands[3]);
    apply_face<2,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
    apply_face<2,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
}

inline bool has_velocity_dependent_scalar_wall (
    const erf_wall_thermodynamics::Boundary& walls) noexcept
{
    for (const auto& wall : walls) {
        if (requires_tangential_speed(wall, ScalarChannel::Heat) ||
            requires_tangential_speed(wall, ScalarChannel::Vapor)) {
            return true;
        }
    }
    return false;
}

AMREX_FORCE_INLINE
amrex::Real
wall_dt_from_max_rate (amrex::Real max_rate,
                       amrex::Real wall_cfl = amrex::Real(0.5)) noexcept
{
    return max_rate > amrex::Real(0.0) ?
        wall_cfl / max_rate : amrex::Real(1.0e30);
}

} // namespace erf_cloud_chamber_wall_flux

#endif
