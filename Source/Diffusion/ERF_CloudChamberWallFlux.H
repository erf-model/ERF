#ifndef ERF_CLOUD_CHAMBER_WALL_FLUX_H_
#define ERF_CLOUD_CHAMBER_WALL_FLUX_H_

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>

#include "ERF_Constants.H"
#include "ERF_DataStruct.H"
#include "ERF_EOS.H"
#include "ERF_IndexDefines.H"
#include "ERF_MicrophysicsUtils.H"
#include "BoundaryConditions/ERF_MOSTUtils.H"
#include "BoundaryConditions/ERF_WallThermodynamics.H"

#include <cmath>
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
wall_uses_most (const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    using namespace erf_wall_thermodynamics;
    return (wall.momentum.model == MomentumModel::BulkAero &&
            wall.momentum.provider == CoefficientProvider::MOST) ||
        (wall.heat.model == ScalarModel::BulkAero &&
         wall.heat.provider == CoefficientProvider::MOST) ||
        (wall.vapor.model == ScalarModel::BulkAero &&
         wall.vapor.provider == CoefficientProvider::MOST);
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
        return wall.thermal.mode == ThermalMode::FixedPhysicalTemperature;
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
    const bool most_momentum =
        wall.momentum.model == MomentumModel::BulkAero &&
        wall.momentum.provider == CoefficientProvider::MOST;
    const bool most_heat =
        wall.heat.model == ScalarModel::BulkAero &&
        wall.heat.provider == CoefficientProvider::MOST;
    const bool most_vapor =
        wall.vapor.model == ScalarModel::BulkAero &&
        wall.vapor.provider == CoefficientProvider::MOST;
    if (!(most_momentum || most_heat || most_vapor)) {
        result.valid = 1;
        return result;
    }

    const amrex::Real theta_wall = wall.thermal.temperature_K *
        std::pow(p_0 / p_hse, rdOcp);
    amrex::Real qv_wall = qv_air;
    if (wall.moisture == MoistureMode::WetEquilibrium) {
        erf_qsatw(wall.thermal.temperature_K, p_hse * amrex::Real(0.01), qv_wall);
    }
    const amrex::Real theta_v_air = theta_air * (amrex::Real(1.0) + epsv * qv_air);
    const amrex::Real theta_v_wall = theta_wall * (amrex::Real(1.0) + epsv * qv_wall);
    const amrex::Real U_most = amrex::max(U_t, amrex::Real(0.01));
    const amrex::Real Ri_b = gravity_sign == 0 ? amrex::Real(0.0) :
        amrex::Real(gravity_sign) * CONST_GRAV * wall_distance / theta_air *
        (theta_v_air - theta_v_wall) / (U_most * U_most);
    const amrex::Real Ri = amrex::max(amrex::Real(-4.0),
        amrex::min(amrex::Real(4.0), Ri_b));

    const amrex::Real z0_m = wall.momentum.z0_m;
    const amrex::Real z0_h = wall.heat.z0 > amrex::Real(0.0) ?
        wall.heat.z0 : z0_m;
    const amrex::Real z0_q = wall.vapor.z0 > amrex::Real(0.0) ?
        wall.vapor.z0 : z0_m;
    const amrex::Real D_m = std::log(wall_distance / z0_m);
    const amrex::Real D_h = std::log(wall_distance / z0_h);
    const amrex::Real D_q = std::log(wall_distance / z0_q);

    similarity_funs sfuns;
    amrex::Real zeta = amrex::Real(0.0);
    bool converged = false;
    constexpr int max_iters = 40;
    for (int iter = 0; iter < max_iters; ++iter) {
        const amrex::Real zeta_old = zeta;
        const amrex::Real psi_m = sfuns.calc_psi_m2(zeta_old);
        const amrex::Real psi_h = sfuns.calc_psi_h2(zeta_old);
        const amrex::Real A_m = amrex::max(D_m - psi_m, amrex::Real(1.0));
        const amrex::Real A_h = amrex::max(D_h - psi_h, amrex::Real(1.0));
        const amrex::Real zeta_new = amrex::Real(0.5) * zeta_old +
            amrex::Real(0.5) * Ri * A_m * A_m / A_h;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            std::isfinite(static_cast<double>(zeta_new)),
            "Cloud Chamber MOST produced a nonfinite stability coordinate");
        zeta = zeta_new;
        if (std::abs(zeta - zeta_old) <= amrex::Real(1.0e-3)) {
            converged = true;
            break;
        }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        converged, "Cloud Chamber MOST fixed-point iteration did not converge");

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
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        std::isfinite(static_cast<double>(result.C_D)) &&
        std::isfinite(static_cast<double>(result.C_H)) &&
        std::isfinite(static_cast<double>(result.C_E)),
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
    // The maximum absolute row sum of the two-component tangential
    // Jacobian of C_D |u_t| u_t is (3 + sqrt(2))/2 * C_D |u_t|.
    return (amrex::Real(3.0) + std::sqrt(amrex::Real(2.0))) / amrex::Real(2.0);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
momentum_row_sum_rate (
    int component,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& low_rate,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& high_rate) noexcept
{
    amrex::Real result = amrex::Real(0.0);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        if (dir != component) {
            result += low_rate[dir] + high_rate[dir];
        }
    }
    return result;
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
    return momentum_infinity_row_sum_factor() * C_D * U_t * dx_inv;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
requires_tangential_speed (
    const erf_wall_thermodynamics::FaceWall& wall,
    ScalarChannel channel) noexcept
{
    using namespace erf_wall_thermodynamics;

    switch (channel) {
    case ScalarChannel::Heat:
        return wall.thermal.mode == ThermalMode::FixedPhysicalTemperature &&
               (wall.heat.model == ScalarModel::BulkAero ||
                wall.heat.model == ScalarModel::NeutralRoughnessLog);
    case ScalarChannel::Vapor:
        return wall.moisture == MoistureMode::WetEquilibrium &&
               (wall.vapor.model == ScalarModel::BulkAero ||
                wall.vapor.model == ScalarModel::NeutralRoughnessLog);
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
    return wall.momentum.model !=
        erf_wall_thermodynamics::MomentumModel::ResolvedNoSlip;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
bool
wall_rate_requires_tangential_speed (
    const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    return requires_tangential_momentum(wall) ||
        requires_tangential_speed(wall, ScalarChannel::Heat) ||
        requires_tangential_speed(wall, ScalarChannel::Vapor);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wall_rate_for_face (
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real U_t,
    amrex::Real dx_inv,
    const FaceRuntimeCoefficients& runtime) noexcept
{
    using namespace erf_wall_thermodynamics;
    if (!wall_rate_requires_tangential_speed(wall)) {
        return amrex::Real(0.0);
    }
    if (wall_uses_most(wall)) {
        AMREX_ALWAYS_ASSERT(runtime.valid != 0);
    }
    const amrex::Real wall_distance = amrex::Real(0.5) / dx_inv;
    amrex::Real rate = momentum_rate_for_face(wall, U_t, dx_inv, runtime);

    if (wall.thermal.mode == ThermalMode::FixedPhysicalTemperature) {
        if (wall.heat.model == ScalarModel::BulkAero) {
            const amrex::Real C_H = wall.heat.provider == CoefficientProvider::MOST ?
                runtime.C_H : wall.heat.coefficient;
            rate = amrex::max(rate, neutral_scalar_rate(C_H, U_t, dx_inv));
        } else if (wall.heat.model == ScalarModel::NeutralRoughnessLog) {
            rate = amrex::max(rate, neutral_scalar_rate(
                neutral_log_scalar_coefficient(
                    wall_distance, wall.momentum.z0_m, wall.heat.z0),
                U_t, dx_inv));
        }
    }
    if (wall.moisture == MoistureMode::WetEquilibrium) {
        if (wall.vapor.model == ScalarModel::BulkAero) {
            const amrex::Real C_E = wall.vapor.provider == CoefficientProvider::MOST ?
                runtime.C_E : wall.vapor.coefficient;
            rate = amrex::max(rate, neutral_scalar_rate(C_E, U_t, dx_inv));
        } else if (wall.vapor.model == ScalarModel::NeutralRoughnessLog) {
            rate = amrex::max(rate, neutral_scalar_rate(
                neutral_log_scalar_coefficient(
                    wall_distance, wall.momentum.z0_m, wall.vapor.z0),
                U_t, dx_inv));
        }
    }
    return rate;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wall_rate_for_face (
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real U_t,
    amrex::Real dx_inv) noexcept
{
    return wall_rate_for_face(wall, U_t, dx_inv, FaceRuntimeCoefficients{});
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
wall_rate_for_face (
    const erf_wall_thermodynamics::FaceWall& wall,
    amrex::Real theta_air,
    amrex::Real qv_air,
    amrex::Real p_hse,
    amrex::Real rdOcp,
    amrex::Real U_t,
    amrex::Real dx_inv,
    int gravity_sign) noexcept
{
    const FaceRuntimeCoefficients runtime = most_wall_coefficients(
        wall, theta_air, qv_air, p_hse, rdOcp, U_t,
        amrex::Real(0.5) / dx_inv, gravity_sign);
    return wall_rate_for_face(wall, U_t, dx_inv, runtime);
}

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
                    amrex::Real rdOcp) noexcept
{
    const amrex::Real theta_wall =
        wall_temperature_K * std::pow(p_0 / p_hse, rdOcp);
    return rho * C_H * U_t * (theta_wall - theta_air);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
bulk_vapor_flux_in (amrex::Real rho,
                    amrex::Real C_E,
                    amrex::Real U_t,
                    amrex::Real wall_temperature_K,
                    amrex::Real qv_air,
                    amrex::Real p_hse) noexcept
{
    amrex::Real qv_wall = amrex::Real(0.0);
    erf_qsatw(wall_temperature_K, p_hse * amrex::Real(0.01), qv_wall);
    return rho * C_E * U_t * (qv_wall - qv_air);
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
    amrex::Real wall_temperature_K,
    amrex::Real qv_air,
    amrex::Real p_hse,
    amrex::Real dx_inv) noexcept
{
    amrex::Real qv_wall = amrex::Real(0.0);
    erf_qsatw(
        wall_temperature_K,
        p_hse * amrex::Real(0.01),
        qv_wall);
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
        const amrex::Real factor = -sample.rho * C_D * sample.U_t;
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
        if (sample.U_t == amrex::Real(0.0) &&
            wall.heat.model != ScalarModel::ResolvedMolecular) {
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
                sample.rdOcp);
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
        if (wall.moisture == MoistureMode::DryImpermeable) {
            return result;
        }
        if (sample.U_t == amrex::Real(0.0) &&
            wall.vapor.model != ScalarModel::ResolvedMolecular) {
            return result;
        }
        switch (wall.vapor.model) {
        case ScalarModel::ResolvedMolecular:
            result.rhoQv_in = resolved_vapor_flux_in(
                sample.rho,
                sample.alpha_C,
                wall.thermal.temperature_K,
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
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse);
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
                wall.thermal.temperature_K,
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
                 const bool cloudy)
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

        const WallFlux physical =
            evaluate_scalar_flux_in(wall, channel, sample);
        if (!owns(physical, channel)) { return; }

        const ScalarFluxReplacement update =
            scalar_flux_replacement<DIR,HIGH>(
                old_flux,
                scalar_flux_in(physical, channel),
                dx_inv[DIR]);

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
                   const bool cloudy = false)
{
    apply_face<0,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
    apply_face<0,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
    apply_face<1,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
    apply_face<1,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp,cloudy);
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
