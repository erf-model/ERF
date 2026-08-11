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
static_assert(std::is_trivially_copyable_v<ScalarFluxReplacement>);
static_assert(std::is_standard_layout_v<ScalarFluxReplacement>);

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
bool
requires_tangential_speed (
    const erf_wall_thermodynamics::FaceWall& wall,
    ScalarChannel channel) noexcept
{
    using namespace erf_wall_thermodynamics;

    switch (channel) {
    case ScalarChannel::Heat:
        return wall.thermal.mode == ThermalMode::FixedPhysicalTemperature &&
               wall.heat.model == ScalarModel::BulkAero;
    case ScalarChannel::Vapor:
        return wall.moisture == MoistureMode::WetEquilibrium &&
               wall.vapor.model == ScalarModel::BulkAero;
    case ScalarChannel::CloudWater:
    case ScalarChannel::None:
        return false;
    }

    return false;
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
    // Physical traction result/orientation support for future momentum
    // closures. This does not map traction into ERF tau_ij arrays and does
    // not activate a production Cloud Chamber momentum wall model.
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
amrex::Real
tangential_speed (int dir,
                  int i, int j, int k,
                  const amrex::Array4<const amrex::Real>& u,
                  const amrex::Array4<const amrex::Real>& v,
                  const amrex::Array4<const amrex::Real>& w,
                  const erf_wall_thermodynamics::FaceWall& wall) noexcept
{
    // x/y/z velocities are staggered. These averages are the established
    // Cartesian face-to-cell convention used by ERF's scalar path.
    const amrex::Real uc =
        amrex::Real(0.5) * (u(i,j,k) + u(i+1,j,k)) - wall.wall_velocity[0];
    const amrex::Real vc =
        amrex::Real(0.5) * (v(i,j,k) + v(i,j+1,k)) - wall.wall_velocity[1];
    const amrex::Real wc =
        amrex::Real(0.5) * (w(i,j,k) + w(i,j,k+1)) - wall.wall_velocity[2];
    const amrex::Real normal = (dir == 0) ? uc : ((dir == 1) ? vc : wc);
    const amrex::Real t0 = uc - ((dir == 0) ? normal : amrex::Real(0.0));
    const amrex::Real t1 = vc - ((dir == 1) ? normal : amrex::Real(0.0));
    const amrex::Real t2 = wc - ((dir == 2) ? normal : amrex::Real(0.0));
    return std::sqrt(t0*t0 + t1*t1 + t2*t2);
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
        if (wall.heat.model == ScalarModel::BulkAero) {
            result.rhoTheta_in = bulk_theta_flux_in(
                sample.rho,
                wall.heat.coefficient,
                sample.U_t,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.rdOcp);
        } else {
            result.rhoTheta_in = resolved_theta_flux_in(
                sample.rho,
                sample.alpha_T,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.dx_inv,
                sample.rdOcp);
        }
        return result;

    case ScalarChannel::Vapor:
        if (wall.moisture == MoistureMode::DryImpermeable) {
            return result;
        }
        if (wall.vapor.model == ScalarModel::BulkAero) {
            result.rhoQv_in = bulk_vapor_flux_in(
                sample.rho,
                wall.vapor.coefficient,
                sample.U_t,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse);
        } else {
            result.rhoQv_in = resolved_vapor_flux_in(
                sample.rho,
                sample.alpha_C,
                wall.thermal.temperature_K,
                sample.scalar_air,
                sample.p_hse,
                sample.dx_inv);
        }
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
                 const amrex::Real rdOcp)
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
            rdOcp
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
                   const amrex::Real rdOcp)
{
    apply_face<0,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<0,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<1,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<1,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<2,false>(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
    apply_face<2,true >(bx,domain,quantity,flux_comp,state,prim,base_state,u,v,w,rhs,xflux,yflux,zflux,dx_inv,walls,alpha_T,alpha_C,rdOcp);
}

inline bool has_bulk_scalar_wall (
    const erf_wall_thermodynamics::Boundary& walls) noexcept
{
    for (const auto& wall : walls) {
        if (wall.heat.model == erf_wall_thermodynamics::ScalarModel::BulkAero ||
            wall.vapor.model == erf_wall_thermodynamics::ScalarModel::BulkAero) {
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
