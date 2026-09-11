#ifndef ERF_CLOUD_CHAMBER_WALL_STRESS_H_
#define ERF_CLOUD_CHAMBER_WALL_STRESS_H_

#include <AMReX_Array4.H>
#include <AMReX_Box.H>
#include <AMReX_GpuContainers.H>

#include "Diffusion/ERF_CloudChamberWallFlux.H"
#include "ERF_IndexDefines.H"

namespace erf_cloud_chamber_wall_stress {

using amrex::Array4;
using amrex::Box;
using amrex::GpuArray;
using amrex::Real;

// DiffusionSrcForMom stores the physical wall traction with the low-side
// coordinate orientation and the opposite of it on the high side.  With the
// operator's outer minus sign, both sides then produce drag for U_t > 0.
template <int DIR, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real stored_tau_from_physical_traction (Real traction) noexcept
{
    amrex::ignore_unused(DIR);
    return HIGH ? -traction : traction;
}

template <int DIR, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
GpuArray<Real, AMREX_SPACEDIM>
wall_traction_at_cell (
    int i, int j, int k,
    const Array4<const Real>& state,
    const Array4<const Real>& base_state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const erf_wall_thermodynamics::FaceWall& wall,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const bool cloudy) noexcept
{
    const auto u_t = erf_cloud_chamber_wall_flux::
        relative_tangential_velocity_cell_centered(DIR, i, j, k, u, v, w, wall);
    const Real U_t = erf_cloud_chamber_wall_flux::tangential_magnitude(u_t);
    const erf_cloud_chamber_wall_flux::MomentumWallSample sample{
        state(i,j,k,Rho_comp), u_t, U_t, Real(0.5) / dx_inv[DIR],
        state(i,j,k,RhoTheta_comp) / state(i,j,k,Rho_comp),
        cloudy ? state(i,j,k,RhoQ1_comp) / state(i,j,k,Rho_comp) : Real(0.0),
        base_state(i,j,k,BaseState::p0_comp),
        DIR == 2 ? (HIGH ? -1 : 1) : 0};
    return erf_cloud_chamber_wall_flux::evaluate_momentum_traction(wall, sample).
        traction_on_fluid;
}

// The closure is evaluated at wall-adjacent cell centers.  Each boundary
// stress node is the midpoint average of its two tangent-adjacent cell-center
// tractions; the normal coordinate is the physical boundary cell.
template <int DIR, int TANGENT, bool HIGH>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Real averaged_physical_traction_component (
    int i, int j, int k,
    const Box& domain,
    const Array4<const Real>& state,
    const Array4<const Real>& base_state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const erf_wall_thermodynamics::FaceWall& wall,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const bool cloudy) noexcept
{
    int cell_lo[AMREX_SPACEDIM] = {i, j, k};
    int cell_hi[AMREX_SPACEDIM] = {i, j, k};
    const int node = cell_lo[TANGENT];
    cell_lo[TANGENT] = amrex::max(domain.smallEnd(TANGENT), node - 1);
    cell_hi[TANGENT] = amrex::min(domain.bigEnd(TANGENT), node);
    cell_lo[DIR] = HIGH ? domain.bigEnd(DIR) : domain.smallEnd(DIR);
    cell_hi[DIR] = cell_lo[DIR];

    const auto traction_lo = wall_traction_at_cell<DIR,HIGH>(
        cell_lo[0], cell_lo[1], cell_lo[2], state, base_state, u, v, w, wall,
        dx_inv, cloudy);
    const auto traction_hi = wall_traction_at_cell<DIR,HIGH>(
        cell_hi[0], cell_hi[1], cell_hi[2], state, base_state, u, v, w, wall,
        dx_inv, cloudy);
    return Real(0.5) * (traction_lo[TANGENT] + traction_hi[TANGENT]);
}

template <int DIR, int TANGENT, bool HIGH>
AMREX_FORCE_INLINE
void apply_tangent_stress (
    const Box& bx,
    const Box& domain,
    const Array4<const Real>& state,
    const Array4<const Real>& base_state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const Array4<Real>& tau,
    const erf_wall_thermodynamics::FaceWall& wall,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const bool cloudy)
{
    Box target = amrex::surroundingNodes(amrex::surroundingNodes(bx, DIR), TANGENT);
    const int face = HIGH ? domain.bigEnd(DIR) + 1 : domain.smallEnd(DIR);
    target.setSmall(DIR, face);
    target.setBig(DIR, face);
    target.setSmall(TANGENT,
                    amrex::max(target.smallEnd(TANGENT), domain.smallEnd(TANGENT)));
    target.setBig(TANGENT,
                  amrex::min(target.bigEnd(TANGENT), domain.bigEnd(TANGENT) + 1));
    if (target.isEmpty()) { return; }

    amrex::ParallelFor(target,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        const Real traction = averaged_physical_traction_component<DIR,TANGENT,HIGH>(
            i, j, k, domain, state, base_state, u, v, w, wall, dx_inv, cloudy);
        tau(i,j,k) = stored_tau_from_physical_traction<DIR,HIGH>(traction);
    });
}

template <int DIR, bool HIGH>
AMREX_FORCE_INLINE
void apply_face (
    const Box& bx,
    const Box& domain,
    const Array4<const Real>& state,
    const Array4<const Real>& base_state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const Array4<Real>& tau12,
    const Array4<Real>& tau13,
    const Array4<Real>& tau23,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const erf_wall_thermodynamics::Boundary& walls,
    const bool cloudy)
{
    if constexpr (HIGH) {
        if (bx.bigEnd(DIR) != domain.bigEnd(DIR)) { return; }
    } else {
        if (bx.smallEnd(DIR) != domain.smallEnd(DIR)) { return; }
    }

    const auto wall = walls[2*DIR + (HIGH ? 1 : 0)];
    if (!erf_cloud_chamber_wall_flux::requires_tangential_momentum(wall)) {
        return;
    }

    if constexpr (DIR == 0) {
        apply_tangent_stress<DIR,1,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau12, wall, dx_inv, cloudy);
        apply_tangent_stress<DIR,2,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau13, wall, dx_inv, cloudy);
    } else if constexpr (DIR == 1) {
        apply_tangent_stress<DIR,0,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau12, wall, dx_inv, cloudy);
        apply_tangent_stress<DIR,2,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau23, wall, dx_inv, cloudy);
    } else {
        apply_tangent_stress<DIR,0,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau13, wall, dx_inv, cloudy);
        apply_tangent_stress<DIR,1,HIGH>(bx, domain, state, base_state, u, v, w,
                                         tau23, wall, dx_inv, cloudy);
    }
}

inline void apply (
    const Box& bx,
    const Box& domain,
    const Array4<const Real>& state,
    const Array4<const Real>& base_state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const Array4<Real>& tau12,
    const Array4<Real>& tau13,
    const Array4<Real>& tau23,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const erf_wall_thermodynamics::Boundary& walls,
    const bool cloudy)
{
    apply_face<0,false>(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
    apply_face<0,true >(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
    apply_face<1,false>(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
    apply_face<1,true >(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
    apply_face<2,false>(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
    apply_face<2,true >(bx, domain, state, base_state, u, v, w, tau12, tau13, tau23, dx_inv, walls, cloudy);
}

inline void apply (
    const Box& bx,
    const Box& domain,
    const Array4<const Real>& state,
    const Array4<const Real>& u,
    const Array4<const Real>& v,
    const Array4<const Real>& w,
    const Array4<Real>& tau12,
    const Array4<Real>& tau13,
    const Array4<Real>& tau23,
    const GpuArray<Real, AMREX_SPACEDIM>& dx_inv,
    const erf_wall_thermodynamics::Boundary& walls)
{
    apply(bx, domain, state, Array4<const Real>{}, u, v, w,
          tau12, tau13, tau23, dx_inv, walls, false);
}

inline bool has_momentum_wall (
    const erf_wall_thermodynamics::Boundary& walls) noexcept
{
    for (const auto& wall : walls) {
        if (erf_cloud_chamber_wall_flux::requires_tangential_momentum(wall)) {
            return true;
        }
    }
    return false;
}

} // namespace erf_cloud_chamber_wall_stress

#endif
