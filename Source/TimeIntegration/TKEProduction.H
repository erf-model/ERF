#ifndef _COMPUTE_TKEPRODUCTION_H_
#define _COMPUTE_TKEPRODUCTION_H_

#include <StressTerm.H>
#include <StrainRate.H>

AMREX_GPU_DEVICE
inline Real
ComputeTKEProduction (const int &i, const int &j, const int &k,
                      const Array4<Real>& u, const Array4<Real>& v, const Array4<Real>& w,
                      const GpuArray<Real, AMREX_SPACEDIM>& cellSizeInv,
                      const Array4<Real>& K_LES,
                      const SolverChoice &solverChoice,
                      bool use_no_slip_stencil_lo,
                      bool use_no_slip_stencil_hi)
{
    amrex::Real TKE_production =
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi);

    TKE_production +=
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::x, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::x, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::x, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi);

    TKE_production +=
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::x, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::x, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::x, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::y, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::y, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::y, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi) *
       (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv,
                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv)) *
        ComputeStressTerm(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv,
                          K_LES, solverChoice, use_no_slip_stencil_lo, use_no_slip_stencil_hi);

    return TKE_production;
}
#endif
