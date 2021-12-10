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
    amrex::Real TKE_production = 0.;
    amrex::Real sig11 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv,
                                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::x, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::x, DiffusionDir::x, K_LES, solverChoice, sig11);

    amrex::Real sig12 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv,
                                          use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::y, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::x, DiffusionDir::y, K_LES, solverChoice, sig12);

    amrex::Real sig13 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::x, DiffusionDir::z, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::x, DiffusionDir::z, K_LES, solverChoice, sig13);

    amrex::Real sig21 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::x, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::x, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::y, DiffusionDir::x, K_LES, solverChoice, sig21);

    amrex::Real sig22 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::y, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::y, DiffusionDir::y, K_LES, solverChoice, sig22);

    amrex::Real sig23 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::y, DiffusionDir::z, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::y, DiffusionDir::z, K_LES, solverChoice, sig23);

    amrex::Real sig31 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::x, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::x, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::z, DiffusionDir::x, K_LES, solverChoice, sig31);

    amrex::Real sig32 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::y, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                        ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::y, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::z, DiffusionDir::y, K_LES, solverChoice, sig32);

    amrex::Real sig33 = (ComputeStrainRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv,
                                           use_no_slip_stencil_lo, use_no_slip_stencil_hi) -
                         ComputeExpansionRate(i, j, k, u, v, w, MomentumEqn::z, DiffusionDir::z, cellSizeInv));
    TKE_production += ComputeStressGivenRates(i, j, k, MomentumEqn::z, DiffusionDir::z, K_LES, solverChoice, sig33);

    return TKE_production;
}
#endif
