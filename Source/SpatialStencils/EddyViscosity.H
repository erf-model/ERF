/** \file EddyViscosity.H */

#ifndef _EDDY_VISCOSITY_H_
#define _EDDY_VISCOSITY_H_

#include <ABLMost.H>
#include <DataStruct.H>
#include <StrainRate.H>

void ComputeTurbulentViscosity(const amrex::MultiFab& xvel, const amrex::MultiFab& yvel, const amrex::MultiFab& zvel,
                               const amrex::MultiFab& cons_in, amrex::MultiFab& eddyViscosity,
                               const amrex::Geometry& geom,
                               const SolverChoice& solverChoice,
                               std::unique_ptr<ABLMost>& most,
                               const amrex::Gpu::DeviceVector<amrex::BCRec> domain_bcs_type_d,
                               bool vert_only = false);

AMREX_GPU_DEVICE
inline
amrex::Real ComputeSmnSmnWithTerrain(const int& i,const int& j,const int& k,
                          const amrex::Array4<amrex::Real const>& u,
                          const amrex::Array4<amrex::Real const>& v,
                          const amrex::Array4<amrex::Real const>& w,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                          const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
  amrex::Real met_h_xi,met_h_eta,met_h_zeta;

  // TO-DO :
  //         1. Calls to TerrainMetrics.H functions to get metrics
  //         2. Pass metrics at each ijk
  //         3. Pass z_nd to ComputeSmnSmn
  //         4. Anisotropy of stress tensor
  met_h_xi   = 0.0;
  met_h_eta  = 0.0;
  met_h_zeta = 1.0;

  amrex::Real S11 = ComputeS11WithTerrain(i+1, j, k, u, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta);
  amrex::Real S22 = ComputeS22WithTerrain(i, j+1, k, v, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta);
  amrex::Real S33 = ComputeS33WithTerrain(i, j, k+1, w, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta);

    amrex::Real S12 = 0.25* (
              ComputeS12WithTerrain(i  , j  , k, u, v, w, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS12WithTerrain(i  , j+1, k, u, v, w, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS12WithTerrain(i+1, j  , k, u, v, w, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS12WithTerrain(i+1, j+1, k, u, v, w, cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            );

    amrex::Real S13 = 0.25* (
              ComputeS13WithTerrain(i  , j, k  , u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS13WithTerrain(i  , j, k+1, u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS13WithTerrain(i+1, j, k  , u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS13WithTerrain(i+1, j, k+1, u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            );

    amrex::Real S23 = 0.25* (
              ComputeS23WithTerrain(i, j  , k  , u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS23WithTerrain(i, j  , k+1, u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS23WithTerrain(i, j+1, k  , u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            + ComputeS23WithTerrain(i, j+1, k+1, u, v, w,  cellSizeInv, met_h_xi, met_h_eta, met_h_zeta, domain, bc_ptr)
            );

    amrex::Real SmnSmn = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;

    return SmnSmn;
}

AMREX_GPU_DEVICE
inline
amrex::Real ComputeSmnSmn(const int& i,const int& j,const int& k,
                          const amrex::Array4<amrex::Real const>& u,
                          const amrex::Array4<amrex::Real const>& v,
                          const amrex::Array4<amrex::Real const>& w,
                          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                          const amrex::Box& domain, const amrex::BCRec* bc_ptr)
{
    amrex::Real S11 = ComputeS11(i+1, j, k, u, cellSizeInv);
    amrex::Real S22 = ComputeS22(i, j+1, k, v, cellSizeInv);
    amrex::Real S33 = ComputeS33(i, j, k+1, w, cellSizeInv);

    amrex::Real S12 = 0.25* (
              ComputeS12(i  , j  , k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i  , j+1, k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i+1, j  , k, u, v, w, cellSizeInv, domain, bc_ptr)
            + ComputeS12(i+1, j+1, k, u, v, w, cellSizeInv, domain, bc_ptr)
            );

    amrex::Real S13 = 0.25* (
              ComputeS13(i  , j, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i  , j, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i+1, j, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS13(i+1, j, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            );

    amrex::Real S23 = 0.25* (
              ComputeS23(i, j  , k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j  , k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j+1, k  , u, v, w,  cellSizeInv, domain, bc_ptr)
            + ComputeS23(i, j+1, k+1, u, v, w,  cellSizeInv, domain, bc_ptr)
            );

    amrex::Real SmnSmn = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;

    return SmnSmn;
}

AMREX_GPU_DEVICE
inline amrex::Real
ComputeTKEProduction (const int &i, const int &j, const int &k,
                      const amrex::Array4<const amrex::Real>& u,
                      const amrex::Array4<const amrex::Real>& v,
                      const amrex::Array4<const amrex::Real>& w,
                      const amrex::Array4<amrex::Real>& K_turb,
                      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                      const amrex::Box& domain, const amrex::BCRec* bc_ptr, const int& use_terrain)
{
    amrex::Real TKE_production;
    if (use_terrain) 
        TKE_production = K_turb(i,j,k,EddyDiff::Mom_h) * ComputeSmnSmnWithTerrain(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);
    else
        TKE_production = K_turb(i,j,k,EddyDiff::Mom_h) * ComputeSmnSmn(i,j,k,u,v,w,cellSizeInv,domain,bc_ptr);

    return TKE_production;
}

/// Compute K_turb (i-1/2, j+1/2, k) etc given K_turb (i, j, k) is known
// Note: This should be at edges for momEqnDir != diffDir, cell centers otherwise
AMREX_GPU_DEVICE
inline amrex::Real
InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                              const enum MomentumEqn &momentumEqn,
                              const enum DiffusionDir &diffDir,
                              const amrex::Array4<amrex::Real>& K_turb) {
  // Assuming we already have 'K_turb' computed for all (i, j, k)
  amrex::Real turbViscInterpolated = 1.0;

  switch (momentumEqn) {
  case MomentumEqn::x: // Reference face is x-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = K_turb(i-1, j, k, EddyDiff::Mom_h);
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_turb(i-1, j, k, EddyDiff::Mom_h) + K_turb(i, j, k, EddyDiff::Mom_h)
                    + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i, j-1, k, EddyDiff::Mom_h) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_turb(i-1, j, k, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::y: // Reference face is y-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_turb(i, j-1, k, EddyDiff::Mom_h) + K_turb(i, j, k, EddyDiff::Mom_h)
                    + K_turb(i-1, j-1, k, EddyDiff::Mom_h) + K_turb(i-1, j, k, EddyDiff::Mom_h) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = K_turb(i, j-1, k, EddyDiff::Mom_h);
      break;
    case DiffusionDir::z:
      turbViscInterpolated = 0.25*( K_turb(i, j-1, k, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j, k-1, EddyDiff::Mom_v) );
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  case MomentumEqn::z: // Reference face is z-face index (i, j, k)
    switch (diffDir) {
    case DiffusionDir::x:
      turbViscInterpolated = 0.25*( K_turb(i, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i-1, j, k-1, EddyDiff::Mom_v) + K_turb(i-1, j, k, EddyDiff::Mom_v) );
      break;
    case DiffusionDir::y:
      turbViscInterpolated = 0.25*( K_turb(i, j, k-1, EddyDiff::Mom_v) + K_turb(i, j, k, EddyDiff::Mom_v)
                    + K_turb(i, j-1, k-1, EddyDiff::Mom_v) + K_turb(i, j-1, k, EddyDiff::Mom_v) );
      break;
    case DiffusionDir::z:
      turbViscInterpolated = K_turb(i, j, k-1, EddyDiff::Mom_v);
      break;
    default:
      amrex::Abort("Error: Diffusion direction is unrecognized");
    }
    break;
  default:
    amrex::Abort("Error: Momentum equation is unrecognized");
  }

  return turbViscInterpolated;
}

#endif
