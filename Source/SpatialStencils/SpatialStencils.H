#ifndef _STENCILS_H_
#define _STENCILS_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_InterpFaceRegister.H>
#include "DataStruct.H"
#include "IndexDefines.H"

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real OmegaFromW(int i, int j, int k, amrex::Real w,
                       const amrex::Array4<const amrex::Real> u,
                       const amrex::Array4<const amrex::Real> v,
                       const amrex::Array4<const amrex::Real> z_nd,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at z-face (i,j,k-1/2)
    amrex::Real met_zlo_xi   = 0.5 * dxInv[0] *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i+1,j  ,k  )    // hi i, lo k
                               -z_nd(i  ,j+1,k  ) - z_nd(i  ,j  ,k  ) ); // lo i, lo k

    // This is dh/deta at z-face (i,j,k-1/2)
    amrex::Real met_zlo_eta  = 0.5 * dxInv[1] *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k  )    // hi j, lo k
                               -z_nd(i+1,j  ,k  ) - z_nd(i  ,j  ,k  ) ); // lo j, lo k

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real vec_zlo_xi   = (k == 0) ? 1.5 * (0.5*(u(i,j,k)+u(i+1,j,k))) - 0.5*(0.5*(u(i,j,k+1)+u(i+1,j,k+1))) :
        0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
    amrex::Real vec_zlo_eta  = (k == 0) ? 1.5 * (0.5*( v(i,j,k)+v(i,j+1,k))) - 0.5*(0.5*(v(i,j,k+1)+v(i,j+1,k+1))) :
        0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );

    amrex::Real omega = -met_zlo_xi * vec_zlo_xi - met_zlo_eta * vec_zlo_eta + w;
    return omega;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real WFromOmega(int i, int j, int k, amrex::Real omega,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // This is dh/dxi at z-face (i,j,k-1/2)
    amrex::Real met_zlo_xi   = 0.5 * dxInv[0] *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i+1,j  ,k  )    // hi i, lo k
                               -z_nd(i  ,j+1,k  ) - z_nd(i  ,j  ,k  ) ); // lo i, lo k

    // This is dh/deta at z-face (i,j,k-1/2)
    amrex::Real met_zlo_eta  = 0.5 * dxInv[1] *
                              ( z_nd(i+1,j+1,k  ) + z_nd(i  ,j+1,k  )    // hi j, lo k
                               -z_nd(i+1,j  ,k  ) - z_nd(i  ,j  ,k  ) ); // lo j, lo k

    // Use extrapolation instead of interpolation if at the bottom boundary
    amrex::Real vec_zlo_xi   = (k == 0) ? 1.5 * (0.5*(u(i,j,k)+u(i+1,j,k))) - 0.5*(0.5*(u(i,j,k+1)+u(i+1,j,k+1))) :
        0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
    amrex::Real vec_zlo_eta  = (k == 0) ? 1.5 * (0.5*( v(i,j,k)+v(i,j+1,k))) - 0.5*(0.5*(v(i,j,k+1)+v(i,j+1,k+1))) :
        0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );

    amrex::Real w = met_zlo_xi * vec_zlo_xi + met_zlo_eta * vec_zlo_eta + omega;
    return w;
}

// This routine follows WFromOmega with clipping at the boundary:
// Terrain metrics at the boundary are used for ghost cells (0 gradient).
// U & V are extrapolated to the top & bottom k faces.
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::Real WFromOmegaBC(int i, int j, int k, amrex::Real omega,
                         const amrex::Array4<const amrex::Real>& u,
                         const amrex::Array4<const amrex::Real>& v,
                         const amrex::Array4<const amrex::Real>& z_nd,
                         const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv)
{
    // Clip metric indices for ghost cell BCs
    amrex::Dim3 hi = ubound(z_nd);
    amrex::Dim3 lo = lbound(z_nd);
    int  im = (i  >= hi.x) ? (hi.x-1) : i;
         im = (im <= lo.x) ? (lo.x+1) : im;
    int  jm = (j  >= hi.y) ? (hi.y-1) : j;
         jm = (jm <= lo.y) ? (lo.y+1) : jm;
    int  km = (k  >= hi.z) ? (hi.z-1) : k;
         km = (km <= lo.z) ? (lo.z+1) : km;

    // This is dh/dxi at z-face (i,j,k-1/2)
    amrex::Real met_zlo_xi   = 0.5 * dxInv[0] *
                              ( z_nd(im+1,jm+1,km  ) + z_nd(im+1,jm  ,km  )    // hi i, lo k
                               -z_nd(im  ,jm+1,km  ) - z_nd(im  ,jm  ,km  ) ); // lo i, lo k

    // This is dh/deta at z-face (i,j,k-1/2)
    amrex::Real met_zlo_eta  = 0.5 * dxInv[1] *
                              ( z_nd(im+1,jm+1,km  ) + z_nd(im  ,jm+1,km  )    // hi j, lo k
                               -z_nd(im+1,jm  ,km  ) - z_nd(im  ,jm  ,km  ) ); // lo j, lo k

    // Clip velocity indices for ghost cell BCs
    amrex::Real vec_zlo_xi;
    amrex::Real vec_zlo_eta;
    hi = ubound(u);
    lo = lbound(u);
    hi.z = hi.z + 1; // nodal for zvelocity

    if( k == hi.z) {
      // foextrap to highest k-face
      vec_zlo_xi  = 0.75*( u(i,j,k-1)+u(i+1,j,k-1) ) - 0.25*( u(i,j,k-2)+u(i+1,j,k-2) );
      vec_zlo_eta = 0.75*( v(i,j,k-1)+v(i,j+1,k-1) ) - 0.25*( v(i,j,k-2)+v(i,j+1,k-2) );
    } else if (k == lo.z) {
      // foextrap to lowest k-face
      vec_zlo_xi  = 0.75*( u(i,j,k)+u(i+1,j,k) ) - 0.25*( u(i,j,k+1)+u(i+1,j,k+1) );
      vec_zlo_eta = 0.75*( v(i,j,k)+v(i,j+1,k) ) - 0.25*( v(i,j,k+1)+v(i,j+1,k+1) );
    } else {
      // Valid ghost cell for u & v
      vec_zlo_xi  = 0.25 * ( u(i,j,k-1) + u(i+1,j,k-1) + u(i,j,k) + u(i+1,j,k) );
      vec_zlo_eta = 0.25 * ( v(i,j,k-1) + v(i,j+1,k-1) + v(i,j,k) + v(i,j+1,k) );
    }

    amrex::Real w = met_zlo_xi * vec_zlo_xi + met_zlo_eta * vec_zlo_eta + omega;
    return w;
}

void MomentumToVelocity(amrex::MultiFab& xvel_out,
                        amrex::MultiFab& yvel_out,
                        amrex::MultiFab& zvel_out,
                        const amrex::MultiFab& cons_in,
                        const amrex::MultiFab& xmom_in,
                        const amrex::MultiFab& ymom_in,
                        const amrex::MultiFab& zmom_in,
                        const amrex::IntVect& ngrow);

void VelocityToMomentum(const amrex::MultiFab& xvel_in, 
                        const amrex::MultiFab& yvel_in,
                        const amrex::MultiFab& zvel_in,
                        const amrex::MultiFab& cons_in,
                        amrex::MultiFab& xmom_out, 
                        amrex::MultiFab& ymom_out, 
                        amrex::MultiFab& zmom_out, 
                        const amrex::IntVect& ngrow);


AMREX_GPU_DEVICE
amrex::Real InterpolateFromCellOrFace(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& qty, const int & qty_index,
                       const amrex::Real& uadv,
                       const enum Coord& coordDir,
                       const int& spatial_order);

AMREX_GPU_DEVICE
amrex::Real InterpolatePertFromCell(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& qty, const int & qty_index,
                       const amrex::Real& uadv,
                       const enum Coord& coordDir,
                       const int& spatial_order,
#ifdef ERF_USE_TERRAIN
                       const Array4<const Real>& r0_arr);
#else
                       const amrex::Real* dptr_hse);
#endif

AMREX_GPU_DEVICE
amrex::Real InterpolateDensityPertFromCellToFace(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& cons_in,
                       const amrex::Real& uadv,
                       const enum Coord& coordDir,
                       const int& spatial_order,
#ifdef ERF_USE_TERRAIN
                       const amrex::Array4<const amrex::Real>& r0_arr);
#else
                       const amrex::Real* dptr_hse);
#endif

/** Meant for {x, y, z}- momentum equations */
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForXMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& u,
#ifdef ERF_USE_TERRAIN
                       const amrex::Array4<const amrex::Real>& z_nd, 
                       const amrex::Array4<const amrex::Real>& detJ,
#endif
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order);
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForYMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& v,
#ifdef ERF_USE_TERRAIN
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
#endif
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order);
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForZMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& w,
#ifdef ERF_USE_TERRAIN
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
#endif
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order);

/** Meant for the continuity, energy, and scalar equations */
AMREX_GPU_DEVICE
amrex::Real AdvectionSrcForState(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& rho_u,
                       const amrex::Array4<const amrex::Real>& rho_v,
                       const amrex::Array4<const amrex::Real>& rho_w,
                       const amrex::Array4<const amrex::Real>& cell_prim, const int &qty_index,
                       const amrex::Array4<amrex::Real>& flux_x,
                       const amrex::Array4<amrex::Real>& flux_y,
                       const amrex::Array4<amrex::Real>& flux_z,
#ifdef ERF_USE_TERRAIN
                       const amrex::Array4<const amrex::Real>& z_nd,
                       const amrex::Array4<const amrex::Real>& detJ,
#endif
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const int& spatial_order);

/** Compute K (i-1/2, j+1/2, k) etc given K(i, j, k) or nut (i, j, k) is known */
AMREX_GPU_DEVICE
amrex::Real InterpolateTurbulentViscosity(const int &i, const int &j, const int &k,
                            const enum MomentumEqn &momentumEqn,
                            const enum DiffusionDir &diffDir,
                            const amrex::Array4<amrex::Real>& Ksmag);

/** Compute tau_ij (m + 1/2), tau_ij (m - 1/2) where m = {i, j, k} for DNS or Smagorinsky */
AMREX_GPU_DEVICE
amrex::Real ComputeStressTerm (
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& w,
                       const enum MomentumEqn &momentumEqn,
                       const enum DiffusionDir &diffDir,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr);

AMREX_GPU_DEVICE
amrex::Real DiffusionSrcForMom(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& u,
                       const amrex::Array4<const amrex::Real>& v,
                       const amrex::Array4<const amrex::Real>& w,
                       const amrex::Array4<const amrex::Real>& cons,
                       const enum MomentumEqn &momentumEqn,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const amrex::Box& domain, const amrex::BCRec* bc_ptr);

AMREX_GPU_DEVICE
amrex::Real ComputeDiffusionFluxForState(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim,
                       const int & prim_index,
                       const amrex::Real invCellWidth,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice,
                       const enum Coord& coordDir);

AMREX_GPU_DEVICE
amrex::Real DiffusionSrcForState(
                       const int &i, const int &j, const int &k,
                       const amrex::Array4<const amrex::Real>& cell_data,
                       const amrex::Array4<const amrex::Real>& cell_prim, const int & qty_index,
                       const amrex::Array4<amrex::Real>& flux_x,
                       const amrex::Array4<amrex::Real>& flux_y,
                       const amrex::Array4<amrex::Real>& flux_z,
                       const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSize,
                       const amrex::Array4<amrex::Real>& Ksmag,
                       const SolverChoice &solverChoice);
#endif
