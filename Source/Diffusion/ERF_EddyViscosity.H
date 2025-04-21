/** \file ERF_EddyViscosity.H */

#ifndef ERF_EDDY_VISCOSITY_H_
#define ERF_EDDY_VISCOSITY_H_

#include "AMReX_BCRec.H"

#include "ERF_SurfaceLayer.H"
#include "ERF_DataStruct.H"

void
ComputeTurbulentViscosity (const amrex::MultiFab& xvel ,
                           const amrex::MultiFab& yvel ,
                           const amrex::MultiFab& Tau11,
                           const amrex::MultiFab& Tau22,
                           const amrex::MultiFab& Tau33,
                           const amrex::MultiFab& Tau12,
                           const amrex::MultiFab& Tau13,
                           const amrex::MultiFab& Tau23,
                           const amrex::MultiFab& cons_in,
                           const amrex::MultiFab& wdist,
                           amrex::MultiFab& eddyViscosity,
                           amrex::MultiFab& Hfx1,
                           amrex::MultiFab& Hfx2,
                           amrex::MultiFab& Hfx3,
                           amrex::MultiFab& Diss,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& mapfac_u,
                           const amrex::MultiFab& mapfac_v,
                           const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                           const SolverChoice& solverChoice,
                           std::unique_ptr<SurfaceLayer>& SurfLayer,
                           const amrex::FArrayBox* z_0,
                           const bool& use_terrain_fitted_coords,
                           const bool& use_moisture,
                           int level,
                           const amrex::BCRec* bc_ptr,
                           bool vert_only = false);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeSmnSmn (int& i, int& j, int& k,
               const amrex::Array4<amrex::Real const>& tau11,
               const amrex::Array4<amrex::Real const>& tau22,
               const amrex::Array4<amrex::Real const>& tau33,
               const amrex::Array4<amrex::Real const>& tau12,
               const amrex::Array4<amrex::Real const>& tau13,
               const amrex::Array4<amrex::Real const>& tau23)
{
    amrex::Real s11bar = tau11(i,j,k);
    amrex::Real s22bar = tau22(i,j,k);
    amrex::Real s33bar = tau33(i,j,k);
    amrex::Real s12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                                + tau12(i+1, j  , k  ) + tau12(i+1, j+1, k  ) );
    amrex::Real s13bar = 0.25 * ( tau13(i  , j  , k  ) + tau13(i  , j  , k+1)
                                + tau13(i+1, j  , k  ) + tau13(i+1, j  , k+1) );
    amrex::Real s23bar = 0.25 * ( tau23(i  , j  , k  ) + tau23(i  , j  , k+1)
                                + tau23(i  , j+1, k  ) + tau23(i  , j+1, k+1) );

    amrex::Real SmnSmn = s11bar*s11bar + s22bar*s22bar + s33bar*s33bar
                       + 2.0*s12bar*s12bar + 2.0*s13bar*s13bar + 2.0*s23bar*s23bar;

    return SmnSmn;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ComputeSmnSmn2D (int& i, int& j, int& k,
                 const amrex::Array4<amrex::Real const>& tau11,
                 const amrex::Array4<amrex::Real const>& tau22,
                 const amrex::Array4<amrex::Real const>& tau12)
{
    amrex::Real sdiff = tau11(i,j,k) - tau22(i,j,k);
    amrex::Real s12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                                + tau12(i+1, j  , k  ) + tau12(i+1, j+1, k  ) );
    return 0.5 * sdiff*sdiff + 2.0*s12bar*s12bar;
}
#endif
