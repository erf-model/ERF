/** \file ERF_EddyViscosity.H */

#ifndef ERF_EDDY_VISCOSITY_H_
#define ERF_EDDY_VISCOSITY_H_

#include "AMReX_BCRec.H"
#include "ERF_SurfaceLayer.H"
#include "ERF_DataStruct.H"
#include "ERF_IndexDefines.H"
#include "ERF_Constants.H"
#include "ERF_EOS.H"

#ifdef ERF_USE_SHOC
#include "ERF_ShocInterface.H"
#endif

void
ComputeTurbulentViscosity (amrex::Real dt,
                           const amrex::MultiFab& xvel ,
                           const amrex::MultiFab& yvel ,
                           amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                           amrex::MultiFab& cons_in,
                           const amrex::MultiFab& wdist,
                           amrex::MultiFab& eddyViscosity,
                           amrex::MultiFab& Hfx1,
                           amrex::MultiFab& Hfx2,
                           amrex::MultiFab& Hfx3,
                           amrex::MultiFab& Diss,
                           const amrex::Geometry& geom,
                           amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                           const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                           const SolverChoice& solverChoice,
                           std::unique_ptr<SurfaceLayer>& SurfLayer,
                           const amrex::MultiFab* z_0,
                           const bool& use_terrain_fitted_coords,
                           const bool& use_moisture,
                           int level,
                           const amrex::BCRec* bc_ptr,
                           bool vert_only = false);

void
ComputeTurbulentViscosityLES (amrex::Vector<std::unique_ptr<amrex::MultiFab>>& Tau_lev,
                              const amrex::MultiFab& cons_in,
                              amrex::MultiFab& eddyViscosity,
                              amrex::MultiFab& Hfx1,
                              amrex::MultiFab& Hfx2,
                              amrex::MultiFab& Hfx3,
                              amrex::MultiFab& Diss,
                              const amrex::Geometry& geom,
                              bool use_terrain_fitted_coords,
                              amrex::Vector<std::unique_ptr<amrex::MultiFab>>& mapfac,
                              const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                              const TurbChoice& turbChoice,
                              const amrex::Real const_grav,
                              std::unique_ptr<SurfaceLayer>& SurfLayer,
                              const MoistureComponentIndices& moisture_indices,
                              const amrex::MultiFab* xvel = nullptr,
                              const amrex::MultiFab* yvel = nullptr);

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
