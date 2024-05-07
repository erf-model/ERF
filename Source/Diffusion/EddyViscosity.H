/** \file EddyViscosity.H */

#ifndef _EDDY_VISCOSITY_H_
#define _EDDY_VISCOSITY_H_

#include <ABLMost.H>
#include <DataStruct.H>
#include <AMReX_BCRec.H>

void
ComputeTurbulentViscosity (const amrex::MultiFab& xvel , const amrex::MultiFab& yvel ,
                           const amrex::MultiFab& Tau11, const amrex::MultiFab& Tau22, const amrex::MultiFab& Tau33,
                           const amrex::MultiFab& Tau12, const amrex::MultiFab& Tau13, const amrex::MultiFab& Tau23,
                           const amrex::MultiFab& cons_in,
                           amrex::MultiFab& eddyViscosity,
                           amrex::MultiFab& Hfx1, amrex::MultiFab& Hfx2, amrex::MultiFab& Hfx3, amrex::MultiFab& Diss,
                           const amrex::Geometry& geom,
                           const amrex::MultiFab& mapfac_u, const amrex::MultiFab& mapfac_v,
                           const std::unique_ptr<amrex::MultiFab>& z_phys_nd,
                           const TurbChoice& turbChoice, const amrex::Real const_grav,
                           std::unique_ptr<ABLMost>& most,
                           const bool& exp_most,
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
               const amrex::Array4<amrex::Real const>& tau23,
               const int& klo,
               const bool& use_most,
               const bool& exp_most)
{
    amrex::Real s11bar = tau11(i,j,k);
    amrex::Real s22bar = tau22(i,j,k);
    amrex::Real s33bar = tau33(i,j,k);
    amrex::Real s12bar = 0.25 * ( tau12(i  , j  , k  ) + tau12(i  , j+1, k  )
                                + tau12(i+1, j  , k  ) + tau12(i+1, j+1, k  ) );

    // NOTES:
    // - If ERF_EXPLICIT_MOST_STRESS is not used, then we do not use the
    //   strains lying on the bottom boundary with MOST. These values are
    //   corrupted with the reflect odd BC used to fill the ghost cells.
    // - Neglecting the strains on the bottom boundary and using a one-
    //   sided average implies that the strains are equal at klo and klo+1.
    //   du/dz and dv/dz can be assumed to be equal through the first cell.
    //   However, whereas dw/dx = dw/dy = 0 on the surface, this is not true at
    //   klo+1. Therefore, the strains are not equal in general.
    // - If ERF_EXPLICIT_MOST_STRESS _is_ used, then strains may be correctly
    //   evaluated at the surface and no assumptions are needed.
    amrex::Real s13bar;
    if (use_most && k==klo && !exp_most) {
        s13bar = 0.5  * ( tau13(i  , j  , k+1) + tau13(i+1, j  , k+1) );
    }
    else
    {
        s13bar = 0.25 * ( tau13(i  , j  , k  ) + tau13(i  , j  , k+1)
                        + tau13(i+1, j  , k  ) + tau13(i+1, j  , k+1) );
    }

    amrex::Real s23bar;
    if (use_most && k==klo && !exp_most) {
        s23bar = 0.5  * ( tau23(i  , j  , k+1) + tau23(i  , j+1, k+1) );
    }
    else
    {
        s23bar = 0.25 * ( tau23(i  , j  , k  ) + tau23(i  , j  , k+1)
                        + tau23(i  , j+1, k  ) + tau23(i  , j+1, k+1) );
    }

    amrex::Real SmnSmn = s11bar*s11bar + s22bar*s22bar + s33bar*s33bar
                       + 2.0*s12bar*s12bar + 2.0*s13bar*s13bar + 2.0*s23bar*s23bar;

    return SmnSmn;
}
#endif
