#ifndef ERF_EBADVECTION_H_
#define ERF_EBADVECTION_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_YAFluxRegister.H>

#include <ERF_DataStruct.H>
#include <ERF_IndexDefines.H>
#include <ERF_ABLMost.H>

/** Compute advection tendency for all scalars other than density and potential temperature */
void EBAdvectionSrcForScalars (const amrex::Box& bx,
                            const int icomp, const int ncomp,
                            const amrex::Array4<const amrex::Real>& avg_xmom,
                            const amrex::Array4<const amrex::Real>& avg_ymom,
                            const amrex::Array4<const amrex::Real>& avg_zmom,
                            const amrex::Array4<const amrex::Real>& cell_prim,
                            const amrex::Array4<amrex::Real>& src,
                            const amrex::Array4<const amrex::EBCellFlag>& cfg_arr,
                            const amrex::Array4<const amrex::Real>& ax_arr,
                            const amrex::Array4<const amrex::Real>& ay_arr,
                            const amrex::Array4<const amrex::Real>& az_arr,
                            const amrex::Array4<const amrex::Real>& vf_arr,
                            const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& cellSizeInv,
                            const amrex::Array4<const amrex::Real>& mf_m,
                            const AdvType horiz_adv_type, const AdvType vert_adv_type,
                            const amrex::Real horiz_upw_frac, const amrex::Real vert_upw_frac,
                            const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx_arr,
                            const amrex::Box& domain,
                            const amrex::BCRec* bc_ptr_h);

#endif
