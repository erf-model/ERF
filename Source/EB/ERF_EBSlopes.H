#ifndef ERF_EBSLOPES_H_
#define ERF_EBSLOPES_H_

#include <AMReX_EBCellFlag.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real,AMREX_SPACEDIM>
erf_calc_slopes_eb (int i, int j, int k, 
                    const amrex::Real ccent_x,
                    const amrex::Real ccent_y,
                    const amrex::Real ccent_z,
                    amrex::Array4<amrex::Real const> const& state,
                    amrex::Array4<amrex::Real const> const& ccent,
                    amrex::Array4<amrex::EBCellFlag const> const& flag);

#endif