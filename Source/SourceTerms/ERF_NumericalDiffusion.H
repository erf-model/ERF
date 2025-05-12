#ifndef ERF_NUMERICALDIFFUSION_H_
#define ERF_NUMERICALDIFFUSION_H_

#include <AMReX.H>
#include <ERF_DataStruct.H>
#include <AMReX_MultiFab.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calc_fifth_order_deriv (const amrex::Real& dnp2,
                        const amrex::Real& dnp1,
                        const amrex::Real& dn,
                        const amrex::Real& dnm1,
                        const amrex::Real& dnm2,
                        const amrex::Real& dnm3)
{
    amrex::Real interp = 10. * (dn   - dnm1)
                        - 5. * (dnp1 - dnm2)
                             + (dnp2 - dnm3);
    return interp;
}

void NumericalDiffusion_Scal (const amrex::Box& bx,
                              const int start_comp,
                              const int num_comp,
                              const amrex::Real dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);

void NumericalDiffusion_Xmom (const amrex::Box& bx,
                              const amrex::Real dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);

void NumericalDiffusion_Ymom (const amrex::Box& bx,
                              const amrex::Real dt,
                              const amrex::Real num_diff_coeff,
                              const amrex::Array4<const amrex::Real>& prim_data,
                              const amrex::Array4<const amrex::Real>& cell_data,
                              const amrex::Array4<      amrex::Real>& rhs,
                              const amrex::Array4<const amrex::Real>& mfx,
                              const amrex::Array4<const amrex::Real>& mfy);
#endif
