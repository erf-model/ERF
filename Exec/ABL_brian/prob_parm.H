#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real Theta_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real A_0;

extern AMREX_GPU_DEVICE_MANAGED amrex::Real U0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real V0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real W0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real U0_Pert_Mag;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real V0_Pert_Mag;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real W0_Pert_Mag;

} // namespace ProbParm

#endif
