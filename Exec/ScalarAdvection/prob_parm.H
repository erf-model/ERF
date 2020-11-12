#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rho_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real rhoe_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real T_0;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real A_0;
} // namespace ProbParm

#endif
