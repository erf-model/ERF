#ifndef _PROB_PARM_H_
#define _PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>

namespace ProbParm {
extern AMREX_GPU_DEVICE_MANAGED amrex::Real M_inf;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real alpha;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real gamma;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real beta;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real sigma;
extern AMREX_GPU_DEVICE_MANAGED amrex::Real R;
} // namespace ProbParm

#endif
