#ifndef _FORCING_H_
#define _FORCING_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include "IndexDefines.H"

struct forcing_params {
  amrex::Real u0 = 0.0;
  amrex::Real v0 = 0.0;
  amrex::Real w0 = 0.0;
  amrex::Real forcing = 0.0;
}; // namespace forcing_params

extern forcing_params fparms;

#endif
