#ifndef _TRANSPORT_PARAMS_H_
#define _TRANSPORT_PARAMS_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

namespace transport_params {

struct transport_param_values {
  amrex::Real const_viscosity = 0.0;
  amrex::Real const_bulk_viscosity = 0.0;
  amrex::Real const_diffusivity = 0.0;
  amrex::Real const_conductivity = 0.0;
};

extern transport_param_values trparms;

void init();
void finalize();

} // namespace transport_params

#endif
