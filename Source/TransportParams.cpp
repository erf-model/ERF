#include "TransportParams.H"

namespace transport_params {

transport_param_values trparms;

void
init()
{
  amrex::ParmParse pp("transport");
  pp.query("const_viscosity", trparms.const_viscosity);
  pp.query("const_bulk_viscosity", trparms.const_bulk_viscosity);
  pp.query("const_conductivity", trparms.const_conductivity);
  pp.query("const_diffusivity", trparms.const_diffusivity);
}

void
finalize()
{
}

} // namespace transport_params
