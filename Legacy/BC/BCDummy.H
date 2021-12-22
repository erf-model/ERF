#ifndef ERF_BCDUMMY_H
#define ERF_BCDUMMY_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCDummy : public BCBase {
  public:
  explicit BCDummy(ERF& erf) : BCBase(erf) {}
  void applyBC (const amrex::Geometry /*geom*/, amrex::Vector<MultiFab*>& /*vars*/, const SolverChoice& /*solverChoice*/, bool is_derived) override {}
};

}

#endif
