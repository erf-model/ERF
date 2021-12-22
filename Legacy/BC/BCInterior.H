#ifndef ERF_BCINTERIOR_H
#define ERF_BCINTERIOR_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCInterior : public BCBase {
  public:
  explicit BCInterior(ERF& erf) : BCBase(erf) {}
  bool isInterior() override { return true; }
  void applyBC (const amrex::Geometry /*geom*/, amrex::Vector<MultiFab*>& /*vars*/, const SolverChoice& /*solverChoice*/, bool is_derived) override {}
};

}

#endif
