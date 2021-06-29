#ifndef ERF_BCINTERIOR_H
#define ERF_BCINTERIOR_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCDummy : public BCBase {
  public:
  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab>& vars, int dim) override {}
  bool isInterior() override { return true; }
};

}

#endif
