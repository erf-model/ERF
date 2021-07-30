#ifndef ERF_BCINTERIOR_H
#define ERF_BCINTERIOR_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCInterior : public BCBase {
  public:
  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) override {}
  bool isInterior() override { return true; }
};

}

#endif
