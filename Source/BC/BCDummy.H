#ifndef ERF_BCDUMMY_H
#define ERF_BCDUMMY_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCDummy : public BCBase {
  public:
  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) override {}
};

}

#endif
