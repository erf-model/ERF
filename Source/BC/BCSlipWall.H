#ifndef ERF_BCSLIPWALL_H
#define ERF_BCSLIPWALL_H

#include "ERF.H"
using namespace amrex;

namespace phys_bcs {

class BCSlipWall : public BCBase {
  public:
  BCSlipWall() = default;

  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars, int dim) override {
    if ((geom.isAllPeriodic()) || (vars.size() == 0)) return;

    // Physical Domain
    amrex::Box dom(geom.Domain());

    int lo = dom.smallEnd(dim);
    int hi = dom.bigEnd(dim);
        
    if(vars[0]->nGrow() > 0) {  // cons_old
      int nghost = vars[0]->nGrow();
      for (amrex::MFIter mfi(*vars[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.growntilebox(nghost);
        const amrex::Array4<amrex::Real> data = vars[0]->array(mfi);
        if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
        if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
      }
    }

    if(vars[0]->nGrow() > 0) {  // xvel
      int nghost = vars[0]->nGrow();
      for (amrex::MFIter mfi(*vars[1], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.growntilebox(nghost);
        const amrex::Array4<amrex::Real> data = vars[0]->array(mfi);
        if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
        if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
      }
    }

    //TODO: Add code for other MultiFabs
  }
};

}

#endif
