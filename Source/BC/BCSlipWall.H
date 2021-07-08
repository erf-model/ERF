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
        
    // setup boundary conditions
    for( auto i = 0; i < vars.size(); ++i) {
      int nghost = vars[i]->nGrow();
      if(nghost > 0) { 
         for (amrex::MFIter mfi(*vars[i], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // get the index box
            amrex::Box b = vars[i]->get(mfi).box();
            auto ncomp = vars[i]->nComp();
            const amrex::Array4<amrex::Real> data = vars[i]->array(mfi);

            // index box type
            amrex::IntVect type = b.type();
 
            if (type == amrex::IntVect::Zero) {  // for state type variable
               if (dim == 0) {
                  // low x 
                  if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data, ncomp);
                  // high x 
                  if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data, ncomp);
               } else if (dim == 1) {
                  // low y 
                  if(b.smallEnd(1) < lo) math_bcs::foextrap_lo(geom, b, lo, data, ncomp);
                  // high y
                  if(b.bigEnd(1) > hi) math_bcs::foextrap_hi(geom, b, hi, data, ncomp);
               } else if (dim == 2) {
                  // low z 
                  if(b.smallEnd(2) < lo) math_bcs::foextrap_lo(geom, b, lo, data, ncomp);
                  // high z
                  if(b.bigEnd(2) > hi) math_bcs::foextrap_hi(geom, b, hi, data, ncomp);
               }
            } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
               if (dim == 0) {
                  // TODO: we need a special extrapolation function for edge based X_VEL
                  // low x 
                  if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // TODO: same here
                  // high x 
                  if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 1) {
                  // low y 
                  if(b.smallEnd(1) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high y 
                  if(b.bigEnd(1) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 2) {
                  // low z
                  if(b.smallEnd(2) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high z 
                  if(b.bigEnd(2) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               }
            
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               if (dim == 0) {
                  // low x
                  if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high x 
                  if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 1) {
                  // low y
                  if(b.smallEnd(1) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high y
                  if(b.bigEnd(1) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 2) {
                  // low z
                  if(b.smallEnd(2) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high z
                  if(b.bigEnd(2) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               }
 
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               if (dim == 0) {
                  // low x
                  if(b.smallEnd(0) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high x
                  if(b.bigEnd(0) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 1) {
                  // low y
                  if(b.smallEnd(1) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high y
                  if(b.bigEnd(1) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               } else if (dim == 2) {
                  // low z
                  if(b.smallEnd(2) < lo) math_bcs::foextrap_lo(geom, b, lo, data);
                  // high z
                  if(b.bigEnd(2) > hi) math_bcs::foextrap_hi(geom, b, hi, data);
               }
            }
         }
     }
  }

  }
};

}

#endif
