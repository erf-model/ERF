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

    amrex::IntVect lo = dom.smallEnd();
    amrex::IntVect hi = dom.bigEnd();
        
    // setup boundary conditions
    for( auto i = 0; i < vars.size(); ++i) {
      int nghost = vars[i]->nGrow();
      if(nghost > 0) { 
         for (amrex::MFIter mfi(*vars[i], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // get the index box
            amrex::Box b = vars[i]->get(mfi).box();
            const amrex::Array4<amrex::Real> data = vars[i]->array(mfi);

            // index box type
            amrex::IntVect type = b.type();
 
            if (type == amrex::IntVect::Zero) {  // for state type variable
               // low x boundary
               if(b.smallEnd(0) < lo[0]) math_bcs::foextrap_lo(geom, b, lo[0], data);
               // high x boundary
               if(b.bigEnd(0) > hi[0]) math_bcs::foextrap_hi(geom, b, hi[0], data);
               // low y boundary
               if(b.smallEnd(1) < lo[1]) math_bcs::foextrap_lo(geom, b, lo[1], data);
               // high y boundary
               if(b.bigEnd(1) > hi[1]) math_bcs::foextrap_hi(geom, b, hi[1], data);
               // low z boundary
               if(b.smallEnd(2) < lo[2]) math_bcs::foextrap_lo(geom, b, lo[2], data);
               // high z boundary 
               if(b.bigEnd(2) > hi[2]) math_bcs::foextrap_hi(geom, b, hi[2], data);
            } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 

               // TODO: we need a special extrapolation function for edge based X_VEL
               // low x 
               if(b.smallEnd(0) < lo[0]) math_bcs::foextrap_lo(geom, b, lo[0], data);
               // TODO: same here
               // high x 
               if(b.bigEnd(0) > hi[0]) math_bcs::foextrap_hi(geom, b, hi[0], data);
               // low y 
               if(b.smallEnd(1) < lo[1]) math_bcs::foextrap_lo(geom, b, lo[1], data);
               // high y 
               if(b.bigEnd(1) > hi[1]) math_bcs::foextrap_hi(geom, b, hi[1], data);
               // low z
               if(b.smallEnd(2) < lo[2]) math_bcs::foextrap_lo(geom, b, lo[2], data);
               // high z 
               if(b.bigEnd(2) > hi[2]) math_bcs::foextrap_hi(geom, b, hi[2], data);
            
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               // low x
               if(b.smallEnd(0) < lo[0]) math_bcs::foextrap_lo(geom, b, lo[0], data);
               // high x 
               if(b.bigEnd(0) > hi[0]) math_bcs::foextrap_hi(geom, b, hi[0], data);
               // low y
               if(b.smallEnd(1) < lo[1]) math_bcs::foextrap_lo(geom, b, lo[1], data);
               // high y
               if(b.bigEnd(1) > hi[1]) math_bcs::foextrap_hi(geom, b, hi[1], data);
               // low z
               if(b.smallEnd(2) < lo[2]) math_bcs::foextrap_lo(geom, b, lo[2], data);
               // high z
               if(b.bigEnd(2) > hi[2]) math_bcs::foextrap_hi(geom, b, hi[2], data);
 
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               // low x
               if(b.smallEnd(0) < lo[0]) math_bcs::foextrap_lo(geom, b, lo[0], data);
               // high x
               if(b.bigEnd(0) > hi[0]) math_bcs::foextrap_hi(geom, b, hi[0], data);
               // low y
               if(b.smallEnd(1) < lo[1]) math_bcs::foextrap_lo(geom, b, lo[1], data);
               // high y
               if(b.bigEnd(1) > hi[1]) math_bcs::foextrap_hi(geom, b, hi[1], data);
               // low z
               if(b.smallEnd(2) < lo[2]) math_bcs::foextrap_lo(geom, b, lo[2], data);
                // high z
               if(b.bigEnd(2) > hi[2]) math_bcs::foextrap_hi(geom, b, hi[2], data);
            }
         }
     }
  }

  }
};

}

#endif
