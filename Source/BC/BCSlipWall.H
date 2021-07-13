#ifndef ERF_BCSLIPWALL_H
#define ERF_BCSLIPWALL_H

#include "ERF.H"
#include "BCBase.H"

using namespace amrex;

namespace phys_bcs {

template <int IDIR, math_bcs::BCBound Bound>
class BCSlipWall : public BCBase {
  public:
  BCSlipWall() = default;

  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars) override {
    if ((geom.isAllPeriodic()) || (vars.size() == 0)) return;

    // setup boundary conditions
    for( auto i = 0; i < vars.size(); ++i) {
      int nghost = vars[i]->nGrow();
      if(nghost > 0) { 
         for (amrex::MFIter mfi(*vars[i], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // get the index box
            amrex::Box b = vars[i]->get(mfi).box();
            auto ncomp = vars[i]->nComp();
            amrex::Array4<amrex::Real> data = vars[i]->array(mfi);

            // index box type
            amrex::IntVect type = b.type();
 
            if (type == amrex::IntVect::Zero) {  // for state type variable
               math_bcs::foextrap<IDIR, Bound>::apply_cell_based(geom, b, data, ncomp);
            } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
               math_bcs::foextrap<0, Bound>::template apply_face_based<(IDIR==0)>(geom, b, data, ncomp);
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               math_bcs::foextrap<1, Bound>::template apply_face_based<(IDIR==1)>(geom, b, data, ncomp);
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               math_bcs::foextrap<2, Bound>::template apply_face_based<(IDIR==2)>(geom, b, data, ncomp);
            }
         }
     }
  }

  }
};

}

#endif
