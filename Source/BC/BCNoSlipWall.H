#ifndef ERF_BCNOSLIPWALL_H
#define ERF_BCNOSLIPWALL_H

#include "ERF.H"
#include "BCBase.H"

using namespace amrex;

namespace phys_bcs {

template <math_bcs::BCBound Bound>
class BCNoSlipWall : public BCBase {
  public:
  BCNoSlipWall() = default;

  void applyBC (const amrex::Geometry geom, amrex::Vector<MultiFab*>& vars, int dim) override {
    if ((geom.isAllPeriodic()) || (vars.size() == 0)) return;

    // Physical Domain
    amrex::Box dom(geom.Domain());

    int extent = (Bound == math_bcs::BCBound::lower) ? dom.smallEnd(dim) : dom.bigEnd(dim);
        
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
               math_bcs::wall_scalars<Bound>(geom, b, extent, data, dim, ncomp);
            } else if (type == amrex::IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL 
               if (dim == 0) {
                  math_bcs::reflect_odd_face_based<Bound>(geom, b, extent, data, dim, 1);
               } else {
                  math_bcs::reflect_odd<Bound>(geom, b, extent, data, dim, 1);
               }
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               if (dim == 1) {
                  math_bcs::reflect_odd_face_based<Bound>(geom, b, extent, data, dim, 1);
               } else {
                  math_bcs::reflect_odd<Bound>(geom, b, extent, data, dim, 1);
               }
            } else if (type == amrex::IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               if (dim == 2) {
                  math_bcs::reflect_odd_face_based<Bound>(geom, b, extent, data, dim, 1);
               } else {
                  math_bcs::reflect_odd<Bound>(geom, b, extent, data, dim, 1);
               }
            }
         }
     }
  }

  }
};

}

#endif
