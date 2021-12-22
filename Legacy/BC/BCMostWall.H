#ifndef ERF_BCMOSTWALL_H
#define ERF_BCMOSTWALL_H

#include "ERF.H"
#include "BCBase.H"
#include "ABLMost.H"
#include "EddyViscosity.H"
#include "PlaneAverage.H"
#include "VelPlaneAverage.H"

using namespace amrex;

struct getViscosity {
   static void get(const MultiFab& velx, const MultiFab& vely, const MultiFab& velz,
                   const MultiFab& states, MultiFab& eta,
                   const Geometry& geom,
                   const SolverChoice& solverChoice) {
     int klo = geom.Domain().smallEnd()[2];
     int khi = geom.Domain().bigEnd()[2];
     const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
     ComputeTurbulentViscosity(velx, vely, velz, states, eta, dx, solverChoice, false, klo, false, khi);
     eta.FillBoundary(geom.periodicity());
   }
};

namespace phys_bcs {

template <int IDIR, math_bcs::BCBound Bound>
class BCMostWall : public BCBase {
  public:
  explicit BCMostWall(ERF& erf) : BCBase(erf) {}
  using Dir = DirectionSelector<IDIR>;

  // we need all velocity variables and temperature/theta as input
  void applyBC (const Geometry geom, Vector<MultiFab*>& vars, const SolverChoice& solverChoice, bool is_derived=false) override {
    if ((geom.isPeriodic(IDIR)) || (vars.size() == 0)) return;

    MultiFab& S_old = erf.get_old_data(0);
    MultiFab& U_old = erf.get_old_data(1);
    MultiFab& V_old = erf.get_old_data(2);
    MultiFab& W_old = erf.get_old_data(3);
    MultiFab eta(S_old.boxArray(),S_old.DistributionMap(),1,S_old.nGrow());

    auto most = erf.get_most();

    // setup the viscosity for most use
    getViscosity::get(U_old, V_old, W_old, S_old, eta, geom, solverChoice);

    // setup boundary conditions
    for( auto i = 0; i < vars.size(); ++i) {
      int nghost = vars[i]->nGrow();

      if(nghost > 0) {
         for (MFIter mfi(*vars[i], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // get the index box
            Box b = vars[i]->get(mfi).box();
            auto ncomp = vars[i]->nComp();
            // index box type
            IntVect type = b.type();
            Array4<Real> data = vars[i]->array(mfi);

            if (type == IntVect::Zero) {  // for state type variable
               Array4<Real> velx_array   = U_old.array(mfi);
               Array4<Real> vely_array   = V_old.array(mfi);
               Array4<Real> velz_array   = W_old.array(mfi);
               Array4<Real> states_array = S_old.array(mfi);
               Array4<Real> eta_array    = eta.array(mfi);
               Vector<Array4<Real>*> sarray{&velx_array, &vely_array, &velz_array, &eta_array, &states_array};
               math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, is_derived, ncomp);

            } else if (type == IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL
               if (IDIR == 0) { // x-vel on x-faces
                   math_bcs::reflect_odd<IDIR, Bound>::template apply_face_based<(IDIR==2)>(geom, b, data, ncomp);
               } else { // x-vel on y- or z-faces
                   Array4<Real> velx_array   = U_old.array(mfi);
                   Array4<Real> vely_array   = V_old.array(mfi);
                   Array4<Real> velz_array   = W_old.array(mfi);
                   Array4<Real> states_array = S_old.array(mfi);
                   Array4<Real> eta_array    = eta.array(mfi);
                   Vector<Array4<Real>*> sarray{&velx_array, &vely_array, &velz_array, &eta_array, &states_array};
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, is_derived, ncomp);
               }

            } else if (type == IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               if (IDIR == 1) { // y-vel on y-faces
                   math_bcs::reflect_odd<IDIR, Bound>::template apply_face_based<(IDIR==2)>(geom, b, data, ncomp);
               } else { // y-vel on x- or z-faces
                   Array4<Real> velx_array   = U_old.array(mfi);
                   Array4<Real> vely_array   = V_old.array(mfi);
                   Array4<Real> velz_array   = W_old.array(mfi);
                   Array4<Real> states_array = S_old.array(mfi);
                   Array4<Real> eta_array    = eta.array(mfi);
                   Vector<Array4<Real>*> sarray{&velx_array, &vely_array, &velz_array, &eta_array, &states_array};
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, is_derived, ncomp);
               }

            } else if (type == IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               if (IDIR == 2) { // z-vel on z-faces
                   math_bcs::reflect_odd<IDIR, Bound>::template apply_face_based<(IDIR==2)>(geom, b, data, ncomp);
               } else { // z-vel on x- or y-faces
                   Array4<Real> states_array = S_old.array(mfi);
                   Array4<Real> eta_array    = eta.array(mfi);
                   Vector<Array4<Real>*> sarray{&states_array, &eta_array};
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, is_derived, ncomp);
               }
           }
         }
     }
  }
 }
};
}

#endif
