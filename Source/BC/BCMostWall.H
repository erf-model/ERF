#ifndef ERF_BCMOSTWALL_H
#define ERF_BCMOSTWALL_H

#include "ERF.H"
#include "BCBase.H"
#include "ABLMost.H"
#include "EddyViscosity.H"
#include "PlaneAverage.H"
#include "VelPlaneAverage.H"

using namespace amrex;

//namespace phys_bcs {

struct GetViscosity {
   static void get(MultiFab* velx, MultiFab* vely, MultiFab* velz,
                   MultiFab* states, MultiFab* eta,
                   const GpuArray<Real, AMREX_SPACEDIM>& cellSize,
                   const SolverChoice& solverChoice,
                   const bool& no_slip_lo_k, const int& klo,
                   const bool& no_slip_hi_k, const int& khi) {
     ComputeTurbulentViscosity(*velx, *vely, *velz, *states, *eta, cellSize, solverChoice, no_slip_lo_k, klo, no_slip_hi_k, khi);
   }
};

namespace phys_bcs {

template <int IDIR, math_bcs::BCBound Bound>
class BCMostWall : public BCBase {
  public:
  BCMostWall() = default;
  using Dir = DirectionSelector<IDIR>;

  // we need all velocity variables and temperature/theta as input
  void applyBC (const Geometry geom, Vector<MultiFab*>& vars, const SolverChoice& solverChoice) override {
    if ((geom.isPeriodic(IDIR)) || (vars.size() == 0)) return;

    auto most = ABLMost();
    int idir0, idir1;

    if (IDIR == 0) {
      idir0 = 1;
      idir1 = 2;
    } else if (IDIR == 1) {
      idir0 = 0;
      idir1 = 2;
    } else if (IDIR == 2) {
      idir0 = 0;
      idir1 = 1;
    }

    MultiFab* velx   = getVelx (vars);
    MultiFab* vely   = getVely (vars);
    MultiFab* velz   = getVelz (vars);
    MultiFab* states = getState(vars);

    bool no_slip_at_lo_k = true;
    bool no_slip_at_hi_k = true;

    int klo = geom.Domain().smallEnd()[2];
    int khi = geom.Domain().bigEnd()[2];

    MultiFab eta(states->boxArray(),states->DistributionMap(),1,1);
    const GpuArray<Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
    GetViscosity::get(velx, vely, velz, states, &eta, dx, solverChoice,
                      no_slip_at_lo_k, klo, no_slip_at_hi_k, khi);

    PlaneAverage save (states, geom, IDIR, true);
    PlaneAverage vxave(velx,   geom, IDIR, true);
    PlaneAverage vyave(vely,   geom, IDIR, true);
    PlaneAverage vzave(velz,   geom, IDIR, true);
    VelPlaneAverage vmagave({velx,vely,velz}, geom, IDIR, true);

    save. compute_averages(Dir(), save.field());
    vxave.compute_averages(Dir(), vxave.field());
    vyave.compute_averages(Dir(), vyave.field());
    vzave.compute_averages(Dir(), vzave.field());
    vmagave.compute_hvelmag_averages(Dir(), idir0, idir1, vmagave.field());

    most.vel_mean[0] = vxave.line_average_interpolated(most.zref, 0);
    most.vel_mean[1] = vyave.line_average_interpolated(most.zref, 0);
    most.vel_mean[2] = vzave.line_average_interpolated(most.zref, 0);
    most.vmag_mean   = vmagave.line_hvelmag_average_interpolated(most.zref);
    most.theta_mean  = save.line_average_interpolated(most.zref, 0);

    most.update_fluxes();

    // setup boundary conditions
    for( auto i = 0; i < vars.size(); ++i) {
      int nghost = vars[i]->nGrow();
      if(nghost > 0) {
         for (MFIter mfi(*vars[i], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            // get the index box
            Box b = vars[i]->get(mfi).box();
            auto ncomp = vars[i]->nComp();
            Array4<Real> data         = vars[i]->array(mfi);
            Array4<Real> velx_array   = velx->array(mfi);
            Array4<Real> vely_array   = vely->array(mfi);
            Array4<Real> velz_array   = velz->array(mfi);
            Array4<Real> states_array = states->array(mfi);
            Array4<Real> eta_array    = eta.array(mfi);
            Vector<Array4<Real>*> sarray{&velx_array, &vely_array, &velz_array, &states_array, &eta_array};

            // index box type
            IntVect type = b.type();

            if (type == IntVect::Zero) {  // for state type variable
               math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, 0, ncomp);
            } else if (type == IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL
               if (IDIR == 0) { // x-vel on x-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_face_based(geom, b, data, ncomp);
               } else { // x-vel on y- or z-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, 1, ncomp);
               }
            } else if (type == IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
               if (IDIR == 1) { // y-vel on y-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_face_based(geom, b, data, ncomp);
               } else { // y-vel on x- or z-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, 2, ncomp);
               }
            } else if (type == IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
               if (IDIR == 2) { // z-vel on z-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_face_based(geom, b, data, ncomp);
               } else { // z-vel on x- or y-faces
                   math_bcs::ablmostwall<IDIR, Bound>::apply_cell_based(geom, b, sarray, data, most, 3, ncomp);
               }
            }
         }
     }
  }
 }

MultiFab* getVelx (const Vector<MultiFab*>& vars)
{
    MultiFab* velx{nullptr};
    for( auto i = 0; i < vars.size(); ++i) {
       IntVect type = vars[i]->boxArray()[0].type();
       if (type == IntVect(AMREX_D_DECL(1, 0, 0))) { // for X_VEL
           velx = vars[i];
           break;
       }
   }
   return velx;
}

MultiFab* getVely (const Vector<MultiFab*>& vars)
{
    MultiFab* vely{nullptr};
    for( auto i = 0; i < vars.size(); ++i) {
       IntVect type = vars[i]->boxArray()[0].type();
       if (type == IntVect(AMREX_D_DECL(0, 1, 0))) { // for Y_VEL
           vely = vars[i];
           break;
       }
   }
   return vely;
 }

MultiFab* getVelz (const Vector<MultiFab*>& vars)
{
    MultiFab* velz{nullptr};
    for( auto i = 0; i < vars.size(); ++i) {
       IntVect type = vars[i]->boxArray()[0].type();
       if (type == IntVect(AMREX_D_DECL(0, 0, 1))) { // for Z_VEL
           velz = vars[i];
           break;
       }
   }
   return velz;
}

MultiFab* getState (const Vector<MultiFab*>& vars)
{
    MultiFab* state{nullptr};
    for( auto i = 0; i < vars.size(); ++i) {
       IntVect type = vars[i]->boxArray()[0].type();
       if (type == IntVect::Zero) {
           state = vars[i];
           break;
       }
   }
   return state;
}

};
}

#endif
