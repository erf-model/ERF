#include <cmath>

#include "ERF.H"
#include "RK3.H"
#include "IndexDefines.H"

using namespace amrex;

Real
ERF::advance(Real time, Real dt, int amr_iteration, int amr_ncycle)
{
  /** the main driver for a single level implementing the time advance.

         @param time the current simulation time
         @param dt the timestep to advance (e.g., go from time to time + dt)
         @param amr_iteration where we are in the current AMR subcycle.  Each
                         level will take a number of steps to reach the
                         final time of the coarser level below it.  This
                         counter starts at 1
         @param amr_ncycle  the number of subcycles at this level
  */

  BL_PROFILE("ERF::advance()");

  int finest_level = parent->finestLevel();

  if (level < finest_level && do_reflux) {
    getFluxReg(level + 1).reset();
  }

  Real dt_new = dt;

  BL_PROFILE("ERF::do_rk3_advance()");

  // Check that we are not asking to advance stuff we don't know to
  // if (src_list.size() > 0) amrex::Abort("Have not integrated other sources
  // into MOL advance yet");

  for (int i = 0; i < num_state_type; ++i) {
    bool skip = false;
    if (!skip) {
      state[i].allocOldData();
      state[i].swapTimeLevels(dt);
    }
  }

  MultiFab& S_old = get_old_data(State_Type);
  MultiFab& S_new = get_new_data(State_Type);

  MultiFab& U_new = get_new_data(X_Vel_Type);
  MultiFab& V_new = get_new_data(Y_Vel_Type);
  MultiFab& W_new = get_new_data(Z_Vel_Type);

  MultiFab& Xmom_old = get_old_data(X_Mom_Type);
  MultiFab& Ymom_old = get_old_data(Y_Mom_Type);
  MultiFab& Zmom_old = get_old_data(Z_Mom_Type);

  MultiFab& Xmom_new = get_new_data(X_Mom_Type);
  MultiFab& Ymom_new = get_new_data(Y_Mom_Type);
  MultiFab& Zmom_new = get_new_data(Z_Mom_Type);

  // Fill level 0 ghost cells (including at periodic boundaries)
  S_old.FillBoundary(geom.periodicity());
  Xmom_old.FillBoundary(geom.periodicity());
  Ymom_old.FillBoundary(geom.periodicity());
  Zmom_old.FillBoundary(geom.periodicity());

  const Real* dx = geom.CellSize();
  const BoxArray&            ba = S_old.boxArray();
  const DistributionMapping& dm = S_old.DistributionMap();

  int nvars = S_old.nComp();

  // Place-holder for source array -- for now just set to 0
  MultiFab source(ba,dm,nvars,1); 
  source.setVal(0.0);

  // Place-holder for eta array -- shear viscosity -- for now just set to 0
  MultiFab eta(ba,dm,1,1); 
  eta.setVal(0.0);

  // Place-holder for zeta array -- bulk viscosity -- for now just set to 0
  MultiFab zeta(ba,dm,1,1); 
  zeta.setVal(0.0);

  // Place-holder for kappa array -- thermal conducitivity --for now just set to 0
  MultiFab kappa(ba,dm,1,1); 
  kappa.setVal(0.0);

  //fluxes (except momentum) at faces
  std::array< MultiFab, AMREX_SPACEDIM > faceflux;
  faceflux[0].define(convert(ba,IntVect(1,0,0)), dmap, nvars, 0);
  faceflux[1].define(convert(ba,IntVect(0,1,0)), dmap, nvars, 0);
  faceflux[2].define(convert(ba,IntVect(0,0,1)), dmap, nvars, 0);

  std::array< MultiFab, 2 > edgeflux_x; // v, w
  std::array< MultiFab, 2 > edgeflux_y; // u, w
  std::array< MultiFab, 2 > edgeflux_z; // u, v

  edgeflux_x[0].define(convert(ba,IntVect(1,1,0)), dmap, 1, 0); // v
  edgeflux_x[1].define(convert(ba,IntVect(1,0,1)), dmap, 1, 0); // w

  edgeflux_y[0].define(convert(ba,IntVect(1,1,0)), dmap, 1, 0); // u
  edgeflux_y[1].define(convert(ba,IntVect(0,1,1)), dmap, 1, 0); // w

  edgeflux_z[0].define(convert(ba,IntVect(1,0,1)), dmap, 1, 0); // u
  edgeflux_z[1].define(convert(ba,IntVect(0,1,1)), dmap, 1, 0); // v

  std::array< MultiFab, AMREX_SPACEDIM > cenflux;
  cenflux[0].define(ba,dmap,1,1); // 0-2: rhoU, rhoV, rhoW
  cenflux[1].define(ba,dmap,1,1);
  cenflux[2].define(ba,dmap,1,1);

  // *****************************************************************
  // Update the cell-centered state and face-based momenta using RK3
  // Inputs:  
  //          S_old    (state on cell centers)
  //          Xmom_old (x-momentum on x-faces)
  //          Ymom_old (y-momentum on y-faces)
  //          Zmom_old (z-momentum on z-faces)
  //          source   (source term on cell centers)
  // Outputs:  
  //          S_new    (state on cell centers)
  //          Xmom_new (x-momentum on x-faces)
  //          Ymom_new (y-momentum on y-faces)
  //          Zmom_new (z-momentum on z-faces)
  //          U_new    (Xmom_new / density on x-faces)
  //          V_new    (Ymom_new / density on y-faces)
  //          W_new    (Zmom_new / density on z-faces)
  // *****************************************************************

  RK3_advance(S_old, S_new, 
              Xmom_old, Ymom_old, Zmom_old,
              Xmom_new, Ymom_new, Zmom_new, 
              U_new, V_new, W_new, 
              source, 
              eta, zeta, kappa, 
              faceflux, 
              edgeflux_x, edgeflux_y, edgeflux_z, 
              cenflux, geom, dx, dt);

  return dt;
}
