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

//  Real dt_new = dt;

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

  MultiFab& U_old = get_old_data(X_Vel_Type);
  MultiFab& V_old = get_old_data(Y_Vel_Type);
  MultiFab& W_old = get_old_data(Z_Vel_Type);

  MultiFab& U_new = get_new_data(X_Vel_Type);
  MultiFab& V_new = get_new_data(Y_Vel_Type);
  MultiFab& W_new = get_new_data(Z_Vel_Type);

  // Fill level 0 ghost cells (including at periodic boundaries)
  //TODO: Check if we should consider the number of ghost cells as a function of spatial order here itself
  S_old.FillBoundary(geom.periodicity());
  U_old.FillBoundary(geom.periodicity());
  V_old.FillBoundary(geom.periodicity());
  W_old.FillBoundary(geom.periodicity());

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

  // TODO: We won't need faceflux, edgeflux, and centflux when using the new code architecture. Remove them.
  // Fluxes (except momentum) at faces. This should comprise of advective as well as diffusive fluxes.
  // There are separate variables to handle the momentum at the faces
  std::array< MultiFab, AMREX_SPACEDIM > faceflux;
  //faceflux[0] is of size (ncells_x + 1, ncells_y    , ncells_z    )
  faceflux[0].define(convert(ba,IntVect(1,0,0)), dmap, nvars, 0);
  //faceflux[1] is of size (ncells_x    , ncells_y + 1, ncells_z    )
  faceflux[1].define(convert(ba,IntVect(0,1,0)), dmap, nvars, 0);
  //faceflux[2] is of size (ncells_x    , ncells_y    , ncells_z + 1)
  faceflux[2].define(convert(ba,IntVect(0,0,1)), dmap, nvars, 0);

  // Edge fluxes for {x, y, z}-momentum equations
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

  // TODO: Better make it a member of the ERF class. Need to deal with static stuff.
   SolverChoice solverChoice(use_advection, use_thermal_diffusion, alpha_T,
                             use_scalar_diffusion, alpha_S,
                             use_momentum_diffusion, kinematicViscosity,
                             use_smagorinsky, use_gravity, spatial_order);
  //solverChoice.display();

  // *****************************************************************
  // Update the cell-centered state and face-based velocity using RK3
  // Inputs:  
  //          S_old    (state on cell centers)
  //          U_old    (x-velocity on x-faces)
  //          V_old    (y-velocity on y-faces)
  //          W_old    (z-velocity on z-faces)
  //          source   (source term on cell centers)
  // Outputs:  
  //          S_new    (state on cell centers)
  //          U_new    (x-velocity on x-faces)
  //          V_new    (y-velocity on y-faces)
  //          W_new    (z-velocity on z-faces)
  // *****************************************************************
  RK3_advance(
              S_old, S_new,
              U_old, V_old, W_old,
              U_new, V_new, W_new,
              source,
              eta, zeta,kappa,
              faceflux,
              edgeflux_x, edgeflux_y, edgeflux_z,
              cenflux, geom, dx, dt,
              solverChoice);

  return dt;
}
