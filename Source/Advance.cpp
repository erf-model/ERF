#include <cmath>

#include "ERF.H"
#include "TimeIntegration.H"
#include "IndexDefines.H"

#include "AMReX_InterpFaceRegister.H"

using namespace amrex;

/** The main driver for a single level time advance.
 *
 *     @param time the current simulation time
 *     @param dt the timestep to advance (e.g., go from time to time + dt)
 *     @param amr_iteration where we are in the current AMR subcycle.  Each
 *                     level will take a number of steps to reach the
 *                     final time of the coarser level below it.  This
 *                     counter starts at 1
 *     @param amr_ncycle  the number of subcycles at this level
 */

Real
ERF::advance(Real time, Real dt, int /*amr_iteration*/, int /*amr_ncycle*/)
{

    BL_PROFILE("ERF::advance()");

    int finest_level = parent->finestLevel();

    if (level < finest_level && do_reflux) {
      getFluxReg(level + 1).reset();
    }

    BL_PROFILE("ERF::do_rk3_advance()");

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

  MultiFab* S_crse;
  MultiFab rU_crse, rV_crse, rW_crse;

  if (level > 0)
  {
      S_crse = &getLevel(level-1).get_old_data(State_Type);

      MultiFab& U_crse = getLevel(level-1).get_old_data(X_Vel_Type);
      MultiFab& V_crse = getLevel(level-1).get_old_data(Y_Vel_Type);
      MultiFab& W_crse = getLevel(level-1).get_old_data(Z_Vel_Type);

      rU_crse.define(U_crse.boxArray(), U_crse.DistributionMap(), 1, 1);
      rV_crse.define(V_crse.boxArray(), V_crse.DistributionMap(), 1, 1);
      rW_crse.define(W_crse.boxArray(), W_crse.DistributionMap(), 1, 1);

      VelocityToMomentum(U_crse,V_crse,W_crse,*S_crse,rU_crse,rV_crse,rW_crse,
                        solverChoice.spatial_order);
  }

  // Fill level 0 ghost cells (including at periodic boundaries)
  //TODO: Check if we should consider the number of ghost cells as a function of spatial order here itself
  S_old.FillBoundary(geom.periodicity());
  U_old.FillBoundary(geom.periodicity());
  V_old.FillBoundary(geom.periodicity());
  W_old.FillBoundary(geom.periodicity());

  const auto& ref_ratio = (level > 0) ? parent->refRatio(level-1) : IntVect(1,1,1);

  InterpFaceRegister ifr;
  if (level > 0)
  {
      ifr.define(S_old.boxArray(), S_old.DistributionMap(), geom, ref_ratio);
  }

  const Real* dx = geom.CellSize();
  const BoxArray&            ba = S_old.boxArray();
  const DistributionMapping& dm = S_old.DistributionMap();

  int nvars = S_old.nComp();

  // Place-holder for source array -- for now just set to 0
  MultiFab source(ba,dm,nvars,1);
  source.setVal(0.0);

  // We keep faceflux so that we can re-fluxing operations for two-way coupling
  std::array< MultiFab, AMREX_SPACEDIM > faceflux;
  //faceflux[0] is of size (ncells_x + 1, ncells_y    , ncells_z    )
  faceflux[0].define(convert(ba,IntVect(1,0,0)), dmap, nvars, 0);
  //faceflux[1] is of size (ncells_x    , ncells_y + 1, ncells_z    )
  faceflux[1].define(convert(ba,IntVect(0,1,0)), dmap, nvars, 0);
  //faceflux[2] is of size (ncells_x    , ncells_y    , ncells_z + 1)
  faceflux[2].define(convert(ba,IntVect(0,0,1)), dmap, nvars, 0);

  // Make sure to fill the ghost cells
  MultiFab state_mf(grids,dmap,nvars,S_old.nGrow());
  FillPatch(*this,state_mf,S_old.nGrow(),time,0,0,nvars);

  // *****************************************************************
  // Update the cell-centered state and face-based velocity using
  // a time integrator.
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
  erf_advance(level,
              state_mf, S_new,
              U_old, V_old, W_old,
              U_new, V_new, W_new,
              rU_crse, rV_crse, rW_crse,
              source,
              faceflux,
              (level > 0) ? parent->Geom(level-1) : geom,
              geom,
              ref_ratio,
              dx, dt, time, &ifr,
              solverChoice);

  return dt;
}
