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

  MultiFab& Sx_old = get_old_data(X_State_Type);
  MultiFab& Sy_old = get_old_data(Y_State_Type);
  MultiFab& Sz_old = get_old_data(Z_State_Type);

  MultiFab& U_old = get_old_data(X_Vel_Type);
  MultiFab& V_old = get_old_data(Y_Vel_Type);
  MultiFab& W_old = get_old_data(Z_Vel_Type);

  const Real* dx = geom.CellSize();
  const BoxArray&            ba = S_old.boxArray();
  const DistributionMapping& dm = S_old.DistributionMap();

  int nvars = S_old.nComp();

  // Temporary MultiFab to hold the primitive variables
  MultiFab prim(ba,dm,nvars,2); 

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

  RK3_advance(S_old, S_new, Sx_old, Sy_old, Sz_old, prim, 
              U_old, W_old, W_old, source, 
              eta, zeta, kappa, 
              faceflux, 
              edgeflux_x, edgeflux_y, edgeflux_z, 
              cenflux, geom, dx, dt);

/*
  void RK3_advance(MultiFab& cu,
                   std::array< MultiFab, AMREX_SPACEDIM >& cumom,
                   MultiFab& prim, 
                   std::array< MultiFab, AMREX_SPACEDIM >& vel,
                   MultiFab& source,
                   MultiFab& eta, MultiFab& zeta, MultiFab& kappa,
                   std::array< MultiFab, AMREX_SPACEDIM>& faceflux,
                   std::array< MultiFab, 2 >& edgeflux_x,
                   std::array< MultiFab, 2 >& edgeflux_y,
                   std::array< MultiFab, 2 >& edgeflux_z,
                   std::array< MultiFab, AMREX_SPACEDIM>& cenflux,
                   const Geometry geom, const Real* dxp, const Real dt)
*/

  // 
  // The code in brackets below is the old do_mol_advance
  // dt_new = do_mol_advance(time, dt, amr_iteration, amr_ncycle);
  // 
  {
  if (do_mol_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  MultiFab& U_old = get_old_data(State_Type);
  MultiFab& U_new = get_new_data(State_Type);
  MultiFab S(grids, dmap, NVAR, 0, MFInfo(), Factory());

  MultiFab S_old, S_new;
  if (mol_iters > 1) {
    S_old.define(grids, dmap, NVAR, 0, MFInfo(), Factory());
    S_new.define(grids, dmap, NVAR, 0, MFInfo(), Factory());
  }

  // Compute S^{n} = MOLRhs(U^{n})
  if (verbose) {
    amrex::Print() << "... Computing MOL source term at t^{n} " << std::endl;
  }
  FillPatch(*this, Sborder, nGrowTr, time, State_Type, 0, NVAR);
  Real flux_factor = 0;
  getMOLSrcTerm(Sborder, S, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_old
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
    ) {
      construct_old_source(
        src_list[n], time, dt, amr_iteration, amr_ncycle, 0, 0);
      MultiFab::Saxpy(S, 1.0, *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  // U^* = U^n + dt*S^n
  MultiFab::LinComb(U_new, 1.0, Sborder, 0, dt, S, 0, 0, NVAR, 0);

  computeTemp(U_new, 0);

  // Compute S^{n+1} = MOLRhs(U^{n+1,*})
  if (verbose) {
    amrex::Print() << "... Computing MOL source term at t^{n+1} " << std::endl;
  }
  FillPatch(*this, Sborder, nGrowTr, time + dt, State_Type, 0, NVAR);
  flux_factor = mol_iters > 1 ? 0 : 1;
  getMOLSrcTerm(Sborder, S, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_new
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
    ) {
      construct_new_source(
        src_list[n], time + dt, dt, amr_iteration, amr_ncycle, 0, 0);
      amrex::MultiFab::Saxpy(S, 1.0, *new_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  // U^{n+1.**} = 0.5*(U^n + U^{n+1,*}) + 0.5*dt*S^{n+1} = U^n + 0.5*dt*S^n +
  // 0.5*dt*S^{n+1} + 0.5*dt*I_R
  amrex::MultiFab::LinComb(U_new, 0.5, Sborder, 0, 0.5, U_old, 0, 0, NVAR, 0);
  amrex::MultiFab::Saxpy(
    U_new, 0.5 * dt, S, 0, 0, NVAR,
    0); //  NOTE: If I_R=0, we are done and U_new is the final new-time state

  computeTemp(U_new, 0);
  }

  return dt;
}

void
ERF::construct_Snew(
  amrex::MultiFab& S_new, const amrex::MultiFab& S_old, amrex::Real dt)
{
  int ng = 0;

  amrex::MultiFab::Copy(S_new, S_old, 0, 0, NVAR, ng);
  for (int n = 0; n < src_list.size(); ++n) {
    amrex::MultiFab::Saxpy(
      S_new, 0.5 * dt, *new_sources[src_list[n]], 0, 0, NVAR, ng);
    amrex::MultiFab::Saxpy(
      S_new, 0.5 * dt, *old_sources[src_list[n]], 0, 0, NVAR, ng);
  }
  if (do_hydro) {
    amrex::MultiFab::Saxpy(S_new, dt, hydro_source, 0, 0, NVAR, ng);
  }
}
