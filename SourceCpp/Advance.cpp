#include <cmath>

#include "ERF.H"
#include "IndexDefines.H"

amrex::Real
ERF::advance(
  amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle)
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

  amrex::Real dt_new = dt;

  // 
  // The code in brackets below is the old do_mol_advance
  // dt_new = do_mol_advance(time, dt, amr_iteration, amr_ncycle);
  // 
  {
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

  if (do_mol_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(0.0);
  }

  amrex::MultiFab& U_old = get_old_data(State_Type);
  amrex::MultiFab& U_new = get_new_data(State_Type);
  amrex::MultiFab S(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());

  amrex::MultiFab S_old, S_new;
  if (mol_iters > 1) {
    S_old.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
    S_new.define(grids, dmap, NVAR, 0, amrex::MFInfo(), Factory());
  }

  // Compute S^{n} = MOLRhs(U^{n})
  if (verbose) {
    amrex::Print() << "... Computing MOL source term at t^{n} " << std::endl;
  }
  FillPatch(*this, Sborder, nGrowTr, time, State_Type, 0, NVAR);
  amrex::Real flux_factor = 0;
  getMOLSrcTerm(Sborder, S, time, dt, flux_factor);

  // Build other (neither spray nor diffusion) sources at t_old
  for (int n = 0; n < src_list.size(); ++n) {
    if (
      src_list[n] != diff_src
    ) {
      construct_old_source(
        src_list[n], time, dt, amr_iteration, amr_ncycle, 0, 0);
      amrex::MultiFab::Saxpy(S, 1.0, *old_sources[src_list[n]], 0, 0, NVAR, 0);
    }
  }

  // U^* = U^n + dt*S^n
  amrex::MultiFab::LinComb(U_new, 1.0, Sborder, 0, dt, S, 0, 0, NVAR, 0);

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

  return dt;
  }
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
