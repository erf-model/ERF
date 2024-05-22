#include <WindFarm.H>
#include <Fitch.H>
#include <EWP.H>
using namespace amrex;


void advance_windfarm (int lev,
                       const Geometry& geom,
                       const Real& dt_advance,
                       MultiFab& cons_in,
                       MultiFab& U_old, MultiFab& V_old, MultiFab& W_old,
                       MultiFab& mf_vars_windfarm, const amrex::MultiFab& mf_Nturb,
                       SolverChoice& solver_choice)
{
    if (solver_choice.windfarm_type == WindFarmType::Fitch) {
        fitch_advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                       mf_vars_windfarm, mf_Nturb);
    } else if (solver_choice.windfarm_type == WindFarmType::EWP) {
          ewp_advance(lev, geom, dt_advance, cons_in, U_old, V_old, W_old,
                      mf_vars_windfarm, mf_Nturb);
    }
}
