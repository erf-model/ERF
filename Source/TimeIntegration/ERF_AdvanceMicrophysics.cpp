#include <ERF.H>

using namespace amrex;

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const Real& dt_advance,
                                const int& iteration,
                                const Real& time )
{
    if (solverChoice.moisture_type != MoistureType::None) {
        micro->Set_RealWidth(lev, real_width);
        cons.FillBoundary(geom[lev].periodicity());
        micro->Update_Micro_Vars_Lev(lev, cons);
        if (plot_micro_src) {
            // Copy current state into MF for microphysics source term
            MultiFab::Copy(*micro_src[lev], cons, 0, 0, cons.nComp(), cons.nGrowVect());
        }
        micro->Advance(lev, dt_advance, iteration, time, solverChoice, vars_new, z_phys_nd, phys_bc_type);
        micro->Update_State_Vars_Lev(lev, cons, *z_phys_nd[lev]);
        if (plot_micro_src) {
            // Subtract updated state variables from original state to compute delta
            MultiFab::Xpay(*micro_src[lev], -1.0, cons, 0, 0, cons.nComp(), cons.nGrowVect());
        }
    }
}
