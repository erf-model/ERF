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
        if (lev > 0) {
            MultiFab& U_new = vars_new[lev][Vars::xvel];
            MultiFab& V_new = vars_new[lev][Vars::yvel];
            MultiFab& W_new = vars_new[lev][Vars::zvel];
            FillPatchFineLevel(lev, time + dt_advance,
                               {&cons, &U_new, &V_new, &W_new},
                               {&cons, &rU_new[lev], &rV_new[lev], &rW_new[lev]},
                               base_state[lev], base_state[lev]);
        } else {
            cons.FillBoundary(geom[lev].periodicity());
        }
        micro->Update_Micro_Vars_Lev(lev, cons);
        micro->Advance(lev, dt_advance, iteration, time, solverChoice, vars_new, z_phys_nd, phys_bc_type);
        micro->Update_State_Vars_Lev(lev, cons, *z_phys_nd[lev]);

        // Sync cons[lev-1] covered cells with the moist state just written
        // to cons[lev].  Without this, the next sub-cycle's FillPatchFineLevel
        // for level lev would pull stale (latent-heat-less) values from
        // cons[lev-1]'s coarse cells via cell-conservative interpolation,
        // causing an artificial outward heat/q flux across the lev/(lev-1)
        // boundary that drains the bubble interior.
        if (lev > 0 && solverChoice.coupling_type == CouplingType::TwoWay &&
            Microphysics::modelType(solverChoice.moisture_type) == MoistureModelType::Lagrangian) {
            AverageDownMoistStateTo(lev - 1);
        }
    }
}
