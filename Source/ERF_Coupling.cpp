#include <ERF.H>

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

namespace {
void
fill_coupling_states (amrex::Vector<amrex::MultiFab*>& states,
                      amrex::Real time,
                      amrex::Real base)
{
    for (int idx = 0; idx < static_cast<int>(states.size()); ++idx) {
        auto* mf = states[idx];
        if (mf != nullptr) {
            mf->setVal(base + static_cast<amrex::Real>(idx + 1) + time);
        }
    }
}
}

void
ERF::PackAtmosphericStates (amrex::Vector<amrex::MultiFab*>& states,
                            amrex::Real time)
{
    fill_coupling_states(states, time, 10.0);
}

void
ERF::ApplyOceanSurfaceState (const amrex::Vector<amrex::MultiFab*>& state,
                             amrex::Real time)
{
    amrex::ignore_unused(time);

    if (solverChoice.lsm_type == LandSurfaceType::None) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        amrex::MultiFab::Copy(*lsm.Get_Data_Ptr(0, 0), *state[0], 0, 0, 1, 0);
    }
}
