#include <ERF.H>

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>

/*
  Coupling reference context (implementation-side):

  1) Legacy state-passing contract:
     Warner et al. (2010), COAWST, Fig. 5 / Block B.
     Atmosphere -> ocean fields are passed as states
     (Uwind, Vwind, Patm, RH, Tair, cloud, rain, SWrad, LWrad), while
     ocean -> atmosphere provides SST.

  2) Future direct flux-passing roadmap:
     COAWST's ATM2OCN_FLUXES pathway (documented in COAWST manuals/workshops
     and exercised in Zambon et al., 2014) motivates a direct
     stress/heat/moisture flux mode.
     COAWST code anchors:
       - Master/mct_roms_wrf.h
       - ROMS/Nonlinear/atm2ocn_flux.F
       - ROMS/Nonlinear/bulk_flux.F
     This file currently implements the legacy state-passing test path only.
*/

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
    // Initial-step-testing example values (deterministic cache verification):
    // at time=t, this writes
    //   states[0] (Uwind) = 11 + t
    //   states[1] (Vwind) = 12 + t
    //   ...
    //   states[8] (LWrad) = 19 + t
    // This is not physical forcing; it is a reproducible exchange check.
    fill_coupling_states(states, time, 10.0);
}

void
ERF::ApplyOceanSurfaceState (const amrex::Vector<amrex::MultiFab*>& state,
                             amrex::Real time)
{
    amrex::ignore_unused(time);

    // Example (legacy state-passing): state[0] carries SST [K].
    // If SST is a constant slab (e.g., 290.5 K), we copy that slab into
    // ERF's active LSM data pointer for ocean-coupled surface application.
    if (solverChoice.lsm_type == LandSurfaceType::None) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        amrex::MultiFab::Copy(*lsm.Get_Data_Ptr(0, 0), *state[0], 0, 0, 1, 0);
    }
}
