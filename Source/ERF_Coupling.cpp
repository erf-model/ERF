#include <ERF.H>

#include <AMReX_Box.H>
#include <AMReX_MFIter.H>
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

int
bottom_cell_k (const amrex::MultiFab& mf)
{
    return mf.boxArray().minimalBox().smallEnd(2);
}

int
top_cell_k (const amrex::MultiFab& mf)
{
    return mf.boxArray().minimalBox().bigEnd(2);
}

void
copy_plane_to_plane_xy (amrex::MultiFab& dst,
                        int dst_k,
                        const amrex::MultiFab& src,
                        int src_k)
{
    AMREX_ALWAYS_ASSERT(dst.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(src.nComp() >= 1);
    AMREX_ALWAYS_ASSERT(dst.boxArray() == src.boxArray());
    AMREX_ALWAYS_ASSERT(dst.DistributionMap() == src.DistributionMap());

    for (amrex::MFIter mfi(dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = amrex::makeSlab(mfi.validbox(), 2, dst_k);

        auto const& dst_arr = dst.array(mfi);
        auto const& src_arr = src.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            dst_arr(i,j,dst_k) = src_arr(i,j,src_k);
        });
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
    // We intentionally copy only one horizontal slab, derived from the
    // interface face convention:
    // - ocean/atmos interface face index on source: k_face = src.bigEnd(2) + 1
    // - source cell below that face (ocean top cell): src_k = k_face - 1
    // - destination LSM surface-facing slab uses its bottom-most k index
    //   (for current level-0 matched-grid tests this is the interface slab).
    // This encodes the physical alignment note from coupled discussions:
    // ERF k=0 atmospheric cell lies above REMORA k=nz-1 ocean cell.
    // For initial matched-grid tests, require identical horizontal BoxArray/DM.
    if (solverChoice.lsm_type == LandSurfaceType::None) {
        return;
    }

    if (!state.empty() && state[0] != nullptr && lsm.Get_Data_Ptr(0, 0) != nullptr) {
        auto* dst = lsm.Get_Data_Ptr(0, 0);
        const int dst_k = bottom_cell_k(*dst);
        const int src_face_k = top_cell_k(*state[0]) + 1;
        const int src_k = src_face_k - 1;
        copy_plane_to_plane_xy(*dst, dst_k, *state[0], src_k);
    }
}
