#include <ERF.H>
#include "ERF_SBMStateManager.H"

using namespace amrex;

namespace {

/**
 * The inflow through the top of a fine level's patches for one Advance of length dt_fine
 * (see ERF_SedimentationFlux.H), from the sedimentation flux crse its parent recorded over
 * its own, enclosing, step. A patch-top face is a z-face at the top of a fine box with no
 * fine cell above it and below the domain top; it gets the coarse flux through the coarse
 * face it lies in, piecewise constant in x and y, in proportion dt_fine / crse.dt.
 *
 * Returns false, leaving inflow alone, when every box reaches the domain top.
 */
bool
make_sed_inflow (const BoxArray& ba, const DistributionMapping& dm, const Box& domain,
                 const IntVect& rr, const SedimentationFlux& crse, Real dt_fine,
                 MultiFab& inflow)
{
    bool any_short = false;
    for (int ib = 0; ib < ba.size() && !any_short; ++ib) {
        any_short = (ba[ib].bigEnd(2) < domain.bigEnd(2));
    }
    if (!any_short) { return false; }

    AMREX_ALWAYS_ASSERT(crse.down && crse.dt > Real(0));
    const int ncomp = crse.down->nComp();

    // The coarse faces over the top face of each fine box, on the fine distribution
    BoxList bl(IndexType(IntVect(0,0,1)));
    for (int ib = 0; ib < ba.size(); ++ib) {
        const Box& b = ba[ib];
        const Box top(IntVect(b.smallEnd(0), b.smallEnd(1), b.bigEnd(2)+1),
                      IntVect(b.bigEnd(0),   b.bigEnd(1),   b.bigEnd(2)+1),
                      IndexType(IntVect(0,0,1)));
        bl.push_back(coarsen(top, rr));
    }
    MultiFab crse_top(BoxArray(std::move(bl)), dm, ncomp, 0);
    crse_top.ParallelCopy(*crse.down);

    inflow.define(convert(ba, IntVect(0,0,1)), dm, ncomp+1, 0);
    inflow.setVal(0.0);

    const Real frac = dt_fine / crse.dt;
    for (MFIter mfi(inflow); mfi.isValid(); ++mfi) {
        const Box& b = ba[mfi.index()];
        if (b.bigEnd(2) >= domain.bigEnd(2)) { continue; }

        Box above(b);
        above.setRange(2, b.bigEnd(2)+1);
        const auto ct = crse_top.const_array(mfi);
        const auto in = inflow.array(mfi);
        // Cells above this box that no fine box covers: their bottom faces are patch tops
        for (const Box& ub : ba.complementIn(above)) {
            ParallelFor(ub, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int ic = amrex::coarsen(i, rr[0]);
                const int jc = amrex::coarsen(j, rr[1]);
                const int kc = amrex::coarsen(k, rr[2]);
                for (int n = 0; n < ncomp; ++n) {
                    in(i,j,k,n) = ct(ic,jc,kc,n) * frac;
                }
                in(i,j,k,ncomp) = Real(1);
            });
        }
    }
    return true;
}

} // namespace

void ERF::advance_microphysics (int lev,
                                MultiFab& cons,
                                const double& dt_advance,
                                const int& iteration,
                                const double& time )
{
    if (solverChoice.moisture_type != MoistureType::None) {
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
        const auto micro_state_before = cloud_chamber_budget ?
            cloud_chamber_budget->state_integrals(cons, geom[lev]) :
            std::array<Real, CloudChamberBudget::NumScalars>{};
        const amrex::MultiFab* base = solverChoice.anelastic[lev] ? &base_state[lev] : nullptr;
        micro->Update_Micro_Vars_Lev(lev, cons, base);

        // A fine patch short of the domain top receives its parent's sedimentation flux
        // through its top face. The parent has already advanced over the step that
        // encloses this one, so its recorded flux is the one to use.
        MultiFab sed_inflow;
        const SedimentationFlux* crse_sed = (lev > 0) ? micro->Get_Sed_Flux(lev-1) : nullptr;
        const bool has_inflow = crse_sed && crse_sed->down &&
            make_sed_inflow(grids[lev], dmap[lev], geom[lev].Domain(), ref_ratio[lev-1],
                            *crse_sed, static_cast<Real>(dt_advance), sed_inflow);
        if (has_inflow) {
            // Inside the RK stages the calls overlap in time and do not add up to one step
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!solverChoice.moisture_tight_coupling,
                "erf.moisture_tight_coupling is not supported with refined levels that stop short of the domain top");
            micro->Set_Sed_Inflow(lev, &sed_inflow);
        }

        micro->Advance(lev, static_cast<Real>(dt_advance), iteration,
                       static_cast<Real>(time), solverChoice, vars_new, z_phys_nd, phys_bc_type);

        if (has_inflow) { micro->Set_Sed_Inflow(lev, nullptr); }
        micro->Update_State_Vars_Lev(lev, cons, *z_phys_nd[lev]);
        if (sbm_state_manager) {
            // The adapter is intentionally a no-op for compact liquid state;
            // refresh the compatibility lanes only from authoritative bins.
            sbm_state_manager->project_to_core(lev, cons,
                                               solverChoice.moisture_indices.qc,
                                               solverChoice.moisture_indices.qr);
        }
        if (cloud_chamber_budget) {
            cloud_chamber_budget->record_internal_source(
                micro_state_before, cloud_chamber_budget->state_integrals(cons, geom[lev]));
        }

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
