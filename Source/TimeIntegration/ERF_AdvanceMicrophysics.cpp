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

/**
 * Add a level's sedimentation fluxes to a flux register, as the coarse level (as_crse) or
 * the fine one, so that Reflux replaces the flux a coarse cell next to a fine patch saw
 * through their shared z-face with the fine level's. Below a patch that stops short of the
 * ground, that hands the patch's outflow to the coarse cell under it; at a patch top the
 * fine flux is the coarse one (make_sed_inflow), so nothing changes there.
 *
 * The register works with the advective form d(rho q)/dt = -(m^2/detJ) dF_z/dz on the
 * computational grid, so the downward amount D per unit area over sf.dt is registered as
 * F_z = -D / (sf.dt m_x m_y). Sedimentation has no x or y flux.
 */
void
register_sed_fluxes (YAFluxRegister& fr, bool as_crse, const SedimentationFlux& sf,
                     const MultiFab& cons, const MultiFab& mfx, const MultiFab& mfy,
                     const Geometry& geom)
{
    const int ncomp = sf.down->nComp();
    const Real inv_dt = Real(1) / sf.dt;
    for (MFIter mfi(cons); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        FArrayBox fx(surroundingNodes(bx,0), ncomp, The_Async_Arena());
        FArrayBox fy(surroundingNodes(bx,1), ncomp, The_Async_Arena());
        FArrayBox fz(surroundingNodes(bx,2), ncomp, The_Async_Arena());
        fx.setVal<RunOn::Device>(0.0);
        fy.setVal<RunOn::Device>(0.0);

        const auto down = sf.down->const_array(mfi);
        const auto mx = mfx.const_array(mfi);
        const auto my = mfy.const_array(mfi);
        const auto f = fz.array();
        ParallelFor(fz.box(), ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            f(i,j,k,n) = -down(i,j,k,n) * inv_dt / (mx(i,j,0) * my(i,j,0));
        });

        const std::array<FArrayBox const*, AMREX_SPACEDIM> flux{{&fx, &fy, &fz}};
        for (int n = 0; n < ncomp; ++n) {
            if (as_crse) {
                fr.CrseAdd(mfi, flux, geom.CellSize(), sf.dt, n, sf.cons_comp[n], 1, RunOn::Device);
            } else {
                fr.FineAdd(mfi, flux, geom.CellSize(), sf.dt, n, sf.cons_comp[n], 1, RunOn::Device);
            }
        }
        // The fluxes live in temporaries: finish with them before the next box
        Gpu::streamSynchronize();
    }
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

        // Under TwoWay coupling, hand the sedimentation fluxes to the flux registers, so
        // that the coarse level gets the outflow from the bottom of a fine patch that stops
        // short of the ground (see register_sed_fluxes)
        const SedimentationFlux* sed = micro->Get_Sed_Flux(lev);
        if (sed && sed->down && solverChoice.coupling_type == CouplingType::TwoWay && finest_level > 0) {
            if (lev > 0) {
                // Inside the RK stages the calls overlap in time and do not add up to one step
                bool any_above_ground = false;
                for (int ib = 0; ib < grids[lev].size(); ++ib) {
                    any_above_ground = any_above_ground ||
                        (grids[lev][ib].smallEnd(2) > geom[lev].Domain().smallEnd(2));
                }
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!(solverChoice.moisture_tight_coupling && any_above_ground),
                    "erf.moisture_tight_coupling is not supported with refined levels that stop short of the ground");
            }
            if (lev < finest_level) {
                register_sed_fluxes(*getAdvFluxReg(lev+1), true, *sed, cons,
                                    *mapfac[lev][MapFacType::m_x], *mapfac[lev][MapFacType::m_y], geom[lev]);
            }
            if (lev > 0) {
                register_sed_fluxes(*getAdvFluxReg(lev), false, *sed, cons,
                                    *mapfac[lev][MapFacType::m_x], *mapfac[lev][MapFacType::m_y], geom[lev]);
            }
        }
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

void
ERF::average_down_precip_accumulators (int lev)
{
    // Only the sedimenting schemes; their Qmoist variables are exactly their surface
    // precipitation accumulators (Kessler: rain; SAM, Morrison, WSM6, WDM6: rain, snow,
    // graupel and the precipitation total), each holding its running total in the cell
    // on the ground.
    if (!micro || !micro->Get_Sed_Flux(lev+1)) { return; }

    const BoxArray& fba = grids[lev+1];
    const int kf_ground = geom[lev+1].Domain().smallEnd(2);
    const int kc_ground = geom[lev].Domain().smallEnd(2);
    const IntVect rr = ref_ratio[lev];

    // The ground cells under each fine box that stands on the ground, as coarse cells
    BoxList bl;
    Vector<int> fine_box, pmap;
    for (int ib = 0; ib < fba.size(); ++ib) {
        if (fba[ib].smallEnd(2) != kf_ground) { continue; }
        Box s = coarsen(fba[ib], rr);
        s.setRange(2, kc_ground);
        bl.push_back(s);
        fine_box.push_back(ib);
        pmap.push_back(dmap[lev+1][ib]);
    }
    if (bl.isEmpty()) { return; }
    const BoxArray cba(std::move(bl));
    const DistributionMapping cdm(std::move(pmap));
    MultiFab avg(cba, cdm, 1, 0);

    const Real inv_area = Real(1) / static_cast<Real>(rr[0] * rr[1]);
    for (int n = 0; n < micro->Get_Qmoist_Size(lev+1); ++n) {
        const MultiFab* fine = micro->Get_Qmoist_Ptr(lev+1, n);
        MultiFab* crse = micro->Get_Qmoist_Ptr(lev, n);
        if (!fine || !crse) { continue; }
        for (MFIter mfi(avg); mfi.isValid(); ++mfi) {
            const auto f = fine->const_array(fine_box[mfi.index()]);
            const auto a = avg.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int ic, int jc, int kc) noexcept
            {
                Real sum = 0.0;
                for (int j = jc*rr[1]; j < (jc+1)*rr[1]; ++j) {
                    for (int i = ic*rr[0]; i < (ic+1)*rr[0]; ++i) {
                        sum += f(i,j,kf_ground);
                    }
                }
                a(ic,jc,kc) = sum * inv_area;
            });
        }
        crse->ParallelCopy(avg, 0, 0, 1);
    }
}
