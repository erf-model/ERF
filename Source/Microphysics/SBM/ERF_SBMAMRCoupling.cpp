#include "ERF.H"
#include "ERF_SBMAMRFluxRegister.H"
#include "ERF_SBMAMRTime.H"
#include "ERF_SBMAMRTransfer.H"
#include "ERF_SBMRestart.H"
#include "ERF_SBMStateManager.H"
#include "ERF_SBMTransport.H"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

using namespace amrex;

namespace {

void discard_sbm_sync_interval (erf_sbm::SBMAMRFluxRegister& flux_register,
                                const std::string& reason)
{
    std::string discard_diagnostic;
    if (!flux_register.discard_interval(discard_diagnostic)) {
        Abort("SBM rejected synchronization could not discard pending YA amounts: " +
              discard_diagnostic + "; original failure: " + reason);
    }
}

} // namespace

void
ERF::make_sbm_flux_register (int lev)
{
    if (static_cast<int>(sbm_amr_flux_regs.size()) <= lev) {
        sbm_amr_flux_regs.resize(lev + 1);
    }
    if (sbm_amr_flux_regs[lev] && sbm_amr_flux_regs[lev]->has_outstanding_amount()) {
        Abort("SBM interface register was rebuilt with a pending reflux amount");
    }
    if (lev > 0) {
        if (solverChoice.coupling_type == CouplingType::TwoWay &&
            solverChoice.moisture_type == MoistureType::SBM) {
            AMREX_ALWAYS_ASSERT(sbm_state_manager != nullptr);
            // Each interface interval is consumed once per level-0 step.
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev == 1,
                "SBM coarse/fine synchronization runs once per level-0 step and "
                "supports only the level-0/level-1 interface");
            sbm_amr_flux_regs[lev] = std::make_unique<erf_sbm::SBMAMRFluxRegister>();
            sbm_amr_flux_regs[lev]->define(
                grids[lev], grids[lev-1], dmap[lev], dmap[lev-1], geom[lev],
                geom[lev-1], ref_ratio[lev-1], lev,
                sbm_state_manager->layout().ncomp());
        } else {
            sbm_amr_flux_regs[lev].reset();
        }
        bind_sbm_amr_flux_registers(lev - 1);
    }
    bind_sbm_amr_flux_registers(lev);
}

void
ERF::bind_sbm_amr_flux_registers (int lev)
{
    if (!sbm_transport || lev < 0 || lev >= sbm_transport->nlevels() ||
        !sbm_transport->is_defined(lev)) {
        return;
    }
    erf_sbm::SBMAMRFluxRegister* as_fine = nullptr;
    erf_sbm::SBMAMRFluxRegister* as_coarse = nullptr;
    if (lev < static_cast<int>(sbm_amr_flux_regs.size())) {
        as_fine = sbm_amr_flux_regs[lev].get();
    }
    if (lev + 1 < static_cast<int>(sbm_amr_flux_regs.size())) {
        as_coarse = sbm_amr_flux_regs[lev + 1].get();
    }
    sbm_transport->set_amr_flux_registers(lev, as_coarse, as_fine);
}

// Build the fine spectral current state only after ERF has finalized the new
// host density, geometry, map factors and Jacobian. The manager's existing
// new-state buffer is the target-time view; the private candidate is
// temporary and is published only after carrier-relative admission.
void
ERF::init_sbm_level_from_coarse (int lev, Real time, const BoxArray& ba,
                                 const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(sbm_state_manager != nullptr &&
                        sbm_transport != nullptr && lev > 0);
    // AmrCore passes the creation time as Real; the parent level's time is exact.
    const double fine_time = t_new[lev - 1];
    AMREX_ALWAYS_ASSERT(std::abs(static_cast<double>(time) - fine_time) <=
                        16.0 * std::numeric_limits<Real>::epsilon() *
                            std::max(1.0, std::abs(fine_time)));
    sbm_state_manager->define(lev, ba, dm, fine_time);
    sbm_transport->define(lev, ba, dm);
    bind_sbm_amr_flux_registers(lev);
    bind_sbm_amr_flux_registers(lev - 1);
    AMREX_ALWAYS_ASSERT(detJ_cc[lev] != nullptr &&
                        mapfac[lev][MapFacType::m_x] != nullptr &&
                        mapfac[lev][MapFacType::m_y] != nullptr);
    std::string diagnostic;
    if (!sbm_transport->rebuild_static_measure(
            lev, *detJ_cc[lev], *mapfac[lev][MapFacType::m_x],
            *mapfac[lev][MapFacType::m_y], diagnostic)) {
        Abort("SBM fine-level mapped measure: " + diagnostic);
    }
    AMREX_ALWAYS_ASSERT(sbm_transport->measure_is_ready(lev - 1));
    const double coarse_time = sbm_state_manager->new_time(lev - 1);
    if (!erf_sbm::semantic_time_matches(coarse_time, fine_time, dt[lev - 1])) {
        Abort("SBM fine-level initialization: the coarse spectrum and host "
              "state times differ");
    }
    const auto& coarse_measure = sbm_transport->static_measure(lev - 1);
    const auto& fine_measure = sbm_transport->static_measure(lev);
    const auto& coarse_spectrum = sbm_state_manager->new_state(lev - 1);
    const auto& fine_target = sbm_state_manager->new_state(lev);
    erf_sbm::SBMAMRStateView coarse_view{
        &coarse_spectrum, coarse_time,
        &vars_new[lev - 1][Vars::cons], Rho_comp, coarse_time,
        &coarse_measure, 0, coarse_time};
    erf_sbm::SBMAMRStateView fine_view{
        &fine_target, fine_time,
        &vars_new[lev][Vars::cons], Rho_comp, fine_time,
        &fine_measure, 0, fine_time};
    MultiFab candidate(ba, dm, sbm_state_manager->layout().ncomp(), 0);
    if (!erf_sbm::ProlongCarrierRelativeSpectrumPrivate(
            sbm_state_manager->layout(), coarse_view, fine_view,
            geom[lev - 1], geom[lev], ref_ratio[lev - 1], lev,
            candidate, diagnostic, sbm_transport->max_groups_per_chunk(), false,
            &sbm_transport->admission())) {
        Abort("SBM fine-level carrier-relative initialization: " + diagnostic);
    }
    if (!sbm_state_manager->accept_synchronized_state(
            lev, candidate, fine_time, diagnostic)) {
        Abort("SBM fine-level candidate publication: " + diagnostic);
    }
    sbm_state_manager->project_to_core(
        lev, vars_new[lev][Vars::cons], solverChoice.moisture_indices.qc,
        solverChoice.moisture_indices.qr);
}

// Remake a level's spectrum on new grids at the synchronized regrid instant,
// after ERF has finalized the level's host density, geometry, map factors
// and Jacobian there. Retained cells keep their accepted values; newly
// covered cells are prolonged from the same-time coarse spectrum. Released
// cells need nothing: the last synchronization already restricted them onto
// their parents. The level is replaced only after its whole candidate is
// admitted; a rejection aborts before any spectrum or projection changes.
void
ERF::remake_sbm_level (int lev, Real time, const BoxArray& ba,
                       const DistributionMapping& dm)
{
    AMREX_ALWAYS_ASSERT(sbm_state_manager != nullptr && sbm_transport != nullptr);
    // AmrCore remakes only refined levels; level 0 is re-gridded on restart by reading it.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(lev > 0,
        "SBM regrid supports remaking refined levels only");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        sbm_state_manager->is_defined(lev) && sbm_state_manager->new_valid(lev) &&
        !sbm_state_manager->step_active(lev),
        "SBM regrid requires a synchronized accepted spectrum on the remade level");
    AMREX_ALWAYS_ASSERT(detJ_cc[lev] != nullptr &&
                        mapfac[lev][MapFacType::m_x] != nullptr &&
                        mapfac[lev][MapFacType::m_y] != nullptr);
    const auto& layout = sbm_state_manager->layout();
    const int ncomp = layout.ncomp();
    // AmrCore passes the regrid time as Real; the level's own time is exact.
    const double regrid_time = t_new[lev];
    AMREX_ALWAYS_ASSERT(std::abs(static_cast<double>(time) - regrid_time) <=
                        16.0 * std::numeric_limits<Real>::epsilon() *
                            std::max(1.0, std::abs(regrid_time)));
    if (!erf_sbm::level_time_matches(sbm_state_manager->new_time(lev), regrid_time,
                                     dt[lev])) {
        Abort("SBM regrid: the remade level's accepted spectrum is not at the regrid instant");
    }
    const double coarse_time = sbm_state_manager->new_time(lev - 1);
    if (!erf_sbm::level_time_matches(coarse_time, regrid_time, dt[lev - 1]) ||
        !sbm_transport->measure_is_ready(lev - 1)) {
        Abort("SBM regrid: the coarse spectrum is not synchronized at the regrid instant");
    }

    const MultiFab& retained = sbm_state_manager->new_state(lev);
#ifdef ERF_SBM_AMR_QUALIFICATION
    const auto inventory_before = sbm_qualification_composite_inventory();
    const auto counts = erf_sbm::CountRegridCells(retained.boxArray(), ba);
    std::vector<std::uint64_t> phase_bytes(4, 0);
    phase_bytes[0] = static_cast<std::uint64_t>(TotalBytesAllocatedInFabs());
#endif
    // The transport level holds only scratch and the measure, so it moves to
    // the new grids first and its measure and admission check the candidate.
    std::string diagnostic;
    sbm_transport->destroy(lev);
    sbm_transport->define(lev, ba, dm);
    if (!sbm_transport->rebuild_static_measure(
            lev, *detJ_cc[lev], *mapfac[lev][MapFacType::m_x],
            *mapfac[lev][MapFacType::m_y], diagnostic)) {
        Abort("SBM regrid mapped measure: " + diagnostic);
    }
    MultiFab candidate(ba, dm, ncomp, retained.nGrowVect());
    candidate.setVal(Real(0.0));
    {
        erf_sbm::SBMTransferAllocationCensus census;

        const MultiFab* retained_source = &retained;
#ifdef ERF_SBM_AMR_QUALIFICATION
        MultiFab faulted_retained;
        if (solverChoice.sbm_amr_qualification_regrid_reject_fault &&
            counts.retained > 0) {
            // Make one cell that survives the regrid inadmissible.
            faulted_retained.define(retained.boxArray(), retained.DistributionMap(),
                                    ncomp, 0);
            MultiFab::Copy(faulted_retained, retained, 0, 0, ncomp, 0);
            for (int box = 0; box < ba.size(); ++box) {
                const auto overlaps = retained.boxArray().intersections(ba[box]);
                if (!overlaps.empty()) {
                    const IntVect cell = overlaps.front().second.smallEnd();
                    faulted_retained.setVal(Real(-1.0), Box(cell, cell), 0, 1, 0);
                    break;
                }
            }
            retained_source = &faulted_retained;
        }
#endif
        const erf_sbm::SBMAMRStateView coarse_view{
            &sbm_state_manager->new_state(lev - 1), regrid_time,
            &vars_new[lev - 1][Vars::cons], Rho_comp, regrid_time,
            &sbm_transport->static_measure(lev - 1), 0, regrid_time};
        if (!erf_sbm::RemapRegriddedSpectrumPrivate(
                layout, coarse_view, *retained_source,
                vars_new[lev][Vars::cons], Rho_comp,
                sbm_transport->static_measure(lev), 0,
                regrid_time, geom[lev - 1], geom[lev], ref_ratio[lev - 1],
                lev, candidate, diagnostic,
                sbm_transport->max_groups_per_chunk(),
                &sbm_transport->admission(), &sbm_transport->admission(), &census)) {
            Abort("SBM regrid candidate rejected before publication: " + diagnostic);
        }
#ifdef ERF_SBM_AMR_QUALIFICATION
        phase_bytes[1] = census.peak_bytes;
#endif
        if (!sbm_state_manager->can_replace_level(lev, candidate, regrid_time,
                                                  diagnostic)) {
            Abort("SBM regrid candidate could not replace the level: " + diagnostic);
        }
    }
#ifdef ERF_SBM_AMR_QUALIFICATION
    phase_bytes[2] = static_cast<std::uint64_t>(TotalBytesAllocatedInFabs());
#endif

    sbm_state_manager->replace_level(lev, std::move(candidate), regrid_time);
    bind_sbm_amr_flux_registers(lev);
    bind_sbm_amr_flux_registers(lev - 1);
    if (lev < finest_level) { bind_sbm_amr_flux_registers(lev + 1); }
    sbm_state_manager->project_to_core(
        lev, vars_new[lev][Vars::cons], solverChoice.moisture_indices.qc,
        solverChoice.moisture_indices.qr);
#ifdef ERF_SBM_AMR_QUALIFICATION
    phase_bytes[3] = static_cast<std::uint64_t>(TotalBytesAllocatedInFabs());
    sbm_qualification_check_regrid(lev, counts.cells, counts.retained,
                                   counts.prolonged, counts.released,
                                   inventory_before, phase_bytes);
#endif
}

// Restrict, reflux and publish each coarse level after its fine level has
// completed the coarse interval. The synchronization candidate is private
// until restriction, YA reflux, and full canonical admission all succeed. It
// remains in U representation except at uncovered interface cells with a
// nonzero YA correction. Direct stage-rate additions remain pending until
// admission.
void
ERF::sbm_synchronize_levels (int nstep)
{
    AMREX_ALWAYS_ASSERT(sbm_state_manager != nullptr && sbm_transport != nullptr);
#ifdef ERF_SBM_AMR_QUALIFICATION
    Real qualification_correction_linf = Real(0.0);
    std::uint64_t qualification_reflux_calls = 0;
    std::uint64_t qualification_reflux_syncs = 0;
#else
    amrex::ignore_unused(nstep);
#endif
    const auto& layout = sbm_state_manager->layout();
    for (int lev = finest_level - 1; lev >= 0; --lev) {
        const int interface_index = lev + 1;
        if (interface_index >= static_cast<int>(sbm_amr_flux_regs.size()) ||
            !sbm_amr_flux_regs[interface_index] ||
            !sbm_transport->measure_is_ready(lev) ||
            !sbm_transport->measure_is_ready(lev + 1)) {
            Abort("SBM AMR synchronization is missing its interface register or static measure");
        }
        auto& flux_register = *sbm_amr_flux_regs[interface_index];
        const double coarse_time = sbm_state_manager->new_time(lev);
        const double fine_time = sbm_state_manager->new_time(lev + 1);
        const double interval_dt_min = flux_register.minimum_fine_dt();
        if (!erf_sbm::semantic_time_matches(coarse_time, t_new[lev], interval_dt_min) ||
            !erf_sbm::semantic_time_matches(fine_time, t_new[lev + 1], interval_dt_min)) {
            discard_sbm_sync_interval(flux_register,
                "SBM spectral and host state times differ");
            Abort("SBM synchronization: the spectral and host state times differ");
        }
        const auto& coarse_measure = sbm_transport->static_measure(lev);
        const auto& fine_measure = sbm_transport->static_measure(lev + 1);
        erf_sbm::SBMAMRStateView coarse_view{
            &sbm_state_manager->new_state(lev), coarse_time,
            &vars_new[lev][Vars::cons], Rho_comp, coarse_time,
            &coarse_measure, 0, coarse_time};
        erf_sbm::SBMAMRStateView fine_view{
            &sbm_state_manager->new_state(lev + 1), fine_time,
            &vars_new[lev + 1][Vars::cons], Rho_comp, fine_time,
            &fine_measure, 0, fine_time};
        MultiFab candidate(grids[lev], dmap[lev], layout.ncomp(), 0);
        const MultiFab& covered = flux_register.covered_parent_mask();
        std::string diagnostic;
        if (!erf_sbm::RestrictMappedSpectrumPrivate(
                layout, fine_view, coarse_view, ref_ratio[lev], lev,
                candidate, diagnostic, sbm_transport->max_groups_per_chunk(),
                false, &sbm_transport->admission(), nullptr,
                interval_dt_min, &covered)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM mapped coarse restriction rejected: " + diagnostic);
        }

        const auto& chunks = sbm_transport->chunks();
        const int max_chunk_components = sbm_transport->max_chunk_components();
        if (max_chunk_components <= 0) {
            discard_sbm_sync_interval(flux_register,
                "SBM interface has no complete constraint-group chunks");
            Abort("SBM interface has no complete constraint-group chunks");
        }
        MultiFab correction(grids[lev], dmap[lev], max_chunk_components, 0);
        MultiFab correction_invalid(grids[lev], dmap[lev], 1, 0);
        correction_invalid.setVal(Real(0.0));
        for (const auto& chunk : chunks) {
            // Sort only this complete semantic chunk. This preserves the
            // group admission boundary while coalescing YA's contiguous
            // global-component source runs.
            std::vector<int> reflux_components = chunk.components;
            std::sort(reflux_components.begin(), reflux_components.end());
            correction.setVal(Real(0.0));
            if (!flux_register.reflux_components(
                    correction, reflux_components, diagnostic)) {
                discard_sbm_sync_interval(flux_register, diagnostic);
                Abort("SBM YA mapped correction rejected: " + diagnostic);
            }
#ifdef ERF_SBM_AMR_QUALIFICATION
            qualification_correction_linf = std::max(
                qualification_correction_linf,
                correction.norm0(0, max_chunk_components, IntVect(0), false));
#endif
            if (!erf_sbm::ApplyMappedRefluxCorrectionChunk(
                candidate, correction, coarse_measure, covered,
                reflux_components, correction_invalid, diagnostic)) {
                discard_sbm_sync_interval(flux_register, diagnostic);
                Abort("SBM YA correction layout rejected: " + diagnostic);
            }
        }
        if (!erf_sbm::MappedRefluxCorrectionMaskIsValid(
                correction_invalid, diagnostic)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM localized YA correction rejected: " + diagnostic);
        }
        if (!sbm_transport->admission().admissible(
                candidate, lev, diagnostic)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM refluxed coarse candidate rejected before publication: " +
                  diagnostic);
        }
        if (!sbm_state_manager->can_accept_synchronized_state(
                lev, candidate, coarse_time, diagnostic)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM refluxed coarse state could not be published: " + diagnostic);
        }
        if (!flux_register.can_accept_reflux(diagnostic)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM refluxed coarse candidate could not consume its pending amount: " +
                  diagnostic);
        }
#ifdef ERF_SBM_AMR_QUALIFICATION
        qualification_reflux_calls += flux_register.reflux_call_count();
        qualification_reflux_syncs += flux_register.reflux_sync_count();
#endif
        if (!flux_register.accept_reflux(diagnostic)) {
            discard_sbm_sync_interval(flux_register, diagnostic);
            Abort("SBM admitted coarse candidate could not consume its pending amount: " +
                  diagnostic);
        }
        sbm_state_manager->publish_synchronized_state(lev, candidate);
        sbm_state_manager->project_to_core(
            lev, vars_new[lev][Vars::cons], solverChoice.moisture_indices.qc,
            solverChoice.moisture_indices.qr);
    }
#ifdef ERF_SBM_AMR_QUALIFICATION
    sbm_qualification_check_step(nstep, qualification_correction_linf,
                                 qualification_reflux_calls,
                                 qualification_reflux_syncs);
#endif
}
