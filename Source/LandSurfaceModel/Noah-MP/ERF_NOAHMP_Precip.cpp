/*
 * Precip bookkeeping for the ERF Noah-MP interface: collects the typed accumulator
 * views, maintains the per-level previous-cumulative snapshots (so interval precip
 * is a clean delta, robust across restarts), and reports guard actions. The actual
 * staging lives in stage_forcing() (ERF_NOAHMP_Advance.cpp).
 */

#include <memory>

#include <AMReX_Print.H>

#include <ERF_NOAHMP.H>
#include <ERF_Constants.H>

using namespace amrex;

// Collect the typed slots (NoahmpPrecipSlot order) into a fixed-index bundle so the
// bookkeeping and delta kernel can loop generically. Missing slots keep a null accum.
erf_noahmp::PrecipSlots
NOAHMP::collect_precip_sources (const SurfacePrecipAccumulationSources& precip_sources)
{
    erf_noahmp::PrecipSlots p;
    const SurfacePrecipAccumulationSource slot_src[NoahmpPrecipSlot::NumSlots] = {
        precip_sources.total, precip_sources.rain, precip_sources.snow,
        precip_sources.graupel, precip_sources.hail };
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        p.accum[s]  = slot_src[s].accum;
        p.factor[s] = slot_src[s].native_to_kg_m2;
        p.has[s]    = (slot_src[s].accum != nullptr);
    }
    p.have_any = surface_precip_has_any_source(precip_sources);
    return p;
}

// Lazily allocate the per-level, per-slot previous-cumulative snapshot. Restart seeds
// from the checkpoint (m_precip_accum_restored); cold start seeds from current (delta 0).
void
NOAHMP::prepare_precip_snapshots (const int& lev,
                                  const erf_noahmp::PrecipSlots& precip)
{
    if (!precip.have_any) { return; }

    if (int(m_precip_accum_prev.size()) <= lev) { m_precip_accum_prev.resize(lev+1); }
    auto& prev = m_precip_accum_prev[lev];
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        if (precip.has[s] && prev[s] == nullptr) {
            // Zero ghost cells + setVal(0) so no cell is uninitialized (garbage prev
            // -> ~1e25 deltas -> Noah-MP water-balance abort); only the slab is read.
            prev[s] = std::make_unique<MultiFab>(
                precip.accum[s]->boxArray(), precip.accum[s]->DistributionMap(), 1, 0);
            prev[s]->setVal(0.0);
            if (m_precip_accum_restored[s]) {
                // Restart: ParallelCopy handles a changed decomposition; consume once.
                prev[s]->ParallelCopy(*m_precip_accum_restored[s], 0, 0, 1);
                m_precip_accum_restored[s].reset();
            } else {
                MultiFab::Copy(*prev[s], *precip.accum[s], 0, 0, 1, 0);
            }
        }
    }
}

// Advance the previous-cumulative snapshot to current values now that this call's
// RAINBL/SR are consumed, so the next call differences against this step.
void
NOAHMP::advance_precip_snapshots (const int& lev,
                                  const erf_noahmp::PrecipSlots& precip)
{
    if (!precip.have_any) { return; }
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        if (precip.has[s]) {
            MultiFab::Copy(*m_precip_accum_prev[lev][s], *precip.accum[s], 0, 0, 1, 0);
        }
    }
}

// Report (per rank, never silent) cells the guard clamped or whose frozen components
// were rescaled to restore MP_SNOW+MP_GRAUP <= MP_RAINNC.
void
NOAHMP::report_precip_diagnostics (const int& lev,
                                   const Vector<erf_noahmp::ClampedPrecipCell>& clamped_cells,
                                   const Vector<erf_noahmp::InvariantPrecipCell>& invariant_cells)
{
    if (!clamped_cells.empty()) {
        const int nshow = amrex::min(int(clamped_cells.size()), 50);
        amrex::AllPrint aprint;
        aprint << "WARNING: NOAHMP precip guard (level " << lev << "): clamped "
               << clamped_cells.size() << " non-physical interval-precip cell(s) (> "
               << lsm_max_precip_interval << " mm):\n";
        for (int c(0); c < nshow; ++c) {
            aprint << "    (i,j)=(" << clamped_cells[c].i << "," << clamped_cells[c].j
                   << ") raw=" << clamped_cells[c].raw_mm << " mm\n";
        }
        if (int(clamped_cells.size()) > nshow) {
            aprint << "    ... and " << (int(clamped_cells.size()) - nshow) << " more\n";
        }
    }

    if (!invariant_cells.empty()) {
        const int nshow = amrex::min(int(invariant_cells.size()), 50);
        amrex::AllPrint aprint;
        aprint << "WARNING: NOAHMP precip invariant (level " << lev << "): rescaled frozen "
               << "components in " << invariant_cells.size() << " cell(s) where "
               << "MP_SNOW+MP_GRAUP exceeded MP_RAINNC:\n";
        for (int c(0); c < nshow; ++c) {
            aprint << "    (i,j)=(" << invariant_cells[c].i << "," << invariant_cells[c].j
                   << ") frozen=" << invariant_cells[c].froz_mm
                   << " mm > total=" << invariant_cells[c].total_mm << " mm\n";
        }
        if (int(invariant_cells.size()) > nshow) {
            aprint << "    ... and " << (int(invariant_cells.size()) - nshow) << " more\n";
        }
    }
}
