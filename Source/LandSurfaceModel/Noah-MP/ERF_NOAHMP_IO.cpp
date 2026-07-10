/*
 * Driver-owned I/O for the ERF Noah-MP interface: the per-step land plotfile
 * and checkpoint/restart of the full Noah-MP prognostic state (plus the precip
 * snapshots). See dev/spec-noahmp-io.md.
 */

#include <string>
#include <memory>

#include <AMReX_VisMF.H>
#include <AMReX_Utility.H>

#include <ERF_NOAHMP.H>

using namespace amrex;

// Per-step NetCDF land output. The cadence is owned by the caller (ERF's time
// loop); this simply fans out over the local boxes. See dev/spec-noahmp-io.md §2.
void
NOAHMP::Plot_Landfile (const int& nstep)
{
    for (NoahmpIO_type &noahmpio : noahmpio_vect) {
        noahmpio.WriteLand(nstep);
    }
}

// Write/read the full NoahMP prognostic state (soil, snow incl. layers, canopy,
// aquifer, albedo history, ...) to/from a NetCDF restart directory so a restart
// reproduces a cold-start trajectory bitwise (issue #3255). Each local block
// writes its tile into the global-domain file collectively.
// See dev/spec-noahmp-io.md §1.
void
NOAHMP::Write_Lsm_Restart (const std::string& dir) const
{
    for (const NoahmpIO_type& noahmpio : noahmpio_vect) {
        const_cast<NoahmpIO_type&>(noahmpio).WriteRestart(dir);
    }

    // Also save the precip snapshots (the cumulative precip Noah-MP last consumed).
    // Restoring these on restart lets the next land call difference against the same
    // baseline, so no precip is dropped or double-counted across the restart.
    if (m_lev < int(m_precip_accum_prev.size())) {
        for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
            if (m_precip_accum_prev[m_lev][s]) {
                amrex::VisMF::Write(*m_precip_accum_prev[m_lev][s],
                                    m_precip_snapshot_name(dir, m_lev, s));
            }
        }
    }
}

void
NOAHMP::Read_Lsm_Restart (const std::string& dir)
{
    for (NoahmpIO_type& noahmpio : noahmpio_vect) {
        noahmpio.ReadRestart(dir);
    }

    // Read the saved precip snapshots into a staging store (the live snapshots are
    // not allocated until the first Advance, which then copies these in). Legacy
    // checkpoints have no such files -> store stays empty -> Advance cold-seeds.
    for (int s(0); s < NoahmpPrecipSlot::NumSlots; ++s) {
        const std::string name = m_precip_snapshot_name(dir, m_lev, s);
        if (amrex::FileExists(name + "_H")) {
            auto mf = std::make_unique<amrex::MultiFab>();
            amrex::VisMF::Read(*mf, name);
            m_precip_accum_restored[s] = std::move(mf);
        }
    }
}
