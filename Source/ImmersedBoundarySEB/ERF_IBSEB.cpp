/**
 * \file ERF_IBSEB.cpp
 * \brief ERF-side hooks of the immersed-boundary surface energy balance
 */
#include <ERF.H>
#include <AMReX_VisMF.H>

using namespace amrex;

/**
 * Build the face set of every level from the blanking, after the immersed
 * forcing has been initialised (fresh start or restart), and restore the
 * face state from the checkpoint when restarting.
 */
void
ERF::init_ibseb ()
{
    if (!ibseb_params.enable) { return; }
    if (solverChoice.buildings_type != BuildingsType::ImmersedForcing &&
        solverChoice.terrain_type   != TerrainType::ImmersedForcing) {
        Abort("erf.ibseb.enable needs erf.buildings_type = ImmersedForcing (or erf.terrain_type = ImmersedForcing)");
    }
    m_ibseb.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        m_ibseb[lev] = std::make_unique<IBFaceSet>(ibseb_params, lev);
        m_ibseb[lev]->build(*terrain_blanking[lev], geom[lev]);
        if (!restart_chkfile.empty()) {
            const std::string name = MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "IBSEBState");
            if (FileExists(name + "_H")) {
                MultiFab state(grids[lev], dmap[lev], m_ibseb[lev]->state_ncomp(), 0);
                VisMF::Read(state, name);
                m_ibseb[lev]->load_state(state);
                Print() << "[IBSEB] Face state restored from " << restart_chkfile << "\n";
            } else {
                Print() << "[IBSEB] Checkpoint has no IBSEBState; keeping the initial face state.\n";
            }
        }
        m_ibseb[lev]->report(t_new[lev], istep[lev], ibseb_params.csv_int > 0);
    }
}

/// Write the face state of one level into the checkpoint.
void
ERF::ibseb_write_checkpoint (const std::string& checkpointname, int lev) const
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    MultiFab state(grids[lev], dmap[lev], m_ibseb[lev]->state_ncomp(), 0);
    m_ibseb[lev]->save_state(state);
    VisMF::Write(state, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "IBSEBState"));
}

/// Periodic report, called from post_timestep.
void
ERF::ibseb_report (int nstep, Real time)
{
    if (!ibseb_params.enable || ibseb_params.csv_int <= 0) { return; }
    if (nstep % ibseb_params.csv_int != 0) { return; }
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (m_ibseb[lev]) { m_ibseb[lev]->report(time, nstep, true); }
    }
}
