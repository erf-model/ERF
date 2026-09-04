/**
 * \file ERF_IBSEB.cpp
 * \brief ERF-side hooks of the immersed-boundary surface energy balance.
 *
 * The balance lives in IBFaceSet (one per level, ``ERF::m_ibseb``); this file
 * holds the three places ERF calls into it:
 *  - init_ibseb() from ERF::InitData_post(), after the immersed forcing has
 *    built the blanking on a fresh start or a restart;
 *  - ibseb_write_checkpoint() from ERF::WriteCheckpointFile(), per level;
 *  - ibseb_report() from ERF::post_timestep().
 * The inputs are parsed in ERF::ReadParameters() into ``ERF::ibseb_params``.
 * Every function is a no-op unless ``erf.ibseb.enable`` is set.
 */
#include <ERF.H>
#include <AMReX_VisMF.H>

using namespace amrex;

/**
 * Build the face set of every level from the blanking and, on a restart,
 * refill its state from the checkpoint.
 *
 * Called from ERF::InitData_post() after restart(), which is the first point
 * where both paths (fresh start and restart) have the blanking of every
 * level built and ghost-filled. The face list is always rebuilt from the
 * blanking rather than read back, so the checkpoint carries only the state
 * (``IBSEBState``, see IBFaceSet::state_ncomp()) and a restart on a different
 * number of ranks works. A checkpoint from a run without the balance has no
 * such field; the initial state is kept and a note is printed.
 */
void
ERF::init_ibseb ()
{
    if (!ibseb_params.enable) { return; }
    if (solverChoice.buildings_type != BuildingsType::ImmersedForcing &&
        solverChoice.terrain_type   != TerrainType::ImmersedForcing) {
        Abort("erf.ibseb.enable needs erf.buildings_type = ImmersedForcing (or erf.terrain_type = ImmersedForcing)");
    }
    // The immersed forcing's own surface-temperature conditions would fight
    // the face balance for the same cells.
    if (solverChoice.if_init_surf_temp > 0.0 ||
        solverChoice.if_surf_temp_flux != Real(1.e-8) ||
        solverChoice.if_Olen_in != Real(1.e-8)) {
        Abort("erf.ibseb.enable: remove erf.if_init_surf_temp, erf.if_surf_temp_flux and erf.if_Olen; "
              "the face balance sets the temperature condition at the buildings");
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
        m_ibseb[lev]->compute_view_fractions();
        m_ibseb[lev]->compute_shortwave(t_new[lev]);
        m_ibseb[lev]->compute_longwave(vars_new[lev][Vars::cons]);
        m_ibseb[lev]->compute_sensible(vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel],
                                       vars_new[lev][Vars::yvel], vars_new[lev][Vars::zvel], solverChoice.c_p);
        m_ibseb[lev]->report(t_new[lev], istep[lev], ibseb_params.csv_int > 0);
    }
}

/**
 * Per-step update of one level, called at the start of ERF::Advance() with
 * the state at the start of the step: shortwave, longwave and the wall
 * function on the faces. The sensible flux computed here is what
 * add_heat_flux_to_source() deposits at every slow stage of the step.
 */
void
ERF::ibseb_advance (int lev, Real time, const MultiFab& cons,
                    const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel)
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    m_ibseb[lev]->compute_shortwave(time);
    m_ibseb[lev]->compute_longwave(cons);
    m_ibseb[lev]->compute_sensible(cons, xvel, yvel, zvel, solverChoice.c_p);
}

/**
 * Write the face state of one level into the checkpoint as ``IBSEBState``.
 * Called inside the level loop of ERF::WriteCheckpointFile(); a no-op unless
 * the balance is on and the level has a face set.
 */
void
ERF::ibseb_write_checkpoint (const std::string& checkpointname, int lev) const
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    MultiFab state(grids[lev], dmap[lev], m_ibseb[lev]->state_ncomp(), 0);
    m_ibseb[lev]->save_state(state);
    VisMF::Write(state, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "IBSEBState"));
}

/**
 * Periodic report from ERF::post_timestep(): every ``erf.ibseb.csv_int``
 * steps, print the summary of each level and append its CSV rows. A
 * non-positive interval disables both.
 */
void
ERF::ibseb_report (int nstep, Real time)
{
    if (!ibseb_params.enable) { return; }
    // With debug on the summary is printed every step, as the fire module
    // does; the CSV rows keep their interval.
    const bool csv_now = (ibseb_params.csv_int > 0) && (nstep % ibseb_params.csv_int == 0);
    if (!csv_now && !ibseb_params.debug) { return; }
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (m_ibseb[lev]) { m_ibseb[lev]->report(time, nstep, csv_now); }
    }
}
