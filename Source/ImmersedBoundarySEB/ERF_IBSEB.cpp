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
#include <ERF_PlaneAverage.H>
#include <ERF_DirectionSelector.H>
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
        const double t_init0 = ParallelDescriptor::second();
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
        m_ibseb[lev]->assign_materials();
        m_ibseb[lev]->compute_view_fractions();
        m_ibseb[lev]->set_init_cost(ParallelDescriptor::second() - t_init0);
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
 * function on the faces, then either the prognostic balance (which finds
 * the skin temperature at the end of the step and advances the slab with
 * it) or, with ``erf.ibseb.prognostic = false``, the slab alone under the
 * fixed skin. The sensible flux left in the set is what
 * add_heat_flux_to_source() deposits at every slow stage of the step, so
 * the air receives exactly the H of the closed balance. The atmosphere is
 * seen at the start of the step and the skin is implicit within it, the
 * usual coupling of a land-surface model.
 */
void
ERF::ibseb_advance (int lev, Real time, Real dt, const MultiFab& cons,
                    const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel)
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    const double t_wall0 = ParallelDescriptor::second();
    m_ibseb[lev]->compute_shortwave(time);
    m_ibseb[lev]->compute_longwave(cons);
    // Phase 8: the ground surface layer's fields and the mixed-layer depth
    // for the wall function beyond neutral (all null / zero unless asked).
    const MultiFab* olen2d = nullptr;
    const MultiFab* pblh2d = nullptr;
    Real z_i_bulk = 0.0;
    if (m_SurfaceLayer && ibseb_params.stability_correction) { olen2d = m_SurfaceLayer->get_olen(lev); }
    if (ibseb_params.convective_velocity == "deardorff") {
        if (m_SurfaceLayer && m_SurfaceLayer->computes_pblh() && ibseb_params.z_i_mode == "pblh") {
            pblh2d = m_SurfaceLayer->get_pblh(lev);
        }
        z_i_bulk = (ibseb_params.z_i_mode == "fixed") ? ibseb_params.z_i
                                                     : ibseb_bulk_richardson_height(lev, cons, xvel, yvel);
        if (ibseb_params.debug) {
            Print() << "[IBSEB DEBUG] lev=" << lev << " mixed-layer depth for w*: " << z_i_bulk << " m ("
                    << ibseb_params.z_i_mode << (pblh2d ? ", pblh per column" : "") << ")\n";
        }
    }
    m_ibseb[lev]->compute_sensible(cons, xvel, yvel, zvel, solverChoice.c_p, olen2d, pblh2d, z_i_bulk);
    if (ibseb_params.prognostic) {
        m_ibseb[lev]->solve_balance(dt);
    } else {
        m_ibseb[lev]->compute_ground(dt);
    }
    m_ibseb[lev]->add_cost(ParallelDescriptor::second() - t_wall0);
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

/**
 * Mixed-layer depth of a level by the bulk Richardson method on the
 * horizontal-mean profile (Troen and Mahrt; Vogelezang and Holtslag).
 *
 * With the first level as the reference, ``Ri_b(z) = g (z - z_1) (theta(z)
 * - theta_1) / (theta_1 (|U(z) - U_1|^2 + 100 u*^2))`` with u* = 0.1 m/s,
 * and the depth is the first cell centre where it exceeds
 * ``erf.ibseb.ri_crit``, or the domain top when it never does (a neutral
 * profile). The profile is the plane average of the conserved state and the
 * face velocities, uniform vertical spacing assumed as elsewhere in the
 * balance; called once per step and level when the convective velocity
 * scale is on and z_i is not fixed, also as the fallback of the pblh mode.
 */
Real
ERF::ibseb_bulk_richardson_height (int lev, const MultiFab& cons, const MultiFab& xvel, const MultiFab& yvel)
{
    MultiFab c2(cons, make_alias, Rho_comp, 2);   // rho and rho theta are the first two components
    PlaneAverage r_ave(&c2, geom[lev], 2);
    r_ave.compute_averages(ZDir(), r_ave.field());
    PlaneAverage u_ave(&xvel, geom[lev], 2);
    u_ave.compute_averages(ZDir(), u_ave.field());
    PlaneAverage v_ave(&yvel, geom[lev], 2);
    v_ave.compute_averages(ZDir(), v_ave.field());
    const int nz = r_ave.ncell_line();
    Gpu::HostVector<Real> rho(nz), rth(nz), uu(u_ave.ncell_line()), vv(v_ave.ncell_line());
    r_ave.line_average(0, rho);
    r_ave.line_average(1, rth);
    u_ave.line_average(0, uu);
    v_ave.line_average(0, vv);
    const Real dz = geom[lev].CellSize(2);
    const Real z_top = geom[lev].ProbHi(2);
    const Real th1 = rth[0] / rho[0];
    const Real U1  = std::sqrt(uu[0] * uu[0] + vv[0] * vv[0]);
    const Real ustar_floor2 = 100.0 * 0.1 * 0.1;
    Real z_i = z_top;
    for (int k = 1; k < nz; ++k) {
        const Real th = rth[k] / rho[k];
        const Real U  = std::sqrt(uu[k] * uu[k] + vv[k] * vv[k]);
        const Real dU = U - U1;
        const Real rib = CONST_GRAV * (k * dz) * (th - th1) / (th1 * (dU * dU + ustar_floor2));
        if (rib > ibseb_params.ri_crit) { z_i = geom[lev].ProbLo(2) + (k + 0.5) * dz - geom[lev].ProbLo(2); break; }
    }
    return z_i;
}
