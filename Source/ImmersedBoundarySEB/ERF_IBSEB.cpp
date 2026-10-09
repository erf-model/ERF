/**
 * \file ERF_IBSEB.cpp
 * \brief ERF-side hooks of the immersed-boundary surface energy balance.
 *
 * The balance lives in IBFaceSet (one per level, ``ERF::m_ibseb``); this file
 * holds the places ERF calls into it:
 *  - init_ibseb() from ERF::InitData_post(), after the immersed forcing has
 *    built the blanking on a fresh start or a restart, with its start-up
 *    checks ibseb_check_refined_levels(), ibseb_check_sun_matches_two_stream()
 *    and ibseb_check_two_stream_provider();
 *  - ibseb_advance() from ERF::Advance(), per level and step, which with
 *    erf.ibseb.sun_mode = two_stream hands the faces the two-stream sun of
 *    the step first (ibseb_set_two_stream_sun(), as init_ibseb() does). It
 *    runs at the start of the step, or with erf.ibseb.radiation = two_stream
 *    after the step's radiation (ibseb_after_radiation()), whose column sweep
 *    the faces then take;
 *  - ibseb_write_checkpoint() from ERF::WriteCheckpointFile(), per level;
 *  - ibseb_report() from ERF::post_timestep().
 * The inputs are parsed in ERF::ReadParameters() into ``ERF::ibseb_params``.
 * Every function is a no-op unless ``erf.ibseb.enable`` is set.
 */
#include <ERF.H>
#include <ERF_PlaneAverage.H>
#include <ERF_DirectionSelector.H>
#include <AMReX_VisMF.H>
#include <AMReX_MultiFabUtil.H>
#include <ERF_IBSEBSolar.H>
#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

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
    if (solverChoice.buildings_type != BuildingsType::ImmersedForcing) {
        Abort("erf.ibseb.enable needs erf.buildings_type = ImmersedForcing");
    }
    // The face detection takes every solid column of the blanking for a
    // building, so terrain by immersed forcing would be put under the
    // balance as well; it is not supported.
    if (solverChoice.terrain_type == TerrainType::ImmersedForcing) {
        Abort("erf.ibseb.enable does not support erf.terrain_type = ImmersedForcing: "
              "the balance runs on building faces only");
    }
    // The face areas, the face heights, the wall-function distance and the
    // ray cast all take the level's constant cell sizes.
    if (solverChoice.mesh_type != MeshType::ConstantDz) {
        Abort("erf.ibseb.enable needs a uniform vertical grid (no erf.terrain_z_levels or stretched mesh): "
              "the face geometry and the ray cast assume constant dz");
    }
    // The face list is built once here from the blanking; a regrid would
    // leave it indexing the old boxes.
    if (regrid_int > 0) {
        Abort("erf.ibseb.enable does not support regridding (erf.regrid_int > 0): the face list is built once at initialisation");
    }
    // The immersed forcing's own surface-temperature conditions would fight
    // the face balance for the same cells.
    if (solverChoice.if_init_surf_temp > 0.0 ||
        solverChoice.if_surf_temp_flux != Real(1.e-8) ||
        solverChoice.if_Olen_in != Real(1.e-8)) {
        Abort("erf.ibseb.enable: remove erf.if_init_surf_temp, erf.if_surf_temp_flux and erf.if_Olen; "
              "the face balance sets the temperature condition at the buildings");
    }
    ibseb_check_refined_levels();
    ibseb_check_sun_matches_two_stream();
    ibseb_check_two_stream_provider();
    m_ibseb.resize(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        m_ibseb[lev] = std::make_unique<IBFaceSet>(ibseb_params, lev);
        const double t_init0 = ParallelDescriptor::second();
        m_ibseb[lev]->build(*terrain_blanking[lev], geom[lev]);
        if (lev > 0) {
            m_ibseb[lev]->add_outside_occluders(*m_ibseb[lev-1], ref_ratio[lev-1],
                                                amrex::coarsen(grids[lev], ref_ratio[lev-1]));
            m_ibseb[lev]->map_buildings_to_level0(*m_ibseb[lev-1], ref_ratio[lev-1]);
            m_ibseb[lev-1]->release_labels();
        }
        // No level above to map: the labels are no longer needed.
        if (lev == finest_level) { m_ibseb[lev]->release_labels(); }
        // The faces take the two-stream columns at their own heights, so from now on
        // every sweep of the level also keeps the direct beam up to the highest interface
        // a face reads. A level without faces needs none. A level with faces must sweep
        // its own columns: one whose grids do not span the domain in z takes its
        // radiation from the level below by interpolation, with no sweep of its own.
        if (ibseb_params.radiation == "two_stream" && m_ibseb[lev]->has_faces() && rad_level_needs_interpolation(lev)) {
            Abort("erf.ibseb.radiation = two_stream: the grids of level " + std::to_string(lev) + " do not span the"
                  " domain in z, so its two-stream radiation comes from level " + std::to_string(lev - 1) + " by"
                  " interpolation and no column sweep of its own supplies its faces; refine the whole height"
                  " (amr.refine_whole_domain_dir) or use erf.ibseb.radiation = prescribed");
        }
        if (ibseb_params.radiation == "two_stream" && m_ibseb[lev]->has_faces()) {
            two_stream_rad.supply_canopy_forcing(lev, m_ibseb[lev]->top_sample_interface());
        }
        std::unique_ptr<MultiFab> restored;
        if (!restart_chkfile.empty()) {
            const std::string name = MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "IBSEBState");
            if (FileExists(name + "_H")) {
                // The field's width is n_slots x (2 + n_slab_layers) and its
                // boxes follow the buildings; a checkpoint written with
                // another layer count or another building set cannot be
                // unpacked.
                const VisMF header(name);
                const int ncomp_chk = header.nComp();
                if (ncomp_chk != m_ibseb[lev]->state_ncomp()) {
                    Abort("erf.ibseb: IBSEBState in " + restart_chkfile + " has " + std::to_string(ncomp_chk)
                          + " components; the deck sets erf.ibseb.n_slab_layers = " + std::to_string(ibseb_params.n_slab_layers)
                          + ", which with the " + std::to_string(m_ibseb[lev]->n_slots())
                          + " face slots per cell of this blanking needs " + std::to_string(m_ibseb[lev]->state_ncomp())
                          + ". The checkpoint was written with another erf.ibseb.n_slab_layers (restart with the"
                            " checkpoint's value) or for another building set.");
                }
                if (header.boxArray() != m_ibseb[lev]->state_boxarray()) {
                    Abort("erf.ibseb: IBSEBState in " + restart_chkfile + " was written for a different building layout ("
                          + std::to_string(header.boxArray().size()) + " boxes against the "
                          + std::to_string(m_ibseb[lev]->state_boxarray().size()) + " the blanking gives); restart from a checkpoint of the same buildings");
                }
                restored = std::make_unique<MultiFab>(m_ibseb[lev]->make_state());
                VisMF::Read(*restored, name);
                m_ibseb[lev]->load_state(*restored);
                Print() << "[IBSEB] Face state restored from " << restart_chkfile << "\n";
            } else {
                Print() << "[IBSEB] Checkpoint has no IBSEBState; keeping the initial face state.\n";
            }
        }
        m_ibseb[lev]->assign_materials();
        m_ibseb[lev]->compute_view_fractions();
        m_ibseb[lev]->set_init_cost(ParallelDescriptor::second() - t_init0);
        // Initial diagnostics for the first report. On a restart they
        // overwrite the sensible flux with a diagnostic value, so the
        // checkpointed flux (which the convective velocity scale of the next
        // step reads as the previous step's) is put back afterwards. With
        // radiation = two_stream no sweep has run yet, so the initial report
        // carries no shortwave on the faces, nor longwave with lw_mode = two_stream.
        const bool two_stream_faces = (ibseb_params.radiation == "two_stream");
        if (two_stream_faces && m_ibseb[lev]->has_faces()) {
            Print() << "[IBSEB] Level " << lev << ": the faces take their radiation from the two-stream columns,"
                       " each at its own height (up to interface " << m_ibseb[lev]->top_sample_interface() << ")."
                       " The first sweep runs in the first step, so the initial report shows no "
                    << (ibseb_params.lw_mode == "two_stream" ? "radiation" : "shortwave") << " on the faces.\n";
        }
        ibseb_set_two_stream_sun(lev, t_new[lev]);
        m_ibseb[lev]->compute_shortwave(t_new[lev]);
        m_ibseb[lev]->compute_longwave(vars_new[lev][Vars::cons]);
        m_ibseb[lev]->compute_sensible(vars_new[lev][Vars::cons], vars_new[lev][Vars::xvel],
                                       vars_new[lev][Vars::yvel], vars_new[lev][Vars::zvel], solverChoice.c_p);
        if (restored) { m_ibseb[lev]->load_state(*restored); }
        const bool restarting = !restart_chkfile.empty();
        // A restart writes a 3D plotfile for its first step (plot_file_on_restart). With
        // two_stream no sweep has run yet, so the face fields that hold radiation are
        // empty there: the absorbed shortwave and the shadow, and the net longwave with
        // lw_mode = two_stream (the skin temperature and fluxes come from the
        // checkpoint, the face counts and view fractions from the geometry). Say so when
        // that plotfile holds one of them.
        if (two_stream_faces && restarting && plot_file_on_restart && lev == 0) {
            const bool lw_two_stream = (ibseb_params.lw_mode == "two_stream");
            auto has_radiation_fields = [lw_two_stream] (const Vector<std::string>& names) {
                return std::any_of(names.begin(), names.end(), [lw_two_stream] (const std::string& n) {
                    return n == "ibseb_sw_abs" || n == "ibseb_shadow" || (lw_two_stream && n == "ibseb_lw_net");
                });
            };
            const bool plot1 = (m_plot3d_int_1 > 0 || m_plot3d_per_1 > 0.0) && has_radiation_fields(plot3d_var_names_1);
            const bool plot2 = (m_plot3d_int_2 > 0 || m_plot3d_per_2 > 0.0) && has_radiation_fields(plot3d_var_names_2);
            if (plot1 || plot2) {
                Print() << "[IBSEB] Restart: the plotfile written at this step shows no "
                        << (lw_two_stream ? "radiation" : "shortwave") << " on the faces (no sweep has run yet)."
                           " If the run before wrote a plotfile at this step, AMReX kept it, renamed with an"
                           " .old suffix.\n";
            }
        }
        // The step a restart starts from was reported, and its faces dumped, by the run
        // before. Reporting it again would add a second row for the step (and replace its
        // dump) with values that differ: the sun here is that of the step's end, not its
        // start, and with two_stream there is no radiation yet. So a restart prints it only.
        m_ibseb[lev]->report(t_new[lev], istep[lev], ibseb_params.csv_int > 0 && !restarting);
    }
}

/**
 * Abort when a building of the level below crosses the edge of a refined
 * level, or comes within one cell of the level below of it.
 *
 * Each level builds its face list, its building ids and its column map from
 * its own blanking (IBFaceSet::build()), so a building the refined level
 * covers in part would get a partial face list there. A building must
 * therefore lie wholly inside the refined level, where the level resolves it,
 * or wholly outside, where only the levels below hold its faces and the
 * refined level's rays see it through the coarser column map
 * (IBFaceSet::add_outside_occluders()). The check marks the cells of level
 * lev-1 that the coarsened grids of level lev cover (amrex::makeFineMask,
 * periodic images included) and counts the solid cells of level lev-1 whose
 * 3 x 3 x 3 block, inside the domain, holds both covered and uncovered cells:
 * a building crossing the edge, or one so close to it that the fluid cells
 * next to its walls or above its roof straddle it. A building only the
 * refined level resolves is checked by IBFaceSet::build(). A no-op on a
 * single level.
 */
void
ERF::ibseb_check_refined_levels () const
{
    for (int lev = 1; lev <= finest_level; ++lev) {
        const MultiFab& crse = *terrain_blanking[lev-1];
        const iMultiFab covered = makeFineMask(crse.boxArray(), crse.DistributionMap(), IntVect(1),
                                               grids[lev], ref_ratio[lev-1], geom[lev-1].periodicity(), 0, 1);
        const Box& domain = geom[lev-1].Domain();
        const Dim3 dlo = lbound(domain);
        const Dim3 dhi = ubound(domain);
        const bool per_x = geom[lev-1].isPeriodic(0);
        const bool per_y = geom[lev-1].isPeriodic(1);
        // Count the cells and keep the lowest one (as a linear index over the
        // domain) for the message.
        const Long ny_d = domain.length(1), nz_d = domain.length(2);
        ReduceOps<ReduceOpSum, ReduceOpMin> reduce_op;
        ReduceData<int, Long> reduce_data(reduce_op);
        for (MFIter mfi(crse); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto const& blank = crse.const_array(mfi);
            auto const& fine  = covered.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GpuTuple<int, Long>
            {
                constexpr Long none = std::numeric_limits<Long>::max();
                if (blank(i, j, k) < 0.5) { return {0, none}; }
                int n_in = 0, n_out = 0;
                for (int dk = -1; dk <= 1; ++dk) {
                for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    const int ii = i + di, jj = j + dj, kk = k + dk;
                    if (kk < dlo.z || kk > dhi.z) { continue; }
                    if (!per_x && (ii < dlo.x || ii > dhi.x)) { continue; }
                    if (!per_y && (jj < dlo.y || jj > dhi.y)) { continue; }
                    if (fine(ii, jj, kk) == 0) { ++n_out; } else { ++n_in; }
                }}}
                if (n_in == 0 || n_out == 0) { return {0, none}; }
                return {1, (Long(i - dlo.x) * ny_d + (j - dlo.y)) * nz_d + (k - dlo.z)};
            });
        }
        const auto rv = reduce_data.value();
        int n_edge = amrex::get<0>(rv);
        Long first = amrex::get<1>(rv);
        ParallelDescriptor::ReduceIntSum(n_edge);
        ParallelDescriptor::ReduceLongMin(first);
        if (n_edge > 0) {
            const int fk = static_cast<int>(first % nz_d);
            const int fj = static_cast<int>((first / nz_d) % ny_d);
            const int fi = static_cast<int>(first / (nz_d * ny_d));
            const Real* dxc = geom[lev-1].CellSize();
            const Real* plo = geom[lev-1].ProbLo();
            std::ostringstream where;
            where << " (the first at x = " << plo[0] + (fi + 0.5) * dxc[0] << " m, y = "
                  << plo[1] + (fj + 0.5) * dxc[1] << " m, z = " << plo[2] + (fk + 0.5) * dxc[2] << " m)";
            Abort("erf.ibseb: " + std::to_string(n_edge) + " solid cells of the buildings on level "
                  + std::to_string(lev-1) + where.str() + " lie on the edge of level " + std::to_string(lev)
                  + " or within one level-" + std::to_string(lev-1) + " cell of it. Each level builds its"
                    " faces from its own cells, so a building must lie wholly inside a refined level, with at"
                    " least one cell of the level below around it, or wholly outside it. To fix the refined box, run"
                    " `python3 Exec/CanonicalTests/SEB/ibseb_refinement_box.py <inputs>`: it checks the deck's"
                    " erf.<name>.in_box_lo / in_box_hi against the height map and prints a box that no building"
                    " crosses (or `... <inputs> --all --fit tight|relaxed` for a box around every building).");
        }
    }
}

/**
 * Abort unless the faces see the sun the two-stream columns see.
 *
 * Called by init_ibseb(). The prescribed provider places the faces' sun, and
 * nothing ties its own solar formulas (Spencer's, with the equation of time)
 * to the orbital formula of the radiation models (no equation of time), which
 * are up to about 4 degrees apart in hour angle through the year. So with
 * erf.radiation_model = TwoStream following start_datetime the faces must take
 * the two-stream sun, erf.ibseb.sun_mode = two_stream (ibseb_set_two_stream_sun()),
 * and with the two-stream sun fixed at erf.fixed_solar_zenith_angle (a cosine,
 * no azimuth) erf.ibseb.sun_mode = fixed at the same zenith. The two-stream sun
 * needs the start date and one site (erf.rad_cons_lat / lon): a grid with
 * per-column latitude and longitude, which the columns follow in a NetCDF
 * build, is not supported. sun_mode = two_stream without the two-stream
 * radiation stops too. With the two-stream shortwave off there is no column
 * sun to match, and only a two_stream request is checked.
 */
void
ERF::ibseb_check_sun_matches_two_stream () const
{
    const IBSEBParams& p = ibseb_params;
    if (solverChoice.rad_type != RadiationType::TwoStream) {
        if (p.sun_mode == "two_stream") {
            Abort("erf.ibseb.sun_mode = two_stream needs erf.radiation_model = TwoStream; "
                  "use sun_mode = solar (or fixed) for the faces' own sun");
        }
        return;
    }
    const RadChoice& rc = solverChoice.radChoice;
    // Without two-stream shortwave there is no column sun to match; the faces
    // may still take its calendar sun, which needs the date checked below.
    if (!rc.sw_enabled && p.sun_mode != "two_stream") { return; }
    if (rc.sw_enabled && rc.fixed_solar_zenith_angle > 0.0) {
        const Real mu_faces = std::cos(p.sun_zenith_deg * PI / Real(180.0));
        if (p.sun_mode != "fixed" || std::abs(mu_faces - rc.fixed_solar_zenith_angle) > Real(1.e-5)) {
            Abort("erf.ibseb: the two-stream sun is fixed at erf.fixed_solar_zenith_angle = "
                  + std::to_string(rc.fixed_solar_zenith_angle) + " (a cosine), so the faces need"
                    " erf.ibseb.sun_mode = fixed with erf.ibseb.sun_zenith_deg = "
                  + std::to_string(std::acos(rc.fixed_solar_zenith_angle) * Real(180.0) / PI)
                  + "; the deck has sun_mode = " + p.sun_mode + " and sun_zenith_deg = "
                  + std::to_string(p.sun_zenith_deg));
        }
        return;
    }
    if (p.sun_mode != "two_stream") {
        Abort("erf.ibseb: the two-stream sun follows start_datetime, so the faces need erf.ibseb.sun_mode = two_stream"
              " (the two-stream sun has no equation of time, so even a matching sun_mode = solar would sit up to about"
              " 4 degrees off it); the deck has sun_mode = " + p.sun_mode);
    }
    if (!use_datetime) {
        Abort("erf.ibseb.sun_mode = two_stream: the two-stream sun follows the calendar and no start date is known;"
              " set start_datetime = \"YYYY-MM-DD HH:MM:SS\" (UTC)");
    }
    // The sweep follows per-column latitude and longitude only in a NetCDF
    // build (advance_radiation passes lat_m / lon_m there and null otherwise).
#ifdef ERF_USE_NETCDF
    const bool has_latlon = !lat_m.empty() && lat_m[0] && !lon_m.empty() && lon_m[0];
#else
    const bool has_latlon = false;
#endif
    if (has_latlon) {
        Abort("erf.ibseb.sun_mode = two_stream takes one site, erf.rad_cons_lat / lon, but this grid carries a"
              " latitude and longitude per column, which the two-stream columns follow instead");
    }
}

/**
 * Abort unless the two-stream columns can supply the faces of every level
 * (erf.ibseb.radiation = two_stream; a no-op otherwise).
 *
 * Called by init_ibseb(). The faces of a level take the canopy forcing its own
 * column sweep writes in the same step (TwoStreamRadiation::canopy_forcing()),
 * so the radiation model must be the two-stream one with its shortwave on (off,
 * the sweep places no sun, and the cosine it writes is not the faces'), and with
 * its longwave on under lw_mode = two_stream (off, the faces would see a 0 K sky
 * and ground). Every level with faces must also sweep its own columns; init_ibseb()
 * checks that once the faces are built (a refined level whose grids do not span
 * the domain in z takes its radiation from the level below by interpolation,
 * rad_level_needs_interpolation(), with no sweep of its own). (A level built from
 * the coarse atmosphere also skips its first sweep; that is
 * erf.interp_atmos_from_coarse with a WRF input, which needs a terrain-fitted grid
 * the balance refuses, or a regrid, which it refuses too.)
 *
 * It also warns, without stopping, when the faces' light and sky depend on the
 * number of layers above them: whenever some of the sky's optical depth is set
 * per layer. That is the clear-sky depth with erf.radiation.tau_model = per_layer
 * (for the longwave, also without erf.radiation.lw_mass_absorption_enable), and,
 * with either model, a cloud layer that is used, the moisture terms and the
 * aerosol. The attenuation above a face then counts the layers, so it changes with
 * the domain's depth and, on a level refined in z, with the refinement, as the
 * columns' own heating rates do.
 */
void
ERF::ibseb_check_two_stream_provider () const
{
    if (ibseb_params.radiation != "two_stream") { return; }
    if (solverChoice.rad_type != RadiationType::TwoStream) {
        Abort("erf.ibseb.radiation = two_stream needs erf.radiation_model = TwoStream; "
              "use erf.ibseb.radiation = prescribed for the faces' own clear-sky radiation");
    }
    const RadChoice& rc = solverChoice.radChoice;
    if (!rc.sw_enabled) {
        Abort("erf.ibseb.radiation = two_stream needs the two-stream shortwave (erf.radiation.sw_enabled = true):"
              " the faces take their direct and diffuse light from it");
    }
    if (ibseb_params.lw_mode == "two_stream" && !rc.lw_enabled) {
        Abort("erf.ibseb.lw_mode = two_stream needs the two-stream longwave (erf.radiation.lw_enabled = true);"
              " without it the faces would see no sky or ground longwave. Use erf.ibseb.lw_mode = gray or fixed");
    }
    // Some of the sky's optical depth is set per layer: the clear-sky depth under
    // per-layer optics, and with either model a cloud layer that is used (a cloud
    // fraction above zero), the moisture terms and the aerosol (every profile gives a
    // depth per layer). A transparent sky is the same whatever the layer count.
    const bool per_layer   = (rc.tau_model == TauModel::PerLayer);
    const bool cloud_layer = (rc.tau_profile_type == TauProfileType::CloudLayer)
                          && (rc.cloud_fraction > 0.0 || rc.cloud_fraction_prog_enable);
    const bool sw_per_layer = (per_layer && rc.tau_per_layer > 0.0)
                           || cloud_layer || rc.aerosol_enable || rc.tau_sw_dynamic_enable;
    const bool lw_per_layer = (ibseb_params.lw_mode == "two_stream")
                           && ((per_layer && !rc.lw_mass_absorption_enable && rc.tau_lw_per_layer > 0.0)
                               || cloud_layer || rc.aerosol_enable || rc.tau_lw_dynamic_enable);
    bool z_refined = false;
    for (int lev = 1; lev <= finest_level; ++lev) { z_refined = z_refined || ref_ratio[lev-1][2] > 1; }
    if (sw_per_layer || lw_per_layer) {
        const std::string bands = (sw_per_layer && lw_per_layer) ? "shortwave and longwave"
                                : (sw_per_layer ? "shortwave" : "longwave");
        Print() << "WARNING: erf.ibseb.radiation = two_stream with optical depth set per layer (" << bands << "):"
                << " the light and sky the faces get depend on the number of layers above them, so on the domain's"
                << " depth" << (z_refined ? " and on the refinement in z, which differs between the levels here" : "")
                << ". erf.radiation.tau_model = mass makes the clear-sky depth independent of the layers; a cloud"
                   " layer, the moisture terms and the aerosol stay per layer.\n";
    }
}

/**
 * Hand one level's face set the two-stream sun at a time, for its next
 * shortwave (erf.ibseb.sun_mode = two_stream; a no-op otherwise): the
 * declination, the distance factor and the calendar day from the start date
 * through two_stream_sun_date(), the shared routine of the two-stream sweep;
 * the hour angle of the same orbital formula (ibseb::orbital_hour_angle()) at
 * erf.rad_cons_lon; erf.rad_cons_lat; and the irradiance the sweep takes,
 * erf.fixed_total_solar_irradiance or the reference scaled by the distance
 * factor.
 *
 * @param[in] lev   Level whose face set takes the sun.
 * @param[in] time  Simulation time [s] of the shortwave that follows.
 */
void
ERF::ibseb_set_two_stream_sun (int lev, Real time)
{
    if (ibseb_params.sun_mode != "two_stream" || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    const RadChoice& rc = solverChoice.radChoice;
    const TwoStreamSunDate sun = two_stream_sun_date(rc, ibseb_orbit, start_time + static_cast<double>(time));
    const Real ha = ibseb::orbital_hour_angle(sun.calday, rc.rad_cons_lon * PI / Real(180.0));
    const Real S0 = (rc.fixed_total_solar_irradiance >= 0.0)
                  ? rc.fixed_total_solar_irradiance
                  : two_stream_tsi_reference * static_cast<Real>(sun.eccf);
    m_ibseb[lev]->set_two_stream_sun(static_cast<Real>(sun.declin), ha, rc.rad_cons_lat, S0);
}

/**
 * Per-step update of one level, called from ERF::Advance() with the state at
 * the start of the step: at the start of the step, or with
 * erf.ibseb.radiation = two_stream after advance_radiation(), whose sweep of
 * this step the faces take (ibseb_after_radiation()). Shortwave, longwave and the wall
 * function on the faces, then either the prognostic balance (which finds
 * the skin temperature at the end of the step and advances the slab with
 * it) or, with ``erf.ibseb.prognostic = false``, the slab alone under the
 * fixed skin. The sensible flux left in the set is what
 * add_heat_flux_to_source() deposits at every slow stage of the step, so
 * the air receives exactly the H of the closed balance. The atmosphere is
 * seen at the start of the step and the skin is implicit within it, the
 * usual coupling of a land-surface model.
 *
 * @param[in] lev     AMR level to update; a no-op if it has no face set.
 * @param[in] time    Time at the start of the step [s], at which the sun
 *                    position of the shortwave is taken.
 * @param[in] dt_lev  Length of the level's step [s], over which the slab is
 *                    advanced and within which the skin is implicit.
 * @param[in] cons    Conserved state at the start of the step: the air
 *                    temperature of the longwave and the wall function, and
 *                    the profile of the bulk Richardson depth.
 * @param[in] xvel    Face-centred x velocity at the start of the step.
 * @param[in] yvel    Face-centred y velocity at the start of the step.
 * @param[in] zvel    Face-centred z velocity at the start of the step; the
 *                    three together drive the wall function of the sensible
 *                    flux.
 */
void
ERF::ibseb_advance (int lev, Real time, Real dt_lev, const MultiFab& cons,
                    const MultiFab& xvel, const MultiFab& yvel, const MultiFab& zvel)
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    const double t_wall0 = ParallelDescriptor::second();
    // With radiation = two_stream the faces take this step's sweep of their level,
    // which advance_radiation() has just run (ERF::Advance() calls this after it).
    TwoStreamCanopyView canopy;
    if (ibseb_params.radiation == "two_stream" && m_ibseb[lev]->has_faces()) {
        canopy = two_stream_rad.canopy_forcing(lev, istep[lev]);
        canopy.fluxes = rad_fluxes[lev].get();
        if (!canopy.valid()) {
            Abort("erf.ibseb.radiation = two_stream: level " + std::to_string(lev) + " has no two-stream sweep at"
                  " step " + std::to_string(istep[lev]) + " for its faces; the balance must run after the step's"
                  " radiation (ERF::Advance()), on a level that sweeps its own columns");
        }
        m_ibseb[lev]->check_canopy_layout(cons, canopy);
    }
    ibseb_set_two_stream_sun(lev, time);
    m_ibseb[lev]->compute_shortwave(time, canopy);
    m_ibseb[lev]->compute_longwave(cons, canopy);
    // The ground surface layer's fields and the mixed-layer depth
    // for the wall function beyond neutral (all null / zero unless asked).
    const MultiFab* olen2d = nullptr;
    const MultiFab* pblh2d = nullptr;
    Real z_i_bulk = 0.0;
    // The surface layer now exists per domain face.  What the wall function wants
    // here is the ground beneath the buildings, so take zlo; a surface layer on a
    // lateral or upper wall says nothing about the stability of this column.
    const auto& ground_sl = m_SurfaceLayer[Orientation(Direction::z, Orientation::low)];
    if (ground_sl && ibseb_params.stability_correction) { olen2d = ground_sl->get_olen(lev); }
    if (ibseb_params.convective_velocity == "deardorff") {
        if (ground_sl && ground_sl->computes_pblh() && ibseb_params.z_i_mode == "pblh") {
            pblh2d = ground_sl->get_pblh(lev);
        }
        // The mixed layer is a property of the whole domain: on every level it
        // comes from level 0's profile, which spans the domain, where a
        // refined level's plane average covers only its patch and, on a
        // level that stops below the top, reads its last height above it.
        // Level 0 computes it from the state at the start of its step and
        // keeps it; a refined level, stepping inside that step, reuses it
        // rather than reading level 0's state at the step's end.
        if (ibseb_params.z_i_mode == "fixed") {
            z_i_bulk = ibseb_params.z_i;
        } else if (lev == 0) {
            z_i_bulk = ibseb_bulk_richardson_height(0, cons, xvel, yvel);
            ibseb_z_i_level0 = z_i_bulk;
        } else {
            z_i_bulk = ibseb_z_i_level0;
        }
        if (ibseb_params.debug) {
            Print() << "[IBSEB DEBUG] lev=" << lev << " mixed-layer depth for w*: " << z_i_bulk << " m ("
                    << ibseb_params.z_i_mode << (pblh2d ? ", pblh per column" : "") << ")\n";
        }
    }
    m_ibseb[lev]->compute_sensible(cons, xvel, yvel, zvel, solverChoice.c_p, olen2d, pblh2d, z_i_bulk);
    if (ibseb_params.prognostic) {
        m_ibseb[lev]->solve_balance(dt_lev);
    } else {
        m_ibseb[lev]->compute_ground(dt_lev);
    }
    m_ibseb[lev]->add_cost(ParallelDescriptor::second() - t_wall0);
}

/**
 * Write the face state of one level into the checkpoint as ``IBSEBState``, a
 * field on the blocks around the buildings (IBFaceSet::state_boxarray())
 * clipped to the cells that own faces, so it scales with the shell of the
 * buildings rather than with the level. Called inside the level loop of
 * ERF::WriteCheckpointFile(); a no-op unless the balance is on and the level
 * has faces.
 *
 * @param[in] checkpointname  Path of the checkpoint directory being written.
 * @param[in] lev             AMR level whose face state is written, as the
 *                            ``IBSEBState`` field of its ``Level_`` group.
 */
void
ERF::ibseb_write_checkpoint (const std::string& checkpointname, int lev) const
{
    if (!ibseb_params.enable || lev >= static_cast<int>(m_ibseb.size()) || !m_ibseb[lev]) { return; }
    // A level without faces has no field to write (and nothing to restore).
    if (!m_ibseb[lev]->has_state()) { return; }
    MultiFab state = m_ibseb[lev]->make_state();
    m_ibseb[lev]->save_state(state);
    VisMF::Write(state, MultiFabFileFullPrefix(lev, checkpointname, "Level_", "IBSEBState"));
}

/**
 * Periodic report from ERF::post_timestep(), with ``nstep`` the number of
 * completed steps (the plotfiles' numbering; the initial state is the step-0
 * report of init_ibseb()): after every ``erf.ibseb.csv_int``-th step, print
 * the summary of each level and append its CSV rows. A non-positive interval
 * disables both.
 *
 * @param[in] nstep  Number of completed steps, tested against
 *                   ``erf.ibseb.csv_int`` and written to the CSV rows.
 * @param[in] time   Simulation time at the end of the step [s], written to
 *                   the summary and the CSV rows.
 */
void
ERF::ibseb_report (int nstep, Real time)
{
    if (!ibseb_params.enable) { return; }
    // With debug on the summary is printed every step, as the fire module
    // does; the CSV rows keep their interval.
    const bool csv_now = (ibseb_params.csv_int > 0) && (nstep % ibseb_params.csv_int == 0);
    if (!csv_now && !ibseb_params.debug) { return; }
    for (int lev = 0; lev <= finest_level && lev < static_cast<int>(m_ibseb.size()); ++lev) {
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
 * ``erf.ibseb.ri_crit``, or the domain depth when it never does (a neutral
 * profile); both are heights above the domain bottom. The profile is the plane average of the conserved state and the
 * face velocities, uniform vertical spacing assumed as elsewhere in the
 * balance; called once per step and level when the convective velocity
 * scale is on and z_i is not fixed, also as the fallback of the pblh mode,
 * on level 0 only (whose profile spans the domain), at the start of its
 * step; the faces of the refined levels take that value.
 *
 * @param[in] lev   AMR level whose horizontal-mean profile is taken.
 * @param[in] cons  Conserved state of the level; the ``Rho_comp`` and
 *                  ``RhoTheta_comp`` averages give the mean potential
 *                  temperature.
 * @param[in] xvel  Face-centred x velocity, for the mean wind profile.
 * @param[in] yvel  Face-centred y velocity, for the mean wind profile.
 * @return Mixed-layer depth above the domain bottom [m]; the depth of the
 *         domain when ``Ri_b`` never exceeds ``erf.ibseb.ri_crit``.
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
    const Real z_top = geom[lev].ProbHi(2) - geom[lev].ProbLo(2);   // depth of the domain
    const Real th1 = rth[0] / rho[0];
    const Real ustar_floor2 = 100.0 * 0.1 * 0.1;
    Real z_i = z_top;
    for (int k = 1; k < nz; ++k) {
        const Real th = rth[k] / rho[k];
        // |U(z) - U_1|^2 of the wind vector, so a veering wind of constant speed
        // still counts as shear.
        const Real dU2 = (uu[k] - uu[0]) * (uu[k] - uu[0]) + (vv[k] - vv[0]) * (vv[k] - vv[0]);
        const Real rib = CONST_GRAV * (k * dz) * (th - th1) / (th1 * (dU2 + ustar_floor2));
        if (rib > ibseb_params.ri_crit) { z_i = (k + 0.5) * dz; break; }
    }
    return z_i;
}
