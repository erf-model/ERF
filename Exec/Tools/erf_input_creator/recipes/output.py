"""Plotfiles, checkpoints, logs, profiles, samplers and the post-processing notes.

Every output is written; the ones named in ``--output`` are active and the
others are commented, so a user turns one on by removing the ``#``.
"""


from .models import moist_plot_vars


def plot_vars(s):
    """3-D plot variables that the chosen models can fill."""
    names = ["density", "rhotheta", "x_velocity", "y_velocity", "z_velocity", "pressure", "theta", "temp"]
    if s.has_tke:
        names.append("rhoKE")
    if s.turbulence not in ("none", "dns"):
        names += ["Kmv", "Khv"]
    if s.gravity:
        names += ["pert_pres", "pert_dens"]
    if s.method == "fitted" or (s.method == "none" and s.stretched):
        names.append("z_phys")      # node heights of the terrain-fitted or stretched grid
    if s.method == "immersed" or s.buildings_on:
        names.append("terrain_IB_mask")
    if s.method == "eb":
        names.append("volfrac")
    names += moist_plot_vars(s)
    if s.real and "z_phys" not in names:
        names.append("z_phys")
    return names


def plot2d_vars(s):
    names = ["z_surf", "surf_pres"]
    if s.zlo_most:
        names += ["u_star", "t_star", "Olen", "z0", "u_10m", "v_10m"]
    return names


def output(deck, s) -> None:
    steps = 100
    mid = [0.5 * (lo + hi) for lo, hi in zip(s.lo, s.hi)]

    sec = deck.section("OUTPUT: SCREEN AND LOG FILES")
    logs = "logs" in s.outputs or "profiles" in s.outputs
    sec.set("erf.v", 1, "ERF verbosity; the log files below stay empty with 0")
    sec.set("amr.v", 1, "AMReX verbosity")
    sec.block(logs, [("erf.sum_interval", 1, "steps between the mass and energy sums and the surface log")])
    sec.alt("erf.sum_period", 60.0, "[s] the same cadence in time; either may trigger")
    if "profiles" in s.outputs:
        sec.set("erf.data_log", ("surf.log", "mean.log", "flux.log", "sgs.log"),
                "1: time u* t* L (surface layer only); 2: mean profiles; 3: resolved fluxes; 4: subgrid fluxes")
        sec.set("erf.profile_int", steps, "steps between rows of the profile logs 2 to 4")
    else:
        sec.block("logs" in s.outputs, [("erf.data_log", "surf.log", "time u* t* L; rows only with a surface layer")])
        if s.turbulence == "none":
            sec.note("Profile logs (erf.profile_int) are left out: without a closure ERF crashes at the first profile write.")
        else:
            sec.alt("erf.data_log", ("surf.log", "mean.log", "flux.log", "sgs.log"),
                    "adds horizontally averaged profiles; needs erf.profile_int")
            sec.alt("erf.profile_int", steps, "steps between profile rows")
    sec.alt("erf.destag_profiles", "false", "profiles on their staggered faces (default true: cell centres)")
    sec.block("logs" in s.outputs, [("erf.der_data_log", "der.log", "time, kinetic energy, enstrophy, internal energy"),
                                    ("erf.energy_data_log", "energy.log", "time, total mass, total energy")])
    sec.alt("erf.check_for_nans", 1, "stop at the first NaN (default 0)")

    sec = deck.section("OUTPUT: PLOTFILES")
    on = "plt" in s.outputs
    sec.note("Stream 1; a second stream uses the same keys ending in _2. Without plot_vars no file is written.")
    sec.block(on, [("erf.plot_file_1", "plt", "prefix; files are plt00000, plt00100, ..."),
                   ("erf.plot_int_1", steps, "steps between plotfiles")])
    sec.alt("erf.plot_per_1", 600.0, "[s] cadence in time instead; plot_int_1 and plot_per_1 together stop the run")
    sec.block(on, [("erf.plot_vars_1", plot_vars(s), "")])
    sec.note("More names: KE scalar soundspeed vorticity_x vorticity_y vorticity_z magvel divU pres_hse dens_hse")
    sec.note("theta_hse buoyancy dpdx dpdy dpdz Tau11..Tau33 hfx1 hfx2 hfx3 nut Kmh Khh Lturb diss walldist;")
    sec.note("with erf.time_avg_vel = true: u_t_avg v_t_avg w_t_avg umag_t_avg; with erf.compute_mean_vars = true:")
    sec.note("u_mean ... wtheta_mean, uu_fluct ... tke_resolved. A name the run cannot fill is dropped with a warning.")
    if "netcdf" in s.outputs:
        sec.set("erf.plotfile_type", "netcdf", "one file per level, <prefix><step>_d01.nc; needs a NetCDF build")
        sec.alt("erf.plotfile_type", "amrex", "native AMReX plotfile (default)")
    else:
        sec.alt("erf.plotfile_type", "netcdf", "NetCDF instead of the native format; needs a NetCDF build")
    sec.alt("erf.plot_face_vels", "true", "also write the staggered velocities (native format only)")
    sec.alt("erf.file_name_digits", 5, "digits of the step in file names (default)")
    sec.blank()
    on = "plt2d" in s.outputs
    sec.note("2-D plotfiles of surface and column fields%s." % ("" if on else " (not active)"))
    sec.block(on, [("erf.plot2d_file_1", "plt2d", "prefix"),
                   ("erf.plot2d_int_1", steps, "steps between 2-D plotfiles"),
                   ("erf.plot2d_vars_1", plot2d_vars(s), "")])
    sec.note("More 2-D names: w_star q_star pblh t_surf sens_flux laten_flux sea_level_pressure lat_m lon_m.")

    sec = deck.section("OUTPUT: CHECKPOINT AND RESTART")
    on = "chk" in s.outputs
    sec.block(on, [("erf.check_file", "chk", "prefix; files are chk00000, chk01000, ..."),
                   ("erf.check_int", 10 * steps, "steps between checkpoints")])
    sec.alt("erf.check_per", 3600.0, "[s] cadence in time; ignored when check_int is also set")
    sec.alt("erf.restart", "chk01000", "restart from this checkpoint")

    sec = deck.section("OUTPUT: SAMPLING")
    on = "line" in s.outputs
    sec.note("Line sampling: one vertical column at the domain centre%s." % ("" if on else " (not active)"))
    ic, jc = s.n[0] // 2, s.n[1] // 2
    sec.block(on, [("erf.do_line_sampling", "true", ""),
                   ("erf.line_sampling_interval", steps, "steps between samples (or erf.line_sampling_per [s])"),
                   ("erf.sample_line_lo", (ic, jc, 0), "level-0 cell index of the first point"),
                   ("erf.sample_line_hi", (ic, jc, s.n[2] - 1), "level-0 cell index of the last point"),
                   ("erf.sample_line_dir", 2, "direction of the line: 0 x, 1 y, 2 z"),
                   ("erf.sample_line_name", "column", "one name per line"),
                   ("erf.line_sampling_vars", ("x_velocity", "y_velocity", "z_velocity", "theta"),
                    "also: density magvel pressure qv qc sgs_tke sgs_tau13 sgs_hfx3"),
                   ("erf.line_sampling_text_output", "true", "text files <name>.<var>; false writes plotfiles")])
    sec.alt("erf.sample_line_lo_real", (mid[0], mid[1], s.lo[2]), "[m] physical coordinates instead of indices")
    sec.alt("erf.sample_line_hi_real", (mid[0], mid[1], s.hi[2]), "[m]")
    sec.blank()
    on = "plane" in s.outputs
    sec.note("Plane sampling: one horizontal plane%s." % ("" if on else " (not active)"))
    if s.levels > 0:
        sec.note("Not available here: plane sampling over refined levels needs amr.ref_ratio_vect = 2 2 2.")
    else:
        sec.block(on, [("erf.do_plane_sampling", "true", ""),
                     ("erf.plane_sampling_interval", steps, "steps between samples (or erf.plane_sampling_per [s])"),
                   ("erf.sample_plane_lo", (s.lo[0], s.lo[1], mid[2]), "[m] lower corner of the plane"),
                   ("erf.sample_plane_hi", (s.hi[0], s.hi[1], mid[2]), "[m] upper corner of the plane"),
                   ("erf.sample_plane_dir", 2, "normal of the plane: 0 x, 1 y, 2 z"),
                   ("erf.sample_plane_name", "plane_z", "prefix of the plotfiles"),
                   ("erf.plane_sampling_vars", ("x_velocity", "y_velocity", "z_velocity", "theta"),
                    "also: density magvel pressure qv qc")])
    sec.blank()
    on = "stations" in s.outputs
    sec.note("Station time series, written to Output_Stations/<name>.dat%s." % ("" if on else " (not active)"))
    if on:
        sec.set("erf.station_names", "mast1", "one name per station")
        sec.set("erf.station_sampling_interval", 10, "steps between samples (or erf.station_sampling_per [s])")
        sec.user("erf.mast1.x", mid[0], "[m] station position (or .lat and .long with real-data cases)")
        sec.user("erf.mast1.y", mid[1], "[m]")
        sec.user("erf.mast1.height_agl", (0.1 * (s.hi[2] - s.lo[2]),), "[m] heights above ground, one or more")
        sec.set("erf.mast1.field", ("x_velocity", "y_velocity", "theta"), "any 3-D or 2-D plot variable name")
    else:
        sec.alt("erf.station_names", "mast1", "one name per station")
        sec.alt("erf.station_sampling_interval", 10, "steps between samples (or erf.station_sampling_per [s])")
        sec.alt("erf.mast1.x", mid[0], "[m]")
        sec.alt("erf.mast1.y", mid[1], "[m]")
        sec.alt("erf.mast1.height_agl", (0.1 * (s.hi[2] - s.lo[2]),), "[m] heights above ground")
        sec.alt("erf.mast1.field", ("x_velocity", "y_velocity", "theta"), "any plot variable name")
    sec.blank()
    on = "subvol" in s.outputs
    nsub = tuple(max(1, n // 2) for n in s.n)
    origin = tuple(lo + (n // 4) * d for lo, n, d in zip(s.lo, s.n, s.dx))
    sec.note("Sub-volume plotfiles: a box of level-0 cells written more often%s." % ("" if on else " (not active)"))
    sec.block(on, [("erf.subvol_file", "subvol", "prefix (default)"),
                   ("erf.subvol_int", 10, "steps between sub-volume files (or erf.subvol_per [s])"),
                   ("erf.subvol.origin", origin, "[m] lower corner; must lie on grid nodes"),
                   ("erf.subvol.nxnynz", nsub, "cells in the box"),
                   ("erf.subvol.dxdydz", s.dx, "[m] must equal the cell size of one level"),
                   ("erf.subvol_sampling_vars", ("x_velocity", "y_velocity", "z_velocity"), "(default)")])

    sec = deck.section("POST-PROCESSING")
    sec.note("Native plotfiles open in ParaView (AMReX/BoxLib Grid Reader) and VisIt.")
    sec.note("Exec/Tools/erf_plotfile.py reads a single-level plotfile in plain Python (planar averages, heights).")
    sec.note("Source/Tools/ERF_PlotfileToNetCDF.cpp converts a native plotfile to NetCDF (CMake target erf_plt_to_nc).")
    sec.note("The log and sampler files are plain text, one row per sample; numpy.loadtxt reads them.")
    sec.note("Check this file before a run:  python3 Exec/Tools/check_erf_inputs.py %s" % s.name)
