"""Problem name, geometry, boundary faces, time step, equation set, advection, levels."""

import math

from ..enums import enum_values

SOUND_SPEED = 347.0   # [m/s] at 300 K, for the acoustic substep estimate


def problem(deck, s) -> None:
    sec = deck.section("PROBLEM")
    sec.set("erf.prob_name", '"%s"' % s.prob_name,
            "selects the initial perturbations and custom sources; an unknown name adds none")
    sec.note("Other names are listed in Exec/ERF_Prob.cpp (for example \"Bubble\", \"Density Current\",")
    sec.note("\"Flow over Witch of Agnesi Hill\", \"BOMEX\"). The prob.* keys below belong to this name;")
    sec.note("ERF does not report a prob.* key that the chosen problem does not read.")


def geometry(deck, s) -> None:
    sec = deck.section("DOMAIN AND GRID")
    sec.set("geometry.prob_lo", s.lo, "[m] lower corner")
    sec.set("geometry.prob_hi", s.hi, "[m] upper corner")
    sec.set("amr.n_cell", s.n, "cells: dx = %g, dy = %g, dz = %g m" % s.dx)
    sec.set("geometry.is_periodic", s.periodic, "x y z; a periodic face must not be given a type")
    if s.real:
        sec.note("amr.n_cell and the extent must match the NetCDF file: nx = west_east, ny = south_north,")
        sec.note("extent = DX * nx by DY * ny. Edit the three lines above to the file's values.")
    sec.blank()
    if s.levels > 0:
        _refinement(sec, s)
    else:
        sec.set("amr.max_level", 0, "0 = one level; refinement needs erf.refinement_indicators")
        sec.alt("erf.refinement_indicators", "box1", "names of the refinement criteria; keys erf.<name>.*")
    sec.alt("amr.max_grid_size", 64, "largest box edge in cells (default 2048)")
    sec.alt("amr.blocking_factor", 8, "box edges are multiples of this (default 1)")
    sec.note("amr.no_box_split_dir defaults to 2: boxes are never split in z.")


def _refinement(sec, s) -> None:
    lx, ly = s.hi[0] - s.lo[0], s.hi[1] - s.lo[1]
    sec.set("amr.max_level", s.levels, "refined levels above level 0")
    sec.set("amr.ref_ratio_vect", (2, 2, 1) * s.levels, "refinement ratio x y z of each level; z is not refined")
    names = ["box%d" % (lev + 1) for lev in range(s.levels)]
    sec.set("erf.refinement_indicators", names, "one criterion per level; keys erf.<name>.*")
    for lev, name in enumerate(names):
        frac = 0.25 + 0.125 * lev          # nested boxes around the domain centre
        sec.set("erf.%s.max_level" % name, lev + 1, "refine up to this level inside the box")
        sec.user("erf.%s.in_box_lo" % name, (s.lo[0] + frac * lx, s.lo[1] + frac * ly), "[m] x y of the lower corner")
        sec.user("erf.%s.in_box_hi" % name, (s.hi[0] - frac * lx, s.hi[1] - frac * ly), "[m] x y of the upper corner")
    sec.alt("erf.box1.field_name", "theta", "refine by a field instead of a box: with value_greater or value_less")
    sec.set("erf.regrid_int", -1, "steps between regrids; -1 keeps the initial boxes")
    sec.set("erf.coupling_type", "TwoWay", "fine levels feed back to coarse levels")
    sec.alt("erf.coupling_type", "OneWay", "coarse to fine only")


_WALL_NOTE = "wall faces also take theta, theta_grad (not both) and, for compressible runs, density"


def boundaries(deck, s) -> None:
    sec = deck.section("BOUNDARY CONDITIONS")
    types = "Outflow, Open (x/y only), HO_Outflow, Inflow, Inflow_Outflow, SlipWall, NoSlipWall, " \
            "Symmetry, surface_layer"
    sec.note("Face types: " + types + ".")
    if s.lateral == "real":
        sec.note("Real-data lateral boundaries: ERF sets the four lateral faces itself and relaxes the solution to")
        sec.note("the boundary data over erf.real_width cells, so xlo, xhi, ylo, yhi need no type here.")
    elif s.lateral == "periodic":
        sec.note("x and y are periodic, so xlo, xhi, ylo, yhi take no type.")
    elif s.lateral == "slipwall":
        for face in ("xlo", "xhi", "ylo", "yhi"):
            sec.set(face + ".type", "SlipWall")
    elif s.lateral in ("outflow", "open"):
        kind = "Outflow" if s.lateral == "outflow" else "Open"
        other = "Open" if kind == "Outflow" else "Outflow"
        for face in ("xlo", "xhi"):
            sec.set(face + ".type", kind)
            sec.alt(face + ".type", other)
        sec.note("y is periodic.")
    else:  # inflow_outflow
        sec.set("xlo.type", "Inflow")
        sec.user("xlo.velocity", (s.wind[0], s.wind[1], 0.0), "[m/s] inflow velocity")
        sec.alt("xlo.density", 1.16, "[kg/m3]; when given, xlo.theta is required")
        sec.alt("xlo.theta", s.theta, "[K]")
        sec.alt("xlo.dirichlet_file", "inflow_profile.txt",
                "rows 'z u v w' or 'z u v w theta', z increasing; replaces xlo.velocity")
        sec.set("xhi.type", "Outflow")
        sec.alt("xhi.type", "Open")
        sec.note("y is periodic.")
    sec.blank()
    if s.method == "immersed":
        sec.set("zlo.type", "NoSlipWall", "below the immersed surface; the wall law is set under TERRAIN (erf.if_*)")
    elif s.surface == "most":
        sec.set("zlo.type", "surface_layer", "Monin-Obukhov surface layer; see SURFACE LAYER below")
        sec.alt("zlo.type", "SlipWall")
        sec.alt("zlo.type", "NoSlipWall")
    elif s.surface == "noslip":
        sec.set("zlo.type", "NoSlipWall")
        sec.alt("zlo.velocity", (0.0, 0.0, 0.0), "[m/s] wall velocity (the wall-normal part is ignored)")
        sec.alt("zlo.type", "SlipWall")
        sec.alt("zlo.type", "surface_layer", "needs a turbulence closure")
    else:
        sec.set("zlo.type", "SlipWall")
        sec.alt("zlo.type", "NoSlipWall")
        sec.alt("zlo.type", "surface_layer", "needs a turbulence closure")
    sec.set("zhi.type", "SlipWall")
    sec.alt("zhi.theta_grad", 0.003, "[K/m] set to the lapse rate at the top of a stratified sounding")
    sec.note(_WALL_NOTE + ".")


def time(deck, s) -> None:
    sec = deck.section("TIME")
    sec.set("stop_time", s.stop_time, "[s] run length")
    if s.real:
        sec.note("The start time is read from the NetCDF file; setting start_time stops the run.")
        sec.alt("start_datetime", '"2024-08-05 00:00:00"', "optional; must equal the file's start time")
        sec.alt("stop_datetime", '"2024-08-05 06:00:00"', "end of the run on the same clock, instead of stop_time")
    sec.alt("max_step", 1000, "stop after this many steps instead (default: no limit)")
    umax, guess = s.umax, s.dt_guess
    comp = s.dycore == "compressible"
    cfl = 0.8 if comp else 0.5
    adaptive = ("adaptive step from the wind and the sound speed (small steps, no substeps needed)" if comp
                else "adaptive step from the advective CFL number")
    if s.dt is not None:
        sec.set("erf.fixed_dt", s.dt, "[s] fixed time step")
        sec.alt("erf.cfl", cfl, adaptive + "; remove erf.fixed_dt")
    elif s.at_rest:
        # the anelastic adaptive step has no velocity to limit it at the start of a run at rest
        sec.user("erf.fixed_dt", guess, "[s] 0.5 dx / %g m/s, the largest speed expected; the flow starts at rest" % umax)
        sec.alt("erf.cfl", cfl, adaptive + ("" if comp else "; unbounded while the flow is at rest"))
    else:
        sec.set("erf.cfl", cfl, adaptive)
        if s.adaptive_anelastic and s.damping != "none":
            sec.set("erf.dt_max", s.dt_est, "[s] largest step; keeps the explicit damping rate x step at or below 0.25")
        sec.alt("erf.fixed_dt", guess, "[s] fixed step instead: 0.5 dx / %g m/s%s"
                % (umax, "; ERF then takes acoustic substeps within it" if comp else ""))
    if comp:
        # the acoustic substep must keep c dt / dx below about 0.7 over the two horizontal directions
        horizontal = [d for d in s.dx[:2] if d > 0.0]
        fast_max = 0.7 * min(horizontal) / (SOUND_SPEED * 2.0 ** 0.5)
        step = s.dt if s.dt is not None else guess
        nsub = max(2, int(math.ceil(step / fast_max)))
        nsub += nsub % 2
        sec.note("With a fixed step ERF chooses the number of acoustic substeps from the sound speed. To fix it:")
        sec.alt("erf.fixed_mri_dt_ratio", nsub, "even number of substeps per step; %d keeps c dt / dx near 0.7 for a %g s step"
                % (nsub, step))
        sec.alt("erf.substepping_cfl", 0.7, "acoustic CFL number used for that choice (default 0.7)")


def dycore(deck, s) -> None:
    sec = deck.section("EQUATION SET")
    chosen = "FullyCompressible" if s.dycore == "compressible" else "Anelastic"
    sec.choice("erf.compressibility_type", chosen,
               enum_values("CompressibilityType", ["FullyCompressible", "Anelastic", "PseudoIncompressible"]),
               describe={"FullyCompressible": "acoustic waves resolved by substeps",
                         "Anelastic": "no acoustic waves; pressure from a Poisson solve",
                         "PseudoIncompressible": "Poisson solve; prototype"})
    if s.dycore == "compressible":
        sec.alt("erf.timestepping_type", "RK3_ImplicitSubstep", "default: RK3 with vertically implicit acoustic substeps")
        sec.alt("erf.timestepping_type", "RK3_NoSubstepping", "RK3 limited by the acoustic CFL number")
    else:
        sec.alt("erf.timestepping_type", "MidPoint", "default for Anelastic")
        sec.alt("erf.timestepping_type", "RK2")
        # FFT is the default of this tool for every anelastic case: at ERF abce396eb the multigrid solve did
        # not converge on cells with aspect ratios of 4:1 and more (coarse flat grids, immersed forcing, two
        # refined levels with z unrefined), and erf.terrain_poisson_solver = gmres_mlmg is read on fitted meshes only
        sec.set("erf.use_fft", "true", "FFT Poisson solve; needs an FFT build (ERF_ENABLE_FFT)")
        sec.alt("erf.use_fft", "false", "multigrid, no FFT build needed; did not converge on cells with aspect ratio 4:1 or more")
        sec.alt("erf.mg_v", 1, "print the Poisson solver's convergence")
        sec.alt("erf.poisson_reltol", 1e-8, "relative tolerance of the Poisson solve")
        sec.alt("erf.poisson_abstol", 1e-8, "absolute tolerance of the Poisson solve")
    sec.blank()
    if s.real:
        sec.set("erf.use_gravity", "true", "ERF forces it on for WRFInput and Metgrid")
    elif s.gravity:
        sec.set("erf.use_gravity", "true", "the default is false")
    else:
        sec.set("erf.use_gravity", "false",
                "init_type = Uniform has no hydrostatic base state; with gravity a compressible run stops")
    if s.dycore == "compressible":
        sec.alt("erf.buoyancy_type", 1, "1 (default) to 4: form of the buoyancy term")


def advection(deck, s) -> None:
    sec = deck.section("ADVECTION")
    sec.note("Values: Centered_2nd, Upwind_3rd, Blended_3rd4th, Centered_4th, Upwind_5th, Blended_5th6th,")
    sec.note("Centered_6th, WENO3, WENOZ3, WENOMZQ3, WENO5, WENOZ5, WENO7, WENOZ7; scalars also Upwind_3rd_SL.")
    sec.set("erf.dycore_horiz_adv_type", "Upwind_3rd", "momentum, density, theta (default)")
    sec.set("erf.dycore_vert_adv_type", "Upwind_3rd", "(default)")
    sec.set("erf.dryscal_horiz_adv_type", "Upwind_3rd", "TKE and passive scalars (default)")
    sec.set("erf.dryscal_vert_adv_type", "Upwind_3rd", "(default)")
    sec.alt("erf.dycore_horiz_adv_type", "Upwind_5th", "less numerical diffusion for LES")
    sec.alt("erf.dycore_horiz_upw_frac", 1.0, "only with a Blended_* type: 1 = upwind, 0 = centred")
    sec.alt("erf.use_efficient_advection", "true", "cheaper stencil evaluation (default false)")
