"""Terrain type and source, vertical stretching, immersed-forcing wall keys, buildings."""

from ..enums import enum_values

TERRAIN_FILE = "terrain.txt"
TERRAIN_USGS_FILE = "terrain_usgs.txt"
BUILDINGS_FILE = "buildings.txt"

_TYPE = {"none": "None", "fitted": "StaticFittedMesh", "immersed": "ImmersedForcing", "eb": "EB"}
_DESCRIBE = {"None": "flat surface; a terrain file or prob.custom_terrain_type is ignored",
             "StaticFittedMesh": "terrain-following grid",
             "MovingFittedMesh": "moving terrain; one level, compressible only",
             "EB": "embedded boundary: cut cells on a uniform grid",
             "ImmersedForcing": "immersed forcing on a uniform grid"}


def terrain(deck, s) -> None:
    sec = deck.section("TERRAIN AND VERTICAL GRID")
    if s.real:
        sec.set("erf.terrain_type", "StaticFittedMesh", "required for WRFInput and Metgrid; terrain and levels come from the file")
        sec.set("erf.terrain_smoothing", 1, "0 basic, 1 smoothed, 2 Sullivan terrain-following")
        sec.alt("erf.terrain_gamma_m", 0.5, "with terrain_smoothing = 1 (default)")
        return
    chosen = _TYPE[s.method]
    if s.method == "none" and s.stretched:
        sec.note("A stretched grid over a flat surface: ERF switches terrain_type None to StaticFittedMesh itself.")
    sec.choice("erf.terrain_type", chosen,
               enum_values("TerrainType", ["None", "StaticFittedMesh", "MovingFittedMesh", "EB", "ImmersedForcing"]),
               describe=_DESCRIBE)
    sec.blank()

    sec.note("Terrain height. ERF takes the first that is set: terrain_file_name_nc, terrain_file_name,")
    sec.note("terrain_file_name_USGS, prob.custom_terrain_type. Hills from custom_terrain_type sit at the domain centre.")
    if s.terrain == "file":
        sec.user("erf.terrain_file_name", TERRAIN_FILE, "table of x, y and heights [m]; replace the dummy file")
    else:
        sec.alt("erf.terrain_file_name", TERRAIN_FILE, "nx, ny, the x values, the y values, then heights with y fastest")
    if s.terrain == "file_usgs":
        sec.user("erf.terrain_file_name_USGS", TERRAIN_USGS_FILE, "'lon lat', 'nx ny', then rows 'x y height'; replace the dummy file")
    else:
        sec.alt("erf.terrain_file_name_USGS", TERRAIN_USGS_FILE, "'lon lat', 'nx ny', then rows 'x y height', x fastest")
    sec.alt("erf.terrain_file_name_nc", "terrain.nc", "NetCDF: height(y, x) with 1-D x and y [m]; needs a NetCDF build")
    if s.terrain == "analytic":
        _analytic(sec, s)
    else:
        sec.alt("prob.custom_terrain_type", "WoA", "analytic hill: WoA, ScharMountain, Cos4Hill, Hemisphere, HalfCylinder, RaisedFlat")
        sec.alt("prob.hmax", s.hill_height, "[m] height (WoA, Cos4Hill)")
        sec.alt("prob.L", s.hill_width, "[m] half-width (WoA, Cos4Hill)")
    sec.blank()

    if s.method == "fitted" or (s.method == "none" and s.stretched):
        sec.note("Terrain-fitted mesh: how the coordinate surfaces relax to flat with height.")
        sec.set("erf.terrain_smoothing", 0, "0 basic terrain-following")
        sec.alt("erf.terrain_smoothing", 1, "smoothed terrain-following; with erf.terrain_gamma_m (default 0.5)")
        sec.alt("erf.terrain_smoothing", 2, "Sullivan terrain-following")
        if s.dycore == "anelastic" and s.method == "fitted":      # ERF ignores the key off fitted terrain
            sec.set("erf.terrain_poisson_solver", "gmres_fft", "default; needs an FFT build (ERF_ENABLE_FFT)")
            sec.alt("erf.terrain_poisson_solver", "mlmg", "multigrid; no FFT build needed")
            sec.alt("erf.terrain_poisson_solver", "gmres_mlmg", "GMRES with a multigrid preconditioner")
    if s.stretched:
        top = s.dz0 * s.stretch ** (s.n[2] - 1)
        sec.set("erf.initial_dz", s.dz0, "[m] first-cell height; ERF solves the ratio (%.4g here, top cell %.4g m)"
                % (s.stretch, top))
        sec.alt("erf.grid_stretching_ratio", 1.05, "fixed ratio between cells (>= 1); the top then need not meet prob_hi z")
    elif s.method in ("none", "fitted"):
        sec.alt("erf.initial_dz", 0.25 * s.dx[2], "[m] vertical stretching from this first-cell height (uniform dz when unset)")
        sec.alt("erf.grid_stretching_ratio", 1.05, "fixed ratio between cells; needs erf.initial_dz")
    if s.method in ("none", "fitted"):
        sec.note("erf.terrain_z_levels = z0 z1 ... lists all %d level heights instead (first = prob_lo z)." % (s.n[2] + 1))

    if s.method == "eb":
        sec.set("eb2.geometry", "terrain", "the EB surface is the terrain height (default); also plane, box, sphere")
        wall = {"slip": "SlipWall", "noslip": "NoSlipWall", "most": "SurfaceLayer"}[s.surface]
        sec.choice("erf.eb_boundary_type", wall, ["SlipWall", "NoSlipWall", "SurfaceLayer"],
                   describe={"NoSlipWall": "(default)",
                             "SurfaceLayer": "required with zlo.type = surface_layer; the run crashes without it"})
        sec.alt("eb2.small_volfrac", 1.0e-4, "cut cells below this volume fraction are covered")
    if s.method == "immersed" or s.buildings_on:
        _immersed(sec, s)

    sec.blank()
    if s.buildings_on:
        sec.set("erf.buildings_type", "ImmersedForcing", "buildings as immersed bodies")
        sec.user("erf.buildings_file_name", BUILDINGS_FILE,
                 "same table layout as erf.terrain_file_name, heights [m] above the ground; replace the dummy file")
    elif s.method != "eb":
        sec.alt("erf.buildings_type", "ImmersedForcing", "buildings as immersed bodies (default None)")
        sec.alt("erf.buildings_file_name", BUILDINGS_FILE, "table of building heights, same layout as erf.terrain_file_name")


def _analytic(sec, s) -> None:
    if s.hill == "woa":
        sec.set("prob.custom_terrain_type", "WoA", "h = hmax / (1 + (r/L)^2)")
        sec.user("prob.hmax", s.hill_height, "[m] hill height")
        sec.user("prob.L", s.hill_width, "[m] half-width")
        sec.set("prob.dir", 0 if s.ridge else 2, "0 ridge across x, 1 ridge across y, 2 round hill")
        sec.alt("prob.z_offset", 0.0, "[m] height of the surrounding ground (default 0)")
        sec.alt("prob.custom_terrain_type", "Cos4Hill", "compact round hill; prob.hmax, prob.L")
        sec.alt("prob.custom_terrain_type", "ScharMountain", "prob.a, prob.hm, prob.lambda")
    elif s.hill == "cos4":
        sec.set("prob.custom_terrain_type", "Cos4Hill", "round hill, zero beyond r = 4 L")
        sec.user("prob.hmax", s.hill_height, "[m] hill height")
        sec.user("prob.L", s.hill_width, "[m] width scale; the hill ends at 4 L")
        sec.alt("prob.z_offset", 0.0, "[m] height of the surrounding ground (default 0)")
        sec.alt("prob.custom_terrain_type", "WoA", "prob.hmax, prob.L, prob.dir")
    else:
        sec.set("prob.custom_terrain_type", "ScharMountain", "wavy ridge across x")
        sec.user("prob.hm", s.hill_height, "[m] height (default 250)")
        sec.user("prob.a", s.hill_width, "[m] half-width of the envelope (default 5000)")
        sec.user("prob.lambda", 0.8 * s.hill_width, "[m] wavelength of the ripples (default 4000)")
        sec.alt("prob.custom_terrain_type", "WoA", "prob.hmax, prob.L, prob.dir")


def _immersed(sec, s) -> None:
    sec.blank()
    sec.note("Immersed forcing: wall treatment on the immersed surface.")
    wall = s.surface == "most"
    sec.block(wall, [("erf.if_use_most", "true", "log-law wall stress on the immersed surface (default false)"),
                     ("erf.if_z0", s.z0, "[m] roughness length of the immersed surface (default 0.1)")])
    heat = s.surface_heat if (wall and s.method == "immersed") else "adiabatic"
    sec.note("Give at most one of if_surf_temp_flux, if_init_surf_temp, if_Olen.")
    if heat == "flux":
        sec.user("erf.if_surf_temp_flux", s.heat_flux, "[K m/s] kinematic heat flux at the immersed surface")
    else:
        sec.alt("erf.if_surf_temp_flux", s.heat_flux, "[K m/s] kinematic heat flux at the immersed surface")
    if heat == "temp":
        sec.user("erf.if_init_surf_temp", s.theta, "[K] temperature of the immersed surface")
    else:
        sec.alt("erf.if_init_surf_temp", s.theta, "[K] temperature of the immersed surface")
    sec.alt("erf.if_surf_heating_rate", -0.25, "[K/h] with if_init_surf_temp")
    sec.alt("eb2.small_volfrac", 0.005, "cells below this fluid volume fraction are treated as solid (default)")
    sec.alt("erf.if_Cd_momentum", 500.0, "drag coefficient inside the body (default 500; 50 with an anelastic level)")
    sec.alt("erf.if_Cd_scalar", 50.0, "(default 50; 5 with an anelastic level)")
