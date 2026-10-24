"""Turbulence closures, forcing, moisture, terrain and damping.

Which spellings ERF accepts, and whether case matters, follows how each key is
read: ``les_type`` and ``rans_type`` go through ``ParmParse::query`` on an
``AMREX_ENUM`` (exact names), ``molec_diff_type`` and the advection types are
compared as strings (exact), while ``pbl_type``, ``moisture_model``,
``terrain_type``, ``init_type`` and ``abl_driver_type`` use
``query_enum_case_insensitive``.  The enum name lists are read from
``Source/DataStructs`` when the repository is at hand, so new members are picked
up; the lists below are the fallback.
"""

import os
import re
from typing import Dict, List

from ..checker import Checker, register
from .keys import _find_docs

LES_TYPES = ["None", "Smagorinsky", "Smagorinsky2D", "Deardorff"]
RANS_TYPES = ["None", "kEqn"]
PBL_TYPES = ["None", "MYJ", "MYNN25", "MYNNEDMF", "YSU", "YSUNew", "MRF", "SHOC", "EAMXX_SHOC", "NATIVE_SHOC"]
MOLEC_TYPES = ["None", "Constant", "ConstantAlpha"]
# ERF_AdvStruct.H string list
ADV_TYPES = ["Centered_2nd", "Upwind_3rd", "Upwind_3rd_SL", "Centered_4th", "Upwind_5th", "Centered_6th",
             "WENO3", "WENO5", "WENO7", "WENOZ3", "WENOZ5", "WENOZ7", "WENOMZQ3",
             "Blended_3rd4th", "Blended_5th6th"]
TERRAIN_TYPES = ["None", "StaticFittedMesh", "MovingFittedMesh", "EB", "ImmersedForcing"]
MOISTURE_MODELS = ["SAM", "SAM_NoIce", "SAM_NoPrecip_NoIce", "Kessler", "Kessler_NoRain", "SatAdj", "Morrison",
                   "Morrison_NoIce", "WSM6", "WDM6", "SuperDroplets", "MoistNoCondensation", "SBM", "None"]
DRIVERS = ["None", "PressureGradient", "GeostrophicWind"]

_ENUM_RE = re.compile(r"AMREX_ENUM\(\s*(\w+)\s*,([^)]*)\)", re.S)
_ENUMS: Dict[str, List[str]] = {}


def source_enum(name: str, fallback: List[str]) -> List[str]:
    """Members of ``AMREX_ENUM(name, ...)`` under Source/DataStructs, or the fallback."""
    if not _ENUMS:
        _ENUMS["_loaded"] = []
        docs = _find_docs()
        if docs:
            root = os.path.abspath(os.path.join(os.path.dirname(docs), "..", ".."))
            d = os.path.join(root, "Source", "DataStructs")
            if os.path.isdir(d):
                for fn in os.listdir(d):
                    if fn.endswith(".H"):
                        try:
                            text = open(os.path.join(d, fn), errors="replace").read()
                        except OSError:
                            continue
                        for m in _ENUM_RE.finditer(text):
                            _ENUMS[m.group(1)] = [x.strip() for x in m.group(2).split(",") if x.strip()]
    return _ENUMS.get(name, fallback)


@register
class TurbulenceChecker(Checker):
    name = "physics.turbulence"

    def run(self) -> None:
        c, d = self.ctx, self.deck
        les_raw = c.per_level_str("erf.les_type", "None")
        rans_raw = c.per_level_str("erf.rans_type", "None")
        pbl_raw = c.per_level_str("erf.pbl_type", "None")
        molec_raw = self.str_("erf.molec_diff_type", "None") or "None"
        self._enum("erf.les_type", les_raw, source_enum("LESType", LES_TYPES), exact=True)
        self._enum("erf.rans_type", rans_raw, source_enum("RANSType", RANS_TYPES), exact=True)
        self._enum("erf.pbl_type", pbl_raw, source_enum("PBLType", PBL_TYPES), exact=False)
        self._enum("erf.molec_diff_type", [molec_raw], MOLEC_TYPES, exact=True)
        les = [v.lower() for v in les_raw]
        rans = [v.lower() for v in rans_raw]
        pbl = [v.lower() for v in pbl_raw]
        molec = molec_raw.lower()
        for lev in range(c.nlev):
            on = [n for n, v in (("les_type", les[lev]), ("rans_type", rans[lev]), ("pbl_type", pbl[lev])) if v != "none"]
            if len(on) > 1 and not (set(on) == {"les_type", "pbl_type"}):
                self.warn("erf.%s" % on[0], "level %d combines %s" % (lev, " and ".join(on)),
                          "Use one closure per level (LES + PBL is the only pair ERF combines).")
        if all(v == "none" for v in les + rans + pbl) and molec == "none":
            self.info("erf.les_type", "no turbulence closure and no molecular diffusion: the run is inviscid apart from the advection scheme")
        if molec == "constant" and not d.has("erf.dynamic_viscosity"):
            self.advise("erf.dynamic_viscosity", "molec_diff_type = Constant uses the default viscosity; set erf.dynamic_viscosity (and alpha_T / alpha_C with ConstantAlpha)")
        if molec == "none":
            for k in ("erf.dynamic_viscosity", "erf.alpha_T", "erf.alpha_C", "erf.rho0_trans"):
                if d.has(k):
                    self.warn(k, "ignored: erf.molec_diff_type is None")
        if not any(v.startswith("smagorinsky") for v in les):
            for k in ("erf.Cs", "erf.smag2d_slope_limiter", "erf.smag2d_kh_cap"):
                if d.has(k):
                    self.warn(k, "ignored without a Smagorinsky LES")
        if "deardorff" not in les:
            for k in ("erf.Ck", "erf.Ce", "erf.Ce_wall", "erf.sigma_k"):
                if d.has(k) and "keqn" not in rans:
                    self.warn(k, "ignored without erf.les_type = Deardorff or erf.rans_type = kEqn")
        if all(v == "none" for v in pbl):
            for k in d.keys_with_prefix("erf.pbl_"):
                if k != "erf.pbl_type":
                    self.warn(k, "ignored: no PBL scheme is selected")
        if "keqn" in rans and self.fixed_dt_none() and not self.bool_("erf.diffusive_dt_limit", False):
            self.advise("erf.diffusive_dt_limit", "k-equation RANS with an adaptive step: the eddy viscosity can grow until the explicit diffusion limits the step",
                        "Set erf.diffusive_dt_limit = true.")
        self._advection()

    def fixed_dt_none(self) -> bool:
        return not self.deck.has("erf.fixed_dt")

    def _enum(self, key, values, allowed, exact) -> None:
        """Report a value ERF would not accept.  ``exact``: the code compares the string
        (or uses ParmParse::query on an enum), so the case must match too."""
        lower = {a.lower(): a for a in allowed}
        for v in values:
            if v in allowed:
                continue
            if v.lower() in lower:
                if exact:
                    self.error(key, "%r differs in case from %s; ERF compares the string exactly and stops" % (v, lower[v.lower()]),
                               "Write %s = %s." % (key, lower[v.lower()]))
                continue
            self.error(key, "%r is not a known value" % v, "Choose from: " + ", ".join(allowed) + ".")
            return

    def _advection(self) -> None:
        d = self.deck
        moist = (self.str_("erf.moisture_model", "None") or "None").lower() != "none"
        for key in ("erf.dycore_horiz_adv_type", "erf.dycore_vert_adv_type", "erf.dryscal_horiz_adv_type",
                    "erf.dryscal_vert_adv_type", "erf.moistscal_horiz_adv_type", "erf.moistscal_vert_adv_type"):
            v = self.str_(key)
            if v is None:
                continue
            if v not in ADV_TYPES:
                lower = {a.lower(): a for a in ADV_TYPES}
                if v.lower() in lower:
                    self.error(key, "%r differs in case from %s; ERF compares the string exactly and stops" % (v, lower[v.lower()]),
                               "Write %s = %s." % (key, lower[v.lower()]))
                else:
                    self.error(key, "%r is not an advection scheme" % v, "Choose from: " + ", ".join(ADV_TYPES) + ".")
            if key.startswith("erf.moistscal") and not moist:
                self.warn(key, "ignored: no moisture model is active")
            if key.startswith("erf.dycore") and v.lower().startswith("weno") and self.ctx.terrain.lower() == "eb":
                self.error(key, "WENO momentum advection is not allowed with erf.terrain_type = EB")
        h = (self.str_("erf.dycore_horiz_adv_type", "Upwind_3rd") or "").lower()
        nd = self.real("erf.num_diff_coeff", 0.0)
        if nd is not None and not (0 <= nd <= 1):
            self.error("erf.num_diff_coeff", "must be in [0, 1]")
        if h.startswith("centered") and (nd or 0) == 0:
            self.advise("erf.dycore_horiz_adv_type", "%s has no upwinding; the docs recommend about 5%% sixth-order numerical diffusion for stability" % h,
                        "Set erf.num_diff_coeff = 0.05, or use an upwind scheme (Upwind_5th horizontal, Upwind_3rd vertical is the WRF-like choice).")
        for key in ("erf.dycore_horiz_upw_frac", "erf.dycore_vert_upw_frac", "erf.dryscal_horiz_upw_frac",
                    "erf.dryscal_vert_upw_frac", "erf.moistscal_horiz_upw_frac", "erf.moistscal_vert_upw_frac"):
            if d.has(key):
                scheme = (self.str_(key.replace("_upw_frac", "_adv_type"), "Upwind_3rd") or "").lower()
                if not scheme.startswith("upwind"):
                    self.warn(key, "the upwind fraction is only used by Upwind_3rd / Upwind_5th schemes")


@register
class ForcingChecker(Checker):
    name = "physics.forcing"

    def run(self) -> None:
        c, d = self.ctx, self.deck
        cor = self.bool_("erf.use_coriolis", False)
        driver = (self.str_("erf.abl_driver_type", "None") or "None").lower()
        if driver not in [x.lower() for x in source_enum("ABLDriverType", DRIVERS)]:
            self.error("erf.abl_driver_type", "%r is not a driver type" % driver, "None, PressureGradient or GeostrophicWind.")
        if cor and not d.has("erf.latitude") and not d.has("erf.coriolis_3d"):
            self.warn("erf.latitude", "use_coriolis = true with the default latitude of 90 degrees (f = 2 Omega)",
                      "Set erf.latitude for the site.")
        if not cor:
            for k in ("erf.latitude", "erf.rotational_time_period", "erf.coriolis_3d", "erf.abl_geo_wind", "erf.abl_geo_wind_table"):
                if d.has(k):
                    self.warn(k, "ignored: erf.use_coriolis is false")
        if driver == "geostrophicwind" and not cor and not d.has("erf.abl_geo_forcing"):
            self.warn("erf.abl_driver_type", "GeostrophicWind with erf.use_coriolis = false: the geostrophic wind only enters through the Coriolis term, so no forcing is applied",
                      "Set erf.use_coriolis = true (with erf.latitude and erf.abl_geo_wind), or erf.abl_driver_type = None.")
        if driver == "geostrophicwind" and not (d.has("erf.abl_geo_wind") or d.has("erf.abl_geo_wind_table") or d.has("erf.abl_geo_forcing")):
            self.warn("erf.abl_geo_wind", "GeostrophicWind driver without erf.abl_geo_wind, abl_geo_wind_table or abl_geo_forcing: zero forcing")
        if driver != "pressuregradient" and d.has("erf.abl_pressure_grad"):
            self.warn("erf.abl_pressure_grad", "only used with erf.abl_driver_type = PressureGradient")
        if driver != "geostrophicwind":
            for k in ("erf.abl_geo_wind", "erf.abl_geo_forcing", "erf.abl_geo_wind_table"):
                if d.has(k):
                    self.warn(k, "only used with erf.abl_driver_type = GeostrophicWind")
        if d.has("erf.abl_geo_wind") and d.has("erf.abl_geo_wind_table"):
            self.warn("erf.abl_geo_wind", "ignored when erf.abl_geo_wind_table is given")
        for k in ("erf.abl_pressure_grad", "erf.abl_geo_forcing", "erf.abl_geo_wind"):
            v = self.reals(k)
            if v and len(v) != 3:
                self.error(k, "needs three values")
        # damping
        rayleigh = any(self.bool_("erf.rayleigh_damp_%s" % s, False) for s in "UVWT")
        zd = self.real("erf.rayleigh_zdamp", 500.0)
        if rayleigh and c.extent and zd is not None and zd >= c.extent[2]:
            self.warn("erf.rayleigh_zdamp", "damping depth %g m is at least the domain height %g m, so the whole column is damped" % (zd, c.extent[2]),
                      "Set erf.rayleigh_zdamp to the depth of the absorbing layer (the default is 500 m).")
        if not rayleigh:
            for k in ("erf.rayleigh_zdamp", "erf.rayleigh_dampcoef", "erf.rayleigh_damping_type"):
                if d.has(k):
                    self.warn(k, "ignored: no erf.rayleigh_damp_U/V/W/T is true")
        elif c.init_mode == "ideal" and (self.str_("zhi.type", "") or "").lower() == "slipwall" and not rayleigh:
            pass
        st = (self.str_("erf.sponge_type", "Standard") or "").lower()
        if self.real("erf.sponge_strength", 0.0) == 0 and d.keys_with_prefix("erf.sponge_"):
            others = [k for k in d.keys_with_prefix("erf.sponge_") if k not in ("erf.sponge_strength", "erf.sponge_type")]
            if others:
                self.warn(others[0], "erf.sponge_* keys are set but erf.sponge_strength is 0, so the sponge is off")
        if st == "none" and d.has("erf.sponge_strength"):
            self.warn("erf.sponge_strength", "ignored: erf.sponge_type = None")
        if c.init_mode == "ideal" and c.extent and self.bool_("erf.use_gravity", False) and not rayleigh \
                and (self.str_("zhi.type", "") or "").lower() == "slipwall" and st != "standard" and c.extent[2] > 2000:
            self.advise("erf.rayleigh_damp_W", "a stratified domain with a rigid lid and no upper damping reflects gravity waves",
                        "Add erf.rayleigh_damp_W = true with erf.rayleigh_zdamp set to the top 20-30 % of the domain.")


@register
class MoistureChecker(Checker):
    name = "physics.moisture"

    def run(self) -> None:
        c, d = self.ctx, self.deck
        m = (self.str_("erf.moisture_model", "None") or "None").lower()
        models = source_enum("MoistureType", MOISTURE_MODELS)
        if m not in [x.lower() for x in models]:
            self.error("erf.moisture_model", "%r is not a moisture model; ERF stops" % m, "Choose from: " + ", ".join(models) + ".")
        if m == "none":
            for k in d.keys():
                if k.startswith("erf.most.surf_moist") or k in ("erf.most.Cq", "erf.moisture_tight_coupling", "erf.implicit_moisture_diffusion"):
                    self.advise(k, "ignored: no moisture model is active", "Remove the key.")
            return
        variable_dz = d.has("erf.grid_stretching_ratio") or d.has("erf.terrain_z_levels") or c.fitted_terrain
        if m.startswith("sam") and variable_dz:
            self.error("erf.moisture_model", "SAM is not correct with variable dz; ERF stops", "Use Kessler or Morrison, or a uniform grid.")
        if not self.bool_("erf.use_gravity", False) and c.init_mode == "ideal":
            self.advise("erf.use_gravity", "a moisture model without gravity gives no buoyancy from latent heating")


@register
class TerrainChecker(Checker):
    name = "physics.terrain"

    def run(self) -> None:
        c, d = self.ctx, self.deck
        t = c.terrain.lower()
        if t not in [x.lower() for x in source_enum("TerrainType", TERRAIN_TYPES)] + ["static", "moving"]:
            self.error("erf.terrain_type", "%r is not a terrain type" % c.terrain, "None, StaticFittedMesh, MovingFittedMesh, EB or ImmersedForcing.")
        src = [k for k in ("erf.terrain_file_name", "erf.terrain_file_name_nc", "erf.terrain_file_name_USGS", "prob.custom_terrain_type") if d.has(k)]
        if len(src) > 1:
            self.warn(src[1], "several terrain sources; the NetCDF file wins over the text file, which wins over USGS and prob.custom_terrain_type")
        if t == "none" and src:
            self.warn(src[0], "a terrain source is given but erf.terrain_type = None, so the mesh stays flat",
                      "Set erf.terrain_type = StaticFittedMesh (or EB / ImmersedForcing).")
        if t in ("staticfittedmesh", "static") and not src and c.init_mode == "ideal" and not d.has("erf.grid_stretching_ratio") \
                and not d.has("erf.terrain_z_levels"):
            self.info("erf.terrain_type", "fitted mesh without a terrain file: the problem's own init_custom_terrain (or flat terrain) is used")
        for k in ("erf.use_rotate_surface_flux", "erf.implicit_terrain_metric", "erf.smag2d_slope_limiter"):
            if d.has(k) and not c.fitted_terrain:
                self.warn(k, "only used on a terrain-fitted mesh")
        gs = self.real("erf.grid_stretching_ratio", 0.0)
        if gs is not None and gs != 0 and gs < 1:
            self.error("erf.grid_stretching_ratio", "must be 0 or >= 1")
        if gs and gs >= 1 and not d.has("erf.initial_dz"):
            self.error("erf.initial_dz", "needed with erf.grid_stretching_ratio")
        if d.has("erf.initial_dz") and not (gs and gs >= 1):
            self.warn("erf.initial_dz", "ignored without erf.grid_stretching_ratio")
        if gs and gs > 1.2:
            self.advise("erf.grid_stretching_ratio", "stretching by %g per cell is steep; the vertical truncation error grows with (ratio - 1)" % gs,
                        "Ratios of 1.02-1.1 are usual.")
        zl = self.reals("erf.terrain_z_levels")
        if zl and c.n_cell and len(zl) != c.n_cell[2] + 1:
            self.error("erf.terrain_z_levels", "%d levels given; n_cell z + 1 = %d are needed" % (len(zl), c.n_cell[2] + 1))
        if zl and any(zl[i + 1] <= zl[i] for i in range(len(zl) - 1)):
            self.error("erf.terrain_z_levels", "levels must increase")
        if zl and gs:
            self.info("erf.grid_stretching_ratio", "ignored when erf.terrain_z_levels is given")
        if (gs or zl) and c.n_cell and c.extent and c.dz_min:
            # the stretched column must end at prob_hi
            if gs and gs > 1 and d.has("erf.initial_dz"):
                dz0 = self.real("erf.initial_dz")
                top = dz0 * (gs ** c.n_cell[2] - 1) / (gs - 1)
                if abs(top - c.extent[2]) > 1e-4:
                    sev = self.warn if abs(top - c.extent[2]) > 0.01 * c.extent[2] else self.advise
                    sev("erf.initial_dz", "initial_dz %g m stretched by %g over %d cells puts the mesh top at %.4g m, not at prob_hi z = %g m; ERF only prints a note and keeps the stretched levels, so the physical domain is %s than the deck says" % (
                        dz0, gs, c.n_cell[2], top, c.extent[2], "taller" if top > c.extent[2] else "shorter"),
                        "Set geometry.prob_hi z = %.6g, or choose initial_dz / ratio so the sum lands on the domain top." % top)
        if t in ("movingfittedmesh", "moving") and c.substepping == "Implicit" and c.max_level and c.max_level > 0:
            self.info("erf.terrain_type", "moving terrain with refinement is lightly tested")
