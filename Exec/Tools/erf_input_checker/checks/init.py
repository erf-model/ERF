"""Initialisation checks: one class for idealised set-ups, one for real data."""

from ..checker import Checker, register

INIT_TYPES = {"none", "input_sounding", "ncfile", "wrfinput", "metgrid", "uniform",
              "constantdensity", "constantdensitylineartheta", "isentropic",
              "moistbasestate", "hindcast"}


class InitChecker(Checker):
    name = "init"

    def common(self) -> None:
        t = self.ctx.init_type
        if t.lower() not in INIT_TYPES:
            self.error("erf.init_type", "%r is not an initialisation type" % t,
                       "Choose from: " + ", ".join(sorted(INIT_TYPES)) + ".")
        if not self.deck.has("erf.init_type"):
            self.warn("erf.init_type", "not set: the problem's own init_custom_pert fills the state and nothing enforces hydrostatic balance",
                      "Set erf.init_type explicitly (Uniform, Isentropic, Input_Sounding, ...).")


@register
class IdealInitChecker(InitChecker):
    name = "init.ideal"

    @classmethod
    def applies(cls, ctx):
        return ctx.init_mode == "ideal"

    def run(self) -> None:
        self.common()
        c, d = self.ctx, self.deck
        t = c.init_type.lower()
        grav = self.bool_("erf.use_gravity", False)
        if t in ("input_sounding", "isentropic", "moistbasestate", "constantdensitylineartheta") and not grav:
            self.warn("erf.use_gravity", "init_type = %s builds a hydrostatic column but erf.use_gravity is false, so the column is not in balance" % c.init_type,
                      "Set erf.use_gravity = true for any stratified atmosphere.")
        if t == "input_sounding":
            files = d.get_strs("erf.input_sounding_file") or ["input_sounding"]
            times = self.reals("erf.input_sounding_time")
            if times and len(times) != len(files):
                self.error("erf.input_sounding_time", "%d times for %d sounding files" % (len(times), len(files)))
            if len(files) > 1 and not self.bool_("erf.nudging_from_input_sounding", False):
                self.advise("erf.input_sounding_file", "several soundings are only used over time with erf.nudging_from_input_sounding = true")
        for key in ("erf.nc_init_file_0", "erf.nc_bdy_file", "erf.nc_low_file", "erf.real_width", "erf.use_real_bcs",
                    "erf.use_wrf_bdy_density", "erf.use_wrf_bdy_qc_qi", "erf.wrfinput_zlevels_from_file"):
            if d.has(key):
                self.warn(key, "only read for init_type = WRFInput / Metgrid / NCFile",
                          "Remove the key, or set erf.init_type to a real-data type.")
        if t in ("uniform", "constantdensity") and grav and not d.has("prob.T_0"):
            self.info("erf.init_type", "a uniform state with gravity relies on the problem's prob.* keys for the base state")


@register
class RealInitChecker(InitChecker):
    name = "init.real"

    @classmethod
    def applies(cls, ctx):
        return ctx.init_mode == "real"

    def run(self) -> None:
        self.common()
        c, d = self.ctx, self.deck
        t = c.init_type.lower()
        if not d.get_strs("erf.nc_init_file_0") and t != "hindcast":
            self.error("erf.nc_init_file_0", "init_type = %s needs a level-0 input file; ERF stops" % c.init_type)
        if t in ("wrfinput", "metgrid"):
            if not c.use_real_bcs and not self.bool_("erf.input_bndry_planes", False):
                self.warn("erf.use_real_bcs", "real-data initial state with idealised lateral boundaries",
                          "Set erf.use_real_bcs = true with erf.nc_bdy_file = wrfbdy_d01 (or metgrid boundary files) to force the edges with the outer model.")
            if t == "wrfinput" and not d.has("erf.nc_bdy_file") and self.bool_("erf.use_real_bcs", False):
                self.error("erf.nc_bdy_file", "use_real_bcs with WRFInput needs the wrfbdy file")
            if c.use_real_bcs:
                rw = self.int_("erf.real_width", 0)
                if rw is not None and rw < 1:
                    self.warn("erf.real_width", "= %d: no relaxation zone, only the outer row is set from the boundary file" % rw,
                              "WRF uses spec_bdy_width = 5; erf.real_width must stay below the width stored in wrfbdy.")
                for i, dim in enumerate("xy"):
                    if c.periodic[i]:
                        self.error("geometry.is_periodic", "%s is periodic while lateral boundaries come from the boundary file" % dim)
            if c.terrain.lower() not in ("staticfittedmesh", "static"):
                self.error("erf.terrain_type", "= %s; ERF stops ('Only terrain_type = StaticFittedMesh are allowed with init_type = WRFInput or Metgrid')" % c.terrain,
                           "Set erf.terrain_type = StaticFittedMesh.")
            if (self.str_("erf.moisture_model", "None") or "None").lower() == "none":
                if t == "wrfinput":
                    self.error("erf.moisture_model", "None with WRFInput; ERF stops ('Can't have moisture_type = None with init_type = WRFInput')",
                               "Set erf.moisture_model (e.g. Kessler, SAM, SatAdj).")
                else:
                    self.advise("erf.moisture_model", "met_em fields carry humidity; without a moisture model it is dropped",
                                "Set erf.moisture_model unless a dry run is intended.")
            if not self.bool_("erf.use_gravity", False):
                self.warn("erf.use_gravity", "real-data runs need gravity for the hydrostatic base state")
        if d.has("start_time") and t in ("wrfinput", "metgrid"):
            self.error("start_time", "the start time is read from the input file; ERF stops if start_time is also set")
        sd = self.str_("start_datetime")
        if sd is not None:
            s = " ".join(d.get_strs("start_datetime"))
            if len(s) not in (16, 19) or s[4] != "-" or s[7] != "-":
                self.error("start_datetime", "%r is not YYYY-MM-DD HH:MM[:SS] (quote it, the value has a space)" % s)
        if d.has("amr.restart") and self.bool_("erf.write_erfbdy", False):
            self.error("erf.write_erfbdy", "cannot be true on a restart; ERF stops")
        if c.max_level and c.max_level > 0:
            has_fine_file = any(d.has("erf.nc_init_file_%d" % lev) for lev in range(1, c.nlev))
            if not has_fine_file:
                self.info("erf.nc_init_file_1", "fine levels are interpolated from level 0 (no fine-level input file)")
