"""Plotfiles, checkpoints and restarts."""

from ..checker import Checker, register


@register
class OutputChecker(Checker):
    name = "output"

    def run(self) -> None:
        d = self.deck
        if d.has("erf.plotfile_type") and (d.has("erf.plotfile_type_1") or d.has("erf.plotfile_type_2")):
            self.error("erf.plotfile_type", "set either plotfile_type or plotfile_type_1/_2, not both; ERF stops")
        for n in ("1", "2"):
            i, p = self.int_("erf.plot_int_%s" % n, -1), self.real("erf.plot_per_%s" % n, -1.0)
            if (i or 0) > 0 and (p or 0) > 0:
                self.error("erf.plot_int_%s" % n, "plot_int and plot_per are both set; ERF stops")
            if ((i or 0) > 0 or (p or 0) > 0) and not d.has("erf.plot_vars_%s" % n):
                self.warn("erf.plot_int_%s" % n, "plotfile %s is written but erf.plot_vars_%s is not set, so it holds no fields" % (n, n),
                          "Set erf.plot_vars_%s = density x_velocity y_velocity z_velocity theta ..." % n)
            if d.has("erf.plot_vars_%s" % n) and (i or 0) <= 0 and (p or 0) <= 0:
                if d.has("erf.plot_int_%s" % n) or d.has("erf.plot_per_%s" % n):
                    self.info("erf.plot_vars_%s" % n, "plotfile %s is switched off (plot_int_%s = %s), so the variable list is unused" % (n, n, i))
                else:
                    self.warn("erf.plot_vars_%s" % n, "variables are listed but neither plot_int_%s nor plot_per_%s is set, so nothing is written" % (n, n),
                              "Set erf.plot_int_%s or erf.plot_per_%s." % (n, n))
        for n in ("1", "2"):
            i, p = self.int_("erf.plot2d_int_%s" % n, -1), self.real("erf.plot2d_per_%s" % n, -1.0)
            if (i or 0) > 0 and (p or 0) > 0:
                self.error("erf.plot2d_int_%s" % n, "plot2d_int and plot2d_per are both set; ERF stops")
        ci = self.int_("amr.check_int", self.int_("erf.check_int", -1))
        cp = self.real("amr.check_per", self.real("erf.check_per", -1.0))
        if (ci or 0) > 0 and (cp or 0) > 0:
            self.warn("amr.check_int", "check_int and check_per are both set; ERF writes a checkpoint whenever either is due",
                      "Keep one of them.")
        max_step = self.int_("max_step", -1)
        fixed_dt = self.ctx.fixed_dt
        stop_time = self.real("stop_time")
        nsteps = None
        if max_step is not None and max_step > 0:
            nsteps = max_step
        if fixed_dt and stop_time and stop_time > 0:
            est = stop_time / fixed_dt
            nsteps = min(nsteps, est) if nsteps else est
        if (ci or 0) <= 0 and (cp or 0) <= 0 and nsteps is not None and nsteps > 2000:
            self.advise("amr.check_int", "about %.0f steps and no checkpoints" % nsteps,
                        "Set amr.check_int (or check_per) so the run can restart after a failure.")
        pi = self.int_("erf.plot_int_1", -1)
        if max_step is not None and 0 < max_step and (pi or 0) > max_step:
            self.info("erf.plot_int_1", "plot_int_1 = %d exceeds max_step = %d; only the initial and the final plotfile are written" % (pi, max_step))
        for key in ("erf.plot_per_1", "erf.plot_per_2", "amr.check_per"):
            per = self.real(key, -1.0)
            if fixed_dt and per and per > 0:
                r = per / fixed_dt
                if abs(r - round(r)) > 1e-6:
                    self.advise(key, "%g s is not a multiple of fixed_dt = %g s; output lands on the step after the target time" % (per, fixed_dt))
        rs = self.str_("amr.restart") or self.str_("erf.restart")
        if rs:
            if self.deck.resolve_path(rs) is None:
                self.warn("amr.restart", "checkpoint %r is not next to the deck or in the current directory; the run stops unless an earlier run writes it first" % rs)
            if d.has("erf.init_type") and self.ctx.init_mode == "ideal":
                self.info("erf.init_type", "ignored on restart; the state comes from the checkpoint")
        for bare in ("erf.plot_file", "erf.plot_int", "erf.plot_per", "erf.plot_vars"):
            if d.has(bare):
                self.warn(bare, "ERF reads only the numbered forms %s_1 and %s_2; the bare key is never read" % (bare, bare),
                          "Rename to %s_1 (or remove it)." % bare)
