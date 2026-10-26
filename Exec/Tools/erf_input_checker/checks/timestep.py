"""Time-step checks, one class per integrator.

:class:`TimeStepChecker` holds what both modes share (fixed_dt against cfl,
diffusive limits, run length).  :class:`CompressibleTimeStep` adds the acoustic
substepping rules and :class:`AnelasticTimeStep` the advective limit and the
keys the anelastic integrator ignores.
"""

import math

from ..checker import Checker, register
from ..context import (DEFAULT_CFL, DEFAULT_CHANGE_MAX, DEFAULT_DT_MAX_INITIAL,
                       DEFAULT_INIT_SHRINK, DEFAULT_SUBSTEPPING_CFL)


class TimeStepChecker(Checker):
    """Shared rules; not registered itself."""

    name = "timestep"

    # -- common pieces ----------------------------------------------------
    def common(self) -> None:
        d = self.deck
        self.fixed_dt = self.real("erf.fixed_dt")
        self.cfl = self.real("erf.cfl", DEFAULT_CFL)
        if self.fixed_dt is not None and self.fixed_dt <= 0:
            self.info("erf.fixed_dt", "= %g means 'not set': the step is adaptive" % self.fixed_dt)
            self.fixed_dt = None
        # ERF only asserts cfl > 0 || fixed_dt > 0 (ParameterSanityChecks); it has no upper bound
        if self.cfl is not None and self.cfl <= 0:
            if self.fixed_dt is None:
                self.error("erf.cfl", "= %g with no erf.fixed_dt; ERF stops ('cfl > 0 || fixed_dt > 0')" % self.cfl)
            self.cfl = None
        elif self.cfl is not None and self.cfl > 1 and self.fixed_dt is None:
            self.warn("erf.cfl", "= %g: ERF accepts it, but the RK3 step is only stable for CFL below about 1 (the docs give 0 < cfl <= 1, default 0.8)" % self.cfl,
                      "Use erf.cfl <= 1.")
        if self.fixed_dt is not None and d.has("erf.cfl"):
            self.warn("erf.cfl", "erf.cfl is ignored while erf.fixed_dt is set",
                      "Remove one of them: fixed_dt for a prescribed step, cfl for an adaptive one.")
        cm = self.real("erf.change_max", DEFAULT_CHANGE_MAX)
        if cm is not None and cm < 1:
            self.error("erf.change_max", "= %g; must be >= 1" % cm)
        sh = self.real("erf.init_shrink", DEFAULT_INIT_SHRINK)
        if sh is not None and not (0 < sh <= 1):
            self.error("erf.init_shrink", "= %g; must be in (0, 1]" % sh)
        if sh is not None and sh != 1 and self.fixed_dt is not None:
            self.info("erf.init_shrink", "the first step is init_shrink * fixed_dt = %g s; later steps use fixed_dt" % (sh * self.fixed_dt))
        for key in ("erf.dt_max", "erf.dt_max_initial"):
            v = self.real(key)
            if v is not None and v <= 0:
                self.error(key, "must be > 0")
            if v is not None and self.fixed_dt is not None:
                self.info(key, "ignored with erf.fixed_dt")
        self._diffusive()
        self._run_length()

    def _diffusive(self) -> None:
        d = self.deck
        check = self.bool_("erf.diffusive_dt_check", False)
        limit = self.bool_("erf.diffusive_dt_limit", False)
        if d.has("erf.diffusive_cfl"):
            v = self.real("erf.diffusive_cfl")
            if v is not None and not (0 < v <= 1):
                self.error("erf.diffusive_cfl", "must be in (0, 1]")
            if not check and not limit:
                self.error("erf.diffusive_cfl", "given with neither erf.diffusive_dt_check nor erf.diffusive_dt_limit; ERF stops")
        if limit:
            if self.fixed_dt is not None:
                self.error("erf.diffusive_dt_limit", "cannot change a fixed dt; ERF stops",
                           "Remove erf.fixed_dt (use erf.cfl), or keep only erf.diffusive_dt_check.")
            if not self._has_eddy_closure():
                self.error("erf.diffusive_dt_limit", "needs an eddy-diffusivity closure (les_type, rans_type or pbl_type); ERF stops")
        elif self._has_eddy_closure() and self.fixed_dt is None and not check:
            self.advise("erf.diffusive_dt_limit", "an adaptive dt with an eddy-viscosity closure is bounded only by the acoustic/advective CFL",
                        "Set erf.diffusive_dt_limit = true so the explicit horizontal diffusion cannot drive the step unstable (k-equation RANS decks need it).")

    def _has_eddy_closure(self) -> bool:
        c = self.ctx
        for key in ("erf.les_type", "erf.rans_type", "erf.pbl_type"):
            if any(v.lower() != "none" for v in c.per_level_str(key, "None")):
                return True
        return False

    def _run_length(self) -> None:
        d = self.deck
        max_step = self.int_("max_step")
        stop_time = self.real("stop_time")
        has_stop_dt = d.has("stop_datetime")
        if max_step is None and stop_time is None and not has_stop_dt:
            self.error("max_step", "neither max_step, stop_time nor stop_datetime is set: the run never ends (unless one is given on the command line)")
        if max_step is not None and max_step < 0 and stop_time is None and not has_stop_dt:
            self.error("max_step", "= %d (unlimited) without stop_time: the run never ends" % max_step)
        if self.fixed_dt and stop_time is not None and stop_time > 0:
            nsteps = stop_time / self.fixed_dt
            self.info("stop_time", "about %.0f steps of %g s to reach stop_time = %g s%s" % (
                nsteps, self.fixed_dt, stop_time,
                "" if max_step is None or max_step < 0 or max_step >= nsteps else "; max_step = %d ends the run first" % max_step))
            if abs(nsteps - round(nsteps)) > 1e-6:
                self.advise("stop_time", "stop_time is not a multiple of fixed_dt; the last step is shortened to %g s" % (
                    stop_time - math.floor(nsteps) * self.fixed_dt))
        if d.has("start_datetime") and not (d.has("stop_datetime") or stop_time is not None or max_step is not None):
            self.error("start_datetime", "a start date needs stop_datetime, stop_time or max_step")

    # -- estimates ----------------------------------------------------------
    def horizontal_dx(self):
        c = self.ctx
        if c.dx is None:
            return None
        dims = [i for i in (0, 1) if c.n_cell[i] > 1]
        if not dims:
            return None
        return min(c.dx[i] for i in dims)


@register
class CompressibleTimeStep(TimeStepChecker):
    name = "timestep.compressible"

    @classmethod
    def applies(cls, ctx):
        return not ctx.is_anelastic

    def run(self) -> None:
        self.common()
        c, d = self.ctx, self.deck
        sub = c.substepping
        ratio = self.int_("erf.fixed_mri_dt_ratio")
        fast = self.real("erf.fixed_fast_dt")
        sub_cfl = self.real("erf.substepping_cfl", DEFAULT_SUBSTEPPING_CFL)
        # ERF aborts only for ratio > 0 and odd; 0 (the default) or negative means "not set"
        if ratio is not None and ratio <= 0:
            ratio = None
        if ratio is not None and ratio % 2:
            self.error("erf.fixed_mri_dt_ratio", "= %d; a positive ratio must be even (ERF stops)" % ratio)
        if fast is not None and fast <= 0:
            fast = None
        if sub == "None":
            for key in ("erf.fixed_fast_dt", "erf.fixed_mri_dt_ratio", "erf.substepping_cfl"):
                if d.has(key):
                    self.warn(key, "ignored: erf.substepping_type = None has a single time step")
            self._no_substep_cfl()
            return
        if self.fixed_dt and fast:
            r = self.fixed_dt / fast
            if abs(r - round(r)) > 1e-9:
                self.error("erf.fixed_fast_dt", "fixed_dt / fixed_fast_dt = %g is not an integer (ERF stops: 'must be an even integer')" % r)
            elif int(round(r)) % 2:
                self.advise("erf.fixed_fast_dt", "fixed_dt / fixed_fast_dt = %d is odd; ERF's message asks for an even ratio but its check only tests for an integer, and the substep scheme pairs steps" % int(round(r)),
                            "Use an even ratio to match the documented rule.")
            if ratio and abs(r - ratio) > 1e-9:
                self.error("erf.fixed_mri_dt_ratio", "fixed_dt / fixed_fast_dt = %g disagrees with fixed_mri_dt_ratio = %d (ERF stops: 'Dt is over-specified')" % (r, ratio))
        if sub_cfl is None or sub_cfl <= 0:
            self.error("erf.substepping_cfl", "must be > 0: the substep count is dt divided by this CFL")
            return
        if sub_cfl > 1:
            self.warn("erf.substepping_cfl", "= %g: the acoustic substep is unstable above a CFL of about 1 (WRF practice is 0.5)" % sub_cfl)
        self._substep_cfl(ratio, fast, sub_cfl)

    def _substep_cfl(self, ratio, fast, sub_cfl) -> None:
        c = self.ctx
        dxh = self.horizontal_dx()
        if dxh is None:
            return
        cs = c.sound_speed()
        u = c.options.umax
        dt_wrf = 6.0 * dxh / 1000.0           # the WRF rule of thumb: dt [s] ~ 6 dx [km]
        dt_acoustic = dxh / (u + cs)           # one fast step at acoustic CFL 1
        if c.horizontal_2d() and c.n_cell[0] == 1 and c.n_cell[1] == 1:
            self.info("erf.fixed_dt", "single-column run: the acoustic step is not limited by dx")
            return
        if self.fixed_dt is None:
            if self.cfl is None:
                return
            self.info("erf.cfl", "adaptive step: ERF sets dt from the acoustic CFL %g (horizontal only with implicit substepping) and chooses the substep count%s" % (
                self.cfl, "" if ratio is None else " (every %d substeps)" % ratio))
            dmi = self.real("erf.dt_max_initial", DEFAULT_DT_MAX_INITIAL)
            if dmi is not None and dt_wrf > 5 * dmi:
                self.advise("erf.dt_max_initial", "the first adaptive step is capped at %g s, then grows by change_max per step toward about %.3g s" % (dmi, dt_wrf),
                            "Set erf.dt_max_initial near the expected step to skip the ramp-up.")
            return
        dt = self.fixed_dt
        if fast or ratio:
            fdt = fast if fast else dt / ratio
            acfl = fdt * (u + cs) / dxh
            n = int(round(dt / fdt))
            key = "erf.fixed_fast_dt" if fast else "erf.fixed_mri_dt_ratio"
            if acfl > 1.0:
                self.warn(key, "fast step %g s gives an acoustic CFL of %.2f for |u| = %g m/s on dx = %g m; the substep is unstable above 1 (0.5 recommended)" % (
                    fdt, acfl, u, dxh), "Use fixed_mri_dt_ratio >= %d, or drop the fast-step keys and let ERF choose." % self._even_ceil(dt / (0.5 * dt_acoustic)))
            elif acfl > 0.6:
                self.advise(key, "acoustic CFL of the fast step is %.2f (WRF practice: <= 0.5)" % acfl,
                            "fixed_mri_dt_ratio = %d would give about 0.5." % self._even_ceil(dt / (0.5 * dt_acoustic)))
            else:
                self.info(key, "fast step %g s, acoustic CFL %.2f with %d substeps" % (fdt, acfl, n))
        else:
            n_est = self._even_ceil(dt / (sub_cfl * dt_acoustic))
            self.info("erf.fixed_dt", "ERF will take about %d acoustic substeps per step (substepping_cfl = %g, |u| = %g m/s, c = %.0f m/s)" % (
                n_est, sub_cfl, u, cs))
            if n_est > 20:
                self.warn("erf.fixed_dt", "about %d substeps per step: the slow step %g s is far above the WRF guideline of %.3g s for dx = %g m" % (
                    n_est, dt, dt_wrf, dxh), "Lower erf.fixed_dt toward %.3g s, or use erf.cfl for an adaptive step." % dt_wrf)
            elif n_est > 10:
                self.advise("erf.fixed_dt", "about %d substeps per step; ERF has run LES with 10 or more but stability is not guaranteed" % n_est,
                            "The WRF guideline for dx = %g m is dt ~ %.3g s (4-6 substeps)." % (dxh, dt_wrf))
        adv = dt * u / dxh
        if adv > 1.0:
            self.warn("erf.fixed_dt", "advective CFL of the slow step is %.2f for |u| = %g m/s on dx = %g m (> 1: the RK3 step will go unstable)" % (adv, u, dxh),
                      "Use erf.fixed_dt <= %.3g s, or erf.cfl so ERF adapts dt to the actual wind." % (0.8 * dxh / u))
        elif adv > 0.7:
            self.advise("erf.fixed_dt", "advective CFL of the slow step is %.2f for |u| = %g m/s; little margin for gusts" % (adv, u),
                        "Pass --umax with the largest wind expected; consider erf.cfl instead of a fixed step.")
        self._fixed_vs_cfl_advice(dt, dxh)

    def _no_substep_cfl(self) -> None:
        c = self.ctx
        dxh = self.horizontal_dx()
        if c.dz_min is None:
            return
        dmin = min([x for x in ([dxh] if dxh else []) + [c.dz_min]])
        cs = c.sound_speed()
        if self.fixed_dt is None:
            if self.cfl is None:
                return
            dt_lim = self.cfl * dmin / (c.options.umax + cs)
            self.info("erf.cfl", "no substepping: adaptive dt is about %.3g s (acoustic CFL %g on the smallest cell %g m)" % (dt_lim, self.cfl, dmin))
            return
        acfl = self.fixed_dt * (c.options.umax + cs) / dmin
        if acfl > 1.0:
            self.warn("erf.fixed_dt", "acoustic CFL %.2f on the smallest cell (%g m) without substepping; the explicit step goes unstable above about 1" % (acfl, dmin),
                      "Use erf.fixed_dt <= %.3g s, keep erf.substepping_type = Implicit, or use erf.cfl." % (0.8 * dmin / (c.options.umax + cs)))
        elif acfl > 0.8:
            self.advise("erf.fixed_dt", "acoustic CFL %.2f without substepping leaves little margin" % acfl)
        else:
            self.info("erf.fixed_dt", "acoustic CFL %.2f without substepping (smallest cell %g m)" % (acfl, dmin))

    def _fixed_vs_cfl_advice(self, dt, dxh) -> None:
        c = self.ctx
        reasons = []
        if c.init_mode == "real":
            reasons.append("the wind comes from %s data" % c.init_type)
        if c.fitted_terrain and self._has_terrain_source():
            reasons.append("terrain accelerates the flow over ridges")
        if c.max_level and c.max_level > 0:
            reasons.append("fine levels take dt / ref_ratio steps")
        if reasons:
            self.advise("erf.fixed_dt", "a fixed step cannot react to the flow: " + "; ".join(reasons),
                        "Prefer erf.cfl (e.g. 0.5-0.8) with erf.fixed_dt removed; keep a fixed step only for regression reproducibility.")

    def _has_terrain_source(self) -> bool:
        d = self.deck
        return any(d.has(k) for k in ("erf.terrain_file_name", "erf.terrain_file_name_nc",
                                      "erf.terrain_file_name_USGS", "prob.custom_terrain_type")) \
            or self.ctx.init_mode == "real"

    @staticmethod
    def _even_ceil(x: float) -> int:
        n = int(math.ceil(x - 1e-12))
        n = max(n, 2)
        return n if n % 2 == 0 else n + 1


@register
class AnelasticTimeStep(TimeStepChecker):
    name = "timestep.anelastic"

    @classmethod
    def applies(cls, ctx):
        return ctx.is_anelastic

    def run(self) -> None:
        self.common()
        c, d = self.ctx, self.deck
        for key in ("erf.fixed_fast_dt", "erf.fixed_mri_dt_ratio", "erf.substepping_cfl", "erf.substepping_type"):
            if d.has(key):
                self.warn(key, "ignored: the anelastic integrator has no acoustic substeps",
                          "Remove the key; it only suggests the deck was copied from a compressible case.")
        if d.has("erf.use_lagged_delta_rt"):
            self.info("erf.use_lagged_delta_rt", "ignored by the anelastic integrator")
        if not self.bool_("erf.use_fft", False):
            self.advise("erf.use_fft", "the anelastic projection uses multigrid; the FFT Poisson solver is much faster on a uniform grid or a terrain-fitted mesh with full columns",
                        "Set erf.use_fft = true (needs ERF_ENABLE_FFT / USE_FFT in the build). Every anelastic regression test uses it.")
        dxh = self.horizontal_dx()
        if dxh is None or c.dz_min is None:
            return
        u = c.options.umax
        if self.fixed_dt is None:
            if self.cfl is None:
                return
            self.info("erf.cfl", "adaptive step from the advective CFL %g: about %.3g s for |u| = %g m/s on dx = %g m" % (
                self.cfl, self.cfl * dxh / u, u, dxh))
            dmi = self.real("erf.dt_max_initial", DEFAULT_DT_MAX_INITIAL)
            if dmi is not None and self.cfl * dxh / u > 5 * dmi:
                self.advise("erf.dt_max_initial", "the first step is capped at %g s and then grows by change_max per step" % dmi,
                            "Set erf.dt_max_initial near %.3g s to skip the ramp-up." % (self.cfl * dxh / u))
            return
        adv = self.fixed_dt * u / dxh
        if adv > 1.0:
            self.warn("erf.fixed_dt", "advective CFL %.2f for |u| = %g m/s on dx = %g m (> 1)" % (adv, u, dxh),
                      "Use erf.fixed_dt <= %.3g s or erf.cfl." % (0.8 * dxh / u))
        elif adv > 0.7:
            self.advise("erf.fixed_dt", "advective CFL %.2f for |u| = %g m/s; little margin" % (adv, u))
        else:
            self.info("erf.fixed_dt", "advective CFL %.2f for |u| = %g m/s on dx = %g m" % (adv, u, dxh))
        wcfl = self.fixed_dt * 0.1 * u / c.dz_min
        if wcfl > 1.0:
            self.warn("erf.fixed_dt", "a vertical velocity of a tenth of |u| gives CFL %.2f on the first cell (%g m)" % (wcfl, c.dz_min),
                      "Thin first cells with a fixed step need erf.cfl, or a smaller step.")
        if c.init_mode == "real" or (c.fitted_terrain and CompressibleTimeStep._has_terrain_source(self)):
            self.advise("erf.fixed_dt", "a fixed step cannot react to %s" % (
                "real-data winds" if c.init_mode == "real" else "flow acceleration over terrain"),
                "Prefer erf.cfl with erf.fixed_dt removed.")
