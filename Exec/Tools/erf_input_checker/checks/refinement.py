"""Mesh refinement: indicators, regridding, coarse-fine coupling."""

from ..checker import Checker, register


@register
class RefinementChecker(Checker):
    name = "refinement"

    def run(self) -> None:
        c, d = self.ctx, self.deck
        if c.max_level is None:
            return
        if c.max_level == 0:
            for key in ("amr.regrid_int", "amr.ref_ratio", "amr.ref_ratio_vect", "erf.refinement_indicators",
                        "erf.coupling_type", "erf.cf_width", "erf.dt_ref_ratio"):
                if d.has(key):
                    self.info(key, "ignored with amr.max_level = 0")
            return
        if d.has("amr.ref_ratio") and d.has("amr.ref_ratio_vect"):
            self.error("amr.ref_ratio_vect", "amr.ref_ratio and amr.ref_ratio_vect may not both be given")
        for lev, rr in enumerate(c.ref_ratio):
            for idim, r in enumerate(rr):
                if r < 1:
                    self.error("amr.ref_ratio", "refinement ratio %d at level %d is < 1" % (r, lev))
                if r > 4:
                    self.advise("amr.ref_ratio", "refinement ratio %d between levels %d and %d is large" % (r, lev, lev + 1),
                                "ERF is normally run with ratios of 2 or 3 per level.")
            if all(r == 1 for r in rr):
                self.warn("amr.ref_ratio", "ratio 1 1 1 between levels %d and %d refines nothing" % (lev, lev + 1))
        indicators = d.get_strs("erf.refinement_indicators")
        has_file = d.has("amr.regrid_file")
        has_nc = any(d.has("erf.nc_init_file_%d" % lev) for lev in range(1, c.nlev))
        if not indicators and not has_file and not has_nc and not d.has("amr.refine_whole_domain_dir"):
            self.warn("erf.refinement_indicators", "amr.max_level = %d but no refinement indicator, regrid file or fine-level input file is given: no fine level is ever created" % c.max_level,
                      "Add erf.refinement_indicators = <name> with erf.<name>.in_box_lo/hi (or a field criterion), or set amr.max_level = 0.")
        for ind in indicators:
            keys = d.keys_with_prefix("erf.%s." % ind)
            if not keys:
                self.error("erf.refinement_indicators", "indicator %r names no erf.%s.* keys" % (ind, ind),
                           "Give erf.%s.in_box_lo / in_box_hi, or field_name with value_greater / value_less / adjacent_difference_greater." % ind)
                continue
            self._check_box(ind)
        regrid = self.int_("amr.regrid_int", -1)
        if regrid is not None and regrid > 0 and indicators and all(d.has("erf.%s.in_box_lo" % i) for i in indicators) \
                and not any(d.has("erf.%s.move_speed_x" % i) or d.has("erf.%s.move_speed_y" % i) for i in indicators):
            self.advise("amr.regrid_int", "regridding every %d steps with only static boxes rebuilds the same grids each time" % regrid,
                        "Set amr.regrid_int = -1 (or leave it unset) for static refinement boxes.")
        coupling = (self.str_("erf.coupling_type", "TwoWay") or "").lower()
        cf = self.int_("erf.cf_width", 0)
        if coupling == "twoway" and cf and cf > 0:
            self.error("erf.cf_width", "two-way coupling needs erf.cf_width = 0; ERF stops")
        if self.int_("erf.cf_set_width", 0):
            self.error("erf.cf_set_width", "must be 0; ERF stops")
        if c.is_anelastic and coupling == "twoway" and c.max_level > 0:
            self.info("erf.coupling_type", "anelastic levels are one-way coupled whatever erf.coupling_type says")

    def _check_box(self, ind: str) -> None:
        c = self.ctx
        lo = self.reals("erf.%s.in_box_lo" % ind)
        hi = self.reals("erf.%s.in_box_hi" % ind)
        if not lo and not hi:
            return
        if len(lo) != len(hi) or len(lo) not in (2, 3):
            self.error("erf.%s.in_box_lo" % ind, "in_box_lo and in_box_hi need the same 2 or 3 values")
            return
        if any(lo[i] >= hi[i] for i in range(len(lo))):
            self.error("erf.%s.in_box_lo" % ind, "in_box_lo must be below in_box_hi in every direction")
        if c.prob_lo and c.prob_hi:
            for i in range(len(lo)):
                if lo[i] < c.prob_lo[i] or hi[i] > c.prob_hi[i]:
                    self.warn("erf.%s.in_box_lo" % ind, "refinement box leaves the domain in direction %d (%g..%g outside %g..%g)" % (
                        i, lo[i], hi[i], c.prob_lo[i], c.prob_hi[i]), "The box is clipped to the domain; say so or shrink it.")
            if c.dx and len(lo) == 3 and c.n_cell:
                # a box narrower than a few coarse cells cannot be refined usefully
                for i in range(3):
                    if (hi[i] - lo[i]) < 2 * c.dx[i] and c.n_cell[i] > 1:
                        self.warn("erf.%s.in_box_hi" % ind, "refinement box is thinner than two coarse cells in direction %d" % i)
