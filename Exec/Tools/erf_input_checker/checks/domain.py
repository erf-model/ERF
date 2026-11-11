"""Domain size, cell counts and cell sizes."""

from ..checker import Checker, register

REQUIRED = ["amr.max_level"]


@register
class DomainChecker(Checker):
    name = "domain"

    def run(self) -> None:
        d, c = self.deck, self.ctx
        if d.has("amr.n_cell") and d.has("erf.n_cell"):
            self.error("erf.n_cell", "both amr.n_cell and erf.n_cell are set; ERF stops",
                       "Keep one of them.")
        if c.n_cell is None:
            self.error("amr.n_cell", "the cell count is missing or not three integers",
                       "Set amr.n_cell = nx ny nz.")
        if c.prob_hi is None:
            self.error("geometry.prob_hi", "the domain size is missing: give geometry.prob_hi or geometry.prob_extent",
                       "Set geometry.prob_hi = x y z (with geometry.prob_lo) or geometry.prob_extent = Lx Ly Lz.")
        if d.has("geometry.prob_hi") and d.has("geometry.prob_extent"):
            self.error("geometry.prob_extent", "both prob_hi and prob_extent are given; AMReX stops (AMReX_Geometry.cpp asserts exactly one)",
                       "Keep one of them.")
        if c.max_level is None:
            self.error("amr.max_level", "amr.max_level must be set (0 for a single level)")
        elif c.max_level < 0:
            self.error("amr.max_level", "must be >= 0")

        if c.n_cell is not None:
            for i, n in enumerate(c.n_cell):
                if n <= 0:
                    self.error("amr.n_cell", "cell count in direction %d is %d" % (i, n))
        if c.extent is not None:
            for i, L in enumerate(c.extent):
                if L <= 0:
                    self.error("geometry.prob_hi", "domain length in direction %d is %g (prob_hi must exceed prob_lo)" % (i, L))

        dx = c.dx
        if dx is None:
            return
        self.info("amr.n_cell", "level-0 cell size dx, dy, dz = %.4g, %.4g, %.4g m%s" % (
            dx[0], dx[1], dx[2], " (dz nominal; first cell %.4g m)" % c.dz_min
            if c.dz_min and abs(c.dz_min - dx[2]) > 1e-9 * dx[2] else ""))

        nx, ny, nz = c.n_cell
        if nx > 1 and ny > 1 and abs(dx[0] - dx[1]) > 1e-6 * dx[0]:
            self.advise("geometry.prob_hi", "dx = %.4g and dy = %.4g differ" % (dx[0], dx[1]),
                        "ERF accepts this, but most advection and LES closures assume dx = dy.")

        # vertical-to-horizontal aspect ratio and the LES closure
        dz0 = c.dz_min or dx[2]
        les = self.str_("erf.les_type", "None").lower()
        if nx > 1 and ny > 1 and les in ("smagorinsky", "deardorff"):
            ar = min(dx[0], dx[1]) / dz0
            if ar > 5:
                self.advise("erf.les_type", "cells are %.0f times wider than they are tall at the surface" % ar,
                            "An isotropic LES length scale (Smagorinsky / Deardorff) is a poor fit for flat cells; "
                            "consider erf.les_type = Smagorinsky2D, erf.mix_isotropic = false, or a PBL scheme.")

        # 2-D decks must say so through the cell count and the boundary types
        if (nx == 1) != (ny == 1) and not c.periodic[1 if ny == 1 else 0]:
            face = "y" if ny == 1 else "x"
            lo, hi = self.str_("%slo.type" % face, ""), self.str_("%shi.type" % face, "")
            if "slipwall" not in (lo or "").lower() or "slipwall" not in (hi or "").lower():
                self.advise("amr.n_cell", "a single cell in %s with non-periodic %s faces" % (face, face),
                            "For a 2-D case set %slo.type = %shi.type = SlipWall or make the direction periodic." % (face, face))
