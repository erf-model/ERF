"""Grid layout: blocking factors, max_grid_size, the number of boxes per rank."""

import math
from typing import List

from ..checker import Checker, register
from ..context import DEFAULT_BLOCKING_FACTOR, DEFAULT_MAX_GRID_SIZE, DEFAULT_NO_BOX_SPLIT_DIR

COLUMN_PBL = {"mynn25", "mynnedmf", "myj", "ysu", "ysunew", "mrf", "native_shoc", "eamxx_shoc", "shoc"}


def _per_level_per_dim(ck: Checker, base: str, default: int, nlev: int) -> List[List[int]]:
    """AMReX reading of amr.<base> and amr.<base>_x/_y/_z: one value or one per level."""
    vals = ck.ints("amr." + base) or [default]
    out = [[(vals[i] if i < len(vals) else vals[-1])] * 3 for i in range(nlev)]
    for idim, suf in enumerate("xyz"):
        v = ck.ints("amr.%s_%s" % (base, suf))
        if v:
            for i in range(nlev):
                out[i][idim] = v[i] if i < len(v) else v[-1]
    return out


def _even_divisors(n: int, lo: int, hi: int) -> List[int]:
    return [m for m in range(lo, hi + 1) if m % 2 == 0 and n % m == 0]


def _is_pow2(k: int) -> bool:
    return k > 0 and (k & (k - 1)) == 0


@register
class GridLayoutChecker(Checker):
    name = "grid"

    def run(self) -> None:
        c = self.ctx
        if c.n_cell is None or c.max_level is None:
            return
        nlev = c.nlev
        bf = _per_level_per_dim(self, "blocking_factor", DEFAULT_BLOCKING_FACTOR, nlev)
        mgs = _per_level_per_dim(self, "max_grid_size", DEFAULT_MAX_GRID_SIZE, nlev)
        nbsd = self.int_("amr.no_box_split_dir", DEFAULT_NO_BOX_SPLIT_DIR)
        self._check_amrex_rules(bf, mgs)
        self._check_odd_and_even(mgs)
        self._check_vertical(mgs, bf, nbsd)
        self._check_boxes(mgs, nbsd)

    # -- what Amr::checkInputExtended enforces (ERF sets amr.no_box_split_dir >= 0, so
    #    AMReX runs the extended check; the legacy check applies only with -1) -------------
    def _check_amrex_rules(self, bf, mgs) -> None:
        c = self.ctx
        nbsd = self.int_("amr.no_box_split_dir", DEFAULT_NO_BOX_SPLIT_DIR)
        legacy = nbsd is None or nbsd < 0
        for lev in range(c.nlev):
            for idim in range(3):
                if not _is_pow2(bf[lev][idim]):
                    self.error("amr.blocking_factor", "blocking factor %d at level %d is not a power of 2; AMReX stops" % (
                        bf[lev][idim], lev), "Use 1, 2, 4, 8, 16, ...")
                skip = (not legacy and idim == nbsd) or bf[lev][idim] > mgs[lev][idim]
                if c.max_level == 0 or skip:
                    continue
                if mgs[lev][idim] % bf[lev][idim] != 0:
                    self.error("amr.max_grid_size", "max_grid_size %d in %s is not a multiple of blocking_factor %d at level %d; AMReX stops" % (
                        mgs[lev][idim], "xyz"[idim], bf[lev][idim], lev))
        for idim in range(3):
            skip = (not legacy and idim == nbsd) or (legacy and bf[0][idim] > mgs[0][idim])
            if not skip and c.n_cell[idim] % bf[0][idim] != 0:
                self.error("amr.n_cell", "n_cell %d in %s is not a multiple of blocking_factor %d; AMReX stops" % (
                    c.n_cell[idim], "xyz"[idim], bf[0][idim]),
                    "Change n_cell to a multiple of the blocking factor, or lower the blocking factor.")
        for lev in range(c.nlev - 1):
            rr = c.ref_ratio[lev]
            if max(rr) < 2:
                self.error("amr.ref_ratio", "refinement ratio %s between levels %d and %d; AMReX stops ('bad ref_ratios')" % (rr, lev, lev + 1))
            if not legacy:
                continue
            for idim in range(3):
                if bf[lev][idim] * rr[idim] < bf[lev + 1][idim]:
                    self.error("amr.blocking_factor", "blocking factors %d (level %d) and %d (level %d) differ by more than the refinement ratio %d; AMReX stops" % (
                        bf[lev][idim], lev, bf[lev + 1][idim], lev + 1, rr[idim]))

    # -- odd numbers --------------------------------------------------------
    def _check_odd_and_even(self, mgs) -> None:
        c = self.ctx
        user_set = any(self.deck.has("amr.max_grid_size" + s) for s in ("", "_x", "_y", "_z"))
        nbsd = self.int_("amr.no_box_split_dir", DEFAULT_NO_BOX_SPLIT_DIR)
        if user_set:
            for lev in range(c.nlev):
                for idim, suf in enumerate("xyz"):
                    m = mgs[lev][idim]
                    if m == DEFAULT_MAX_GRID_SIZE:
                        continue
                    odd_rr = lev > 0 and c.ref_ratio[lev - 1][idim] > 1 and c.ref_ratio[lev - 1][idim] % 2 == 1
                    if m % 2 == 1 and idim != nbsd and not odd_rr:
                        self.advise("amr.max_grid_size", "max_grid_size %d in %s at level %d is odd: the Inputs docs say it must be even, boxes of odd length cannot be halved for load balancing or coarsened by 2 at a coarse-fine interface" % (
                            m, suf, lev), "Use an even value, ideally a multiple of 8 (CPU) or 32 (GPU).")
                    elif m % 2 == 1 and idim == nbsd:
                        self.info("amr.max_grid_size", "odd max_grid_size %d in %s is harmless: boxes are never split in that direction" % (m, suf))
                    elif m < 8 and c.n_cell[idim] > 8:
                        self.advise("amr.max_grid_size", "max_grid_size %d in %s is very small; the ghost-cell overhead per box becomes large" % (m, suf))
        for idim, suf in enumerate("xyz"):
            n = c.n_cell[idim]
            if n > 1 and n % 2 == 1:
                hint = "cannot be split into equal halves" + (
                    " and no fine level can cover it evenly" if c.max_level > 0 else "")
                self.advise("amr.n_cell", "n_cell = %d in %s is odd and %s" % (n, suf, hint),
                            "Prefer a multiple of 8 so boxes divide evenly across ranks (and of 32 on GPUs).")
            elif n > 8 and idim < 2 and n % 8 != 0:
                self.info("amr.n_cell", "n_cell = %d in %s is not a multiple of 8; multiples of 8 (CPU) or 32 (GPU) divide into equal boxes more easily" % (n, suf))

    # -- vertical decomposition and column physics -------------------------------
    def _check_vertical(self, mgs, bf, nbsd) -> None:
        c, d = self.ctx, self.deck
        nz = c.n_cell[2]
        single_value = d.has("amr.max_grid_size") and not d.has("amr.max_grid_size_z")
        z_limited = any(mgs[lev][2] < nz * (c.ref_ratio[lev - 1][2] if lev > 0 else 1) for lev in range(c.nlev))
        needs_columns = self._needs_full_columns()
        if nbsd == 2:
            if d.has("amr.max_grid_size_z") and z_limited:
                self.info("amr.max_grid_size_z", "ignored: amr.no_box_split_dir = 2 (the ERF default) never splits boxes in z")
            if self.int_("amr.refine_grid_layout_z", 0) == 1:
                self.info("amr.refine_grid_layout_z", "ignored while amr.no_box_split_dir = 2")
            return
        if nbsd not in (-1, 0, 1, 2):
            self.error("amr.no_box_split_dir", "must be -1, 0, 1 or 2")
        if c.substepping == "Implicit" and nbsd != 2:
            self.error("amr.no_box_split_dir", "= %d with implicit acoustic substepping; ERF stops (the column solve needs unsplit columns)" % nbsd,
                       "Remove amr.no_box_split_dir, or set erf.substepping_type = None.")
        if z_limited and needs_columns:
            self.error("amr.max_grid_size", "boxes are split in z (max_grid_size_z < n_cell z) while %s needs whole columns; ERF stops" % needs_columns,
                       "Set amr.max_grid_size_x / _y instead of a single amr.max_grid_size, and leave z unlimited.")
        elif z_limited and single_value:
            self.warn("amr.max_grid_size", "a single max_grid_size with amr.no_box_split_dir = %d also chops the domain in z" % nbsd,
                      "Use amr.max_grid_size_x and amr.max_grid_size_y to limit only the horizontal box size.")
        if self.int_("amr.refine_grid_layout_z", 0) == 1 and needs_columns:
            self.warn("amr.refine_grid_layout_z", "= 1 lets load balancing chop boxes in z, which %s does not allow" % needs_columns)

    def _needs_full_columns(self) -> str:
        c, d = self.ctx, self.deck
        pbl = [p.lower() for p in c.per_level_str("erf.pbl_type", "None")]
        if any(p in COLUMN_PBL for p in pbl):
            return "the PBL scheme"
        if (d.get_str("zlo.type", "") or "").lower() == "surface_layer":
            return "the surface layer at zlo"
        if c.substepping == "Implicit":
            return "implicit acoustic substepping"
        vi = self.bool_("erf.vert_implicit", not c.is_anelastic)
        fac = self.reals("erf.vert_implicit_fac")
        if vi and (not fac or any(f > 0 for f in fac)):
            return "the implicit vertical diffusion"
        return ""

    # -- boxes and ranks ----------------------------------------------------------
    def _check_boxes(self, mgs, nbsd) -> None:
        c, d = self.ctx, self.deck
        ranks = c.options.ranks
        n = c.n_cell
        user_set = any(d.has("amr.max_grid_size" + s) for s in ("", "_x", "_y", "_z"))
        split_z = (nbsd != 2) and mgs[0][2] < n[2]
        if not user_set:
            if ranks:
                self.info("amr.max_grid_size", "not set: ERF splits level 0 into about %d boxes in x and y, one per rank" % ranks,
                          "" if ranks <= n[0] * n[1] else "There are more ranks than horizontal cells.")
            else:
                self.info("amr.max_grid_size", "not set: ERF makes one level-0 box per rank (never split in z); pass --ranks N to check the balance")
            return
        dims = [0, 1] + ([2] if split_z else [])
        counts = []
        uneven = []
        for idim in dims:
            m = mgs[0][idim]
            k = int(math.ceil(n[idim] / float(m)))
            counts.append(k)
            if n[idim] % m != 0 and k > 1:
                uneven.append((idim, m, n[idim] % m))
        nbox = 1
        for k in counts:
            nbox *= k
        self.info("amr.max_grid_size", "level 0 has %d box(es) of up to %s cells%s" % (
            nbox, "x".join(str(min(mgs[0][i], n[i])) for i in range(3)),
            "" if nbsd != 2 else " (never split in z)"))
        for idim, m, rem in uneven:
            suf = "xyz"[idim]
            opts = _even_divisors(n[idim], max(8, m // 2), min(n[idim], 2 * m))
            adv = "Even boxes come from a max_grid_size that divides n_cell: " + (
                ", ".join(str(o) for o in opts[:6]) if opts else "none between %d and %d; change n_cell" % (m // 2, 2 * m)) + "."
            self.advise("amr.max_grid_size", "n_cell %d in %s is not a multiple of max_grid_size %d: the last box has only %d cells, so the work per box is uneven" % (
                n[idim], suf, m, rem), adv)
        if ranks:
            if nbox < ranks:
                self.warn("amr.max_grid_size", "%d level-0 boxes for %d ranks; AMReX will chop boxes further (refine_grid_layout) but they are not what the deck says" % (nbox, ranks),
                          "Lower max_grid_size so that the box count is a multiple of the rank count, or leave it unset.")
            elif nbox % ranks != 0:
                per = int(math.ceil(nbox / float(ranks)))
                self.advise("amr.max_grid_size", "%d boxes on %d ranks: %d ranks hold %d boxes and the rest fewer (%.0f%% of the time idle)" % (
                    nbox, ranks, nbox - (per - 1) * ranks, per, 100.0 * (1 - nbox / float(per * ranks))),
                    "Choose max_grid_size so the box count is a multiple of the rank count.")
        gpu = c.options.gpu
        cells = 1
        for i in range(3):
            cells *= min(mgs[0][i], n[i])
        if gpu and cells < 32 * 32 * 32 and n[0] * n[1] * n[2] > 4 * cells:
            self.advise("amr.max_grid_size", "boxes of %d cells are small for a GPU; kernel launch overhead dominates" % cells,
                        "On GPUs prefer boxes of 64-128 cells per side (one or a few boxes per GPU).")
        if not gpu and cells > 4 * 1024 ** 2 and ranks and nbox < ranks:
            self.advise("amr.max_grid_size", "boxes of %d cells are large for CPU ranks" % cells,
                        "On CPUs 32-64 cells per side per box usually balance cache use and ghost-cell overhead.")
        if d.has("fabarray.mfiter_tile_size"):
            ts = self.ints("fabarray.mfiter_tile_size")
            if gpu and ts and min(ts) < 1024:
                self.info("fabarray.mfiter_tile_size", "tiling is ignored on GPUs")
