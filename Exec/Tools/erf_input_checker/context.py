"""Facts derived from the deck that several checkers need.

:class:`RunContext` is built once per deck.  It resolves the ERF conventions
(``amr.n_cell`` or ``erf.n_cell``, ``geometry.prob_hi`` or ``prob_extent``,
per-level lists that may be a single value) and classifies the run:

* ``is_anelastic``        -- ``erf.anelastic`` = 1 at level 0
* ``substepping``         -- "Implicit" (the compressible default) or "None"
* ``init_mode``           -- "real" for WRFInput / Metgrid / NCFile / HindCast,
                             "ideal" otherwise
* ``terrain``             -- the ``erf.terrain_type`` value, normalised
"""

import math
from typing import List, Optional

from .deck import Deck, DeckError

REAL_INIT_TYPES = {"wrfinput", "metgrid", "ncfile", "hindcast"}
FITTED_TERRAIN = {"staticfittedmesh", "movingfittedmesh", "static", "moving"}

# AMReX / ERF defaults (Source/main.cpp add_par and AMReX_AmrMesh)
DEFAULT_MAX_GRID_SIZE = 2048
DEFAULT_BLOCKING_FACTOR = 1
DEFAULT_NO_BOX_SPLIT_DIR = 2
DEFAULT_CFL = 0.8
DEFAULT_SUBSTEPPING_CFL = 1.0
DEFAULT_CHANGE_MAX = 1.1
DEFAULT_INIT_SHRINK = 1.0
DEFAULT_DT_MAX_INITIAL = 1.0


class Options:
    """Command-line knobs that change the advice."""

    def __init__(self, ranks: Optional[int] = None, gpu: bool = False,
                 umax: float = 20.0, theta: float = 300.0):
        self.ranks = ranks
        self.gpu = gpu
        self.umax = umax      # largest wind speed [m/s] expected, for CFL estimates
        self.theta = theta    # potential temperature [K] used for the sound speed


class RunContext:
    def __init__(self, deck: Deck, options: Optional[Options] = None):
        self.deck = deck
        self.options = options or Options()
        self.problems: List[str] = []   # parse problems met while deriving facts

        self.n_cell = self._n_cell()
        self.prob_lo, self.prob_hi = self._extent()
        self.max_level = self._safe(lambda: deck.get_int("amr.max_level", None), None)
        self.nlev = (self.max_level + 1) if self.max_level is not None else 1
        self.ref_ratio = self._ref_ratio()

        self.anelastic = self._per_level_int("erf.anelastic", 0)
        self.is_anelastic = bool(self.anelastic and self.anelastic[0] == 1)
        self.substepping = self._substepping()
        self.init_type = (deck.get_str("erf.init_type", "None") or "None")
        self.init_mode = "real" if self.init_type.lower() in REAL_INIT_TYPES else "ideal"
        self.terrain = (deck.get_str("erf.terrain_type", "None") or "None")
        self.fitted_terrain = self.terrain.lower() in FITTED_TERRAIN
        self.periodic = self._periodic()
        # ERF defaults use_real_bcs to true for WRFInput / Metgrid when a wrfbdy file is
        # named, and lets the deck turn it off; without the file it is false.
        it = self.init_type.lower()
        default_real = it == "metgrid" or (it == "wrfinput" and deck.has("erf.nc_bdy_file"))
        self.use_real_bcs = self._safe(lambda: deck.get_bool("erf.use_real_bcs", default_real), default_real)
        if it == "wrfinput" and not deck.has("erf.nc_bdy_file"):
            self.use_real_bcs = False
        self.fixed_dt = self._safe(lambda: deck.get_real("erf.fixed_dt"), None)
        if self.fixed_dt is not None and self.fixed_dt <= 0:
            self.fixed_dt = None   # ERF's own "not set" sentinel

    # -- helpers ----------------------------------------------------------
    def _safe(self, fn, default):
        try:
            return fn()
        except DeckError as exc:
            self.problems.append(str(exc))
            return default

    def _n_cell(self) -> Optional[List[int]]:
        d = self.deck
        key = "amr.n_cell" if d.has("amr.n_cell") else ("erf.n_cell" if d.has("erf.n_cell") else None)
        if key is None:
            return None
        v = self._safe(lambda: d.get_ints(key), [])
        return v if len(v) == 3 else None

    def _extent(self):
        d = self.deck
        lo = self._safe(lambda: d.get_reals("geometry.prob_lo"), []) or [0.0, 0.0, 0.0]
        hi = self._safe(lambda: d.get_reals("geometry.prob_hi"), [])
        if len(hi) != 3:
            ext = self._safe(lambda: d.get_reals("geometry.prob_extent"), [])
            hi = [lo[i] + ext[i] for i in range(3)] if len(ext) == 3 and len(lo) == 3 else None
        if len(lo) != 3:
            lo = None
        return lo, hi

    def _ref_ratio(self) -> List[List[int]]:
        d = self.deck
        nfine = max(self.nlev - 1, 0)
        if d.has("amr.ref_ratio_vect"):
            v = self._safe(lambda: d.get_ints("amr.ref_ratio_vect"), [])
            out = [v[3 * i:3 * i + 3] for i in range(len(v) // 3)]
        else:
            v = self._safe(lambda: d.get_ints("amr.ref_ratio"), []) or [2]
            out = [[r, r, r] for r in v]
        while len(out) < nfine:
            out.append(out[-1] if out else [2, 2, 2])
        return out[:nfine]

    def _per_level_int(self, key: str, default: int) -> List[int]:
        v = self._safe(lambda: self.deck.get_ints(key), [])
        if not v:
            return [default] * self.nlev
        if len(v) == 1:
            return v * self.nlev
        return (v + [v[-1]] * self.nlev)[:self.nlev]

    def per_level_str(self, key: str, default: str) -> List[str]:
        v = self.deck.get_strs(key)
        if not v:
            return [default] * self.nlev
        if len(v) == 1:
            return v * self.nlev
        return (v + [v[-1]] * self.nlev)[:self.nlev]

    def _substepping(self) -> str:
        v = self.per_level_str("erf.substepping_type", "Implicit")[0]
        if self.is_anelastic:
            return "None"
        return "Implicit" if v.lower() == "implicit" else "None"

    def _periodic(self) -> List[bool]:
        v = self._safe(lambda: self.deck.get_ints("geometry.is_periodic"), [])
        if len(v) != 3:
            return [False, False, False]
        return [x != 0 for x in v]

    # -- derived geometry -------------------------------------------------
    @property
    def extent(self) -> Optional[List[float]]:
        if self.prob_lo is None or self.prob_hi is None:
            return None
        return [self.prob_hi[i] - self.prob_lo[i] for i in range(3)]

    @property
    def dx(self) -> Optional[List[float]]:
        """Level-0 cell size (nominal dz for a stretched grid)."""
        if self.extent is None or self.n_cell is None:
            return None
        if any(n <= 0 for n in self.n_cell):
            return None
        return [self.extent[i] / self.n_cell[i] for i in range(3)]

    @property
    def dz_min(self) -> Optional[float]:
        """Smallest vertical cell: the first cell of a stretched grid, else dz."""
        d = self.deck
        dx = self.dx
        if d.has("erf.terrain_z_levels"):
            z = self._safe(lambda: d.get_reals("erf.terrain_z_levels"), [])
            if len(z) > 1:
                return min(z[i + 1] - z[i] for i in range(len(z) - 1))
        if d.has("erf.initial_dz"):
            return self._safe(lambda: d.get_real("erf.initial_dz"), None)
        return dx[2] if dx else None

    def sound_speed(self) -> float:
        """Dry sound speed at the surface for the option temperature."""
        gamma, r_d = 1.4, 287.0
        return math.sqrt(gamma * r_d * self.options.theta)

    def horizontal_2d(self) -> bool:
        return self.n_cell is not None and (self.n_cell[0] == 1 or self.n_cell[1] == 1)
