#!/usr/bin/env python3
"""Check or propose the refined box of a building deck for the face balance.

    python3 ibseb_refinement_box.py inputs                              # check every in_box of the deck
    python3 ibseb_refinement_box.py inputs --all --fit tight            # the smallest box around every building
    python3 ibseb_refinement_box.py inputs --all --fit relaxed          # the same, padded by 3 coarse cells
    python3 ibseb_refinement_box.py inputs --region 200 440 200 440 --fit relaxed --margin 5 --name city

With erf.ibseb.enable on more than one level, every building must lie wholly
inside a refined level, with at least one cell of the level below around it,
or wholly outside it: each level builds its faces from its own cells, and a
building the level covers in part would get a partial face list. ERF checks
that at start-up and stops (ERF::ibseb_check_refined_levels()); this script
says before the run where the edges of a static refined box may go.

It reads the deck (following FILE includes) for geometry.prob_lo/prob_extent,
amr.n_cell, amr.max_level, amr.ref_ratio(_vect), amr.blocking_factor,
amr.n_error_buf, erf.buildings_file_name and erf.refinement_indicators with
their in_box_lo / in_box_hi, and the height map itself (ERF's terrain text
format: nx, ny, the x nodes, the y nodes, the heights z[ix*ny + iy]).

How a building looks on the coarse level: the map is interpolated to the
nodes of the finest level, and the embedded boundary ramps from a roof to
the ground over a cell, so on the coarse level a building can fill the cell
beyond its footprint (a 20 m cell next to a 10 m ramp, for example). The
script takes a coarse column as built when any node of the finest level on
its closed footprint stands above the ground, which is the most the
building can occupy. That is conservative: a box this script accepts passes
the start-up check; a box it rejects may still pass by a cell.

A proposed box needs --fit. ``tight`` is the smallest box the start-up check
accepts: every building it holds with one coarse cell around its coarse
extent. ``relaxed`` pads that box by --margin coarse cells (3 by default) on
every side before the check, which keeps the edge of the refined level, where
the coarse level fills the fine one, away from the buildings and the flow
around them; the padded box is grown again if the padding brings a building
onto its edge. In check mode (no --all or --region) --fit is not used.

Only refinement in x and y by a static box is handled (level 1 boxes; the
deck's other levels are not checked), spanning the whole depth (in_box with
two values, or three that cover it). A proposed box has its edges on
multiples of amr.blocking_factor fine cells, so the grids ERF builds are
exactly the box; a box given off those blocks is judged as the grids AMReX
grows it to. The deck must carry amr.n_error_buf = 0 (ERF stops on an
explicit box otherwise): checking a deck's boxes stops without it, and a
proposal says so and prints the line.
"""
import argparse
import bisect
import os
import sys


# ---------------------------------------------------------------------------
# Deck reading
# ---------------------------------------------------------------------------

def read_deck(path, keys=None, seen=None):
    """Return {key: [tokens]} from a ParmParse deck, following FILE = includes."""
    keys = {} if keys is None else keys
    seen = set() if seen is None else seen
    path = os.path.abspath(path)
    if path in seen:
        return keys
    seen.add(path)
    with open(path) as f:
        for raw in f:
            line = raw.split("#", 1)[0].strip()
            if "=" not in line:
                continue
            k, v = line.split("=", 1)
            k, toks = k.strip(), v.split()
            toks = [t.strip('"') for t in toks]
            if k == "FILE":
                read_deck(os.path.join(os.path.dirname(path), toks[0]), keys, seen)
            else:
                keys[k] = toks       # a later line overrides, as an include may be overridden
    return keys


def get(keys, name, n=None, typ=float, default=None):
    if name not in keys:
        if default is None:
            raise SystemExit(f"the deck does not set {name}")
        return default
    v = [typ(t) for t in keys[name]]
    if n is not None and len(v) < n:
        raise SystemExit(f"{name} needs {n} values, the deck gives {len(v)}")
    return v


def read_map(path):
    """The height map as node lists x, y and heights z[ix][iy]."""
    with open(path) as f:
        vals = [float(t) for t in f.read().split()]
    nx, ny = int(vals[0]), int(vals[1])
    x = vals[2:2 + nx]
    y = vals[2 + nx:2 + nx + ny]
    flat = vals[2 + nx + ny:2 + nx + ny + nx * ny]
    if len(flat) != nx * ny:
        raise SystemExit(f"{path}: expected {nx * ny} heights, found {len(flat)}")
    for name, nodes in (("x", x), ("y", y)):
        if any(b <= a for a, b in zip(nodes, nodes[1:])):
            raise SystemExit(f"{path}: the {name} nodes must increase strictly")
    z = [flat[i * ny:(i + 1) * ny] for i in range(nx)]
    return x, y, z


def brackets(nodes, targets):
    """For each target, the node index k and weight w with value (1 - w) f[k] + w f[k + 1],
    held at the end values outside the nodes (as numpy.interp)."""
    if len(nodes) == 1:
        return [(0, 0.0)] * len(targets)
    out = []
    for v in targets:
        if v <= nodes[0]:
            out.append((0, 0.0))
        elif v >= nodes[-1]:
            out.append((len(nodes) - 2, 1.0))
        else:
            k = bisect.bisect_right(nodes, v) - 1
            out.append((k, (v - nodes[k]) / (nodes[k + 1] - nodes[k])))
    return out


def interp_map(x, y, z, xs, ys):
    """Bilinear interpolation of the map at the nodes xs (x) and ys (y): along y for
    every map row first, then along x, each bracket found once."""
    by, bx = brackets(y, ys), brackets(x, xs)
    single_y = len(y) == 1
    zy = [[row[k] if single_y else (1.0 - w) * row[k] + w * row[k + 1] for k, w in by] for row in z]
    if len(x) == 1:
        return [list(zy[0]) for _ in xs]
    return [[(1.0 - w) * p + w * q for p, q in zip(zy[k], zy[k + 1])] for k, w in bx]


# ---------------------------------------------------------------------------
# Buildings on the coarse level
# ---------------------------------------------------------------------------

def built_columns(h_fine, ratio):
    """Coarse columns that any finest-level node on their closed footprint stands above,
    and the tallest such node; h_fine are the heights at the finest nodes,
    (nxc*rx + 1) x (nyc*ry + 1)."""
    rx, ry = ratio
    nxc = (len(h_fine) - 1) // rx
    nyc = (len(h_fine[0]) - 1) // ry
    built = [[False] * nyc for _ in range(nxc)]
    top = [[0.0] * nyc for _ in range(nxc)]
    for i in range(nxc):
        for j in range(nyc):
            hmax = max(h_fine[a][b] for a in range(i * rx, (i + 1) * rx + 1) for b in range(j * ry, (j + 1) * ry + 1))
            if hmax > 1e-6:
                built[i][j] = True
                top[i][j] = hmax
    return built, top


def label(built, per_x, per_y):
    """4-connected labels of the built columns (periodic wrap where asked), 0 = open."""
    nx, ny = len(built), len(built[0])
    lab = [[0] * ny for _ in range(nx)]
    n = 0
    for i in range(nx):
        for j in range(ny):
            if not built[i][j] or lab[i][j]:
                continue
            n += 1
            lab[i][j] = n
            stack = [(i, j)]
            while stack:
                a, b = stack.pop()
                for da, db in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    p, q = a + da, b + db
                    if per_x:
                        p %= nx
                    if per_y:
                        q %= ny
                    if 0 <= p < nx and 0 <= q < ny and built[p][q] and not lab[p][q]:
                        lab[p][q] = n
                        stack.append((p, q))
    return lab, n


_COLUMNS = {}


def columns_of(lab, b):
    """Columns of building b, from one pass over the labels per label map (cached):
    the scans below ask for every building at every growth step."""
    key = id(lab)
    if key not in _COLUMNS:
        by_label = {}
        for i, row in enumerate(lab):
            for j, v in enumerate(row):
                if v:
                    by_label.setdefault(v, []).append((i, j))
        _COLUMNS[key] = by_label
    return _COLUMNS[key].get(b, [])


def margin_cells(lab, b, per_x, per_y):
    """Coarse columns of building b grown by one column (the 3 x 3 rule of ERF's check)."""
    nx, ny = len(lab), len(lab[0])
    cells = set()
    for i, j in columns_of(lab, b):
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                p, q = i + di, j + dj
                if per_x:
                    p %= nx
                elif not 0 <= p < nx:
                    continue
                if per_y:
                    q %= ny
                elif not 0 <= q < ny:
                    continue
                cells.add((p, q))
    return cells


def covered(boxes, i, j):
    """Whether coarse column (i, j) is in any of boxes, each (ilo, ihi, jlo, jhi) inclusive."""
    return any(b[0] <= i <= b[1] and b[2] <= j <= b[3] for b in boxes)


def classify(lab, nb, boxes, per_x, per_y):
    """inside / outside / crossing per building for the union of coarse index boxes (the refined level)."""
    out = {}
    for b in range(1, nb + 1):
        cov = [covered(boxes, i, j) for i, j in margin_cells(lab, b, per_x, per_y)]
        out[b] = "inside" if all(cov) else ("outside" if not any(cov) else "crossing")
    return out


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("deck", help="the inputs file")
    ap.add_argument("--buildings", help="height map (default: erf.buildings_file_name of the deck)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--all", action="store_true", help="propose a box around every building")
    g.add_argument("--region", type=float, nargs=4, metavar=("XLO", "XHI", "YLO", "YHI"),
                   help="propose the smallest box holding this region [m] that no building crosses")
    ap.add_argument("--fit", choices=("tight", "relaxed"),
                    help="with --all or --region (required there): the smallest box the check accepts, "
                         "or one padded by --margin coarse cells")
    ap.add_argument("--margin", type=int, default=3,
                    help="padding of a relaxed box in coarse cells on every side (default 3)")
    ap.add_argument("--name", default="city", help="refinement indicator name for the proposed lines (default city)")
    a = ap.parse_args()
    if (a.all or a.region) and a.fit is None:
        ap.error("a proposed box needs --fit tight (the smallest box the check accepts) "
                 "or --fit relaxed (padded by --margin coarse cells)")
    if a.margin < 0:
        ap.error("--margin must not be negative")

    keys = read_deck(a.deck)
    deck_dir = os.path.dirname(os.path.abspath(a.deck))
    plo = get(keys, "geometry.prob_lo", 3, default=[0.0, 0.0, 0.0])
    if "geometry.prob_extent" in keys:
        ext = get(keys, "geometry.prob_extent", 3)
    else:
        phi = get(keys, "geometry.prob_hi", 3)
        ext = [phi[d] - plo[d] for d in range(3)]
    ncell = get(keys, "amr.n_cell", 3, typ=int)
    per = get(keys, "geometry.is_periodic", 3, typ=int, default=[0, 0, 0])
    per_x, per_y = bool(per[0]), bool(per[1])
    if "amr.ref_ratio_vect" in keys:
        rr = get(keys, "amr.ref_ratio_vect", 3, typ=int)[:2]
    else:
        r = get(keys, "amr.ref_ratio", typ=int, default=[2])[0]
        rr = [r, r]
    bf = get(keys, "amr.blocking_factor", typ=int, default=[8])[0]
    nerr = get(keys, "amr.n_error_buf", typ=int, default=[1])[0]
    max_level = get(keys, "amr.max_level", typ=int, default=[0])[0]
    bfile = a.buildings or os.path.join(deck_dir, get(keys, "erf.buildings_file_name", typ=str)[0])

    dxc = [ext[0] / ncell[0], ext[1] / ncell[1]]
    dxf = [dxc[0] / rr[0], dxc[1] / rr[1]]
    # Grid edges stay where they are asked for when they fall on whole blocks
    # of the finest level, which the blocking factor of the refined level sets.
    step = [bf * dxf[0], bf * dxf[1]]

    x, y, z = read_map(bfile)
    xs = [plo[0] + dxf[0] * n for n in range(ncell[0] * rr[0] + 1)]
    ys = [plo[1] + dxf[1] * n for n in range(ncell[1] * rr[1] + 1)]
    h = interp_map(x, y, z, xs, ys)
    built, top = built_columns(h, rr)
    lab, nb = label(built, per_x, per_y)
    print(f"{bfile}: {nb} buildings on the {ncell[0]} x {ncell[1]} coarse columns of "
          f"{dxc[0]:g} x {dxc[1]:g} m (refined by {rr[0]} x {rr[1]})")
    for b in range(1, nb + 1):
        cols = columns_of(lab, b)
        ii = [c[0] for c in cols]; jj = [c[1] for c in cols]
        print(f"  building {b}: x {plo[0] + min(ii) * dxc[0]:g}-{plo[0] + (max(ii) + 1) * dxc[0]:g} m, "
              f"y {plo[1] + min(jj) * dxc[1]:g}-{plo[1] + (max(jj) + 1) * dxc[1]:g} m, "
              f"up to {max(top[i][j] for i, j in cols):g} m (coarse extent, ramps included)")
    nerr_note = (f"amr.n_error_buf = {nerr}{' (the AMReX default)' if 'amr.n_error_buf' not in keys else ''}: "
                 "ERF stops when a box is set explicitly with n_error_buf > 0; set amr.n_error_buf = 0")

    def to_index(lo, hi):
        """Coarse index box of the grids a real in_box becomes (n_error_buf is 0, checked above)."""
        # Clamped to the domain laterally first, as ERF_RefineBox.cpp does.
        lo = [max(lo[0], plo[0]), max(lo[1], plo[1])]
        hi = [min(hi[0], plo[0] + ext[0]), min(hi[1], plo[1] + ext[1])]
        if hi[0] <= lo[0] or hi[1] <= lo[1]:
            raise SystemExit(f"the box {lo} - {hi} m is empty inside the domain")
        ilo = int((lo[0] - plo[0]) / dxf[0]); ihi = int((hi[0] - plo[0]) / dxf[0] - 1)
        jlo = int((lo[1] - plo[1]) / dxf[1]); jhi = int((hi[1] - plo[1]) / dxf[1] - 1)
        # Snapped to the refinement ratio (ERF_RefineBox.cpp), then grown to
        # whole blocks of amr.blocking_factor fine cells, as AMReX makes the
        # grids, then to coarse columns.
        for m in (rr[0], bf):
            ilo -= ilo % m; ihi = ihi + (m - 1 - ihi % m)
        for m in (rr[1], bf):
            jlo -= jlo % m; jhi = jhi + (m - 1 - jhi % m)
        box = [ilo // rr[0], ihi // rr[0], jlo // rr[1], jhi // rr[1]]
        return [max(box[0], 0), min(box[1], ncell[0] - 1), max(box[2], 0), min(box[3], ncell[1] - 1)]

    def report(boxes, title):
        cls = classify(lab, nb, boxes, per_x, per_y)
        bad = [b for b, c in cls.items() if c == "crossing"]
        inside = [b for b, c in cls.items() if c == "inside"]
        where = "; ".join(f"x {plo[0] + b[0] * dxc[0]:g}-{plo[0] + (b[1] + 1) * dxc[0]:g} m, "
                          f"y {plo[1] + b[2] * dxc[1]:g}-{plo[1] + (b[3] + 1) * dxc[1]:g} m" for b in boxes)
        print(f"{title}: {where}: "
              f"{len(inside)} buildings inside, {nb - len(inside) - len(bad)} outside, "
              f"{len(bad)} crossing the edge{(' (' + ', '.join(map(str, bad)) + ')') if bad else ''}")
        return bad

    def propose(box, others=()):
        """Grow box until no building crosses the refined level it makes with the
        others, keeping the edges on whole blocks."""
        sx, sy = int(round(step[0] / dxc[0])) or 1, int(round(step[1] / dxc[1])) or 1
        def snap(bx):
            return [bx[0] - bx[0] % sx, bx[1] + (sx - 1 - bx[1] % sx),
                    bx[2] - bx[2] % sy, bx[3] + (sy - 1 - bx[3] % sy)]
        box = snap(box)
        for _ in range(4 * (ncell[0] + ncell[1])):
            cls = classify(lab, nb, [box] + list(others), per_x, per_y)
            bad = [b for b, c in cls.items() if c == "crossing"]
            if not bad:
                break
            for b in bad:
                for i, j in margin_cells(lab, b, per_x, per_y):
                    box = [min(box[0], i), max(box[1], i), min(box[2], j), max(box[3], j)]
            box = snap([max(box[0], 0), min(box[1], ncell[0] - 1), max(box[2], 0), min(box[3], ncell[1] - 1)])
        else:
            raise SystemExit("no box found: a building crosses a periodic seam; refine the whole domain instead")
        return box

    status = 0
    if a.all or a.region:
        if a.all:
            if nb == 0:
                raise SystemExit("the map has no buildings")
            cols = [(i, j) for i in range(len(lab)) for j in range(len(lab[0])) if lab[i][j] > 0]
            box = [min(c[0] for c in cols), max(c[0] for c in cols), min(c[1] for c in cols), max(c[1] for c in cols)]
        else:
            box = to_index(a.region[0::2], a.region[1::2])
        box = propose(box)
        if a.fit == "relaxed" and a.margin > 0:
            m = a.margin
            box = propose([max(box[0] - m, 0), min(box[1] + m, ncell[0] - 1),
                           max(box[2] - m, 0), min(box[3] + m, ncell[1] - 1)])
        report([box], f"proposed box ({a.fit}{f', {a.margin} coarse cells of padding' if a.fit == 'relaxed' else ''})")
        lo = (plo[0] + box[0] * dxc[0], plo[1] + box[2] * dxc[1])
        hi = (plo[0] + (box[1] + 1) * dxc[0], plo[1] + (box[3] + 1) * dxc[1])
        if nerr != 0:
            print(f"note: {nerr_note} (the lines below do)")
        print(f"\n# refined level for erf.ibseb (ibseb_refinement_box.py --fit {a.fit})")
        print(f"amr.max_level = 1\namr.n_error_buf = 0\namr.refine_whole_domain_dir = 2")
        print(f"erf.refinement_indicators = {a.name}")
        print(f"erf.{a.name}.max_level = 1")
        print(f"erf.{a.name}.in_box_lo = {lo[0]:g} {lo[1]:g}")
        print(f"erf.{a.name}.in_box_hi = {hi[0]:g} {hi[1]:g}")
    else:
        # Checking the deck's own boxes: ERF would stop on them as they stand.
        if nerr != 0:
            raise SystemExit(nerr_note)
        names = keys.get("erf.refinement_indicators", [])
        boxes = [n for n in names if f"erf.{n}.in_box_lo" in keys]
        if max_level < 1 or not boxes:
            print("the deck refines no static box (amr.max_level < 1 or no erf.<name>.in_box_lo); nothing to check")
            return 0
        ibox = {}
        for n in boxes:
            lev = get(keys, f"erf.{n}.max_level", typ=int, default=[max_level])[0]
            lo = get(keys, f"erf.{n}.in_box_lo", 2)
            hi = get(keys, f"erf.{n}.in_box_hi", 2)
            if len(keys[f"erf.{n}.in_box_lo"]) > 2 and (lo[2] > plo[2] or hi[2] < plo[2] + ext[2]):
                print(f"note: erf.{n} does not span the depth; buildings must also end two coarse cells below its top")
            if lev > 1:
                print(f"note: erf.{n} reaches level {lev}; only its level-1 box is checked here")
            for v, nm in ((lo[0], "in_box_lo x"), (lo[1], "in_box_lo y"), (hi[0], "in_box_hi x"), (hi[1], "in_box_hi y")):
                d = 0 if "x" in nm else 1
                if abs(((v - plo[d]) / step[d]) - round((v - plo[d]) / step[d])) > 1e-6:
                    print(f"note: erf.{n}.{nm} = {v:g} is not on a whole block of {step[d]:g} m; "
                          "the grids grow to the blocks, and those grids are what is judged")
            ibox[n] = to_index(lo, hi)
        # Level 1 is the union of the boxes: a building two boxes cover together is inside it.
        bad = report(list(ibox.values()), "level 1 (" + ", ".join(f"erf.{n}" for n in ibox) + ")")
        if bad:
            status = 1
            for n, box in ibox.items():
                others = [b for m, b in ibox.items() if m != n]
                part = [b for b in bad if any(covered([box], i, j) for i, j in margin_cells(lab, b, per_x, per_y))]
                if not part:
                    continue
                fix = propose(box, others)
                report([fix] + others, f"  erf.{n} grown so no building crosses level 1")
                print(f"  erf.{n}.in_box_lo = {plo[0] + fix[0] * dxc[0]:g} {plo[1] + fix[2] * dxc[1]:g}")
                print(f"  erf.{n}.in_box_hi = {plo[0] + (fix[1] + 1) * dxc[0]:g} {plo[1] + (fix[3] + 1) * dxc[1]:g}")
    return status


if __name__ == "__main__":
    sys.exit(main())
