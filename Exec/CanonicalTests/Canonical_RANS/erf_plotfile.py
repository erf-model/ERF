"""Minimal reader for single-level AMReX (ERF) plotfiles.

Standard library only, so the regression check scripts run on a bare CI
runner. Returns planar (x-y) averages of cell-centred fields as a function
of the vertical index, plus the cell-centre heights of a flat mesh.

Limitations: one AMR level, double precision, native byte order as written
in the FAB headers (little endian on every platform ERF is tested on).
"""

import array
import os
import re
import struct


def _read_header(plotfile):
    with open(os.path.join(plotfile, "Header")) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    it = iter(lines)
    next(it)                          # version string
    ncomp = int(next(it))
    names = [next(it) for _ in range(ncomp)]
    ndim = int(next(it))
    time = float(next(it))
    finest = int(next(it))
    prob_lo = [float(v) for v in next(it).split()]
    prob_hi = [float(v) for v in next(it).split()]
    next(it)                          # ref ratios (an empty line when finest == 0)
    dom = next(it)                    # "((0,0,0) (nx-1,ny-1,nz-1) (0,0,0))"
    lo, hi = re.findall(r"\(([-\d,]+)\)", dom)[:2]
    lo = [int(v) for v in lo.split(",")]
    hi = [int(v) for v in hi.split(",")]
    next(it)                          # step numbers
    dx = [float(v) for v in next(it).split()]
    return dict(names=names, ndim=ndim, time=time, finest=finest,
                prob_lo=prob_lo, prob_hi=prob_hi, lo=lo, hi=hi, dx=dx)


def _read_cell_h(level_dir):
    """Return the box list and (file, offset) per box from Cell_H."""
    with open(os.path.join(level_dir, "Cell_H")) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    # version, how, ncomp, nghost, then the BoxArray
    ncomp = int(lines[2])
    nghost_line = lines[3].strip()
    idx = 4
    m = re.match(r"\((\d+)\s+(\d+)", lines[idx])
    nbox = int(m.group(1))
    idx += 1
    boxes = []
    for _ in range(nbox):
        lo, hi = re.findall(r"\(([-\d,]+)\)", lines[idx])[:2]
        boxes.append(([int(v) for v in lo.split(",")],
                      [int(v) for v in hi.split(",")]))
        idx += 1
    idx += 1                          # closing paren of the BoxArray
    nfab = int(lines[idx]); idx += 1
    fabs = []
    for _ in range(nfab):
        m = re.match(r"FabOnDisk:\s+(\S+)\s+(\d+)", lines[idx])
        fabs.append((m.group(1), int(m.group(2))))
        idx += 1
    nghost = [int(v) for v in re.findall(r"-?\d+", nghost_line)]
    ng = nghost[0] if nghost else 0
    return ncomp, ng, boxes, fabs


def _read_fab(path, offset, ncomp):
    with open(path, "rb") as fh:
        fh.seek(offset)
        header = b""
        while not header.endswith(b"\n"):
            header += fh.read(1)
        htxt = header.decode()
        # FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (1 2 3 4 5 6 7 8)))((lo) (hi) (0)) ncomp
        nbytes = int(re.search(r"FAB \(\((\d+)", htxt).group(1))
        # AMReX byte-order descriptor: "(1 2 ... 8)" is big endian,
        # "(8 7 ... 1)" is little endian (native on x86 and Apple silicon).
        order = re.search(r"\(\d+, \(([\d ]+)\)\)\)\(", htxt).group(1).split()
        little = order[0] == "8"
        lo, hi = re.findall(r"\(([-\d,]+)\)", htxt.split(")))")[1])[:2]
        lo = [int(v) for v in lo.split(",")]
        hi = [int(v) for v in hi.split(",")]
        n = 1
        for a, b in zip(lo, hi):
            n *= (b - a + 1)
        fmt = "d" if nbytes == 8 else "f"
        data = array.array(fmt)
        data.frombytes(fh.read(n * ncomp * nbytes))
        if little != (struct.pack("=h", 1) == struct.pack("<h", 1)):
            data.byteswap()
        return lo, hi, data


def planar_averages(plotfile, fields):
    """Planar (x-y) average of each named field per k index.

    Returns (z_cc, {field: [avg per k]}) for a flat single-level mesh.
    Raises KeyError if a field is missing from the plotfile.
    """
    hdr = _read_header(plotfile)
    names = hdr["names"]
    for f in fields:
        if f not in names:
            raise KeyError("field %s not in plotfile (have: %s)" % (f, ", ".join(names)))
    comps = [names.index(f) for f in fields]
    lo, hi, dx = hdr["lo"], hdr["hi"], hdr["dx"]
    nz = hi[2] - lo[2] + 1
    nxy = (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1)
    sums = {f: [0.0] * nz for f in fields}
    level_dir = os.path.join(plotfile, "Level_0")
    ncomp, ng, boxes, fabs = _read_cell_h(level_dir)
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = _read_fab(os.path.join(level_dir, fname), off, ncomp)
        nx = fhi[0] - flo[0] + 1
        ny = fhi[1] - flo[1] + 1
        nzf = fhi[2] - flo[2] + 1
        npts = nx * ny * nzf
        for f, c in zip(fields, comps):
            base = c * npts
            for k in range(blo[2], bhi[2] + 1):
                kk = k - flo[2]
                s = 0.0
                for j in range(blo[1], bhi[1] + 1):
                    jj = j - flo[1]
                    row = base + (kk * ny + jj) * nx + (blo[0] - flo[0])
                    s += sum(data[row: row + (bhi[0] - blo[0] + 1)])
                sums[f][k - lo[2]] += s
    avgs = {f: [v / nxy for v in sums[f]] for f in fields}
    z_cc = [hdr["prob_lo"][2] + (k + 0.5) * dx[2] for k in range(nz)]
    return z_cc, avgs, hdr


def read_surf_hist(path):
    """Return the last row of surf_hist.dat as a dict."""
    with open(path) as fh:
        rows = [ln.split() for ln in fh if ln.strip()]
    header = rows[0]
    last = [float(v) for v in rows[-1]]
    return dict(zip(header, last))


def read_fields(plotfile, fields):
    """Full 3D cell data of the named fields on a single-level plotfile.

    Returns (hdr, data) where data[field][i][j][k] is a nested list over the
    whole domain (0-based from the domain's low corner). Meant for the small
    canonical grids; it is pure Python.
    """
    hdr = _read_header(plotfile)
    names = hdr["names"]
    for f in fields:
        if f not in names:
            raise KeyError("field %s not in plotfile (have: %s)" % (f, ", ".join(names)))
    comps = [names.index(f) for f in fields]
    lo, hi = hdr["lo"], hdr["hi"]
    nx, ny, nz = [hi[d] - lo[d] + 1 for d in range(3)]
    out = {f: [[[0.0] * nz for _ in range(ny)] for _ in range(nx)] for f in fields}
    level_dir = os.path.join(plotfile, "Level_0")
    ncomp, ng, boxes, fabs = _read_cell_h(level_dir)
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = _read_fab(os.path.join(level_dir, fname), off, ncomp)
        fnx = fhi[0] - flo[0] + 1
        fny = fhi[1] - flo[1] + 1
        fnz = fhi[2] - flo[2] + 1
        npts = fnx * fny * fnz
        for f, c in zip(fields, comps):
            base = c * npts
            arr = out[f]
            for k in range(blo[2], bhi[2] + 1):
                for j in range(blo[1], bhi[1] + 1):
                    row = base + ((k - flo[2]) * fny + (j - flo[1])) * fnx
                    for i in range(blo[0], bhi[0] + 1):
                        arr[i - lo[0]][j - lo[1]][k - lo[2]] = data[row + (i - flo[0])]
    return hdr, out

