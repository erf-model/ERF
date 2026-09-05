#!/usr/bin/env python3
"""Check the phase 1 face storage against the blanking mask in the plotfile.

Reads a single-level AMReX plotfile, takes terrain_IB_mask, counts the wall
faces per direction as transitions between fluid (< 0.5) and solid (>= 0.5)
cells, and compares with the counts the [IBSEB] line reports. Also checks
that ibseb_nfaces summed over the domain equals the total face count and
that ibseb_tskin equals T_skin_init wherever a face exists.
"""
import sys, os, re, numpy as np

def read_plotfile(path):
    """Single-level AMReX plotfile into {name: array[nx, ny, nz]}."""
    hdr = open(os.path.join(path, "Header")).read().split("\n")
    nvar = int(hdr[1]); names = hdr[2:2 + nvar]
    i = 2 + nvar
    dim = int(hdr[i]); i += 1
    i += 2  # time, finest level
    i += 3  # prob lo, prob hi, ref ratio
    m = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", hdr[i])[0]
    lo = [int(m[0]), int(m[1]), int(m[2])]
    n = [int(m[3]) - lo[0] + 1, int(m[4]) - lo[1] + 1, int(m[5]) - lo[2] + 1]
    cell_h = open(os.path.join(path, "Level_0", "Cell_H")).read().split("\n")
    k = 0
    while not cell_h[k].startswith("("): k += 1
    nfab = int(cell_h[k].strip("(").split()[0]); k += 1
    boxes = []
    for _ in range(nfab):
        b = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", cell_h[k])[0]
        boxes.append([int(v) for v in b]); k += 1
    fabs = []
    for line in cell_h:
        if line.startswith("FabOnDisk:"):
            parts = line.split(); fabs.append((parts[1], int(parts[2])))
    assert len(fabs) == nfab, (len(fabs), nfab)
    data = {nm: np.zeros(n, dtype=np.float64) for nm in names}
    for b, (fname, off) in zip(boxes, fabs):
        with open(os.path.join(path, "Level_0", fname), "rb") as f:
            f.seek(off)
            f.readline()  # FAB header line
            nx, ny, nz = b[3] - b[0] + 1, b[4] - b[1] + 1, b[5] - b[2] + 1
            arr = np.frombuffer(f.read(8 * nx * ny * nz * nvar), dtype="<f8").reshape(nvar, nz, ny, nx)
        sl = (slice(b[0]-lo[0], b[3]-lo[0]+1), slice(b[1]-lo[1], b[4]-lo[1]+1), slice(b[2]-lo[2], b[5]-lo[2]+1))
        for c, nm in enumerate(names):
            data[nm][sl] = arr[c].transpose(2, 1, 0)
    return data

def expected_faces(mask, periodic=(True, True, False)):
    solid = mask >= 0.5
    counts = []
    for d in range(3):
        a = solid
        b = np.roll(solid, -1, axis=d)
        trans = a != b
        if not periodic[d]:
            sl = [slice(None)] * 3; sl[d] = slice(-1, None); trans[tuple(sl)] = False
        counts.append(int(trans.sum()))
    return counts

if __name__ == "__main__":
    plt, log, tskin_init = sys.argv[1], sys.argv[2], float(sys.argv[3])
    d = read_plotfile(plt)
    exp = expected_faces(d["terrain_IB_mask"])
    line = [l for l in open(log) if "[IBSEB] lev=0" in l][-1]
    got = [int(re.search(r"%s=(\d+)" % ax, line).group(1)) for ax in ("x", "y", "z")]
    total_reported = int(re.search(r"faces=(\d+)", line).group(1))
    nf_sum = int(round(d["ibseb_nfaces"].sum()))
    has = d["ibseb_nfaces"] > 0
    tmax_dev = float(np.abs(d["ibseb_tskin"][has] - tskin_init).max()) if has.any() else 0.0
    ok = (exp == got) and (nf_sum == total_reported) and (tmax_dev < 1e-9)
    print(f"expected x/y/z faces from mask: {exp}  reported: {got}  nfaces_sum: {nf_sum}  tskin_dev: {tmax_dev:.1e}  -> {'PASS' if ok else 'FAIL'}")
    sys.exit(0 if ok else 1)
