# SEB/Phase1_Storage

Phase 1 of the immersed-boundary surface energy balance: the face storage.
No physics yet; the skin temperature stays at its initial value. The test
checks that the wall faces of resolved buildings are found, stored, written
to the plotfile and the per-building CSV, and survive a checkpoint and
restart, and that none of it depends on the number of ranks.

```
./run_storage.sh /path/to/erf_exec        # NP=4 by default for the parallel runs
```

## The scenario

The ImmersedForcingTest skyscraper: 640 m cube domain on a 5 m grid, one
building from `skyscraper_5m_128x128.txt` as immersed forcing, a convective
sounding, four steps. `erf.ibseb.enable` with four slab layers and a 300 K
initial skin temperature; a report row every step.

## What is checked

1. **Face detection against the mask.** `check_storage.py` reads
   `terrain_IB_mask` from the plotfile, counts the fluid-to-solid
   transitions per direction, and compares with the `[IBSEB]` line's face
   counts. It also sums `ibseb_nfaces` over the domain against the total
   and checks `ibseb_tskin` equals the initial value wherever a face exists.
2. **Rank independence.** The same counts and area on one and four ranks.
3. **Checkpoint round trip.** A run to step 2 writes `IBSEBState`; the
   restart to step 4 reports that it restored the face state, and its last
   CSV row equals the straight run's.

## Plotting

```
python3 plot_storage.py plt00004            # writes plots/*.png
```

`plot_storage.py` uses yt to draw a horizontal slice of the mask and of
`ibseb_nfaces` through the building, and vertical slices of `ibseb_tskin`
through its centroid, with the mask outline. `--z`, `--x`, `--y` move the
slices and `--out` the directory. The face-bearing cells form a one-cell
ring around the footprint, two faces at the corners, and a cap over the
roof; in this phase every face is at the initial skin temperature.

## Reference

```
== straight, 1 rank
[IBSEB] lev=0 step=3 t=0.5 faces=2056 (x=900 y=900 z=256) buildings=1 area=51400 m2 T_skin_min=300 T_skin_max=300
expected x/y/z faces from mask: [900, 900, 256]  reported: [900, 900, 256]  nfaces_sum: 2056  tskin_dev: 0.0e+00  -> PASS
== straight, 4 ranks
expected x/y/z faces from mask: [900, 900, 256]  reported: [900, 900, 256]  nfaces_sum: 2056  tskin_dev: 0.0e+00  -> PASS
rank independence: PASS (faces=2056 (x=900 y=900 z=256) buildings=1 area=51400 m2)
== checkpoint at step 2, restart to step 4 (4 ranks)
state restored: yes
restart CSV row matches straight: PASS
ALL PASS
```
