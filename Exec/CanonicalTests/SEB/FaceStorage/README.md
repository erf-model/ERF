# SEB/FaceStorage

Immersed-boundary surface energy balance: the face storage.
The test checks that the wall faces of resolved buildings are found, stored,
written to the plotfile and the per-building CSV, and survive a checkpoint
and restart face by face, and that none of it depends on the number of
ranks. The straight legs run the balance diagnostic, so the skin temperature
stays at its initial value and the counts can be checked on their own; the
checkpoint legs run it prognostic, so the state moves and a lost record
would show.

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
3. **Checkpoint round trip.** A prognostic run to step 2 writes
   `Level_0/IBSEBFaces`; the restart to step 4 reports that it restored the
   face state, its last CSV row equals the straight prognostic run's, and
   `compare_faces.py` finds every column of every face equal, to the 12
   digits the dumps are written with, between the two runs' face dumps (`faces_straight_prognostic.rank*.csv` against
   `faces_restart.rank*.csv`). The checkpoint legs run the balance prognostic
   (`inputs_roundtrip`) so the skin temperature and the slab move between
   steps 2 and 4 and a lost or misplaced record would show.
   `IBSEBFaces` holds one record per face: a header line and a marker word,
   then the 2056 face keys and 2056 x (2 + 4) values, 8 bytes each, 115 176
   bytes on this deck. The
   cell field it replaced (`IBSEBState`, 32 boxes over 7872 cells at 18
   components) carried 141 696 values for the same 12 336, about eleven
   times the data (erf-model/ERF#4019). The script checks the file is the
   records alone from its size. The same restart on one rank gives the same
   last CSV row and the same faces as the restart on the checkpoint's rank
   count, to 1e-8 relative: every face finds its record by its key, whatever
   the decomposition. That leg is not exact because steps 3 and 4 then run on
   another rank count and ERF's own solve differs in the last digit between
   rank counts (a straight run on one rank differs from one on four the same
   way, in the air temperature and wind the faces read); the balance residual,
   |balance| after the Newton solve at about 1e-8 W/m2, is compared within
   1e-6 W/m2 since it carries no digits. A lost or misplaced record would move
   the skin temperature by kelvins and fail both bands.
4. **Mismatched slab layers.** A restart from that checkpoint with
   `erf.ibseb.n_slab_layers = 6` must abort naming the deck's layer count,
   since the records hold four layers each.
5. **Other buildings.** A restart from that checkpoint with the height map
   rotated by 16 rows must abort: the faces of the rotated building have no
   records in the file. The rotation keeps the face count, so it is the
   per-face lookup, not the count check, that stops the run, and the script
   greps for that message.

## Files

| file | role |
| --- | --- |
| `inputs_base` | the deck: domain, building, sounding, `erf.ibseb.*` with the balance diagnostic |
| `inputs_straight` | four steps straight, diagnostic (checks 1 and 2) |
| `inputs_roundtrip` | `inputs_base` with the balance prognostic, the base of the checkpoint legs |
| `inputs_straight_prognostic` | four steps straight, prognostic, with a face dump: the reference of check 3 |
| `inputs_chk`, `inputs_restart` | checkpoint at step 2, restart to step 4 with a face dump |
| `inputs_restart_layers` | the restart with six slab layers (check 4) |
| `check_storage.py` | face counts from the mask against the `[IBSEB]` line (check 1) |
| `compare_faces.py` | two face dumps, face by face and column by column (check 3) |
| `plot_storage.py` | slices of the mask, `ibseb_nfaces` and `ibseb_tskin` |
| `run_storage.sh` | runs everything above |

## Plotting

```
python3 plot_storage.py plt00004            # writes plots/*.png
```

`plot_storage.py` uses yt to draw a horizontal slice of the mask and of
`ibseb_nfaces` through the building, and vertical slices of `ibseb_tskin`
through its centroid, with the mask outline. `--z`, `--x`, `--y` move the
slices and `--out` the directory. The face-bearing cells form a one-cell
ring around the footprint, two faces at the corners, and a cap over the
roof; with the balance diagnostic every face is at the initial skin temperature.

## Reference

```
== straight, 1 rank
[IBSEB] lev=0 step=4 t=0.5 faces=2056 (x=900 y=900 z=256) buildings=1 area=51400 m2 T_skin_min=300 T_skin_max=300 SW_abs_mean=208.888473 SW_abs_max=477.5777772 shadow_frac=0.01410505837 LW_net_mean=-69.87071641 H_mean=229.5998632 G_mean=6.063298012e-12 Q_ext_mean=0 resid_max=0 w_star_max=0 H_total_W=11801432.9706
expected x/y/z faces from mask: [900, 900, 256]  reported: [900, 900, 256]  nfaces_sum: 2056  tskin_dev: 0.0e+00  -> PASS
== straight, 4 ranks
expected x/y/z faces from mask: [900, 900, 256]  reported: [900, 900, 256]  nfaces_sum: 2056  tskin_dev: 0.0e+00  -> PASS
rank independence: PASS (faces=2056 (x=900 y=900 z=256) buildings=1 area=51400 m2)
== straight prognostic run to step 4 (4 ranks): the reference of the restart legs
[IBSEB] lev=0 step=4 t=0.5 faces=2056 (x=900 y=900 z=256) buildings=1 area=51400 m2 T_skin_min=293.1469161 T_skin_max=312.4323691 SW_abs_mean=208.888473 SW_abs_max=477.5777772 shadow_frac=0.01410505837 LW_net_mean=-66.80831905 H_mean=156.2510905 G_mean=-14.17093651 Q_ext_mean=0 resid_max=2.760521056e-08 w_star_max=0 H_total_W=8031306.05179
== checkpoint at step 2, restart to step 4 (4 ranks)
IBSEBFaces: 2056 faces x (2 + 4) values, 115176 bytes (the records alone): PASS
state restored: yes
restart CSV row matches straight (exact): PASS
2056 faces in faces_straight_prognostic, 2056 in faces_restart: all columns equal
restart faces match straight (every column of every face): PASS
== the same restart on 1 rank (every face finds its record in the 4-rank checkpoint)
restart on 1 rank equals the restart on the checkpoint's rank count (to 1e-8, resid_max within 1e-6): PASS
2056 faces in faces_restart, 2056 in faces_restart_np1: all columns equal to rtol 1e-08 (resid within 1e-06)
restart faces on 1 rank match the 4-rank restart (every face, to 1e-8): PASS
== restart with erf.ibseb.n_slab_layers = 6 against the 4-layer checkpoint (must abort)
mismatched slab layers rejected: PASS
== restart against a checkpoint of other buildings (must abort): the height map rotated by 16 rows
checkpoint of other buildings rejected: PASS
ALL PASS
```

The restart comparison on the checkpoint's rank count is exact for every column of the report row and of every face, including the columns read from the atmosphere (net longwave, sensible flux). Until October 2026 a restart rebuilt the immersed forcing's blanking without clearing the almost-fluid cells, and those two differed by up to 5e-5 relative (sensible flux 222.9486 against 222.9373 W/m2).
