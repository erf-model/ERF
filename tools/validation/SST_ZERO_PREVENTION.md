# Preventing SST = 0 / coastal "holes" in ERF runs — full workflow

**Standing rule: never launch an ERF run on un-gated inputs.** Every
wrfinput/wrflowinp pair that drives a NoahMP+MOST run must pass the SST=0 gate
(`check_sst_zero.py`) *before* it is staged. A water cell with `SST <= 0` hands
MOST an impossible surface temperature (`t_surf = 0 K`), which produces an
enormous surface heat/momentum flux and blows the dycore up along the coastline.

This note records the complete story (what the defect is, why it happened, how we
fixed it at the source, and the mandatory pre-launch check) so we do not
re-litigate it. See also `WPS_doctor/METGRID_TBL_SST_FIX.md` (source fix) and
`WPS_doctor/README.md` (post-hoc repair).

---

## 1. The defect

When driving WRF/ERF from **ERA5**, the default `SST` stanza in `METGRID.TBL`
leaves **coastal water cells with `SST = 0`** — the cells ERA5's coarse land/sea
mask calls "land" but the finer geogrid mask calls "water." metgrid finds no
valid source SST for those cells and falls through to `fill_missing = 0.`.

- `SST = 0` **over land** is normal (SST is only defined over water) — leave it.
- The bug is specifically **water cells** (`XLAND == 2` / `LANDMASK == 0`) with
  `SST <= 0`.

ERF reads `SST -> t_surf`. With `sst_update = 1` (hourly `wrflowinp` re-read),
the holes recur every hour, so a checkpoint restart inherits them — **cold start
is required** once inputs are cleaned.

Measured on the Sep-2024 ERA5 case:
- 9 km domain: **609** water SST=0 cells in wrfinput, ~103k across 169 hourly
  `wrflowinp` frames.
- 3 km domain: **6107** water cells (more coastline at finer resolution).

---

## 2. What did NOT fully fix it

1. **`interp_mask = LANDSEA(1)` + `masked = land` alone.** Necessary (stops
   land-contaminated SST bleeding onto the ocean) but **insufficient** — it does
   nothing for coastal water cells the ERA5 source simply doesn't cover. Those
   still fall through to `SST = 0`.
2. **A narrow `interp_option = sixteen_pt+four_pt`.** The basic interpolators
   can't reach isolated coastal water cells.
3. **real.exe.** It does **not** patch SST=0 water cells; it passes them straight
   into wrfinput/wrflowinp.

The hi-freq inputs had the mask applied and *still* shipped 609/6107 holes
because the `interp_option` was too narrow. This is the trap to remember.

---

## 3. The source fix (preferred) — widen the SST `interp_option`

Give SST the **same wide interpolation chain that the clean `SKINTEMP` stanza
uses**, ending in `search` (which scans outward for the nearest valid value):

```
name=SST
        interp_option=sixteen_pt+four_pt+wt_average_4pt+wt_average_16pt+search
        fill_missing=0.
        missing_value=-1.E30
        flag_in_output=FLAG_SST
        interp_mask = LANDSEA(1)
        masked=land
```

Both parts are required:
- `interp_mask = LANDSEA(1)` + `masked = land` — excludes land-contaminated
  source SST from ocean points.
- the full chain **ending in `search`** — fills coastal water cells from real
  neighbouring ocean SST during metgrid.

Then re-run **metgrid -> real.exe**.

**Verified result:** 9 km water SST holes **609 -> 0** across all 169 frames;
min SST over water 284 K (physical). Same fix applied to 3 km: **6107 -> 0**.

This is strictly better than post-hoc patching: the SST is *physically
interpolated* from neighbouring ocean values, not back-filled with a guess.

---

## 4. The fallback — WPS_doctor (inputs you can't re-metgrid)

For legacy inputs already produced without the wide stanza:

```bash
python WPS_doctor/wrf_preflight.py --wrfinput wrfinput_d01 --wrflowinp wrflowinp_d01 --fix
# run a SECOND pass on the .fixed files: pass 1 fills SST before flipping XLAND
# single-cell outliers, leaving the few reclassified land->water cells as new
# holes; pass 2 fills those.
python WPS_doctor/wrf_preflight.py --wrfinput wrfinput_d01.fixed --wrflowinp wrflowinp_d01.fixed --fix
```

Then cold-start from the cleaned files. WPS_doctor fills water holes from water
neighbours (isolated lakes from their own skin temp TSK, flagged suspicious).

---

## 5. The mandatory pre-launch gate (NEVER skip)

Run `check_sst_zero.py` on the exact files (or symlinks) the run will use. It
must report **ALL CLEAN** — 0 water `SST<=0` and 0 suspicious `0<SST<150` — for
both wrfinput and every wrflowinp frame:

```bash
module load frameworks            # provides netCDF4 on Aurora
python tools/validation/check_sst_zero.py wrfinput_d01 wrflowinp_d01 --mask wrfinput_d01
# RESULT: ALL CLEAN   -> ok to launch
# RESULT: HOLES/ERROR -> do NOT run; go to section 3 or 4
```

`wrflowinp` has no land mask, so `--mask` must point at the matching wrfinput
(XLAND). The gate checks every time level. Exit 0 = clean, 1 = holes, 2 = usage
error. Drop it into the job script as a hard gate:

```bash
python tools/validation/check_sst_zero.py wrfinput_d01 wrflowinp_d01 --mask wrfinput_d01 \
  || { echo "SST=0 holes — refusing to launch"; exit 1; }
```

---

## 6. Staging checklist (what we actually do per run)

1. Build inputs with the **wide SST `interp_option`** in METGRID.TBL (section 3).
2. real.exe -> wrfinput/wrfbdy/wrflowinp.
3. **SST=0 gate** on the real.exe output (section 5) -> must be ALL CLEAN.
4. Stage the run dir by **mirroring** a known-good run dir (copy all non-output
   files — NoahmpTable.TBL, RRTMGP coeffs, my_iofields_list.txt, etc.; hand-picked
   symlinks miss hidden Fortran unit-N opens).
5. Point wrfinput/wrfbdy/wrflowinp symlinks at the **gated** files, and
   **re-run the gate through the staging symlinks** (catches a stale/wrong link).
6. Cold start (holes corrupt surface state from t=0; a restart would inherit
   them).
7. Keep the in-job SST gate as a hard guard (section 5).

---

## 7. History / provenance

- Surfaced repairing an ERA5-driven Sep-2024 extratropical-cyclone case for ERF.
- The original GFS/HRRR-driven WFIP3 inputs had complete coastal SST and did
  **not** show the defect — it is specific to ERA5's coastal SST coverage.
- 9 km and 3 km validation runs (Sep 2024) both launched only after the wide
  `interp_option` produced 0 holes and the gate reported ALL CLEAN.
- Reference forum thread: forum.mmm.ucar.edu thread 9879 ("How to process ERA5
  invariant data — issue with land-sea mask").
