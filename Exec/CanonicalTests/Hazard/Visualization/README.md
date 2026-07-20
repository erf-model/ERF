# ERF-Hazard Visualisation Scripts

Two Python scripts for visualising ERF-Hazard canonical test outputs.
No C++ build required — run directly after an ERF-Hazard simulation.

---

## Scripts

| Script | Input | Requires |
|---|---|---|
| `plot_hazard_fields.py` | AMReX plotfile directory (`plt_NNNNN`) | `yt`, `matplotlib`, `numpy` |
| `plot_hazard_timeseries.py` | Plain-text CSV/DAT diagnostic files | `matplotlib`, `pandas`, `numpy` |

Install dependencies:
```bash
pip install yt matplotlib pandas numpy
```

---

## plot_hazard_fields.py — AMReX Plotfile Figures

Reads ERF plotfiles and produces PNG slice plots.

```bash
# Smoke plume + terrain + wind for HaboobFireHill
python plot_hazard_fields.py --plotdir path/to/plt_00020 --case HaboobFireHill

# Dust emission asymmetry on Gaussian hill
python plot_hazard_fields.py --plotdir path/to/plt_00020 --case DustGaussianHill

# Wind recirculation in open pit mine
python plot_hazard_fields.py --plotdir path/to/plt_00020 --case HaboobFirePit
```

### Output files per case

| Case | Outputs |
|---|---|
| `HaboobFireHill` | `_smoke_plan.png`, `_smoke_xz.png`, `_wind_sfc.png`, `_dust_emission.png`, `_dust_xz.png`, `_theta_xz.png` |
| `HaboobFireFlat` | Same as HaboobFireHill |
| `HaboobFirePit` | `_wind_recirculation.png`, `_smoke_xz.png`, `_dust_xz.png`, `_theta_xz.png` |
| `DustGaussianHill` | `_dust_emission.png`, `_wind_sfc.png`, `_dust_xz.png` |
| `DustGaussianPit` | `_wind_recirculation.png`, `_dust_emission.png`, `_dust_xz.png` |

### Required plotfile fields

Make sure the ERF `inputs` file includes these in `erf.plot_vars_1`:
```
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta
```
Smoke (`smoke`) and dust (`rhoadv_dust`) are written automatically
when `ERF_ENABLE_FIRE=ON` and `ERF_USE_DUST=ON`.

---

## plot_hazard_timeseries.py — Diagnostic Time Series

Reads `dust_diag.dat`, `smoke_diag.dat` (plain CSV, no yt needed).

### Terrain amplification comparison
```bash
python plot_hazard_timeseries.py --mode terrain_amplification \
    --flat  HaboobFireFlat/dust_diag.dat \
    --hill  HaboobFireHill/dust_diag.dat \
    --pit   HaboobFirePit/dust_diag.dat
# Output: terrain_amplification.png
```

### Single dust diagnostic
```bash
python plot_hazard_timeseries.py --mode dust_diag \
    --file HaboobFireHill/dust_diag.dat \
    --label "HaboobFireHill"
# Output: dust_diag_HaboobFireHill.png
```

### Smoke diagnostic
```bash
python plot_hazard_timeseries.py --mode smoke_diag \
    --file HaboobFireHill/smoke_diag.dat \
    --label "HaboobFireHill"
# Output: smoke_diag_HaboobFireHill.png
```

### Fire-dust coupling interaction contributions
```bash
python plot_hazard_timeseries.py --mode coupling \
    --baseline     FireDustBaseline/dust_diag.dat \
    --interaction1 FireDustInteraction1/dust_diag.dat \
    --interaction2 FireDustInteraction2/dust_diag.dat \
    --interaction3 FireDustInteraction3/dust_diag.dat \
    --all          FireDustInteractions123/dust_diag.dat
# Output: fire_dust_coupling.png
```

---

## Most Compelling Single Figure

The terrain amplification plot from `HaboobFireFlat` vs `HaboobFireHill` vs
`HaboobFirePit` is the clearest single figure for a paper or presentation —
it shows in one image how terrain geometry controls dust emission, using
only the plain-text diagnostic output with no post-processing tools.
