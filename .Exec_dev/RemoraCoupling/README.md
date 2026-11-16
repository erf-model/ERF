# ERF + REMORA coupling: the Perlin et al. (2007) upwelling problem

This directory holds the ERF side of the idealized coastal upwelling validation
used in the ERF+REMORA coupling paper, following

> N. Perlin, E. D. Skyllingstad, R. M. Samelson and P. L. Barbour, 2007:
> *Numerical simulation of air--sea coupling during coastal upwelling.*
> J. Phys. Oceanogr. **37**, 2081--2093.  (hereafter P07)

`inputs_oneway` reproduces P07's **uncoupled** control, in which the sea surface
temperature seen by the atmosphere is held fixed.  Everything needed is in that
file; no command-line overrides are required.

```
mpiexec -n <N> ../../Exec/ERF3d.<...>.ex inputs_oneway
```

72 h at dt = 3 s, writing a plotfile every 12 h.  The problem is uniform in x
and y, so for parameter studies a 4 x 4 column gives the identical answer (to
seven significant figures) in a few minutes:

```
mpiexec -n 1 ../../Exec/ERF3d.<...>.ex inputs_oneway \
        erf.n_cell="4 4 48" erf.prob_hi="4000 4000 9000."
```

## What is non-default about the configuration

Three groups of settings in `inputs_oneway` depart from ERF's MYNN defaults.
They are commented in place; briefly:

| setting | value | why |
|---|---|---|
| `pbl_mynn_{A1,A2,B1,B2,C1}` | MY82 constants | COAMPS uses Mellor--Yamada (1982), not Nakanishi--Niino |
| `pbl_mynn_{C2,C3,C5}` | 0 | NN09 extensions; zeroing them recovers the MY82 stability functions |
| `pbl_mynn_Lt_alpha` | 0.10 | the MY82 value, not NN09's 0.23 |
| `pbl_mynn_surface_lm` | true | surface-layer length scale alone in the first cell, so the first-cell TKE approaches the MY82 surface value instead of dipping |

`erf.prob_name = Perlin07` selects this case: it sets the same initial TKE
profile as the `WPS` problem and additionally supplies the uniform potential
temperature tendency (`prob.advection_heating_rate`) that stands in for the
large-scale subsidence and advection P07's domain omits, which `WPS` does not.
P07 use +1 K/day to offset their radiation's average cooling; we use
+0.68 K/day, the daily-mean mixed-layer cooling RRTMGP gives here in July
(SW +0.97, LW -1.65 K/day).  With +1 K/day on top of RRTMGP the layer is net
heated, the surface layer goes stable and the boundary layer stalls near 300 m.  Radiation is computed with RRTMGP (`erf.radiation_model = RRTMGP`);
set `erf.rrtmgp_file_path` to the directory holding the RRTMGP data files.
The radiation starts at `2007-07-01 12:00:00` UTC, i.e. 0400 LST as in P07;
P07 give no date, only "summertime" conditions off Oregon.

The SST is an absolute temperature (`erf.most.abs_surf_temp`), as P07 plot it,
and the sounding starts at 1015 hPa.  P07 plot SST in K/degC against potential
temperature, so with a 1000 hPa surface the SST would be about 1 K too warm
in theta and the mixed layer would settle near 287 K instead of P07's 286 K.

Two further points that are easy to get wrong:

* `Kmv` and `Khv` in the plotfile are stored as `rho*K` (kg/m/s).  `nut` is the
  eddy viscosity in m2/s.  Compare `nut` against Fig. 8a of P07.
* With `erf.is_land = 0` the surface mixing ratio is overwritten every step with
  the saturation value at the sea surface temperature, so `most.surf_moist` has
  no effect.

## Known discrepancies against P07

The discrepancies previously listed here were measured with an earlier
configuration and need to be regenerated for the current `inputs_oneway`.

## Scripts

All take plotfile or run-directory arguments and need `yt`, `numpy` and
`matplotlib`.

| script | what it does |
|---|---|
| `make_paper_figure.py` | The six-panel validation figure for the paper: theta, u, v, TKE, K_m and L_m against height, overlaid on values pixel-traced from P07's Figs. 7 and 8.  Writes `perlin_validation.pdf/.png`.  Takes a run directory, default `run_oneway`; run it from this directory. |
| `compare_perlin.py` | Quicker exploratory version of the same comparison for one or more plotfiles, plus a one-line summary table (peak K_m, boundary layer depth, surface fluxes).  `compare_perlin.py <plt> [<plt> ...] --labels a,b` |
| `mynn_diag.py` | Recomputes the MYNN length scales (`l_S`, `l_T`, `l_B`), the stability functions and `K_m` offline from a plotfile column, exactly as `ERF_ComputeDiffusivityMYNN25.cpp` does, and compares against what ERF stored.  Used to verify the scheme is evaluating NN09 faithfully; reproduces ERF to 0.5%.  `mynn_diag.py <plt>` |
| `heat_budget.py` | Mixed-layer heat budget below 900 m: storage, surface sensible flux, flux through the column top, the large-scale heating (read from each run's `job_info`) and the RRTMGP heating (from the `qsrc_sw`/`qsrc_lw` plot variables), plus an entrainment diagnostic.  The radiation and surface flux are time-averaged over every plotfile in the interval, so write plotfiles hourly (`erf.plot_int_1 = 1200`); the residual is then a few W/m2.  `heat_budget.py <run_dir> [...] [--t0 plt57600 --t1 plt86400]` |
| `PlotPlts.py` | Pre-existing general-purpose profile plotter for u, v and `Kmv` from one or more plotfiles.  Note it labels `Kmv` as K_mv although the stored quantity is `rho*K`. |

`mynn_diag.py` and `heat_budget.py` assume the `bulk_coeff` surface layer with
the coefficients and the 287.0 K absolute SST in `inputs_oneway`.  If those
change, `mynn_diag.py` takes `--Cd/--Ch/--Cq/--surf_temp`; in `heat_budget.py`
edit `CD, CH, CQ` and `SST` at the top of the file.

## Other inputs files here

`inputs_abl`, `inputs_couette`, `inputs_shear_ERFtoREMORA`,
`inputs_shear_two_way` and the other soundings belong to separate coupling
tests and are unrelated to the P07 validation.
