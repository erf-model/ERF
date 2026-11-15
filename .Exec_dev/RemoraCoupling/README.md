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

Two groups of settings in `inputs_oneway` depart from ERF's MYNN defaults.
They are commented in place; briefly:

| setting | value | why |
|---|---|---|
| `pbl_mynn_{A1,A2,B1,B2,C1}` | MY82 constants | COAMPS uses Mellor--Yamada (1982), not Nakanishi--Niino |
| `pbl_mynn_{C2,C3,C5}` | 0 | NN09 extensions; zeroing them recovers the MY82 stability functions |
| `pbl_mynn_Lt_alpha` | 0.10 | the MY82 value, not NN09's 0.23 |

`erf.prob_name = Perlin07` selects this case: it sets the same initial TKE
profile as the `WPS` problem and additionally supplies the uniform potential
temperature tendency (`prob.advection_heating_rate`, +1 K/day) that stands in
for the large-scale subsidence and advection P07's domain omits, which `WPS`
does not.  Radiation is computed with RRTMGP (`erf.radiation_model = RRTMGP`);
set `erf.rrtmgp_file_path` to the directory holding the RRTMGP data files.

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
| `heat_budget.py` | Mixed-layer heat budget (storage, surface flux, flux through the control-volume top, radiative sink) plus an entrainment diagnostic.  Closes to 0.2--0.4 W/m2.  Expects run directories named `run_<case>` containing `plt72000` and `plt86400`.  `heat_budget.py <case> [<case> ...]` |
| `PlotPlts.py` | Pre-existing general-purpose profile plotter for u, v and `Kmv` from one or more plotfiles.  Note it labels `Kmv` as K_mv although the stored quantity is `rho*K`. |

`mynn_diag.py` and `heat_budget.py` assume the `bulk_coeff` surface layer with
the coefficients in `inputs_oneway`; both take `--Cd/--Ch/--Cq` or equivalent
if those are changed.

## Other inputs files here

`inputs_abl`, `inputs_couette`, `inputs_shear_ERFtoREMORA`,
`inputs_shear_two_way` and the other soundings belong to separate coupling
tests and are unrelated to the P07 validation.
