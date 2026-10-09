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

Three settings in `inputs_oneway` depart from ERF's MYNN defaults, and all
three are needed to reproduce P07.  They are commented in place; briefly:

| setting | value | why |
|---|---|---|
| `pbl_mynn_{A1,A2,B1,B2,C1}` | MY82 constants | COAMPS uses Mellor--Yamada (1982), not Nakanishi--Niino |
| `pbl_mynn_{C2,C3,C5}` | 0 | NN09 extensions; zeroing them recovers the MY82 stability functions |
| `pbl_mynn_Lt_alpha` | 0.10 | the MY82 value, not NN09's 0.23 |
| `pbl_mynn_Lt_taper_exp` | 1.35 | boundary-layer-depth taper on `l_T` (see below) |

The taper matters most.  Without it the MYNN master length scale saturates once
`kappa*z` exceeds `l_T` and stays flat to the inversion, giving a flat-topped
eddy viscosity; P07's falls steeply above about 150 m.  **The exponent 1.35 was
calibrated against their Fig. 8a and is not independently derived.**

`erf.prob_name = Perlin07` selects this case: it sets the same initial TKE
profile as the `WPS` problem and additionally supplies the uniform potential
temperature tendency used for the radiative cooling, which `WPS` does not.

Two further points that are easy to get wrong:

* `Kmv` and `Khv` in the plotfile are stored as `rho*K` (kg/m/s).  `nut` is the
  eddy viscosity in m2/s.  Compare `nut` against Fig. 8a of P07.
* With `erf.is_land = 0` the surface mixing ratio is overwritten every step with
  the saturation value at the sea surface temperature, so `most.surf_moist` has
  no effect.

## Known discrepancies against P07

Documented in the validation section of the paper, and summarized here so the
numbers are not lost:

* eddy viscosity agrees to within 15% from 50 to 400 m; peak 4.3 m2/s at 95 m
  against their 4.7 at 100 m;
* turbulent kinetic energy is good near the surface (within 9% to 100 m) but
  decays too fast above, reaching -58% at 350 m;
* the mixed layer is 0.45--0.7 K warmer than theirs, traceable to the use of
  constant rather than stability-dependent surface exchange coefficients;
* the boundary layer is shallower (420 m against 475 m) and is still deepening
  at 72 h.

## Scripts

All take plotfile or run-directory arguments and need `yt`, `numpy` and
`matplotlib`.

| script | what it does |
|---|---|
| `make_paper_figure.py` | The six-panel validation figure for the paper: theta, u, v, TKE, K_m and L_m against height, overlaid on values pixel-traced from P07's Figs. 7 and 8.  Writes `perlin_validation.pdf/.png`.  Takes a run directory, default `run_taper1.35`; run it from this directory. |
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
