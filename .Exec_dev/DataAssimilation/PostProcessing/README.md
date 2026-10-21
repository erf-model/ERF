# Post-processing for Data Assimilation

This directory contains a post-processing tool that quantifies how close a simulation is to the "ground truth" (the observations) by computing the RMS error between pairs of plotfiles at a series of output times. It is used to compare, for example, an ensemble/free run without data assimilation against the truth, or a data-assimilation run against the truth.

## What it does

For each pair of plotfiles (ground truth, simulation) listed in the input files, the code:

1. Reads both plotfiles and extracts the variables listed in `vars.txt` as cell-centered data.
2. Copies the simulation data onto the ground-truth `BoxArray`/`DistributionMapping`, so the two plotfiles may use different box decompositions. Both must have the same mesh resolution and domain.
3. Computes the difference (truth - simulation) and the RMS over all cells and all listed variables:

   `rms = sqrt( sum_{cells, vars} (truth - sim)^2 / (Ncells * Nvars) )`

4. Appends the result to `rms_error.txt`.

Note that the RMS is taken over all variables together without normalization, so variables with larger magnitudes (e.g. `density`, `theta`) dominate the value.

## Inputs

All inputs are read from the **current working directory** (file names are hard-coded in `main.cpp`):

| File | Description |
|------|-------------|
| `ground_truth_plt_files.txt` | Paths to the ground-truth plotfiles, one per line. |
| `coarse_simulation_plt_files.txt` | Paths to the simulation plotfiles, one per line. |
| `vars.txt` | Names of the variables to include in the error, one per line (e.g. `density`, `theta`, `x_velocity`, `y_velocity`, `z_velocity`, `qv`, `qc`, `qrain`). |

Requirements:

- The two lists must have the same number of entries; line *i* of one is compared with line *i* of the other (they should correspond to the same time).
- Empty lines are ignored.
- Every variable in `vars.txt` must exist in both plotfiles.
- The ground-truth and simulation plotfiles must be on the same mesh (the simulation plotfiles are on the coarse mesh, so the ground truth should be the coarsened version produced by the pre-processing step).

## Output

`rms_error.txt`, with one line per plotfile pair:

```
# iter rms_error
0 2.451202805974638
1 5.175938235335174
...
```

`iter` is the zero-based index of the pair in the input lists. The values are also printed to stdout.

## Build

Edit `AMREX_INC` and `AMREX_LIB` in the `Makefile` to point to your AMReX installation, then:

```bash
make -j8
```

This produces the executable `out`. Use `make clean` to remove build products.

## Run

From this directory, after preparing the three input files:

```bash
./out
```

or in parallel (e.g. with MPI):

```bash
srun -n 4 ./out
```

No command-line arguments are required. Only the I/O rank writes `rms_error.txt`.
