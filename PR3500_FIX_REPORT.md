# ERF PR3500 Fix Report

Date: 2026-09-09
Branch: `Use_SHOC_cldfrac_with_RRTMGP`

## Merge

`origin/development` was fetched and merged into this branch.

- Development tip merged: `489a2aa64b8b49986096a3445285d74bcea73739`
- Pre-merge branch tip: `2515e8db62903e9aa20f28bc389466e9e0883a71`
- Merge base: `0f7b33d46713f464d4eda60453bbca3b60e13d8b`
- Merge commit: `185ae6713`
- Repair commit: the branch `HEAD` containing this report.

The only textual conflict was `Docs/sphinx_doc/Inputs.rst`. The development
version of the table was retained and the PR3500 `erf.rad_use_shoc_cldfrac`
entry was restored in the merged table.

## Repair disposition

### R1: Simple-radiation interface

`RadiationSimple::Run` now accepts the same optional liquid-cloud-fraction
pointer as the common radiation interface. Simple radiation intentionally
ignores the pointer, preserving its existing physics while restoring the
interface contract.

### R2: Shared cloud-mass mask

Both liquid and ice condensate conversions in `Radiation::run_impl` now use
the shared total cloud mask. The existing liquid-only SHOC diagnostic, binary
positive-ice rule, fraction floor, in-cloud mixing-ratio cap, effective radii,
snow handling, and single kg/m2-to-g/m2 conversion are retained.

The production converter and the regression witness share
`radiation_cloud_mass`. For `q=1e-4`, `rho=1`, and `dz=100`, the witness gives:

- liquid-only fraction `0.25`: in-cloud mass `0.04`, grid-mean mass `0.01`;
- ice-only or mixed total fraction `1`: in-cloud mass `0.01`, grid-mean mass
  `0.01` for each phase;
- fraction floor `1e-4` and `q=1e-8`: in-cloud mass `0.01`, grid-mean mass
  `1e-6`;
- cap case `q=0.01`, total fraction `1`: in-cloud and grid-mean mass `0.5`.

The shared-mask limitation is documented: liquid occupying only part of an
ice-cloud layer cannot be represented by this contract.

### R3: Sampling flag

`rad_do_subcol_sampling` now selects the existing MCICA maximum-random-overlap
path when true and a deterministic binary clear/cloudy band-to-g-point map
when false. Both shortwave and longwave paths use the selection. Deterministic
mapping initializes all shortwave optical quantities (`tau`, `ssa`, `g`) and
longwave optical depth for clear and cloudy states.

Native SHOC plus `rad_use_shoc_cldfrac=true` plus disabled sampling is rejected
at initialization with actionable guidance. Binary cloud fractions remain
available with sampling disabled.

### R4: Native SHOC ordering and pointer contract

The merged development ordering is preserved: native SHOC advances first,
updated state and face quantities are refilled, and Native-SHOC RRTMGP receives
the same-step `shoc_cldfrac` diagnostic afterward. Existing optional-pointer
validation continues to check cell centering, component count, BoxArray, and
DistributionMapping compatibility.

### R5: Documentation

The Inputs, Plotfile3DReference, PBLschemes, and Native SHOC README documents
now describe the liquid-only diagnostic, host `qc` ownership, shared total
mask, binary ice rule, sampling behavior, and the mixed-phase limitation.

### R6: Fallback and cadence

`rad_use_shoc_cldfrac=false` remains a fraction-only binary comparison. It does
not restore pre-coupling radiation ordering. Native SHOC/RRTMGP still runs
radiation after the Native SHOC update and refill.

### R7: Numeric test portability

SHOC radiation tests no longer use `EXPECT_DOUBLE_EQ` for `amrex::Real`.
Exact comparisons cast literals to `amrex::Real`; mass-contract checks use a
precision-scaled tolerance.

### R8: Evidence and review status

This report records the implementation, build, and test evidence below. The
requested external P3 review was not available in this workspace; no external
review claim is made. Full scientific flux/heating-rate and multi-level
runtime validation remain follow-up work requiring usable RRTMGP coefficient
files and a CUDA device.

## Verification

Passed:

- `bash -n MyBuildcldfrac/cmake.sh`
- `git diff --check`
- Clean requested build using `MyBuildcldfrac/cmake.sh` with its configured
  CUDA 12.9, MPI, Kokkos, HDF5/NetCDF, Noah-MP, and RRTMGP options.
- `erf_exec` and `erf_shoc_test` linked successfully.
- `ERF_ShocRadiationCouplingTests.cpp` compiled in the CUDA/unit-test
  configuration, including the new shared-mask and sampling tests.
- The focused executable linked successfully using the generated CMake link
  script.
- `LD_PRELOAD=/nopt/cuda/12.9/lib64/stubs/libcuda.so ./erf_shoc_unit_tests
  --gtest_list_tests` exited 0 and listed the complete SHOC suite, including
  `ShocRadiationCloudFraction.*`.

The focused runtime filter was attempted with:

```text
LD_PRELOAD=/nopt/cuda/12.9/lib64/stubs/libcuda.so \
  ./erf_shoc_unit_tests --gtest_filter=ShocRadiationCloudFraction.*
```

It returned the project skip code 77 because the node reports the CUDA stub
driver and zero devices. Assertions were therefore compiled but not executed
on this host.

The clean `cmake.sh` build log is:
`MyBuildcldfrac/cmake_build_20260909_203757.log`.

## Suggested PR text

Summary: Merge current development into the Native SHOC/RRTMGP branch and
repair PR3500 cloud-fraction coupling. Native SHOC now supplies a same-step
liquid-only PDF diagnostic to RRTMGP after state refill; liquid and ice use one
shared total optical mask; and the sampling flag has explicit MCICA and
deterministic binary behavior with an invalid-combination startup check.

Tests: clean `MyBuildcldfrac/cmake.sh` build passed; the focused SHOC test
executable compiled and linked; runtime SHOC tests were skipped on this host
because no CUDA device/driver is available.
