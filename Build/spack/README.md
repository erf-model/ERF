# ERF's Spack package

`spack_repo/erf` is a [Spack](https://spack.io) package repository with one
package, `erf`, and one override, `amrex`.

## Getting Started

Use Spack 1.1 or newer and a recent `spack/spack-packages` checkout.

Register the repository from a local checkout or directly from GitHub:

```sh
# From a local ERF checkout
spack repo add /path/to/ERF/Build/spack/spack_repo/erf

# Or directly from GitHub, without a local checkout
spack repo add https://github.com/erf-model/ERF.git
```

Then install ERF or inspect its build options:

```sh
spack install erf                                      # latest release, MPI on, double precision
spack install erf@26.11 +netcdf +noahmp +rrtmgp +fft     # release with optional physics and I/O
spack info erf                                         # available variants
```

Keep Spack's builtin repository (`spack/spack-packages`, registered by default)
registered too. It supplies the build systems and dependencies, including
AMReX, MPI, NetCDF, and GoogleTest.

## Build options and limitations

For GPU builds, use `+cuda cuda_arch=<n>`, `+rocm amdgpu_target=<gfx>`, or
`+sycl %oneapi`.

The `+rrtmgp`, `+eamxx_shoc`, and `+noahmp` variants enable RTE-RRTMGP radiation,
EAMxx SHOC turbulence, and Noah-MP, respectively. The recipe fetches their
sources and the required EKAT stack at the commits pinned by the selected ERF
version.

The following platform and toolchain limits apply:

* RRTMGP and EAMxx SHOC are unavailable on Apple Silicon with ERF's EKAT pins.
* For ERF releases through 26.09, SYCL builds with RRTMGP or SHOC need oneAPI
  older than 2026. ERF 26.11 and development use the newer EKAT/Kokkos stack,
  which supports oneAPI 2026.
* OpenMP on macOS needs GCC for C and C++ throughout the dependency tree:
  `spack install erf +openmp %%c,cxx=gcc`.
* CUDA builds with `+tools` require ERF 26.09 or newer.
* GPU builds with RRTMGP or SHOC require a single GPU architecture for their
  bundled Kokkos.

## Documentation

See [Getting Started](https://erf.readthedocs.io/en/latest/GettingStarted.html)
for instructions on running an ERF example. The
[ERF documentation](https://erf.readthedocs.io/en/latest/) describes the model,
build options, and input parameters.

The [Spack documentation](https://spack.readthedocs.io/en/latest/) covers
package installation, compilers, environments, and package development.

## Development model

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for how to contribute to ERF
development.

This repository contains ERF's maintained Spack recipe. Spack's builtin
repository also carries an `erf` package, which may lag behind this version.
Registering this repository puts it ahead of builtin, so `spack install erf`
uses this recipe.

| File | Purpose |
|---|---|
| `spack_repo/erf/packages/erf/package.py` | ERF recipe for releases from 23.10 and `develop`, with variants for supported build options |
| `spack_repo/erf/packages/amrex/package.py` | AMReX override adding `particles_precision` for `ERF_PARTICLES_PRECISION` |
| `spack_repo/erf/repo.yaml` | Repository namespace (`erf`) and package API version |
| `check_pins.py` | Checks dependency pins against ERF's source tree and EKAT's dependencies |
| `../../spack-repo-index.yaml` | Locates this repository for `spack repo add <git URL>` |

When updating the package:

* After changing dependency pins, run `python3 Build/spack/check_pins.py`
  from the ERF root. The RRTMGP/EKAT and Noah-MP commits, and the AMReX commit
  used by `@develop` or a release in `_amrex_release_commits`, must match
  ERF's submodule pins. The E3SM commit must match
  `Build/GNU_Ekat/eamxx_clone.sh`. Use `--fix` to update the recipe's development
  Noah-MP and AMReX pins.
* For a new release, add the version and tarball's SHA-256 checksum
  (`spack checksum erf <version>`), plus its Noah-MP pin
  (`git ls-tree <tag> Submodules/Noah-MP`). Check the AMReX gitlink too; if
  the release needs fixes beyond the AMReX tag, add its exact commit to
  `_amrex_release_commits`.
* To expose a supported CMake option, add a variant with
  `when="@<first release>:"`, its `define_from_variant` in `cmake_args`, and
  any conflicts ERF's CMake enforces.
* Remove the AMReX override once spack-packages carries the
  `particles_precision` variant, which maps to `ERF_PARTICLES_PRECISION`.

ERF releases from 24.10 usually use the AMReX release of the same month.
ERF 24.10 uses AMReX 24.11 for an API added after the AMReX 24.10 tag. ERF
26.11 uses its exact AMReX submodule commit through `amrex@develop commit=...`:
that commit includes fixes beyond AMReX 26.11, including the FillPatch fix
for multilevel tests with Clang floating-point exception trapping. The commit
is fixed; it does not follow the AMReX development branch. Before ERF 25.03,
embedded boundaries are off, matching ERF's CMake settings. ERF 23.10 through
24.09 use their bundled AMReX.

## Container images

`.github/workflows/spack-containers.yml` publishes runnable images built
from this recipe: `ghcr.io/erf-model/erf:development` and the latest
release, each with a dated copy that never moves. They carry MPICH, NetCDF
and HDF5, Noah-MP, RRTMGP, FFT and particles for `x86_64_v3`; `erf_exec`
is on `PATH` and the examples are under `share/erf/examples`.

## License

See [LICENSE.md](../../LICENSE.md) for ERF's license.

## Citation

See the [main README](../../README.rst#citation) for publications to cite
when using ERF.
