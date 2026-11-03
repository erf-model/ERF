# ERF's Spack package

`spack_repo/erf` is a [Spack](https://spack.io) package repository with one
package, `erf`, and one override, `amrex`. It is kept here, next to ERF's
CMake, so that a new release, a new CMake option or a submodule bump and the
recipe that follows it land in the same pull request.

```
spack repo add /path/to/ERF/Build/spack/spack_repo/erf   # Spack 1.1 or newer
spack repo add https://github.com/erf-model/ERF.git        # Spack 1.1+, no checkout needed
spack install erf                                          # latest release, MPI on, double precision
spack install erf@develop +netcdf +noahmp +rrtmgp +fft     # ERF's development branch
spack info erf                                             # available variants
```

Spack's builtin repository (`spack/spack-packages`, registered by default) must
be present too: the recipe uses its build systems and its `amrex`, `mpi`,
`netcdf` and `googletest` packages. GPU builds take `+cuda cuda_arch=<n>`,
`+rocm amdgpu_target=<gfx>` or `+sycl %oneapi`; `+rrtmgp` and `+eamxx_shoc`
build the RTE-RRTMGP/EKAT stack and SHOC from the commits ERF pins, `+noahmp`
Noah-MP from its pin, each at the version the ERF release vendors.

| File | Purpose |
|---|---|
| `spack_repo/erf/packages/erf/package.py` | the recipe: release versions from 23.10 and `develop`, with variants for the supported build options |
| `spack_repo/erf/packages/amrex/package.py` | Spack's `amrex` plus a `particles_precision` variant, needed for `ERF_PARTICLES_PRECISION`; delete once spack-packages carries it |
| `spack_repo/erf/repo.yaml` | repository namespace (`erf`) and package API version |
| `check_pins.py` | the recipe's RRTMGP/EKAT, Noah-MP, E3SM and AMReX pins must equal ERF's submodule gitlinks; `--fix` rewrites the develop Noah-MP and AMReX pins |
| `../../spack-repo-index.yaml` | tells `spack repo add <git URL>` where this repository is |

## Requirements and limits

Use Spack 1.1 or newer and a recent `spack/spack-packages` checkout (the AMReX
recipe must include the tracer-header fix merged on 2026-10-19). The recipe has
been tested with Spack 1.2.2.

* RRTMGP and EAMxx SHOC are unavailable on Apple Silicon with ERF's EKAT pins.
* Through ERF 26.09, SYCL builds with RRTMGP or SHOC need oneAPI older than
  2026; development's newer EKAT/Kokkos stack supports oneAPI 2026.
* OpenMP on macOS needs GCC for C and C++ throughout the dependency tree:
  `spack install erf +openmp %%c,cxx=gcc`.
* CUDA builds cannot enable `+tools`. Builds with RRTMGP or SHOC require a
  single GPU architecture for their bundled Kokkos.

## Policies

* **This recipe is canonical.** Spack's builtin repository also carries an
  `erf` package; it is a snapshot that may lag. Registering this repository
  puts it ahead of builtin, so `spack install erf` uses this one.
* **Pins move with the tree.** The RRTMGP/EKAT, Noah-MP and E3SM commits the
  recipe fetches, and the AMReX commit `@develop` builds, must equal ERF's
  gitlinks; run `python3 Build/spack/check_pins.py` after a submodule bump
  (`--fix` rewrites the develop Noah-MP and AMReX pins). Releases from 24.10
  build the AMReX release of the same month (24.10 builds 24.11: its submodule
  is five commits past the 24.10 tag), with embedded boundaries off before
  25.03 as ERF's CMake sets them; 23.10 to 24.09 build the AMReX they bundle.
* **A release adds one step.** Add the version with the tarball's sha256
  (`spack checksum erf <version>`) and a Noah-MP pin line for it
  (`git ls-tree <tag> Submodules/Noah-MP`).
* **A CMake option adds a variant** with `when="@<first release>:"`, its
  `define_from_variant` in `cmake_args` and any conflict ERF's CMake enforces.
* **Overrides are temporary.** `packages/amrex` exists only for
  `particles_precision`; delete it once spack-packages carries that variant.

Other projects keep their recipe the same way, for example ParaView
(`Utilities/spack/repo`), Axom (`scripts/spack`), COSMA (`spack_repo/cosma`),
Rayleigh (`spack-repo`) and Amanzi (`config/spack`).
