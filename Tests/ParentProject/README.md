# Parent project test

## Purpose

Checks that ERF builds from a parent CMake project with `add_subdirectory`, the way
a repository that holds ERF as a git submodule builds it. The parent compiles in the
example structure model (`Tests/StructureModelExample`) from a folder outside the ERF
tree, through `ERF_STRUCTURE_MODEL_DIRS`, builds `erf_exec` and runs one short case.

This folder is not part of the ERF build. Configure it on its own, as below.

## Files

| File | Role |
|---|---|
| `CMakeLists.txt` | The parent project: copies the example model outside the ERF tree, sets `ERF_ENABLE_MOVING_BODIES` and `ERF_STRUCTURE_MODEL_DIRS` as normal variables, calls `add_subdirectory` on ERF and adds the test `ParentProject_StructureExample`. |

The input deck is the model's own, `Tests/StructureModelExample/inputs_structure_example`.

| Cache variable | Default | Effect |
|---|---|---|
| `ERF_PARENT_ERF_DIR` | the ERF tree two levels up | ERF source tree passed to `add_subdirectory` |
| `ERF_PARENT_MODEL_DIR` | `<build>/OutsideModel` | Where the example model is copied; must be outside `ERF_PARENT_ERF_DIR`, or configure stops |

## How to run

Use a build directory outside the ERF tree, so the model copy is outside it too:

```bash
cmake -S Tests/ParentProject -B /tmp/erf-parent-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/erf-parent-build --target erf_exec fvolumesum -j 4
ctest --test-dir /tmp/erf-parent-build -R ParentProject --output-on-failure
```

Other ERF options can be added on the command line (for example `-DERF_ENABLE_MPI=ON`).

## Expected results

One point pushes on still air with 50 N along x, in a periodic box with no friction
or gravity. At t = 4 s the x-momentum of the air is F t = 200 kg m/s to round-off.
The volume integral of `x_velocity` is F t / rho_0 = 50 N x 4 s / (1 kg/m^3) =
200 m^4/s to about 1e-6 relative; the small difference comes from density changes of
order u/c. The test passes when:

- the plotfile holds the model's field `structure_example_fx`, so the model is compiled in;
- the log holds `t = 4 impulse_x = 200`, so the model's `advance()` ran every step;
- `fvolumesum` gives 200 to 1e-5 relative (it prints 6 significant digits).

ERF's own `Structures_ExampleModel` test checks the same values in the standalone
build, on 2 MPI ranks. The parent test runs on one rank (MPI is off unless the parent
turns it on), so it checks the build, not the parallel force normalisation. Like
`Structures_ExampleModel`, it is registered only with `ERF_PRECISION=DOUBLE`.

## CTest

`ParentProject_StructureExample`, registered in the parent build only (double precision). The
`parent-project` workflow in `.github/workflows` runs it on Linux.
