# Parent project test

## Purpose

Checks that ERF builds from a parent CMake project with `add_subdirectory`, the way
a repository that holds ERF as a git submodule builds it. The parent compiles in the
example structure model (`Tests/StructureModelExample`) from a folder outside the ERF
tree, through `ERF_STRUCTURE_MODEL_DIRS`, builds `erf_exec` and runs one short case.
Three more tests register through ERF's general CTest helpers
(`Tests/ERFTestHelpers.cmake`), the way a parent repository registers its own cases.
With unit tests on, one more runs the model's GoogleTest, which the model folder adds to
`erf_unit_tests` through `erf_structure_model_add_unit_tests`.

This folder is not part of the ERF build. Configure it on its own, as below.

## Files

| File | Role |
|---|---|
| `CMakeLists.txt` | The parent project: copies the example model outside the ERF tree, sets `ERF_ENABLE_MOVING_BODIES` and `ERF_STRUCTURE_MODEL_DIRS` as normal variables, calls `add_subdirectory` on ERF, adds the test `ParentProject_StructureExample`, includes `Tests/ERFTestHelpers.cmake` and registers `ParentProject_HelperAbort`, `ParentProject_HelperRegression` and `ParentProject_HelperRegressionDetects` with it, and adds `ParentProject_ModelUnitTest` when `erf_unit_tests` is built. |
| `test_files/ParentProject_HelperRegression/ParentProject_HelperRegression.i` | The deck of `ParentProject_HelperRegression`: the model's deck with one plotfile at the last step and no checkpoint. Keep its other lines the same as the model's deck. |
| `test_files/ParentProject_HelperRegressionDetects/ParentProject_HelperRegressionDetects.i` | The same deck under the name `add_test_r` reads for `ParentProject_HelperRegressionDetects`. |

`ParentProject_StructureExample` and `ParentProject_HelperAbort` run the model's own deck,
`Tests/StructureModelExample/inputs_structure_example`.

| Cache variable | Default | Effect |
|---|---|---|
| `ERF_PARENT_ERF_DIR` | the ERF tree two levels up | ERF source tree passed to `add_subdirectory` |
| `ERF_PARENT_MODEL_DIR` | `<build>/OutsideModel` | Where the example model is copied; must be outside `ERF_PARENT_ERF_DIR`, or configure stops |

## Data

No external data. The decks are in this folder and in `Tests/StructureModelExample`,
and the gold plotfile is made by the first test of each run.

## How to run

Use a build directory outside the ERF tree, so the model copy is outside it too:

```bash
cmake -S Tests/ParentProject -B /tmp/erf-parent-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/erf-parent-build --target erf_exec fvolumesum fcompare -j 4
ctest --test-dir /tmp/erf-parent-build -R ParentProject --output-on-failure
```

Other ERF options can be added on the command line (for example `-DERF_ENABLE_MPI=ON`).
With `-DERF_ENABLE_UNIT_TESTS=ON` (and the `Submodules/googletest` submodule), the build
of `erf_unit_tests` also takes the example model's unit test
`StructureModelExample.SupportIsThePointAndItsKernelReach`, added through
`erf_structure_model_add_unit_tests` from the copy outside the ERF tree.

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

The two tests registered through the helpers pass when:

- `ParentProject_HelperAbort` (`add_test_abort`): the run with
  `structure.example.epsilon=0` stops with `structure.example.epsilon must be > 0 and
  finite, got 0`. It shows that the helper finds `erf_exec` in ERF's binary folder from
  the parent's directory.
- `ParentProject_HelperRegression` (`add_test_r`, double precision only): `fcompare`
  finds the plotfile `plt00016` equal, to ERF's default tolerances (2e-10), to the gold
  plotfile in the parent's own gold folder, `<build>/GoldFiles`. That gold is the plotfile
  `ParentProject_StructureExample` checked against the closed form (CTest fixture
  `ParentProjectGold`), written by the same executable from the same case, so the
  comparison checks the helper's plumbing (deck folder, gold folder, launcher, `fcompare`)
  and not the physics; the physics is the closed-form check above. Both run on one rank,
  also with `ERF_ENABLE_MPI=ON` (the parent sets `ERF_TEST_NRANKS` to 1).
- `ParentProject_HelperRegressionDetects` (`add_test_r`, double precision only): the same
  case with 51 N instead of 50 N, against the same gold. The flow is linear in the force,
  so `fcompare` must report a relative difference of 1/50 = 0.02 in `x_velocity`
  (measured 0.0200000124); the test passes when its output shows that number. It shows
  that the comparison above would catch a different plotfile.
- `ParentProject_ModelUnitTest` (with `ERF_ENABLE_UNIT_TESTS=ON`): `erf_unit_tests
  --gtest_filter=StructureModelExample.*` runs exactly one test and it passes. A hook that
  adds no file gives "0 tests" and fails it.

## CTest

`ParentProject_StructureExample`, `ParentProject_HelperRegression` and
`ParentProject_HelperRegressionDetects` (double precision), `ParentProject_HelperAbort`,
and `ParentProject_ModelUnitTest` (unit tests on), registered in the parent build only.
The `parent-project` workflow in `.github/workflows` runs the first four on Linux, with
MPI and unit tests off. The model's GoogleTest also runs in ERF's own CI (the
particles-ON leg of `ci.yml`), which checks that `StructureModelExample` is listed.
