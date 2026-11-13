# Parent executable test

## Purpose

Checks that a parent CMake project can use ERF the way ERF uses AMReX. ERF is a library
(`ERF::erf`), and the parent builds its own executable with `erf_setup_executable`, from its
own sources and its own problem file. The parent hands ERF none of its model-directory hooks
(`ERF_STRUCTURE_MODEL_DIRS`, `ERF_FIRE_ENGINE_DIRS`, `ERF_DUST_MODEL_DIRS`), and adds ERF
with `EXCLUDE_FROM_ALL`, so ERF builds only what the parent's executable needs.

Its sources are copies of the example structure model and of ERF's problem file, kept
outside the ERF tree as a repository keeps its own. The model registers itself from the
parent's executable.

`Tests/ParentProject` is the other layout: there, ERF compiles the parent's model into
ERF's own `erf_exec` through `ERF_STRUCTURE_MODEL_DIRS`.

This folder is not part of the ERF build. Configure it on its own, as below.

## Files

| File | Role |
|---|---|
| `CMakeLists.txt` | The parent project. It copies its sources outside the ERF tree, sets `ERF_ENABLE_MOVING_BODIES` as a normal variable, adds ERF with `add_subdirectory(... EXCLUDE_FROM_ALL)`, and builds `parent_model` with `erf_setup_executable`. It registers `ParentExecutable_StructureExample`, then includes `Tests/ERFTestHelpers.cmake` and registers `ParentExecutable_HelperAbort` with `ERF_TEST_EXE` set to `parent_model`. |

| Test | Checks |
|---|---|
| `ParentExecutable_StructureExample` | `parent_model` runs the example model's deck. The model's force is in the plotfile (`structure_example_fx`), its `advance()` ran every step (`t = 4 impulse_x = 200`), and the volume integral of `x_velocity` at 4 s is F t / rho_0 = 200 m^4/s to 1e-5 relative. Needs double precision. |
| `ParentExecutable_HelperAbort` | `add_test_abort` from ERF's helpers runs the parent's executable (`ERF_TEST_EXE`): the model stops on `structure.example.epsilon=0`. |

| Cache variable | Default | Effect |
|---|---|---|
| `ERF_PARENT_ERF_DIR` | the ERF tree two levels up | ERF source tree passed to `add_subdirectory` |
| `ERF_PARENT_SOURCE_DIR` | `<build>/OutsideSources` | Where the parent's sources are copied; it must be outside `ERF_PARENT_ERF_DIR`, or the configure stops |

## How to run

Use a build directory outside the ERF tree, so that the copied sources are outside it too:

```
cmake -S Tests/ParentExecutable -B ../parent-exec-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../parent-exec-build --target parent_model fvolumesum --parallel 4
ctest --test-dir ../parent-exec-build -R ParentExecutable --output-on-failure
```

The `ParentProject` workflow (`.github/workflows/parent-project.yml`) runs it.
