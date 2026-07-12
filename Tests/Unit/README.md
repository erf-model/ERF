# ERF unit tests

ERF unit tests require CMake 3.18 or newer. Configure them with MPI when
parallel tests are needed:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
  -DERF_ENABLE_MPI=ON -DERF_ENABLE_TESTS=ON \
  '-DERF_PARALLEL_TEST_NRANKS=1;2;4'
cmake --build build --parallel
```

The main and SHOC binaries use the same `gtest_discover_tests()`-based
registration helper. Serial cases are discovered after the binaries are built
and appear individually in CTest. MPI cases are listed by the binary at build
time and registered once per case and configured rank count. Useful commands
are:

```bash
ctest --test-dir build -N -L unit
ctest --test-dir build -N -L parallel
ctest --test-dir build -L unit -LE parallel --output-on-failure --no-tests=error
ctest --test-dir build -L parallel --output-on-failure --no-tests=error
ctest --test-dir build -L infrastructure --output-on-failure --no-tests=error
ctest --test-dir build -L regression --output-on-failure --no-tests=error
```

Use `TEST` for one cohesive property, `TEST_P` when cases have meaningful
names and independently useful filters, and an internal loop for a dense
mathematical sweep whose iterations collectively establish one property.
Parameter suffixes should describe the physical or numerical identity, such as
`DryWarm`, `ColdUpperAir`, or `Upwind3Negative`.

Use the ERF test-only assertion helpers when a scaled comparison needs more
diagnostic context:

```cpp
ERF_EXPECT_NEAR(actual, expected, tolerance);
ERF_EXPECT_FINITE(actual);
ERF_EXPECT_NONNEGATIVE(condensate, tolerance);
```

Failures report actual and expected values, absolute and relative errors,
tolerance, and precision mode. Keep tolerance provenance next to the test
contract and use precision-aware values for single-precision builds. Add a
GoogleTest `PrintTo(const T&, std::ostream*)` overload for domain structures
when it makes parameterized failures or assertion messages readable; do not
add production stream operators solely for tests.

Device kernels may compute values, counters, or normalized errors, but all
GoogleTest assertions remain on the host. Synchronize device work before host
code reads results. Keep CPU/device coverage separate when initialization or
allocation requirements differ.

MPI tests must keep all ranks on the same control-flow path. Do not use a
rank-asymmetric `ASSERT_*` before a later collective, skip only a subset of
ranks when collectives remain, or apply different GoogleTest filters on
different ranks. Prefer nonfatal rank-local expectations before collective
operations, and reduce collective preconditions before deciding to skip.
Non-root ranks are quiet by default; set `ERF_GTEST_VERBOSE_RANKS=1` while
debugging if rank-local output is needed.

GoogleTest death and exit tests are prohibited because they are fragile with
MPI, AMReX, and accelerator runtimes. Invalid-domain behavior should use an
explicit status/predicate test or a CMake/CTest subprocess wrapper. The MPI
failure-listener infrastructure test intentionally fails on rank 1 and checks
the wrapper output; it is not a death test.

To rerun one serial case:

```bash
build/Tests/Unit/erf_unit_tests \
  '--gtest_filter=RepresentativeStates/ERFEOSRepresentativeStateTest.*/*DryWarm'
```

To rerun one MPI case at two ranks, select its CTest name:

```bash
ctest --test-dir build -R 'erf_parallel_tests\.SatAdjParallel\.DistributedMultiBoxMatchesScalarReference\.np2' -VV
```

The whole-binary order-dependence check is available with the `stress` label:

```bash
build/Tests/Unit/erf_unit_tests --gtest_shuffle --gtest_repeat=3
ctest --test-dir build -L stress -VV
```

The equivalent whole-binary MPI stress run is manual or scheduled rather than
part of PR CI:

```bash
mpiexec -n 2 build/Tests/Unit/erf_parallel_tests \
  --gtest_shuffle --gtest_repeat=3
```

GoogleTest XML files are written below `build/Tests/Unit/test-results/gtest`.
MPI ranks add deterministic `.npN.rankR.xml` suffixes so concurrent ranks do
not overwrite one another. CI uploads these reports and the CTest temporary
logs even when a test fails. When a shuffled run fails, preserve the reported
`random seed` and rerun with `--gtest_shuffle --gtest_random_seed=<seed>`.
