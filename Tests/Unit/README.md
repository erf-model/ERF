# ERF unit tests

The canonical guide is the Sphinx
[Unit Testing page](../../Docs/sphinx_doc/UnitTests.rst). The rendered guide is
available in the
[ERF documentation](https://erf.readthedocs.io/en/latest/UnitTests.html).

## Quick start

```bash
git submodule update --init --recursive

cmake -S . -B build-tests \
  -DCMAKE_BUILD_TYPE=Debug \
  -DERF_ENABLE_UNIT_TESTS=ON \
  -DERF_ENABLE_MPI=ON \
  '-DERF_PARALLEL_TEST_NRANKS=1;2'

cmake --build build-tests --parallel

ctest --test-dir build-tests -L unit --output-on-failure --no-tests=error
ctest --test-dir build-tests -L parallel --output-on-failure --no-tests=error
ctest --test-dir build-tests -L infrastructure --output-on-failure --no-tests=error
```

## Rules that prevent fragile tests

- Keep GoogleTest assertions on the host.
- Synchronize device work before reading results on the host.
- Keep every MPI rank on the same collective-safe control-flow path.
- Justify numerical tolerances.
- Use named parameterized cases only when independent filtering helps.
- Use `SCOPED_TRACE` for dense internal loops.
- Do not use GoogleTest death or exit tests.
- Do not edit generated CTest files.

Read the Sphinx guide before adding a test. It includes build variants,
test labels, source registration, assertion helpers, parameterized tests,
device and MPI patterns, portability rules, troubleshooting, and a checklist
for contributors and coding agents.
