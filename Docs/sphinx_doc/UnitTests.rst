.. _UnitTests:

Unit Testing
============

ERF uses GoogleTest for focused tests of equations, algorithms, interfaces,
and parallel behavior. CMake discovers and registers the cases with CTest.
CTest selects and runs the registered cases.

Unit tests complement regression tests. A unit test checks a small contract
and reports the failed condition. A regression test runs an ERF problem and
compares its output with a reference result. Use both when a change affects
both a local rule and a complete model configuration. See
:ref:`RegressionTests` for the regression-test workflow.

This guide explains how to build, run, debug, and add unit tests.

.. contents:: On this page
   :local:
   :depth: 2

Agent and contributor quick path
--------------------------------

Use this sequence for each test change.

1. Classify the change as serial, SHOC, device-path, MPI, or regression.
2. Inspect the nearest existing test and its common header.
3. State one observable contract and choose an independent oracle.
4. Add the test to the correct existing source file, or add a new source file
   to the correct CMake target.
5. Build the narrowest configuration that exercises the change.
6. Verify discovery with ``ctest -N``.
7. Run the exact test, then its full CTest label.
8. Run each affected precision, feature, backend, and MPI rank configuration.
9. Record changed files, commands, results, and configurations not run.

Do not widen a tolerance, weaken an assertion, or skip a hardware path merely
to obtain a passing result. First identify whether the failure is in the
implementation, reference, test design, build configuration, or test harness.

Choose the right test
---------------------

Use the smallest test that proves the contract.

Use this decision path.

1. Does the change require a complete ERF input case or output comparison?

   - **Yes:** add or update a regression test.
   - **No:** continue.

2. Does the behavior execute inside an AMReX kernel?

   - **Yes:** add a device-path unit test. Assert the result on the host.
   - **No:** continue.

3. Does the behavior require communication or more than one MPI rank?

   - **Yes:** add a parallel unit test.
   - **No:** add a serial unit test.

SHOC is a domain-specific source location and test binary, not a competing
test category. A SHOC case may still exercise scalar, device, or public-flow
behavior.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Test type
     - Use it for
     - Typical result
   * - Serial unit test
     - A scalar function, local update, public class flow, or host-side rule.
     - A case in ``erf_unit_tests``.
   * - Device-path unit test
     - Code that must compile and execute inside an AMReX kernel.
     - A host assertion on values computed by a device kernel.
   * - Parallel unit test
     - MPI ownership, reductions, ghost exchange, decomposition independence,
       or rank-sensitive behavior.
     - One CTest case per GoogleTest case and MPI rank count.
   * - Infrastructure test
     - The test harness itself.
     - A case with the ``infrastructure`` label.
   * - Stress test
     - Test-order dependence or intermittent shared-state failures.
     - A whole-binary shuffle/repeat run with the ``stress`` label.
   * - Regression test
     - A complete ERF problem, output field, or coupled behavior.
     - A case with the ``regression`` label.

Prerequisites
-------------

ERF unit tests require:

- CMake 3.20 or newer;
- a C++17 compiler;
- the ERF submodules;
- MPI for the parallel suite;
- the libraries and runtime required by the selected ERF backend.

Initialize the submodules after cloning:

.. code-block:: bash

   git submodule update --init --recursive

The GoogleTest source is vendored in ``Submodules/googletest``. ERF builds it
as part of the test build. Do not install another GoogleTest package for this
workflow.

Build serial unit tests
-----------------------

For a unit-only CPU build:

.. code-block:: bash

   cmake -S . -B build-unit \
     -DCMAKE_BUILD_TYPE=Debug \
     -DERF_ENABLE_UNIT_TESTS=ON \
     -DERF_ENABLE_MPI=OFF

   cmake --build build-unit --parallel

   ctest --test-dir build-unit \
     -L unit \
     --output-on-failure \
     --no-tests=error

``ERF_ENABLE_UNIT_TESTS=ON`` builds the GoogleTest targets without requiring
the regression suite. A Debug build gives the clearest local diagnostics.
Use the same compiler and ERF feature options that matter to the code under
test.

Build unit and regression tests
-------------------------------

Use ``ERF_ENABLE_TESTS=ON`` when one build should contain both suites:

.. code-block:: bash

   cmake -S . -B build-tests \
     -DCMAKE_BUILD_TYPE=Debug \
     -DERF_ENABLE_TESTS=ON \
     -DERF_ENABLE_UNIT_TESTS=ON \
     -DERF_ENABLE_REGRESSION_TESTS_ONLY=OFF \
     -DERF_ENABLE_MPI=ON \
     '-DERF_PARALLEL_TEST_NRANKS=1;2;4'

   cmake --build build-tests --parallel

This configuration registers serial unit tests, parallel unit tests at one,
two, and four ranks, infrastructure tests, stress tests, and regression tests.

Use a smaller rank list for a quick local build:

.. code-block:: bash

   cmake -S . -B build-tests \
     -DCMAKE_BUILD_TYPE=Debug \
     -DERF_ENABLE_UNIT_TESTS=ON \
     -DERF_ENABLE_MPI=ON \
     '-DERF_PARALLEL_TEST_NRANKS=1;2'

   cmake --build build-tests --parallel

CMake caches the rank list. Reconfigure the build directory after changing it.

Windows and multi-config generators
-----------------------------------

Visual Studio and other multi-config generators choose the configuration at
build and test time:

.. code-block:: powershell

   cmake -S . -B build-unit `
     -DERF_ENABLE_UNIT_TESTS=ON `
     -DERF_ENABLE_MPI=OFF

   cmake --build build-unit --config Debug --parallel

   ctest --test-dir build-unit `
     -C Debug `
     -L unit `
     --output-on-failure `
     --no-tests=error

Do not assume that a direct test-binary path is the same on a single-config
and a multi-config build. Prefer CTest when writing portable instructions.

Cross-compiling
---------------

ERF discovers GoogleTest cases by running each compiled test binary with
``--gtest_list_tests`` after the binary is built. The listing path does not
initialize MPI, AMReX, or the accelerator runtime, but the build host must
still be able to execute the target binary.

For a true cross-compile, provide a working
``CMAKE_CROSSCOMPILING_EMULATOR`` through the toolchain or perform the test
build on the target system. Do not treat a cross-compiled binary as a
validated test run until it has executed on the target hardware.

When ``add_test`` names an executable target directly, CMake preserves the
target's configuration-specific path and cross-compiling emulator. A wrapper
that launches a test binary itself must propagate ``CROSSCOMPILING_EMULATOR``
explicitly. Do not replace a target-name command with a raw executable path
without preserving that behavior.

Build for a GPU backend
-----------------------

Use the normal ERF backend and machine options. For example, enable CUDA, HIP,
or SYCL in the same build that enables unit tests. Build and run on a machine
that provides the selected compiler, libraries, and device runtime.

A successful GPU compilation proves that the test and production code compile
for that backend. It does not prove runtime behavior. Run the tests on a node
with the required device.

Test configuration options
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Option
     - Effect
   * - ``ERF_ENABLE_TESTS``
     - Enables the regression suite. In a fresh build directory, it also
       defaults ``ERF_ENABLE_UNIT_TESTS`` to ``ON`` unless regression-only mode
       is active.
   * - ``ERF_ENABLE_UNIT_TESTS``
     - Enables the GoogleTest targets and unit-test registration. Set it
       explicitly when reusing a build directory whose cache may contain
       ``OFF``.
   * - ``ERF_ENABLE_REGRESSION_TESTS_ONLY``
     - When used with ``ERF_ENABLE_TESTS=ON``, forces
       ``ERF_ENABLE_UNIT_TESTS=OFF`` so only regression tests are configured.
   * - ``ERF_ENABLE_MPI``
     - Builds the parallel GoogleTest target and enables MPI-aware ERF code.
   * - ``ERF_PARALLEL_TEST_NRANKS``
     - Sets the semicolon-separated MPI rank counts used to register each
       parallel case.
   * - ``ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE``
     - Adds the supported Open MPI or PRRTE oversubscription flag when needed.
   * - ``ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS``
     - Adds launcher-specific flags before the test executable.
   * - ``ERF_PRECISION``
     - Selects ``DOUBLE`` or ``SINGLE`` mesh precision.
   * - ``ERF_ENABLE_PARTICLES``
     - Enables particle code and tests that depend on it.
   * - ``ERF_PARTICLES_PRECISION``
     - Selects particle precision when particles are enabled.

List the registered tests
-------------------------

List tests before running a narrow selection:

.. code-block:: bash

   ctest --test-dir build-tests -N -L unit
   ctest --test-dir build-tests -N -L shoc
   ctest --test-dir build-tests -N -L parallel
   ctest --test-dir build-tests -N -L infrastructure
   ctest --test-dir build-tests -N -L stress
   ctest --test-dir build-tests -N -L regression

``-N`` lists tests without running them. Add ``-V`` to show each test command:

.. code-block:: bash

   ctest --test-dir build-tests -N -V -L parallel

This command is useful when checking an MPI launcher, rank count, environment
variable, or generated filter.

Run test groups
---------------

Run serial and SHOC unit tests:

.. code-block:: bash

   ctest --test-dir build-tests \
     -L unit \
     --output-on-failure \
     --no-tests=error

Run only SHOC tests:

.. code-block:: bash

   ctest --test-dir build-tests \
     -L shoc \
     --output-on-failure \
     --no-tests=error

Run parallel tests:

.. code-block:: bash

   ctest --test-dir build-tests \
     -L parallel \
     --output-on-failure \
     --no-tests=error

Run the test-harness self-test:

.. code-block:: bash

   ctest --test-dir build-tests \
     -L infrastructure \
     --output-on-failure \
     --no-tests=error

Run the serial order-dependence check:

.. code-block:: bash

   ctest --test-dir build-tests \
     -L stress \
     --output-on-failure \
     --no-tests=error

Use ``--no-tests=error`` in scripts. A misspelled label or filter should fail
the command, not produce a false pass.

Run one serial test
-------------------

Use the CTest name when possible:

.. code-block:: bash

   ctest --test-dir build-tests \
     -R '^ERFEOSConstants\.KappaGammaContract$' \
     -VV \
     --no-tests=error

You may also run the serial binary directly on a single-config build:

.. code-block:: bash

   ./build-tests/Tests/Unit/erf_unit_tests \
     --gtest_filter='ERFEOSConstants.KappaGammaContract'

List the binary's GoogleTest cases with:

.. code-block:: bash

   ./build-tests/Tests/Unit/erf_unit_tests --gtest_list_tests

The SHOC tests use a separate binary:

.. code-block:: bash

   ./build-tests/Tests/Unit/Shoc/erf_shoc_unit_tests --gtest_list_tests

Direct paths differ on multi-config builds. CTest avoids that difference.

Run one parameterized test
--------------------------

Parameterized tests include their instance and parameter names. List the tests
first, then copy the exact name or use a narrow GoogleTest filter.

For example:

.. code-block:: bash

   ./build-tests/Tests/Unit/erf_unit_tests \
     --gtest_filter='RepresentativeStates/ERFEOSRepresentativeStateTest.*/*DryWarm'

Give parameters stable names such as ``DryWarm``, ``ColdUpperAir``, or
``Upwind3Negative``. A name should describe the physical or numerical case,
not its position in an array.

Run one MPI test
----------------

MPI CTest names include the suite, case, and rank count:

.. code-block:: text

   erf_parallel_tests.<suite>.<case>.np<ranks>

Run one SatAdj case at two ranks:

.. code-block:: bash

   ctest --test-dir build-tests \
     -R '^erf_parallel_tests\.SatAdjParallel\.DistributedMultiBoxMatchesScalarReference\.np2$' \
     -VV \
     --no-tests=error

CTest supplies the launcher, rank count, filter, XML path, and processor count.
Use the CTest command instead of reconstructing the launcher by hand.

For rank-local output during debugging:

.. code-block:: bash

   ERF_GTEST_VERBOSE_RANKS=1 \
   ctest --test-dir build-tests \
     -R '^erf_parallel_tests\..*\.np2$' \
     -VV \
     --no-tests=error

On PowerShell:

.. code-block:: powershell

   $env:ERF_GTEST_VERBOSE_RANKS = "1"
   ctest --test-dir build-tests -R '\.np2$' -VV --no-tests=error
   Remove-Item Env:ERF_GTEST_VERBOSE_RANKS

Run shuffle and repeat checks
-----------------------------

The registered serial stress test shuffles the complete binary and repeats it
three times:

.. code-block:: bash

   ctest --test-dir build-tests -L stress -VV --no-tests=error

You may run the binary directly:

.. code-block:: bash

   ./build-tests/Tests/Unit/erf_unit_tests \
     --gtest_shuffle \
     --gtest_repeat=3

GoogleTest prints the random seed. Preserve it when reporting a failure. Repeat
the same order with:

.. code-block:: bash

   ./build-tests/Tests/Unit/erf_unit_tests \
     --gtest_shuffle \
     --gtest_random_seed=<seed>

The whole-binary MPI stress run is not registered in pull-request CI. Run it
manually, or add it to a machine-specific scheduled workflow:

.. code-block:: bash

   mpiexec -n 2 ./build-tests/Tests/Unit/erf_parallel_tests \
     --gtest_shuffle \
     --gtest_repeat=3

Use the launcher and flags required by the machine.

Find test output
----------------

GoogleTest XML files appear below:

.. code-block:: text

   build-tests/test-results/gtest

MPI report names include the rank count and rank so concurrent processes do
not overwrite one another.

CTest writes its latest logs below:

.. code-block:: text

   build-tests/Testing/Temporary

Start with:

.. code-block:: text

   LastTest.log
   LastTestsFailed.log

The main ``ERF CI`` workflow uploads the GoogleTest report tree and the CTest
temporary logs even when a test step fails.

Test layout
-----------

The unit-test tree follows the production domains:

.. code-block:: text

   Tests/Unit/
     BoundaryConditions/
     Diagnostics/
     IO/
     Microphysics/
     Shoc/
     Utils/

The main targets are:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - Target
     - Main label
     - Purpose
   * - ``erf_unit_tests``
     - ``unit``
     - Serial, scalar, public-flow, MultiFab, and device-path tests.
   * - ``erf_shoc_unit_tests``
     - ``unit`` and ``shoc``
     - SHOC-specific tests in a separate binary.
   * - ``erf_parallel_tests``
     - ``parallel``
     - MPI tests registered once per case and rank count.
   * - ``erf_mpi_gtest_listener_self_test``
     - ``infrastructure``
     - Verifies that a non-root failure reaches the root-rank report.

Common headers near a test family hold fixtures, representative states,
reference calculations, tolerance rules, and test-only helpers. Keep helpers
close to the tests that define their contract.

Add a serial test
-----------------

Follow these steps.

1. Find the nearest test family and read its common header.
2. State the contract in one sentence.
3. Choose ``TEST``, ``TEST_P``, or a loop.
4. Write a deterministic case with the smallest useful state.
5. Add the source file to the correct CMake target if it is new.
6. Rebuild the target.
7. Confirm that GoogleTest and CTest discover the case.
8. Run the narrow case, then its label.
9. Run any affected precision, particle, MPI, or backend configuration.

A new source file for the main serial binary belongs in the
``target_sources(${erf_exe_name} ...)`` list in
``Tests/Unit/CMakeLists.txt``. A new ``TEST`` in an existing source file needs
no CMake edit.

Use a name that states behavior:

.. code-block:: cpp

   TEST(MyFeatureScalar, ConservesTotalQuantity)
   {
       const amrex::Real total_before = compute_total(state);

       advance_one_step(state);

       const amrex::Real total_after = compute_total(state);
       const amrex::Real tolerance = conservation_tolerance(total_before);
       ERF_EXPECT_NEAR(total_after, total_before, tolerance);
   }

Place a short motivation comment above the test when the contract is not
obvious:

.. code-block:: cpp

   // Motivation: The local source terms only exchange mass among species.
   // Their sum must remain unchanged after one update.
   TEST(MyFeatureScalar, LocalSourcesConserveTotalMass)
   {
       // ...
   }

Do not write ``Works``, ``Basic``, or ``Smoke`` when a precise behavior name
is available.

Add a SHOC test
---------------

Place SHOC tests in ``Tests/Unit/Shoc``. Add a new ``TEST`` to an existing
source file when it belongs to that file's contract. Add a new source file to
the ``erf_shoc_unit_tests`` source list in
``Tests/Unit/Shoc/CMakeLists.txt``.

Run the narrow case first, then the SHOC label:

.. code-block:: bash

   ctest --test-dir build-tests \
     -R '<exact-shoc-test-name>' \
     -VV \
     --no-tests=error

   ctest --test-dir build-tests \
     -L shoc \
     --output-on-failure \
     --no-tests=error

Keep reusable SHOC fixtures under ``Tests/Unit/Shoc/fixtures``. Use
``ERF_SHOC_UNIT_FIXTURE_DIR`` when a test must locate a source fixture. Do not
hard-code an absolute source-tree path.

Choose TEST, TEST_P, or a loop
------------------------------

Use ``TEST`` for one cohesive property.

Use ``TEST_P`` when each case has a meaningful identity and an independent
rerun is useful. Provide stable parameter names:

.. code-block:: cpp

   using StateCase = std::tuple<std::string, State>;

   class RepresentativeStateTest
       : public ::testing::TestWithParam<StateCase> {};

   TEST_P(RepresentativeStateTest, RoundTripRecoversInput)
   {
       const auto& name = std::get<0>(GetParam());
       const auto& state = std::get<1>(GetParam());
       SCOPED_TRACE(name);

       const auto recovered = round_trip(state);
       ERF_EXPECT_NEAR(recovered, state.value, tolerance_for(state));
   }

   INSTANTIATE_TEST_SUITE_P(
       RepresentativeStates,
       RepresentativeStateTest,
       ::testing::Values(
           StateCase{"DryWarm", make_dry_warm_state()},
           StateCase{"ColdUpperAir", make_cold_upper_air_state()}),
       [](const ::testing::TestParamInfo<StateCase>& info) {
           return std::get<0>(info.param);
       });

Use an internal loop for a dense mathematical sweep whose iterations together
prove one property. Add ``SCOPED_TRACE`` so a failure identifies the case:

.. code-block:: cpp

   for (int degree = 0; degree <= maximum_degree; ++degree) {
       SCOPED_TRACE("degree=" + std::to_string(degree));
       // Evaluate and assert the property.
   }

Do not parameterize hundreds of cells or quadrature points merely to create
more CTest entries.

Write strong contracts
----------------------

Prefer tests that enforce a scientific or numerical rule:

- conservation;
- positivity or boundedness;
- inverse or round-trip consistency;
- monotonicity;
- symmetry;
- branch and threshold behavior;
- agreement with an independently derived reference;
- decomposition independence;
- host/device agreement;
- a public-flow regression for a confirmed bug.

A sign-only assertion is weak when the exact relation is known. A frozen
number is weak when a conservation law or round trip gives a better oracle.

Do not duplicate the production implementation and call the copy an
independent reference. When a reference must resemble production logic,
state how it was derived and which implementation choices it intentionally
does not share.

Use numerical assertions
------------------------

Include the test-only assertion header with the relative path used by the
neighboring tests.

ERF provides:

.. code-block:: cpp

   ERF_EXPECT_NEAR(actual, expected, tolerance);
   ERF_ASSERT_NEAR(actual, expected, tolerance);
   ERF_EXPECT_FINITE(value);
   ERF_EXPECT_NONNEGATIVE(value, tolerance);

``ERF_EXPECT_NEAR`` and ``ERF_ASSERT_NEAR`` pass when the absolute error is no
larger than the supplied tolerance. A failure also reports the relative error
and precision mode. The relative error is diagnostic; it does not replace the
supplied tolerance.

Use standard GoogleTest assertions for exact, integer, Boolean, and ordering
contracts:

.. code-block:: cpp

   EXPECT_EQ(actual_count, expected_count);
   EXPECT_TRUE(is_valid);
   EXPECT_LT(lower, upper);
   ASSERT_NE(pointer, nullptr);

Use ``EXPECT_*`` when later checks remain safe and useful. Use ``ASSERT_*``
when the test cannot continue safely. In an MPI test, a rank-local fatal
assertion before a later collective can deadlock the other ranks.

Justify every tolerance
-----------------------

Put tolerance provenance next to the contract. Base it on one or more of:

- machine precision;
- single- versus double-precision behavior;
- the number and order of floating-point operations;
- accumulation length;
- a documented approximation error;
- backend math variation;
- finite-difference truncation error;
- a measured discontinuity at an intentional branch.

Use an absolute tolerance near zero. Use a scaled absolute tolerance when the
magnitude varies. Keep single- and double-precision choices explicit.

Do not widen a tolerance until CI passes. First determine whether the
difference follows from roundoff, a backend implementation, an invalid
reference, nondeterminism, or a real defect.

Print test structures clearly
-----------------------------

Add a test-only ``PrintTo`` overload when a structure appears in a
parameterized test or failure message. Place it in the same namespace as the
type so GoogleTest finds it through argument-dependent lookup.

Use ``erf_gtest::print_field`` to keep field formatting consistent:

.. code-block:: cpp

   inline void PrintTo (const CellState& state, std::ostream* os)
   {
       *os << "{";
       erf_gtest::print_field(os, "temperature", state.temperature, true);
       erf_gtest::print_field(os, "pressure", state.pressure);
       erf_gtest::print_field(os, "qv", state.qv);
       *os << "}";
   }

Do not add a production ``operator<<`` solely for a test.

Test device code
----------------

GoogleTest assertions are host code. Never call ``EXPECT_*``, ``ASSERT_*``,
``GTEST_SKIP()``, stream output, or host-only formatting from an AMReX device
lambda.

A device-path test should:

1. allocate or initialize the smallest useful AMReX object;
2. run the production function inside the intended AMReX kernel;
3. compute values, counters, or normalized errors on the device;
4. synchronize the device work;
5. copy or reduce the result to the host;
6. assert the contract on the host.

Use ``amrex::Real`` and the production host/device annotations. Keep scalar
and kernel tests separate when they require different setup. A scalar test
proves the formula. A kernel test proves that the production path compiles and
executes in the backend.

Do not assume identical last-bit results across CPU, CUDA, HIP, and SYCL.
Assert the contract with a justified tolerance.

Add an MPI test
---------------

Add a new MPI source file to the ``erf_parallel_tests`` target in
``Tests/Unit/CMakeLists.txt``. Each ``TEST`` in that binary is registered at
every rank count in ``ERF_PARALLEL_TEST_NRANKS``.

All ranks run the same GoogleTest filter. Keep their collective sequence
identical.

A small reduction test may follow this pattern:

.. code-block:: cpp

   TEST(MyFeatureParallel, GlobalSumIncludesEveryRank)
   {
       const int rank = amrex::ParallelDescriptor::MyProc();
       const int nprocs = amrex::ParallelDescriptor::NProcs();

       int global_sum = rank + 1;
       amrex::ParallelDescriptor::ReduceIntSum(global_sum);

       const int expected = nprocs * (nprocs + 1) / 2;
       EXPECT_EQ(global_sum, expected);
   }

Do not:

- use a rank-asymmetric ``ASSERT_*`` before a later collective;
- call ``GTEST_SKIP()`` on only some ranks while other ranks continue;
- return early on one rank before a collective;
- use different filters on different ranks;
- rely on root-only state without broadcasting or reducing it;
- assume a particular box-to-rank assignment unless that assignment is the
  contract under test.

The MPI failure listener gathers non-root assertion messages. It improves
diagnostics. It cannot repair a deadlock caused by divergent control flow.

When a precondition may differ by rank, reduce it first. Then make the same
decision on every rank.

Use fixtures with restraint
---------------------------

Use ``TEST_F`` when shared setup makes the contract clearer. Keep fixtures
small and deterministic. Reset mutable state for every case.

CTest normally runs each discovered GoogleTest case in a separate process.
``SetUpTestSuite()`` therefore does not amortize expensive setup across normal
CTest cases. Do not build a large fixture to gain performance that the CTest
execution model cannot provide.

Avoid hidden global state and test-order dependence. Use the ``stress`` test
to expose order-sensitive failures.

Test invalid input without death tests
--------------------------------------

Do not use:

.. code-block:: cpp

   EXPECT_DEATH
   ASSERT_DEATH
   EXPECT_EXIT
   ASSERT_EXIT

Death and exit tests are fragile with MPI, AMReX initialization, accelerator
runtimes, and platform launchers.

Prefer one of these designs:

- expose and test a status or predicate;
- test a validation helper directly;
- test a returned error value;
- use a CMake or CTest subprocess wrapper when process termination is the
  contract.

The MPI listener self-test uses a CTest wrapper around an intentionally failing
binary. It is an infrastructure test, not a model for ordinary numerical
tests.

Keep tests portable
-------------------

A new test must preserve ERF's supported build configurations.

- Use standard C++17.
- Use ``amrex::Real`` for ERF real-valued state.
- Respect single and double precision.
- Keep host-only code out of device paths.
- Synchronize before the host reads device results.
- Avoid backend-specific math assumptions unless the test is backend-specific.
- Do not require MPI in a serial target.
- Do not require particles when ``ERF_ENABLE_PARTICLES=OFF``.
- Guard feature-specific tests in the same way as the production feature.
- Do not assume POSIX paths in C++ test logic.
- Do not assume a single-config build layout.
- Do not depend on wall-clock timing, process scheduling, or test order.
- Seed any random input and print the seed on failure.
- Keep MPI results independent of the tested decomposition unless the
  decomposition itself is the subject.
- Treat successful compilation and successful runtime validation as separate
  claims.

When a required backend or machine is unavailable, state what you did not run.
Do not claim portability from inspection alone.

Coding standards
----------------

Follow the ERF coding guide in ``CONTRIBUTING.md``.

For tests:

- use four spaces, not tabs;
- use braces for every control-flow body;
- keep one clear contract per test;
- name the suite for the subsystem or path;
- name the case for observable behavior;
- place reusable test-only logic in a nearby common header;
- keep production changes separate from test-only cleanup;
- comment the reason for a test, not each line of mechanics;
- keep the state small enough that a failure is easy to inspect;
- prefer public behavior over private implementation details;
- report enough context to reproduce a failure.

Avoid these patterns
--------------------

Do not add:

- death or exit tests;
- unseeded random tests;
- sleeps or wall-clock performance thresholds;
- rank-dependent control flow around MPI collectives;
- GoogleTest assertions inside device code;
- arbitrary golden values when a stronger invariant exists;
- magic tolerances without provenance;
- a copied production formula presented as an independent oracle;
- tests that only print output;
- tests that pass without asserting a contract;
- broad fixtures with hidden mutable state;
- order-dependent tests;
- host-only assumptions in a portable kernel test;
- a production stream operator needed only for diagnostics;
- a new test file that is not added to the correct CMake target;
- edits to generated CTest files.

Debug a failure
---------------

Start with the narrowest failing CTest command:

.. code-block:: bash

   ctest --test-dir build-tests \
     -R '<exact-test-name>' \
     -VV \
     --no-tests=error

Then:

1. read the assertion and ``SCOPED_TRACE`` context;
2. inspect ``build-tests/Testing/Temporary/LastTest.log``;
3. inspect the GoogleTest XML report;
4. rerun the exact GoogleTest filter when direct execution is useful;
5. preserve the shuffle seed for order-dependent failures;
6. enable ``ERF_GTEST_VERBOSE_RANKS=1`` for MPI diagnostics;
7. compare the failing precision, backend, feature flags, and rank count with a
   passing configuration;
8. decide whether the failure is numerical, parallel, configuration-specific,
   or a test-harness error before changing a tolerance.

A hang in an MPI case usually requires inspection of rank control flow and
collective order. Run the case at the smallest failing rank count and enable
rank-local output.

Troubleshooting
---------------

No unit tests are listed
~~~~~~~~~~~~~~~~~~~~~~~~

Check that:

- ``Submodules/googletest`` is initialized;
- ``ERF_ENABLE_UNIT_TESTS=ON`` or ``ERF_ENABLE_TESTS=ON`` was used;
- the new source file is part of the intended target;
- the target was rebuilt after the test was added;
- the selected label or regular expression is correct.

When a reconfigured build omits unit tests, inspect ``ERF_ENABLE_UNIT_TESTS``
in the CMake cache or set it explicitly.

Use ``--no-tests=error`` so an empty selection fails.

No parallel tests are listed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check that:

- ``ERF_ENABLE_MPI=ON``;
- CMake found an MPI launcher;
- ``ERF_PARALLEL_TEST_NRANKS`` contains positive integers;
- ``erf_parallel_tests`` was built;
- ``ctest -N -V -L parallel`` shows valid launcher commands.

The generated registration file is in the build tree. Inspect it when
debugging, but do not edit it.

The MPI launcher rejects local ranks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the scheduler or machine launcher rules. On Open MPI or PRRTE systems,
``ERF_PARALLEL_TEST_ALLOW_OVERSUBSCRIBE=ON`` can add the supported
oversubscription flag. Use ``ERF_PARALLEL_TEST_MPIEXEC_PREFLAGS`` for required
launcher-specific flags.

A GPU test builds but does not run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run on a node with the selected device and runtime. Check the backend, compiler,
driver, visible devices, and scheduler allocation. A login-node build does not
establish runtime support.

A test fails only in single precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check the tolerance derivation and the operation count. Do not replace a
single-precision failure with the double-precision tolerance. Also check for
literal types, unintended promotion, cancellation, and invalid exact-equality
assumptions.

A test fails only at more MPI ranks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check ownership, ghost cells, reductions, collective order, and assumptions
about box distribution. Compare global counts and invariants, not rank-local
layout, unless layout is the contract.

Completion criteria
-------------------

A unit-test change is complete when:

- the test states one observable contract;
- the test contains at least one meaningful assertion;
- its oracle is independent enough to detect the intended defect;
- each numerical tolerance has written provenance where the reason is not
  obvious;
- a new source file is part of the intended CMake target;
- ``ctest -N`` lists the new case;
- the exact test passes;
- its full CTest label passes;
- affected precision, particle, MPI-rank, and backend configurations were run
  or recorded as not run;
- device work is synchronized before host assertions read its result;
- MPI control flow remains collective-safe;
- ``git diff --check`` passes;
- the final report lists changed files, commands, results, and untested
  configurations.

For a bug fix, verify that the new test fails for the intended reason before
the fix when practical. After the fix, verify that it passes and that nearby
tests remain green.

For generated code or agent-written tests, review the contract by hand. A test
that compiles is not necessarily a useful test. A test that passes is not
necessarily correct.

Related documentation
---------------------

- :ref:`GettingStartedWithERFTests` gives a first CPU-only test build and run.
- :ref:`RegressionTests` covers complete ERF regression cases.
- :doc:`buildingConfiguration` covers machine and library configuration.
- The `GoogleTest primer <https://google.github.io/googletest/primer.html>`_
  explains the base assertion and fixture APIs.
- The `GoogleTest advanced guide <https://google.github.io/googletest/advanced.html>`_
  covers parameterized tests, printers, listeners, and filters.
- The `CTest manual <https://cmake.org/cmake/help/latest/manual/ctest.1.html>`_
  documents labels, regular expressions, parallel execution, and test output.
