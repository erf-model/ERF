.. _Testing:
.. _GettingStartedWithERFTests:

Getting Started with ERF Tests
==============================

ERF uses two complementary kinds of tests.

A unit test checks one focused behavior, such as an equation, numerical
property, interface, or parallel rule. A regression test runs a complete ERF
case and compares its output with a reference result.

This page follows one path: build and run the CPU-only serial unit tests in a
Debug configuration. You do not need MPI or a GPU. Use :ref:`UnitTests` when
you are ready to add tests or work with SHOC, MPI, device code, accelerator
backends, multiple precisions, or cross-compiling. Use :ref:`RegressionTests`
when a change affects a complete ERF problem or its output.

What you need
-------------

Run these commands from the root of an ERF checkout. You need:

- CMake 3.20 or newer;
- a C++17 compiler;
- the ERF source tree;
- the repository submodules.

Initialize the submodules after cloning:

.. code-block:: bash

   git submodule update --init --recursive

Configure and build
-------------------

Configure a Debug build with unit tests enabled and MPI disabled:

.. code-block:: bash

   cmake -S . -B build-unit \
     -DCMAKE_BUILD_TYPE=Debug \
     -DERF_ENABLE_UNIT_TESTS=ON \
     -DERF_ENABLE_MPI=OFF

   cmake --build build-unit --parallel

The ``-B build-unit`` argument keeps generated files outside the source tree.
Remove that directory when you need a completely clean configuration.

List the tests
--------------

CMake registers tests in the build directory. CTest is the program that lists,
selects, and runs those registered tests.

List the unit tests without running them:

.. code-block:: bash

   ctest --test-dir build-unit -N -L unit

``-N`` means “show only.” ``-L unit`` selects tests with the ``unit`` label.

Run the unit tests
------------------

Run the unit-test label:

.. code-block:: bash

   ctest --test-dir build-unit \
     -L unit \
     --output-on-failure \
     --no-tests=error

``--output-on-failure`` prints the failed test's output. ``--no-tests=error``
makes an empty selection fail instead of appearing successful.

Run one test
------------

Use ``-R`` with an exact CTest-name regular expression:

.. code-block:: bash

   ctest --test-dir build-unit \
     -R '^ERFEOSConstants\.KappaGammaContract$' \
     -VV \
     --no-tests=error

``-VV`` shows the full test command and output. List the tests first when you
do not know the exact name.

Understand the result
---------------------

A successful run reports that all selected tests passed and returns a zero
exit status.

A failed test returns a nonzero status. Read the assertion message and the
output printed by ``--output-on-failure``.

An empty selection also returns a nonzero status because the commands use
``--no-tests=error``. Check the build options, label, regular expression, and
whether the test target was rebuilt.

CTest writes its latest detailed log to:

.. code-block:: text

   build-unit/Testing/Temporary/LastTest.log

Where to go next
----------------

- **Adding or reviewing a unit test?** Read :ref:`UnitTests`.
- **Need SHOC, MPI, GPU, device-path, precision, or cross-compiling guidance?**
  Read :ref:`UnitTests`.
- **Changing a complete ERF case or reference output?** Read
  :ref:`RegressionTests`.
