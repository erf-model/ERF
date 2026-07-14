.. _Testing:

Testing and Verification
========================

ERF uses two complementary test layers.

.. list-table::
   :header-rows: 1
   :widths: 24 38 38

   * - Layer
     - Purpose
     - Guide
   * - Unit tests
     - Check focused equations, algorithms, interfaces, device paths, and MPI
       behavior with GoogleTest and CTest.
     - :ref:`UnitTests`
   * - Regression tests
     - Run complete ERF problems and compare their output with reference data.
     - :ref:`RegressionTests`

Enable both layers with:

.. code-block:: bash

   git submodule update --init --recursive

   cmake -S . -B build-tests \
     -DCMAKE_BUILD_TYPE=Debug \
     -DERF_ENABLE_TESTS=ON \
     -DERF_ENABLE_MPI=ON \
     '-DERF_PARALLEL_TEST_NRANKS=1;2'

   cmake --build build-tests --parallel

Run the main groups with:

.. code-block:: bash

   ctest --test-dir build-tests -L unit --output-on-failure --no-tests=error
   ctest --test-dir build-tests -L parallel --output-on-failure --no-tests=error
   ctest --test-dir build-tests -L regression --output-on-failure --no-tests=error

Use ``ctest --test-dir build-tests -N`` to list tests and ``-R`` to select a
test by regular expression.

Read :ref:`UnitTests` before adding or modifying a GoogleTest case. It covers
test selection, numerical tolerances, device code, MPI safety, portability,
debugging, and coding standards.

Read :ref:`RegressionTests` before adding a complete ERF problem or changing a
gold file.
