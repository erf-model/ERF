
.. role:: cpp(code)
   :language: c++

.. _sec:FireFramework:

Fire Model - Framework Implementation
======================================

Overview
--------

The fire model framework establishes the foundation for wildfire propagation simulation using the Rothermel fire spread model with FARSITE-based elliptical expansion. This section documents the core implementation of fire layer structure, integration with ERF, and computational kernels.

Core Implementation
-------------------

Fire Model Source Code
~~~~~~~~~~~~~~~~~~~~~~

**Location:** ``Source/Fire/``

**Files:** ``ERF_Fire.H``, ``ERF_Fire.cpp``, ``Make.package``

**Status:** Core implementation with all required variables and computational functions

Fire Model Variables
~~~~~~~~~~~~~~~~~~~~

The implementation includes storage for:

- **Rothermel Parameters:** Fuel moisture, bed depth, particle density, energy content, fuel load, wind speed, slope
- **FARSITE Variables:** Ellipse ratio, eccentricity, major/minor axes
- **Fire State:** Head/flank/back spread rates, fire line intensity, flame length

Core Computational Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All core functions are implemented with full or placeholder logic:

- ``Define()`` - Initialize parameters
- ``Init()`` - Initialize fire variables
- ``Advance()`` - Main time stepping
- ``Update_Fire_Vars()`` - Update fire state
- ``ComputeRothermellSpreadRate()`` - Rothermel calculation
- ``ComputeEllipticalExpansion()`` - Ellipse calculation
- ``ComputeFireIntensity()`` - Intensity calculation

Integration with ERF
~~~~~~~~~~~~~~~~~~~~

- Fire model header included in main ERF class
- Fire member variable added to ERF
- Initialization function implemented and called in constructor
- Proper integration with build system

Test Infrastructure
~~~~~~~~~~~~~~~~~~~

**Location:** ``Exec/CanonicalTests/Fire/``

**Files:** ``inputs_fire_dummy``, ``test_fire_dummy.py``, ``README.md``

**Status:** Complete with executable test script

Documentation
~~~~~~~~~~~~~~

- **Main Documentation:** ``Docs/Fire_Model_Documentation.md``
- **Verification Checklist:** ``PHASE1_CHECKLIST.md``
- **Test Documentation:** ``Exec/CanonicalTests/Fire/README.md``

Running Framework Tests
-----------------------

To run the framework regression test:

.. code-block:: bash

   cd Exec/CanonicalTests/Fire/
   python3 test_fire_dummy.py --erf_exe ./erf --input_file inputs_fire_dummy

Key Files Reference
--------------------

.. list-table::
   :widths: 30 50
   :header-rows: 1

   * - File
     - Purpose
   * - ``Source/Fire/ERF_Fire.H``
     - Fire model class definition
   * - ``Source/Fire/ERF_Fire.cpp``
     - Fire model implementation
   * - ``Source/Fire/Make.package``
     - Build configuration
   * - ``Source/ERF.H``
     - Modified to include Fire
   * - ``Source/ERF.cpp``
     - Fire initialization code
   * - ``Source/ERF_Constructors.cpp``
     - Fire initialization call
   * - ``Exec/CanonicalTests/Fire/inputs_fire_dummy``
     - Test input file
   * - ``Exec/CanonicalTests/Fire/test_fire_dummy.py``
     - Test script
   * - ``Docs/Fire_Model_Documentation.md``
     - Complete documentation
   * - ``PHASE1_CHECKLIST.md``
     - Verification checklist

Implementation Status
---------------------

Framework implementation is complete with:

- Fire layer structure with all required variables
- Computational function calls
- Regression test infrastructure
- Documentation
- Proper integration with ERF
- Build system configuration

Subsequent development extends the computational kernels for full physics representation.

Next Steps
----------

Further development focuses on implementing the complete Rothermel model equations:

- Wind factor calculations
- Slope factor calculations
- Reaction intensity computations
- Heat absorption factors
- Propagating flux ratios

References
----------

- Rothermel, R. C. (1972). A mathematical model for predicting fire spread in wildland fuels. Res. Paper INT-115, USDA Forest Service, Intermountain Forest and Range Experiment Station.
- Finney, M. A. (2004). FARSITE: Fire Area Simulator model development and evaluation. Res. Paper RMRS-RP-4 Revised, USDA Forest Service, Rocky Mountain Research Station.
- Andrews, P. L. (2018). Current status and future needs of the BehavePlus Fire Modeling System. International Journal of Wildland Fire, 27(9), 558-566.

