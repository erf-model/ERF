:orphan:

ERF
---

ERF solves the compressible Navier-Stokes on a Arakawa C-grid for large-scale weather modeling.

ERF is built on `AMReX <https://github.com/AMReX-Codes/amrex>`_,
an adaptive mesh refinement software framework, which provides the underlying software infrastructure for
block structured AMR operations.
The full AMReX documentation can be found `here <https://amrex-codes.github.io/amrex/docs_html/>`_ and the tutorials can be found `here <https://amrex-codes.github.io/amrex/tutorials_html/>`_.

ERF is designed to run on machines from laptops to multicore CPU and hybrid CPU/GPU systems.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   div#userguide.section,
   div#theory.section,
   div#implementation.section,
   div#goals.section,
   </style>

.. toctree::
   :hidden:

   coc

.. toctree::
   :caption: USER GUIDE
   :maxdepth: 1
   :hidden:

   GettingStarted.rst
   Inputs.rst

.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/DryEquations.rst
   theory/WetEquations.rst
   theory/Buoyancy.rst
   theory/Microphysics.rst
   theory/DNSvsLES.rst
   theory/PBLschemes.rst
   theory/Forcings.rst
   theory/UnitsAndConstants.rst

.. toctree::
   :caption: IMPLEMENTATION
   :maxdepth: 1
   :hidden:

   ArakawaCGrid.rst
   MapFactors.rst
   TimeAdvance.rst
   Discretizations.rst
   MeshRefinement.rst
   BoundaryConditions.rst
   Derived.rst
   Checkpoint.rst
   Plotfiles.rst
   Visualization.rst

.. toctree::
   :caption: COUPLING TO AMR-WIND
   :maxdepth: 1
   :hidden:

   CouplingToAMRWind.rst

.. toctree::
   :caption: ERF vs WRF
   :maxdepth: 1
   :hidden:

   ERFvsWRF.rst

.. toctree::
   :caption: TESTING
   :maxdepth: 1
   :hidden:

   RegressionTests.rst

.. toctree::
   :caption: GOALS
   :maxdepth: 1
   :hidden:

   Applications_Requirements.rst
