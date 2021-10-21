:orphan:

ERF
---

ERF solves the compressible Navier-Stokes on a Arakawa C-grid for large-scale weather modeling.

ERF is built on `AMReX <https://github.com/AMReX-Codes/amrex>`_,
an adaptive mesh refinement software framework, which provides the underlying software infrastructure for
block structured AMR operations.
The full AMReX documentation can be found `here <https://amrex-codes.github.io/amrex/docs_html/>`_ and the tutorials can be found `here <https://amrex-codes.github.io/amrex/tutorials_html/>`_.

ERF is designed to run on machines from laptops to multicore CPU and hybrid CPU/GPU systems.

For details on the equations that ERF solves, see the :ref:`theory section <theory>`.

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   div#theory.section,
   div#implementation.section,
   div#goals.section,
   </style>

Contents:

.. toctree::
   :maxdepth: 2

   Introduction.rst

Theory
------
.. toctree::
   :caption: THEORY
   :maxdepth: 1
   :hidden:

   theory/Euler_Equations.rst
   theory/UnitsAndConstants.rst
   theory/Algorithms.rst
   theory/ArakawaCGrid.rst

Implementation
--------------
.. toctree::
   :caption: IMPLEMENTATION
   :maxdepth: 1
   :hidden:

   Euler_Discretization.rst
   GettingStarted.rst
   BoundaryConditions.rst
   Visualization.rst

Goals
-----
.. toctree::
   :caption: GOALS
   :maxdepth: 1
   :hidden:

   Applications_Requirements.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`search`


