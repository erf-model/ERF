.. ERF documentation master file, created by
   sphinx-quickstart on Tue Nov 15 14:07:58 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ERF's documentation!
=================================

.. raw:: html

   <style>
   /* front page: hide chapter titles
    * needed for consistent HTML-PDF-EPUB chapters
    */
   div#theory.section,
   </style>

Contents:

.. toctree::
   :maxdepth: 2

   Introduction.rst
   Applications_Requirements.rst
   Euler_Discretization.rst
   GettingStarted.rst
   BoundaryConditions.rst
   Visualization.rst

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

README
======
.. _README:

.. include:: ../../README.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`search`


