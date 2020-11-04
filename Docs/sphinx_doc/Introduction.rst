.. highlight:: rst

.. Warning:: This documentation is a work in progress. It is reasonably complete and hopefully useful but is likely error prone and in places misleading.


Introduction
============

ERF solves the compressible Navier-Stokes on a Arakawa C-grid for large-scale weather modeling.

Dependencies
------------

ERF is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), 
an adaptive mesh refinement software framework, which provides the underlying software infrastructure for 
block structured AMR operations. 
The full AMReX documentation can be found `here <https://amrex-codes.github.io/AMReXUsersGuide.pdf>`_. 


Development
-----------

A separate developers guide does not yet exist; along with the algorithmic description in this Users' Guide, doxygen documentation exists in place and an input file exists in `ERF/Docs` that can be build using:

::

	doxygen Doxyfile
