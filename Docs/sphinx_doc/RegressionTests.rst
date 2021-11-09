

Regression Tests
================

There are currently 10 tests which are run as part of every PR and also as part
of a nightly regression suite.  The CI tests use cmake and are based on the version
of AMReX in the ERF submodule; the nightly tests use GNUMake and use the current
development branch of AMReX.

The following problems are currently tested:

Couette Flow
------------

Isentropic Vortex
-----------------

This tests advection of an isentropic vortex tith triply periodic boundaries.

Advecting Isentropic Vortex
---------------------------

This tests advection of an isentropic vortex tith triply periodic boundaries.

Ekman Spiral
---------------------------

This tests the Coriolis and geostrophic forcing.

Scalar Advection/Diffusion by Uniform Flow
------------------------------------------------

This tests scalar advection and diffusion with triply periodic boundaries.

Scalar Diffusion
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Scalar Diffusion: SineL
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Scalar Advection Sheared Rigid Rotation
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Scalar Advection by Uniform Flow: U
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Scalar Advection by Uniform Flow: UV
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Scalar Advection by Sheared Flow: U
------------------------------------------------

This tests scalar advection with triply periodic boundaries.

Taylor Green Vortex
------------------------------------------------

This tests advection and diffusion with triply periodic boundaries.

Advecting Taylor Green Vortex
------------------------------------------------

This tests advection and diffusion with triply periodic boundaries.

