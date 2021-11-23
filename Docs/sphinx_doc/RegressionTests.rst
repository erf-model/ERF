

Regression Tests
================
Last update: 2021-11-11.

There are currently 12 tests which are run as part of every PR and also as part
of a nightly regression suite.  The CI tests use cmake and are based on the version
of AMReX in the ERF submodule; the nightly tests use GNUMake and use the current
development branch of AMReX.

Results from the nightly CPU tests can be found here: `CPU tests`_

Results from the nightly GPU tests can be found here: `GPU tests`_

.. _`CPU tests`: https://ccse.lbl.gov/pub/RegressionTesting1/ERF

.. _`GPU tests`: https://ccse.lbl.gov/pub/GpuRegressionTesting/ERF

The following problems are currently tested in the CI:

+-------------------------------+----------+-----+-----+-----+-------+----------------+
| Test                          | nx ny nz | xbc | ybc | zbc | Ext   | Other          |
+===============================+==========+=====+=====+=====+=======+================+
| ScalarAdvectionUniformU       | 16 16 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarAdvectionUniformUV      | 16 16  4 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarAdvectionShearedU       | 16  4 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarAdvectionRigidRotation  | 16 16  4 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarAdvectionDiffusion      | 16 16 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarDiffusionGaussian       | 16 16 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| ScalarDiffusionSine           | 16 16  4 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| IsentropicVortexAdvecting     | 48 48  4 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| IsentropicVortexStationary    | 48 48  4 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| TaylorGreenAdvecting          | 16 16 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| TaylorGreenAdvectingDiffusing | 16 16 16 | per | per | per | None  |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| CouetteFlow                   | 32 16  4 | per | NSW | per | None  |                |
|                               |          |     | Dir |     |       |                |
+-------------------------------+----------+-----+-----+-----+-------+----------------+
| EkmanSpiral                   | 4 4 400  | per | per | NSW | Geo   | +Coriolis      |
|                               |          |     |     | SW  |       | +gravity       |
+-------------------------------+----------+-----+-----+-----+-------+----------------+

Scalar Advection by Uniform Flow in X-Direction
------------------------------------------------

Test Location: `Tests/test_files/ScalarAdvectionUniformU`_

.. _`Tests/test_files/ScalarAdvectionUniformU`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarAdvectionUniformU

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a1| image:: figures/tests/scalar_advec_uniform_u_start.png
        :width: 300

.. |b1| image:: figures/tests/scalar_advec_uniform_u_end.png
        :width: 300

.. _fig:scalar_advection_u

.. table:: Advection of a spherical blob in a uniform velocity field (100,0,0)

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a1|                         |                       |b1|                           |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.0264788).  |
   +-----------------------------------------------------+------------------------------------------------------+

Scalar Advection by Uniform Flow in XY Plane
------------------------------------------------
This tests scalar advection with triply periodic boundaries.

Test Location: `Tests/test_files/ScalarAdvectionUniformUV`_

.. _`Tests/test_files/ScalarAdvectionUniformUV`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarAdvectionUniformUV

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a2| image:: figures/tests/scalar_advec_uniform_uv_start.png
        :width: 300

.. |b2| image:: figures/tests/scalar_advec_uniform_uv_end.png
        :width: 300

.. _fig:scalar_advection_uv

.. table:: Advection of a spherical blob in a uniform velocity field (10,5,0)

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a2|                         |                        |b2|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.6937161).  |
   +-----------------------------------------------------+------------------------------------------------------+

Scalar Advection by Sheared Flow
------------------------------------------------
This tests scalar advection in horizontal flow in the x-direction with triply periodic boundaries.

Test Location: `Tests/test_files/ScalarAdvectionShearedU`_

.. _`Tests/test_files/ScalarAdvectionShearedU`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarAdvectionShearedU

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a3| image:: figures/tests/scalar_advec_sheared_u_start.png
        :width: 300

.. |b3| image:: figures/tests/scalar_advec_sheared_u_end.png
        :width: 300

.. _fig:scalar_advection_sheared_u

.. table:: Advection of a spherical blob in a uniform shearing velocity field (8 log( (z+z0)/z0 ) / log ( (zref+z0)/z0 )
   with z0 = 0.1 and zref = 80 in a domain 8x8x8

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a3|                         |                        |b3|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.9819669.   |
   +-----------------------------------------------------+------------------------------------------------------+

Scalar Advection: Rigid Rotation
----------------------------------
This tests scalar advection in a flow field representing rigid body rotation.

Test Location: `Tests/test_files/ScalarAdvectionRigidRotation`_

.. _`Tests/test_files/ScalarAdvectionRigidRotation`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarAdvectionRigidRotation

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvecAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a4| image:: figures/tests/scalar_advec_rigid_rot_start.png
        :width: 300

.. |b4| image:: figures/tests/scalar_advec_rigid_rot_end.png
        :width: 300

.. _fig:scalar_advection_rigid_rot

.. table::   Advection of a 2D blob in a rotating velocity field (.5-y, x-.5, 0) in a domain 1x1x1

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a4|                         |                        |b4|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 6.283185).   |
   +-----------------------------------------------------+------------------------------------------------------+

See http://ammar-hakim.org/sj/je/je16/je16-ldg.html#rigid-body-rotating-flow

Scalar Diffusion: Sphere of Scalar
------------------------------------------------
This tests scalar diffusion with triply periodic boundaries.

Test Location: `Tests/test_files/ScalarDiffusionGaussian`_

.. _`Tests/test_files/ScalarDiffusionGaussian`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarDiffusionGaussian

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a5| image:: figures/tests/scalar_diff_start.png
        :width: 300

.. |b5| image:: figures/tests/scalar_diff_end.png
        :width: 300

.. _fig:scalar_diffusion_gaussian

.. table:: Diffusion of a spherical blob of scalar

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a5|                         |                        |b5|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.01).       |
   +-----------------------------------------------------+------------------------------------------------------+

Scalar Diffusion: Sinusoidal Variation of Scalar
------------------------------------------------
This tests scalar diffusion with triply periodic boundaries.

Test Location: `Tests/test_files/ScalarDiffusionSine`_

.. _`Tests/test_files/ScalarDiffusionSine`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarDiffusionSine

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a6| image:: figures/tests/scalar_diff_sine_start.png
        :width: 300

.. |b6| image:: figures/tests/scalar_diff_sine_end.png
        :width: 300

.. _fig:scalar_diffusion_sine

.. table:: Diffusion of a scalar initialized as sin(x)

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a6|                         |                        |b6|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.2).        |
   +-----------------------------------------------------+------------------------------------------------------+


Scalar Advection/Diffusion by Uniform Flow
------------------------------------------------
This tests scalar advection and diffusion with triply periodic boundaries.

Test Location: `Tests/test_files/ScalarAdvectionDiffusionUniformU`_

.. _`Tests/test_files/ScalarAdvectionDiffusionUniformU`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/ScalarAdvectionDiffusionUniformU

Problem Location: `Exec/ScalarAdvDiff`_

.. _`Exec/ScalarAdvDiff`: https://github.com/erf-model/ERF/tree/development/Exec/ScalarAdvDiff

.. |a7| image:: figures/tests/scalar_advec_diff_start.png
        :width: 300

.. |b7| image:: figures/tests/scalar_advec_diff_end.png
        :width: 300

.. _fig:scalar_diffusion_sine

.. table:: Advection and diffusion of a spherical blob in a uniform velocity field (100,0,0)

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a7|                         |                        |b7|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Scalar concentration at t=0.                      |   Scalar concentration at 20 steps (t = 0.01).       |
   +-----------------------------------------------------+------------------------------------------------------+

Isentropic Vortex: Stationary
---------------------------------
This tests advection of an isentropic vortex with triply periodic boundaries.

Test Location: `Tests/test_files/IsentropicVortexStationary`_

.. _`Tests/test_files/IsentropicVortexStationary`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/IsentropicVortexStationary

Problem Location: `Exec/IsentropicVortex`_

.. _`Exec/IsentropicVortex`: https://github.com/erf-model/ERF/tree/development/Exec/IsentropicVortex

Isentropic Vortex: Advecting
---------------------------
This tests advection of an isentropic vortex with triply periodic boundaries.

Test Location: `Tests/test_files/IsentropicVortexAdvecting`_

.. _`Tests/test_files/IsentropicVortexAdvecting`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/IsentropicVortexAdvecting

Problem Location: `Exec/IsentropicVortex`_

.. _`Exec/IsentropicVortex`: https://github.com/erf-model/ERF/tree/development/Exec/IsentropicVortex

Taylor Green Vortex: Advection
------------------------------------------------
This tests advection and diffusion with triply periodic boundaries.

Test Location: `Tests/test_files/TaylorGreenAdvecting`_

.. _`Tests/test_files/TaylorGreenAdvecting`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/TaylorGreenAdvecting

Problem Location: `Exec/TaylorGreenVortex`_

.. _`Exec/TaylorGreenVortex`: https://github.com/erf-model/ERF/tree/development/Exec/TaylorGreenVortex

Taylor Green Vortex: Advection and Diffusion
------------------------------------------------
This tests advection and diffusion with triply periodic boundaries.

Test Location: `Tests/test_files/TaylorGreenAdvectingDiffusing`_

.. _`Tests/test_files/TaylorGreenAdvectingDiffusing`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/TaylorGreenAdvectingDiffusing

Problem Location: `Exec/TaylorGreenVortex`_

.. _`Exec/TaylorGreenVortex`: https://github.com/erf-model/ERF/tree/development/Exec/TaylorGreenVortex

.. |a8| image:: figures/tests/TGV_start.png
        :width: 300

.. |b8| image:: figures/tests/TGV_end.png
        :width: 300

.. _fig:taylor_green_vortex

.. table:: Scalar concentration

   +-----------------------------------------------------+------------------------------------------------------+
   |                        |a8|                         |                        |b8|                          |
   +-----------------------------------------------------+------------------------------------------------------+
   |   Flow field at t=0.                                |   Flow field at 10 steps (t = 1.6).                  |
   +-----------------------------------------------------+------------------------------------------------------+

Channel Flow: DNS
------------------------

This tests DNS flow in a channel which is periodic in x and z, and no-slip-wall on both y-faces

Test Location:

Problem Location: `Exec/ChannelFlow`_

.. _`Exec/ChannelFlow`: https://github.com/erf-model/ERF/tree/development/Exec/ChannelFlow

Channel Flow: LES
------------------------

This tests LES flow in a channel which is periodic in x and z, and no-slip-wall on both y-faces

Test Location:

Problem Location: `Exec/ChannelFlow`_

.. _`Exec/ChannelFlow`: https://github.com/erf-model/ERF/tree/development/Exec/ChannelFlow

Couette Flow
------------
Test Location: `Tests/test_files/CouetteFlow`_

.. _`Tests/test_files/CouetteFlow`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/CouetteFlow

Problem Location: `Exec/CouetteFlow`_

.. _`Exec/CouetteFlow`: https://github.com/erf-model/ERF/tree/development/Exec/CouetteFlow

Ekman Spiral
---------------------------
This tests the Coriolis and geostrophic forcing.

Test Location: `Tests/test_files/EkmanSpiral`_

.. _`Tests/test_files/EkmanSpiral`: https://github.com/erf-model/ERF/tree/development/Tests/test_files/EkmanSpiral

Problem Location: `Exec/EkmanSpiral`_

.. _`Exec/EkmanSpiral`: https://github.com/erf-model/ERF/tree/development/Exec/EkmanSpiral
