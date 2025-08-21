
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

ERF has the option to include Lagrangian particles in addition to the mesh-based solution.
Currently there is one example of particle types available in ERF: tracer_particles.

The particle functionality is very simple and meant for demonstration.

The particles are initialized one per mesh cell in a vertical plane at :math:`i = 3` for tracer particles.

The tracer particles are advected by the velocity field with optional sedimentation.

We note that unless the domain is periodic in the vertical direction, any particles that
cross the bottom boundary during the advection step will be moved back into the domain
at a location 1/5 of the way between the bottom boundary and the top of the cell at k = 0.

However, the AMReX particle data structure is very general and particles may take on a number of
different roles in future.

.. figure:: figures/ERFParticles.gif
   :alt: Particles in a squall line
   :align: center
   :width: 100%

   Two-dimensional squall line simulation with particles. The particles and contours of cloud water mixing ratio are shown.

To enable the use of particles, one must set

::

   USE_PARTICLES = TRUE

in the GNUmakefile if using gmake, or add

::

   -DERF_ENABLE_PARTICLES:BOOL=ON \

to the cmake command if using cmake.  (See, e.g., ``Build/cmake_with_particles.sh``)

One must also set

::

   erf.use_tracer_particles = 1

in the inputs file or on the command line at runtime.

The time at which the particles are initialize can be controlled by a parameter in the inputs file.
For tracer particles one would set this as

::

   tracer_particles.start_time = 0.5


Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec/DevTests/ParticlesOverWoA.

To visualize the number of particles per cell as a mesh-based variable, add
``tracer_particle_count`` (if you have set ``erf.use_tracer_particles``) and
to the line in the inputs file that begins

::

   erf.plot_vars_1 =



