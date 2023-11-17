
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

ERF has the option to include Lagrangian particles in addition to the mesh-based solution.  Currently
there are two particle types available in ERF: tracer_particles and hydro_particles.
The particle functionality is very simple and meant for demonstration.
The particles are initialized one per mesh cell in a
vertical plane at :math:`i = 3` for tracer particles and a horizontal plane at :math:`k = 23` for hydro particles.
The tracer particles are advected by the velocity field; the hydro particles fall with a velocity determined by gravity minus drag.

However, the AMReX particle data structure is very general and particles may take on a number of
different roles in future.

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

and / or

::

   erf.use_hydro_particles = 1

in the inputs file or on the command line at runtime.

Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec/DevTests/ParticlesOverWoA.

To visualize the number of particles per cell as a mesh-based variable, add
``tracer_particle_count`` (if you have set ``erf.use_tracer_particles``) and
 ``hydro_particle_count`` (if you have set ``erf.use_tracer_particles``)
to the line in the inputs file that begins

::

   erf.plot_vars_1 =

