
 .. role:: cpp(code)
    :language: c++

 .. _Particles:

Particles
=========

ERF has the option to include Lagrangian particles in addition to the mesh-based solution.  Currently the
particle functionality is very simple -- the particles are initialized randomly, one per mesh cell
in a particular plane, and are advected by the velocity field.

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

   erf.use_tracer_particles = true

in the inputs file or on the command line at runtime.

Currently, by default, the particles are initialized at cell centers, one per cell when the cell index is
(3,j,k), with zero initial velocity.  They are advanced in time every time step.

Caveat: the particle information is currently output when using the AMReX-native plotfile format, but not
when using netcdf.  Writing particles into the netcdf files is a WIP.

To see an example of using the particle functionality, build the executable using gmake in Exec/DevTests/ParticlesOverWoA.

To visualize the number of paritcles per cell as a mesh-based variable, add ``particle_count`` to the line in the inputs file

::

   erf.plot_vars_1 =

