.. _sec:ActuatorCore:

Actuator core for bodies in the flow
====================================

The actuator core (``Source/MovingBodies/Core``, namespace ``erf_actuator``)
is the part of ERF that connects a body in the flow to the mesh: it samples the
flow velocity at a set of points, spreads point forces back onto the momentum
equations, keeps sampling lines in the body's wake and accumulates running
statistics of the body's diagnostics. It knows nothing about what the body is.
A body model supplies the positions of its points every step, receives the
velocity sampled there, and hands back the force each point exerts on the
fluid; the core does the rest on uniform, stretched and terrain-following
meshes alike. The core is built with ``-DERF_ENABLE_MOVING_BODIES=ON`` (CMake)
or ``USE_MOVING_BODIES = TRUE`` (GNU make), which defines
``ERF_USE_MOVING_BODIES``; the default build does not compile it. The sampler
and the spreader take one flat list of point coordinates (``x, y, z`` per
point), whatever the bodies are; a body model owns its points and refreshes
their positions every step.

Velocity sampling
-----------------

``sample_velocity`` interpolates each velocity component to a point from its own
staggered grid: ``u`` from the x faces, ``v`` from the y faces and ``w`` from
the z faces, bilinearly in the horizontal between the four surrounding columns
and linearly in the physical height within each column, using the face heights
of ``z_phys_nd`` on a stretched or terrain-following mesh (``face_height`` in
``ERF_ActuatorGeometry.H``, shared with the spreader so that both see the same
face positions). A field linear in ``x``, ``y`` and the physical height is
therefore reproduced exactly, on every mesh type. Every MPI rank samples the
points whose containing cell lies in one of its boxes, reading neighbours from
the ghost cells (the velocity fields need one filled ghost cell, and
``z_phys_nd`` one filled ghost node), and the
values are summed across ranks, so the result does not depend on the domain
decomposition and every rank holds the same velocities on return. A point
outside the domain, or one found by no rank or by more than one, aborts the
run with its coordinates: the caller owns the points and must keep them inside
the domain.

``sample_cell_scalar`` does the same for one component of a cell-centred field
(bilinear between the four cell centres around the point, linear in the
physical height between the two cell centres that bracket it in the column).
``terrain_heights`` returns the height of the terrain surface under each
point: the ``k = 0`` node plane of ``z_phys_nd``, bilinear between the four
nodes around the point's ``(x, y)``, or the domain's lower ``z`` on a
uniform-dz mesh. A body model uses it to place a body's base at a height above
the ground and to measure clearances.

``points_covered_by`` tells whether every point, with the cells within a given
reach around it, lies inside the union of boxes of a level: the sampler needs
the neighbouring cells and the spreader the kernel's reach, and on a refined
level that is not the whole domain. Periodic directions wrap; the first point
outside is named. Its ``CoverZ`` argument says how the vertical is checked: by
the cells within reach of the point's own height (``Reach``, a flat uniform
mesh), by the boxes' footprints at any height (``Footprint``, a stretched or
terrain-following mesh, where the sampler's cell is found by height), or by
every cell the heights within reach may fall in (``Column``, for the
spreader), from the mesh's ground and spacing extremes given by
``mesh_z_bounds``. Those extremes are taken over the whole level, so on steep
terrain a point over high ground is asked for cells that only the columns over
the lowest ground need. ``column_cell_ranges`` gives, for each point, the cells
the heights within reach fall in over the point's own columns (``ColumnCells``),
and ``points_covered_by`` takes them in place of the mesh-wide extremes. A
caller that tests many sets of points against the same boxes builds the boxes'
footprints once (``CoverageGrids``) and passes them instead of the ``BoxArray``.

``region_touches_grids`` answers the opposite question for a caller that must
stay clear of a level: whether the level's boxes may hold a cell of a region
given in metres. In a periodic direction every image of the region counts. A
region that overlaps a cell of the level is never missed; one that only meets a
cell at a face may not count. On a stretched or terrain-following mesh it takes
the level's table of node heights from ``mesh_z_levels``: for each node level,
the lowest and the highest height over every rank's boxes, and the cells the
level holds. The cells of a height range are then those with some bottom node
not above it and some top node not below it, which holds whether or not the
boxes reach the ground. With the ground and spacing extremes of
``mesh_z_bounds`` alone the answer is wider, and where a box does not stand on
the ground it is the whole column. The cell indices are bounded before they
become integers, so a region far outside the domain, or an infinite one, is
safe; a height range with an infinite end takes the cells of its finite end.
``LevelHeights`` keeps the table between calls and builds it again when the
level's boxes or its domain change, or when the caller asks (a moving terrain).

``sample_velocity`` is collective: every rank passes the same points, and the
answers are summed over the ranks. ``sample_velocity_local`` is the part on one
rank, with no reduction: it returns the velocity at the points this rank's
boxes hold and a count of 1 for each, and zeros elsewhere. A body model whose
ranks hold different points calls it and exchanges the answers itself.

Force spreading and momentum source
-----------------------------------

``spread_forces`` carries the force each point exerts on the fluid onto ERF's
face-centred momentum sources with a 3-D Gaussian kernel of width ``epsilon``
(in metres; one width for every point, or one per point), cut off beyond three
widths. The kernel is normalised discretely
on each staggered grid, with the face volumes (``dx dy dz detJ`` on a
terrain-following mesh, from ``face_volume``), so the source integrates back to
the point force exactly for every point and component, whatever the resolution
and however much of the kernel the ground or the domain top cuts off; the
kernel's shape is Gaussian only where it is resolved. Each point is visited
over the faces within three widths of it, and again as its periodic images, so
a kernel wraps across a periodic boundary and the cost grows with the number
of points, not with the mesh: on a uniform mesh each point visits the
``(6 epsilon / dx)^3`` faces around it; on a stretched or terrain-following
mesh the faces within ``3 epsilon`` in physical height are not known from the
index, so the search runs, within the horizontal reach, over the cells whose
node heights anywhere on the level can lie within that height range. Every rank accumulates the faces its boxes own,
each physical face once however the boxes share their boundary faces, and the
normalisations are summed across ranks, so the result does not depend on the
domain decomposition. The sources are set, not added to, on every valid face;
faces on a non-periodic domain boundary get no source. ``integrate_source``
returns the integral of a face-centred source over the mesh, which the
regression tests of a body model compare with the force it carries.

.. note::

   **Why the normalisation is discrete.** Actuator codes in the AMR-Wind and
   SOWFA lineage divide each point's force by the analytic Gaussian volume,
   ``epsilon^3 pi^(3/2)``, and then sum the kernel over the cells. The sum over
   cell volumes equals that analytic volume only when the kernel lies wholly
   inside the domain, is resolved by the mesh and covers cells of one size.
   Otherwise the momentum the fluid receives is not the force the point
   carries, and nothing reports it:

   * the part of a kernel below the ground or above the domain top is never
     deposited; a point at height ``h`` loses the fraction ``erfc(h/epsilon)/2``
     of its force, up to half for a point near the ground;
   * with ``epsilon`` below about one cell the cell sum of a Gaussian differs
     from its integral (by Poisson summation the worst case in one direction is
     ``2 exp(-pi^2 (epsilon/dx)^2)``: 3e-4 of the force in 3-D at
     ``epsilon = dx``, several per cent at ``0.7 dx``, tens of per cent at
     ``0.5 dx``), and the error changes sign with the point's position in the
     cell; the usual ``epsilon >= 2 dx`` rule of the actuator-line literature
     is about the smoothness of the body force, not this sum;
   * on a stretched or terrain-following mesh the cell volumes vary across the
     kernel, which no constant can account for (about 0.1 per cent of the
     force for a resolved kernel on a mesh stretched by 10 per cent per cell;
     up to 5 per cent if one cell's volume stands for all).

   The ground cut is the large one. Note that the renormalised kernel of a
   point cut by the ground is applied above the point: its centroid moves up
   by ``0.39 epsilon`` for a point ``0.3 epsilon`` above the ground. The
   cut-off at three widths drops 0.044 per cent of the mass, which the
   normalisation puts back inside, narrowing the kernel by 0.13 per cent.

   The symptoms are a wake deficit slightly too weak, a force felt by the flow
   that depends on the grid spacing and on how low the body sits, and momentum
   budgets that do not close. ERF pays a second pass per point to sum the
   kernel over the faces it actually reaches, with their volumes, so the
   integrated source equals the force exactly in every case.

``spread_forces`` works in two passes, which are also available on their own.
``spread_partial_norms`` returns, for each point, the part of its normalisation
(the sum of the kernel weight times the face volume) over the faces this rank's
boxes own, with no reduction. ``add_scaled_kernels`` adds each point's force,
divided by its whole normalisation, times the kernel weight to the faces of
this rank's boxes. ``spread_forces`` sums the parts over the ranks between
the two; a body model whose ranks hold different points sums them over the
ranks a point reaches. ``add_spread_forces`` is ``spread_forces`` without the
zeroing of the sources, for several sets of points spread in turn. On a stretched or terrain-following mesh each kernel
visits only the cells its heights can fall in, from a table of the level's node
heights (``mesh_z_bounds`` with ``over_all_ranks = false``, which is this
rank's table and is not a collective call).

Limits
------

The core takes the mesh as Cartesian in the horizontal: the faces sit at
``prob_lo + i dx`` and the face volumes are ``dx dy dz detJ``, with no map
factors. On a map-projected run (map factors other than 1) the spread force
would no longer integrate to the point force and the sampling positions would
be off, so a body model must refuse such a run until the core carries the map
factors.

A body model works on one level, which it chooses, and keeps its points with
the kernel's reach (three widths) on that level's grids, checked with
``points_covered_by`` at start-up and after every regrid: the normalisation
runs over the faces of one level, so on a refined level whose patch cuts a
kernel the whole force would land in the covered part. The level that samples
is the level that spreads.

The spreader adds its terms in a fixed order, without atomic adds: the
normalisation of a point over a rank's boxes in box order, then over the ranks;
the source at a face point by point. The source therefore does not change from
run to run on a GPU, as it does not on a CPU.

Bodies on terrain
-----------------

The sampler and the spreader work in the physical heights of a fitted mesh
(``z_phys_nd``) and with the cell volumes (``detJ``), so a body's points and
the kernel follow the terrain without further change. A body model that
places its base at a height above the ground adds the terrain height from
``terrain_heights`` once at start-up; each vertical wake line stops at the
ground under its own centre, which the body model takes from
``terrain_heights`` at the line centres (one ground height per downstream
distance); a single ground height is the hub's, and on terrain rising
downstream a line's lowest points would then lie inside the ground, where the
sampler stops the run naming them.

Wake sampling lines
-------------------

``WakeLines`` puts sampling lines behind one body: at each requested distance,
in units of the body's diameter along the horizontal projection of the body's
axis through its centre (a wake follows the wind at the body's height, and a
tilted axis would otherwise carry the far lines into the ground), one
horizontal line normal to that direction and one vertical line, each of
``num_points`` points spanning ``half_width`` diameters either side of the
axis. The vertical line is clipped at the ground height given, so its lowest
point may sit closer to the axis than ``half_width``; each point's offset ``s``
along its line, in diameters, is in the files. The velocity sampled at the
points with the actuator sampler is written as it is to ``<output_root>_wake.csv``
(``time, xD, line, s, x, y, z, u, v, w``) and accumulated into a running time
average that ``<output_root>_wake_avg.csv`` always holds (``samples`` is the
number of samples in it). The running sums, with the sampling geometry, are
written to and read from a checkpoint directory
(``<name>_wake_avg.dat``: a ``key = value`` header with the sample count, the
number of values, the hub, the axis and the diameter, then the offsets, the
point positions and the running sums), so the average continues across a
restart; a checkpoint without the file starts the average afresh and says so.

Running statistics
------------------

``RunningStats`` keeps, for a named body and a list of named scalar
quantities, the mean, the root mean square, the minimum and the maximum of
every sample added since the averaging start, with the sample count and the
first and last sample times (sample averages, not time averages; the root
mean square is that of the signal, not of its fluctuation, and the sums are
accumulated in double precision); ``<output_root>_stats.csv`` always holds the
current values, one row per quantity. The sums are written to and read from a
checkpoint directory (``<name>_stats.dat``: a ``key = value`` header with the
sample count, the number of quantities and the first and last sample times,
then one row per quantity with its sum, sum of squares, minimum and maximum),
so the statistics continue across a restart; a checkpoint without the file
starts them afresh and says so. Quantity names carry no spaces. These are
the numbers a turbulent-inflow run is compared on, where a single instant
means little.

Diagnostics files
-----------------

``open_log`` opens a diagnostics file (CSV or whitespace-separated) for
writing, creating its directory; the run stops if the directory cannot be
created or the file cannot be opened.
A fresh run truncates the file; a restarted run appends to what the run it
continues wrote, and the header is written only when the file does not exist yet or is empty.
