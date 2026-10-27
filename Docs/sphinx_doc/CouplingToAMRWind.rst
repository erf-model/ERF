
 .. role:: cpp(code)
    :language: c++

 .. _CouplingToAMRWind:

Coupling To AMR-Wind
====================

The simplest form of coupling is one-way, file-based coupling. With this approach, an
ERF simulation is run and output files with relevant simulation data are stored at
regular time intervals. Once the ERF simulation is complete, an AMR-Wind simulation
begins, using the ERF output files as time-dependent boundary conditions or forcing.

ERF previously also wrote 1D vertical column output (``erf.output_1d_column``) for
time-height profile assimilation in AMR-Wind. That output has been removed; the 2D
file-based coupling below is the supported route.

2D File-based coupling
----------------------

For 2D file-based coupling to AMR-Wind, vertical planes are saved from an ERF simulation
to define inflow boundary conditions for AMR-Wind.  These files are in native-AMReX format and
contain the AMReX data structure known as a :cpp:`BndryRegister`.
These files are specifically meant to be read by AMR-Wind, which is also AMReX-based, or ERF itself.

To generate the boundary data files for each of the four planes,
the following (or similar) should be added to the input file:

.. code-block:: none

  erf.output_bndry_planes = 1
  erf.bndry_output_planes_interval = 2
  erf.bndry_output_start_time = 0.0
  erf.bndry_output_planes_file = "BndryFiles"
  erf.bndry_output_var_names = temperature velocity density
  erf.bndry_output_box_lo = 256. 256.
  erf.bndry_output_box_hi = 768. 768.

The above inputs will output boundary planes of data at x=xlo, x=xhi, y=ylo and y=yhi, where
xlo and ylo are defined by :cpp:`bndry_output_box_lo` and
xhi and yhi are defined by :cpp:`bndry_output_box_hi`. Note that the selected lo and hi values must
be interior to the first and last cell centers on the coarse grid (i.e., at least one half cell width
from the boundary). Except if the domain is periodic
in that direction, any location up to and including the domain boundary may be chosen.
In this case the variables that are
written are temperature, velocity and density, and they are written every 2 coarse time steps starting at
:cpp:`bndry_output_start_time` which is 0 in this case.

The start time counts from the start of the run (elapsed simulation time), also when
the run's clock starts at a date (:cpp:`start_datetime`, or the date of a WRF input file).

Within :cpp:`BndryFiles`, each output step writes a folder :cpp:`bndry_outputNNNNN` named by
its step number and adds a row ``step time`` to the ascii file :cpp:`time.dat`; the time is
that of the run's clock, so it includes :cpp:`start_datetime` when one is given.
A run restarted from a checkpoint continues the series it finds there. It keeps the rows of
:cpp:`time.dat` up to the restart step and drops the later ones, which a run that went on past
the checkpoint wrote and the restarted run writes again. :cpp:`time.dat` is rewritten only
when the run writes its first plane, so a run that stops before then leaves it as it was.
The plane folders of dropped rows are not removed, and AMReX keeps each plane folder it
writes again as :cpp:`bndry_outputNNNNN.old.<n>`. The restarted run writes a plane at the
restart step itself only if the series lacks it and the step is an output step, so that with
the same inputs and time steps the series comes out as a run without the restart writes it.
Output switched on at a restart starts at the first output step. A run started from scratch
begins a new series and drops any rows an earlier run left in :cpp:`time.dat`. A run cannot
write its output planes into the folder it reads its input planes from.

We also have the functionality in ERF to read in these types of files;
for this one would add the following (or similar) line to the inputs file:

.. code-block:: none

  erf.input_bndry_planes = 1
  erf.bndry_file = "BndryFiles"
  erf.bndry_input_var_names = density temperature velocity

When run with these inputs, ERF will read in the time sequence of files contained in the folder :cpp:`BndryFiles`,
and perform time interpolation as necessary. The only assumption about the times associated with the files
is that the start and end times of the current simulation
lie in the time period covered by the files in :cpp:`BndryFiles`.  Within :cpp:`BndryFiles` there is an
ascii file :cpp:`time.dat` which contains the (originating) timesteps and physical times associated with each of the files.

It is assumed at this point that the physical domain of the simulation reading the files is exactly the physical
domain specified by :cpp:`bndry_output_box_lo` and :cpp:`bndry_output_box_hi` when the files were written.  If not, ERF will
abort with an error message.

We note that the boundary plane data will only be used on faces identified in the inputs file as inflow faces, i.e. if
we specific inflow/outflow in the x-direction, and periodic in the y-direction, as below, then only the "xlo" boundary data
from :cpp:`BndryFiles` will actually be used.

::

    geometry.is_periodic = 0 1 0

    xlo.type = "Inflow"
    xhi.type = "Outflow"

