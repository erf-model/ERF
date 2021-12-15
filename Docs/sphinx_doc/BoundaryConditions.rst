
 .. role:: cpp(code)
    :language: c++

.. _sec:domainBCs:

Domain Boundary Conditions
--------------------------

ERF allows users to specify boundary condition with keywords.
For a more detailed discussion of boundaries, see :ref:`sec:boundaries`

Domain Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify boundary conditions with keywords, we use the following options
preceded by “xlo”, “xhi”, “ylo”, “yhi”, “zlo”, and “zhi”:

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'Outflow'                                                               |             |           |
|                    | * 'Inflow'                                                                |             |           |
|                    | * 'SlipWall'                                                              |             |           |
|                    | * 'NoSlipWall'                                                            |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

To use the same example problem as above, the following:

::

    xlo.type = "Inflow"
    xhi.type = "Outflow"
    zlo.type = "SlipWall"
    zhi.type = "SlipWall"

    geometry.is_periodic = 0 1 0

would define a problem with inflow in the low-\ :math:`x` direction,
outflow in the high-\ :math:`x` direction, periodic in the :math:`y`-direction,
and slip wall on the low and high :math:`y`-faces, and
Note that no keyword is needed for a periodic boundary, here only the
specification in ``geometry.is\_periodic`` is needed.

.. _sec:dirichlet:

Dirichlet Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERF provides the ability to specify constant Dirichlet BCs in the inputs file. We use the following options
preceded by “xlo”, “xhi”, “ylo”, “yhi”, “zlo”, and “zhi”:

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| velocity           | Sets boundary velocity for mass inflows                                   |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| density            | Sets boundary density for mass inflows                                    |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| tracer             | Sets boundary tracer concentration at inflow faces                        |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| theta              | Sets potential temperature at inflow faces                                |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+


As an example,

::

    xlo.type                =   "Inflow"
    xlo.velocity            =   1.  0.  0.
    xlo.density             =   1.
    xlo.tracer              =   0.
    xlo.theta               =   1.

sets the boundary condtion type at the low x face to be an inflow with
xlo.type = “Inflow”.
Then xlo.velocity = 1. 0. 0. sets the inflow velocity,
xlo.density = 1. sets the inflow density,
xlo.tracer = 0. sets the inflow tracer value, and
xlo.theta = 1. sets the inflow potential temperature.

Users can create more complex Dirichlet boundary condtions by writing
their own fill function in ``BCFill.H``, then using that function to create
an ``amrex::StateDescriptor::BndryFunc`` object and specifying which variables
will use it in ``Setup.cpp``. More information on boundary conditions is in
the section :ref:`sec:physicalBCs` below.

.. _sec:physicalBCs:

Physical Boundaries
-------------------

Physical boundaries are special for several reasons.  There are a number of
standard types typical of PDE problems (reflecting, extrapolated, etc),
and a special one that indicates external Dirichlet. In the case of Dirichlet,
the user supplies data to fill grow cells.

ERF provides the ability to specify constant Dirichlet BCs
in the inputs file (see section :ref:`sec:dirichlet`).
Users can create more complex Dirichlet boundary condtions by writing
their own fill function in ``BCFill.H``, then using that function to create
an ``amrex::StateDescriptor::BndryFunc`` object and specifying which variables
will use it in ``Setup.cpp``.

It is important to note that external Dirichlet boundary data is to be specified as
if applied on the face of the cell bounding the domain, even for cell-centered
state data. For cell-centered data, the array passed into the
boundary condition code is filled with cell-centered values in the valid
region and in fine-fine, and coarse-fine grow cells. Additionally, grow cells
for standard extrapolation and reflecting boundaries are pre-filled. The
differential operators throughout ERF are aware of the special boundaries
that are Dirichlet and wall-centered, and the stencils are adjusted accordingly.

For convenience, ERF provides a limited set of mappings from a physics-based boundary condition
specification to a mathematical one that the code can apply. This set
(See ``AMReX/Src/Base/AMReX_BC_TYPES.H`` for more detail):
includes

-  *Outflow*:

   -  velocity: FOEXTRAP

   -  temperature: FOEXTRAP

   -  scalars: FOEXTRAP

-  *No Slip Wall with Adiabatic Temp*:

   -  velocity: EXT_DIR, :math:`u=v=0`

   -  temperature: REFLECT_EVEN, :math:`dT/dt=0`

   -  scalars: HOEXTRAP

-  *Slip Wall with Adiabatic Temp*:

   -  velocity: EXT_DIR, :math:`u_n=0`; HOEXTRAP, :math:`u_t`

   -  temperature: REFLECT_EVEN, :math:`dT/dn=0`

   -  scalars: HOEXTRAP

The keywords used above are defined:

-  INT_DIR: data taken from other grids or interpolated

-  EXT_DIR: data specified on EDGE (FACE) of bndry

-  HOEXTRAP: higher order extrapolation to EDGE of bndry

-  FOEXTRAP: first order extrapolation from last cell in interior

-  REFLECT_EVEN: :math:`F(-n) = F(n)` true reflection from interior cells

-  REFLECT_ODD: :math:`F(-n) = -F(n)` true reflection from interior cells
