
 .. role:: cpp(code)
    :language: c++

.. _sec:domainBCs:

Domain Boundary Conditions
--------------------------

ERF provides two ways to specify boundary condition types: integer keys or keywords.
An inputs file must choose one style or the other, they cannot be mixed.
Here we provide examples of both styles. We then discuss how to provide Dirichlet
boundary conditions. For a more detailed discussion of boundaries, see :ref:`sec:boundaries`

Boundary conditions with integer keys
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify boundary conditions with integer keys we use

-  erf.lo_bc: boundary type of each low face

-  erf.hi_bc: boundary type of each high face

The valid boundary types are:

+---------------------------+--------------------+
| 0 – Interior / Periodic   | 3 – Symmetry       |
+---------------------------+--------------------+
| 1 – Inflow                | 4 – Slip Wall      |
+---------------------------+--------------------+
| 2 – Outflow               | 5 – No Slip Wall   |
+---------------------------+--------------------+

Note: ``erf.lo_bc`` and ``erf.hi_bc`` must be consistent with
``geometry.is_periodic`` —if the domain is periodic in a particular
direction then the low and high bc’s must be set to 0 for that direction.

As an example, the following:

::

    erf.lo_bc = 1 0 4 
    erf.hi_bc = 2 0 4 

    geometry.is_periodic = 0 0 1

would define a problem with inflow (1) in the low- :math:`x` direction,
outflow (2) in the high- :math:`x` direction, periodic in the :math:`y`-direction.
and slip wall (4) on the low and high :math:`z`-faces.

Boundary conditions with keywords
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify boundary conditions with keywords, we use the following options
preceded by “xlo”, “xhi”, “ylo”, “yhi”, “zlo”, and “zhi”:

+--------------------+---------------------------------------------------------------------------+-------------+-----------+
|                    | Description                                                               |   Type      | Default   |
+====================+===========================================================================+=============+===========+
| type               | Used to define boundary type. Available options include:                  |  String     |  None     |
|                    |                                                                           |             |           |
|                    | * 'po'  or 'pressure_outflow'                                             |             |           |
|                    | * 'mi'  or 'mass_inflow'                                                  |             |           |
|                    | * 'sw'  or 'slip_wall'                                                    |             |           |
|                    | * 'nsw' or 'no_slip_wall'                                                 |             |           |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+

To use the same example problem as above, the following:

::

    xlo.type = mi
    xhi.type = po
    zlo.type = sw
    zhi.type = sw

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
| tracer             | Sets boundary tracer concentration for mass inflows                       |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+
| temp               | Sets temperature for mass inflows                                         |    Real     |  None     |
+--------------------+---------------------------------------------------------------------------+-------------+-----------+


As an example,

::

    xlo.type                =   "mass_inflow"
    xlo.velocity            =   1.  0.  0.
    xlo.density             =   1.
    xlo.tracer              =   0.
    xlo.temp                =   1.

sets the boundary condtion type at the low x face to be an inflow with
xlo.type = “mass_inflow”.
Then xlo.velocity = 1. 0. 0. sets the inflow velocity,
xlo.density = 1. sets the inflow density,
xlo.tracer = 0. sets the inflow tracer value, and
xlo.temp = 1. sets the inflow temperature.

Another example, from the Couette flow example,

::

    erf.lo_bc                =  4 4 5
    erf.hi_bc                =  5 5 5

    # 0 = Interior/Periodic  3 = Symmetry
    # 1 = Inflow             4 = SlipWall
    # 2 = Outflow            5 = NoSlipWall

    # Boundary condition
    zhi.velocity            =   2.  0.  0.

Here, erf.hi_bc = 5 5 5 sets the boundary conditions on all high faces to no slip walls.
zhi.velocity = 2. 0. 0. sets the wall at the high-z face to be moving in the x-direction.
Note that ERF allows walls to move tangentially, but not in the normal direction.

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
