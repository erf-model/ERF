
.. _sec:TowerFrame:

Lattice tower frame model
=========================

The frame model (``Source/MovingBodies/Towers/ERF_Frame.H``) is a linear,
three-dimensional finite-element model of a lattice tower: every member is a
two-node beam with six degrees of freedom per node, the base joints are fixed
or stand on six-component springs, and the static response to nodal loads and
gravity, the natural modes and the time response are solved in double
precision (``ERF_Frame.H``, ``ERF_FrameDynamics.H``). Its element stiffness and
mass, local axes, section properties, concentrated masses and gravity loads
follow OpenFAST's SubDyn module, and it reads SubDyn input files, so that a tower
described for SubDyn has the same stiffness and natural frequencies in ERF.
A tower type that gives ``erf.conductors.<type>.frame_file``, or generates its
frame with ``erf.conductors.<type>.frame_panels`` (section "Generated lattice
towers"), stands on this frame in a run (section "Coupling to the conductor
lines"); without either, towers use the one-mode tower (``ERF_TowerDynamics.H``,
``erf.conductors.<type>.frequency``). With the members' design data, every
member is checked against its strength in tension and compression (section
"Member checks", ASCE 10-15), and the steel can be set at a temperature that
lowers its stiffness and strength (section "Steel at temperature", EN 1993-1-2).

Conventions
-----------

This section fixes the frame, units and ordering every other section uses.

- Coordinates are those of the input file, in m, with :math:`z` up; gravity acts
  along :math:`-z`.
- Every node has six degrees of freedom in the order :math:`u_x, u_y, u_z` (m),
  then rotations :math:`\theta_x, \theta_y, \theta_z` (rad). Loads use the same order:
  forces (N), then moments (N m). Node :math:`n` (0-based) holds degrees of
  freedom :math:`6n` to :math:`6n+5`.
- Nodes are the joints in input order, then the interior nodes of members that
  are divided into more than one element (SubDyn's ``NDiv``), member by member,
  from the member's first joint towards its second.

Elements
--------

This section gives the element stiffness and its axes.

An element of length :math:`L` (m) from node :math:`a` to node :math:`b` has
local axes :math:`(x, y, z)` with :math:`z` along the element. As SubDyn's
``GetDirCos`` builds them, with :math:`(\Delta x, \Delta y, \Delta z)` the
element's extent and :math:`\Delta_{xy} = \sqrt{\Delta x^2 + \Delta y^2}`:

.. math::

   \hat x = \left(\frac{\Delta y}{\Delta_{xy}}, -\frac{\Delta x}{\Delta_{xy}}, 0\right), \quad
   \hat y = \left(\frac{\Delta x \Delta z}{L \Delta_{xy}}, \frac{\Delta y \Delta z}{L \Delta_{xy}}, -\frac{\Delta_{xy}}{L}\right), \quad
   \hat z = \frac{(\Delta x, \Delta y, \Delta z)}{L},

so a horizontal element has a horizontal :math:`\hat x`. A vertical element
pointing up has the frame's axes; one pointing down has :math:`\hat y` and
:math:`\hat z` reversed. The member's spin (SubDyn's ``MSpin``, degrees on input)
then turns :math:`\hat x` and :math:`\hat y` about :math:`\hat z`.

In its local axes the element is SubDyn's ``ElemK_Beam``: an Euler-Bernoulli
beam (``FEMMod`` 1), or a Timoshenko beam (``FEMMod`` 3) whose bending terms carry
the shear factors

.. math::

   \Phi_x = \frac{12 E I_{yy}}{G \kappa_x A L^2}, \qquad
   \Phi_y = \frac{12 E I_{xx}}{G \kappa_y A L^2},

where :math:`E` and :math:`G` are the Young's and shear moduli (Pa), :math:`A` the
area (m\ :sup:`2`), :math:`I_{xx}` and :math:`I_{yy}` the second moments about the local
:math:`x` and :math:`y` axes (m\ :sup:`4`) and :math:`\kappa_x A`, :math:`\kappa_y A` the
shear areas along local :math:`x` and :math:`y`. :math:`\Phi_x` softens bending in the
local :math:`x`-:math:`z` plane, :math:`\Phi_y` in the :math:`y`-:math:`z` plane;
both are 0 for Euler-Bernoulli beams. The axial stiffness is :math:`EA/L` and the
torsional stiffness :math:`G J_t / L`, with :math:`J_t` the torsion constant
(m\ :sup:`4`). The element's stiffness in the frame's axes is
:math:`T K_\mathrm{local} T^T`, with :math:`T` four copies of the direction cosine
matrix on its diagonal.

Sections
--------

This section gives the beam properties of each of SubDyn's three cross-section
tables, as SubDyn's ``SetElementProperties`` derives them.

.. list-table::
   :header-rows: 1
   :widths: 15 55 30

   * - Shape
     - Properties
     - Timoshenko shear coefficient
   * - Circular (1c)
     - :math:`D, t` (m): :math:`A`, :math:`I_{xx} = I_{yy}`, :math:`J_t = J_0 = 2 I_{xx}`; :math:`t = 0` is solid
     - Steinboeck et al. (SubDyn theory manual, eq. 13) with the diameter ratio :math:`(D-2t)/D`, so 1 for a
       solid section, as SubDyn takes it
   * - Rectangular (1r)
     - :math:`S_a` along local :math:`x`, :math:`S_b` along local :math:`y`, :math:`t` (m):
       :math:`I_{xx} = S_a S_b^3/12`, :math:`I_{yy} = S_a^3 S_b/12`, less the hollow; :math:`J_t` for a solid or
       a thin-walled section
     - solid: :math:`10(1+\nu)/(12+11\nu)`; hollow: SubDyn's thin-wall formula; :math:`\nu = E/(2G) - 1`
   * - Arbitrary (4)
     - :math:`A, A_{sx}, A_{sy}, I_{xx}, I_{yy}, J_0, J_t` as given
     - :math:`\kappa_x = A_{sx}/A`, :math:`\kappa_y = A_{sy}/A`

:math:`J_0` (the polar moment) enters only the rotary inertia; the stiffness uses
:math:`J_t`. A lattice tower's angle members are best given as arbitrary sections,
with their principal axes set by the spin.

Supports and gravity
--------------------

This section describes the boundary conditions and the loads from weight.

A support (SubDyn's base reaction joint) fixes any of its joint's six degrees
of freedom (flag 1) and leaves the others free (flag 0). The free ones may stand
on a symmetric 6 x 6 spring read from the row's SSI file: the 21 entries of its
upper triangle, column by column, named ``Kxx Kxy Kyy Kxz Kyz Kzz Kxtx Kytx Kztx
Ktxtx Kxty Kyty Kzty Ktxty Ktyty Kxtz Kytz Kztz Ktxtz Ktytz Ktztz`` (N/m, N/rad,
N m/rad), each line holding a value then its name.

The weight of an element, :math:`w = \rho A g` per length (N/m), is applied as
SubDyn's ``ElemG`` applies it: :math:`-wL/2` along :math:`z` at each node and the end
moments of a uniformly loaded fixed-fixed beam, :math:`\pm (L^2/12)\,\hat e \times \mathbf q`
with :math:`\hat e` the element's axis and :math:`\mathbf q = (0, 0, -w)`: of magnitude
:math:`w L^2 \sin\alpha / 12` for an element at :math:`\alpha` to the vertical (the
part of the weight across it), about a horizontal axis. With these consistent loads the nodal displacements of a
uniformly loaded beam are exact. A concentrated mass adds its weight at its joint
and the moment of that weight about the joint when its centre is offset.

Solution
--------

This section states what a static solve returns and in which frame.

The stiffness of the free degrees of freedom, springs included, is factored once
by a dense Cholesky factorisation. A degree of freedom that nothing holds stops
the factorisation; the frame is refused with a message naming that degree of
freedom and its joint (for example "the rotation about z of joint 2").
A solve returns:

- the displacement of every degree of freedom (m, rad), in the frame's axes;
- each support's reaction, the force and moment on the frame by the ground at
  its joint (N, N m), in the frame's axes, springs included;
- each element's end forces, the forces and moments on the element at its first
  node and then at its second (N, N m), in the element's local axes: the axial
  force (tension positive) is the local :math:`z` force at the second node.

The model is linear: no geometric (stress) stiffness and no large rotations,
as in SubDyn's beams.

Mass
----

This section gives the mass matrix the modes and the time response use.

An element's mass is SubDyn's ``ElemM_Beam``: the consistent mass of the cubic
beam, :math:`\rho A L` in translation, the rotary inertia of the cross-section
(:math:`\rho I_{xx}`, :math:`\rho I_{yy}`) and the torsional inertia :math:`\rho J_0 L`,
rotated to the frame's axes as the stiffness is. A concentrated mass adds
SubDyn's rigid-body matrix at its joint: the mass :math:`m` at the offset
:math:`\mathbf r = (x, y, z)` of its centre from the joint, with the inertia tensor
entries :math:`J_{xx}, J_{yy}, J_{zz}, J_{xy}, J_{xz}, J_{yz}` about its centre,

.. math::

   M_{66} = \begin{pmatrix} m I & -m [\mathbf r]_\times \\ m [\mathbf r]_\times & J + m (|\mathbf r|^2 I - \mathbf r \mathbf r^T) \end{pmatrix},

with :math:`[\mathbf r]_\times` the cross-product matrix of :math:`\mathbf r`. A support's SSI
file adds its 6 x 6 mass (entries ``Mxx`` ... ``Mtztz``, in the order of the
stiffness entries) to the support's joint.

Natural modes
-------------

This section describes how the lowest natural modes are found and checked.

``frame_modes`` finds the lowest :math:`p` solutions of
:math:`K \boldsymbol\phi = \omega^2 M \boldsymbol\phi` on the free degrees of freedom by
subspace iteration with :math:`q = \max(2p, p + 8)` trial vectors (Bathe): each
iteration solves :math:`K Y = M X` with the factored stiffness and solves the
projected :math:`q \times q` problem by Jacobi rotations, through the projected
stiffness so that a degree of freedom without mass gives an infinite frequency
instead of a singular matrix. The iteration stops when the first :math:`p + 1`
eigenvalues change by less than :math:`10^{-12}` of themselves. A Sturm sequence
count then checks that no mode was missed: the number of negative pivots of
:math:`K - \sigma M`, :math:`\sigma` between the :math:`p`-th and the next
eigenvalue, must be :math:`p`. The modes are returned with their frequencies (Hz)
and their shapes, mass-normalised (:math:`\boldsymbol\phi^T M \boldsymbol\phi = 1`) on
every degree of freedom (0 on the fixed ones).

Time response
-------------

This section gives the time integration and its properties.

``FrameDynamics`` advances :math:`M \ddot{\mathbf u} + C \dot{\mathbf u} + K \mathbf u = \mathbf f(t)`
by Newmark's average-acceleration method (:math:`\beta = 1/4`, :math:`\gamma = 1/2`)
with Rayleigh damping :math:`C = a_0 M + a_1 K`. Over a step :math:`h` (s), with the
loads at the end of the step,

.. math::

   \left(K + \tfrac{2}{h} C + \tfrac{4}{h^2} M\right) \mathbf u_{n+1} =
   \mathbf f_{n+1} + M\left(\tfrac{4}{h^2} \mathbf u_n + \tfrac{4}{h} \dot{\mathbf u}_n + \ddot{\mathbf u}_n\right)
   + C\left(\tfrac{2}{h} \mathbf u_n + \dot{\mathbf u}_n\right),

   \ddot{\mathbf u}_{n+1} = \tfrac{4}{h^2} (\mathbf u_{n+1} - \mathbf u_n) - \tfrac{4}{h} \dot{\mathbf u}_n - \ddot{\mathbf u}_n, \qquad
   \dot{\mathbf u}_{n+1} = \dot{\mathbf u}_n + \tfrac{h}{2} (\ddot{\mathbf u}_n + \ddot{\mathbf u}_{n+1}).

The effective stiffness on the left is factored once per step size. The method
is unconditionally stable and second-order accurate and adds no numerical
damping: an undamped free vibration keeps its energy
:math:`\tfrac12 \dot{\mathbf u}^T M \dot{\mathbf u} + \tfrac12 \mathbf u^T K \mathbf u` exactly, and a
mode of angular frequency :math:`\omega` advances by the phase
:math:`2 \arctan(\omega h / 2)` per step (a period lengthened by about
:math:`(\omega h)^2/12`). ``rayleigh_coefficients`` gives :math:`a_0` (1/s) and
:math:`a_1` (s) for damping ratios at two frequencies; a mode of angular frequency
:math:`\omega` then has the ratio :math:`(a_0/\omega + a_1 \omega)/2`. A run starts at rest
in static equilibrium, or from a given displacement and velocity (the
acceleration then follows from the equation of motion), and its state (the
time, displacement, velocity and acceleration) can be saved and restored.

Coupling to the conductor lines
-------------------------------

This section describes how a tower in a run stands on the frame
(``ERF_FrameTower.H``), behind the same interface as the one-mode tower.

- Axes: the frame file's coordinates are tower-local, origin at the centre of
  the tower's base on the terrain, :math:`x` along the line
  (:math:`(a_y, -a_x, 0)` for the cross-arm direction :math:`\mathbf a`), :math:`y`
  along the cross-arm, :math:`z` up. Every load and motion is turned between
  these axes and ERF's.
- Links: each drag node of the tower (the equivalent lattice that takes the
  wind, as for any tower) and each line attachment is tied rigidly to its four
  nearest frame nodes, or to one node it lies on, never to the crossings of
  diagonals, so a panel's wind acts at its leg joints: on a generated frame
  every joint but the crossings, on a frame from a file every joint a member
  within 20 degrees of the vertical or the horizontal meets (a leg, a chord or
  a strut). The frame has one interface joint, its cross-arm's centre, which must
  stand within half the cross-arm's face depth of the tower's cross-arm height:
  a frame built for another cross-arm height is refused (a line sharing the
  tower may hang at any height on it). A force :math:`\mathbf F` at the
  point :math:`\mathbf p` is shared as
  :math:`\mathbf F_i = \mathbf F/n + \mathbf w \times \mathbf d_i`, with :math:`\mathbf d_i` the
  node's offset from the nodes' centroid :math:`\mathbf c`,
  :math:`\mathbf w = I^+ ((\mathbf p - \mathbf c) \times \mathbf F)`,
  :math:`I = \sum (|\mathbf d_i|^2 1 - \mathbf d_i \mathbf d_i^T)` and :math:`^+` the
  pseudo-inverse: the force and its moment are kept. The point moves as the
  adjoint, :math:`\mathbf u_p = \bar{\mathbf u} + \boldsymbol\theta \times (\mathbf p - \mathbf c)`,
  :math:`\boldsymbol\theta = I^+ \sum \mathbf d_i \times (\mathbf u_i - \bar{\mathbf u})`, so
  the links do no work. Every drag node and attachment must lie within the
  type's ``base_width`` of a frame node, and the frame must stand on four
  supports at :math:`z = 0`, one per quadrant of the base, which give the four
  legs.
- Motion: the frame moves about its static equilibrium under its own weight,
  so the cross-arm starts where the lines were built. Each coupling step
  advances it under the drag and the lines' mean pull over that step, both
  held over it, by the implicit midpoint rule (Newmark's method started from
  the acceleration the held loads give). The frame then takes from the lines
  the work they give up as MoorDyn moves the cross-arm at a constant velocity
  (to the coupling tolerance), so the exchange gains no energy from motions the
  coupling step does not resolve, such as the frame's stiff local modes under
  the insulator strings' axial modes (about 45 Hz). The footings' reactions are those under the held loads. The coupling resolves the frame's
  first natural frequency, and Newmark's method integrates the higher modes
  stably. The members' drag is held from the start of each coupling step, so the frame's
  ``damping_ratio`` must be positive (its modes above about 6 times the first
  frequency would otherwise grow under the drag). The type's ``damping_ratio`` sets
  the Rayleigh damping at the first natural frequency and at ten times it;
  between the two the damping ratio is lower, down to 0.575 times it at
  :math:`\sqrt{10}` times the first frequency, where the second sway, torsion
  and cross-arm modes often lie.
- Footings: the reactions at the four supports are those under the weight
  plus :math:`K\mathbf u + a_1 K \dot{\mathbf u} + M\ddot{\mathbf u} - \mathbf f` at their
  nodes, with :math:`K` and :math:`M` those of the members and concentrated masses only
  (the mass-proportional damping excluded). A support spring's force is part of
  its support's reaction, as in the static solve; at a support on springs the
  reaction so found also holds the springs' share of the stiffness-proportional
  damping, the support mass's inertia and the mass-proportional damping at that
  node. Their resultant gives the base shear and
  the overturning moment, and each support's upward reaction is its leg's
  compression. Before the first step the footings take the static reactions
  under the tower's present loads.
- Restart: a frame's state for a checkpoint (``state()`` and ``set_state()``)
  is its Newmark state, its last loads and its members' temperatures;
  ``restart_mismatch()`` refuses a restart whose ``steel_temperature`` gives
  other temperatures, naming the tower. A frame checkpointed before its first
  step restarts with the static footings, as the run that wrote it had them.

Generated lattice towers
------------------------

This section describes the frame ERF builds from a tower type's dimensions when
the type gives ``erf.conductors.<type>.frame_panels`` (``ERF_LatticeFrame.H``), in
the tower-local axes of the previous section. Every member is an equal-leg
angle, one Euler-Bernoulli element per member, steel with
:math:`E = 200` GPa, :math:`G = 77` GPa, :math:`\rho = 7850` kg/m\ :sup:`3` and the
type's ``yield_strength`` (345 MPa by default).

- Shaft: four legs from the base corners, :math:`(\pm b/2, \pm b/2, 0)` with
  :math:`b` the ``base_width``, tapering linearly to ``top_width`` at the
  cross-arm's height :math:`H` (the height the tower's line hangs at); the peak,
  when ``peak`` > 0, continues at ``top_width`` to :math:`H` + ``peak``. Levels:
  ``frame_panels`` equal panels from the ground to the cross-arm's bottom,
  :math:`H` minus its depth (``arm_depth``, or ``top_width`` when that is 0), then
  the cross-arm's top, then the peak's, in equal panels no taller than the
  shaft's. Horizontal struts join the corners of every level above the ground.
- Faces: ``bracing = crossed`` (the default) puts two diagonals on every panel
  face, bolted where they cross, so each half is a member; the struts are then
  redundant members. ``bracing = single`` puts one diagonal per face,
  alternating in direction from panel to panel.
- Cross-arm: on each side a truss of four chords from the shaft's corners at the
  cross-arm's top and bottom to its tip, :math:`(0, \pm L_a/2, H)` with :math:`L_a` the
  ``arm_length``, in panels about as long as the cross-arm is deep, with a frame of four struts and
  a diagonal at each inner station and a diagonal on each face of every panel
  but the last. The tower's drag model takes the same shape: its cross-arm's
  face is ``arm_depth`` deep over the shaft and tapers to the tips
  (``ERF_Tower.H``), so the drag nodes tied to this truss carry the
  area of its side view.
- Loads: the tower's drag and its lines' pull are tied to every joint but the
  crossings of diagonals, which only the diagonals' bending holds out of their
  face; a panel's wind acts at its leg joints, as lattice tower analyses apply it.
- Hanger: a joint at the cross-arm's centre, :math:`(0, 0, H)`, where a line that
  hangs at the tower's centre is attached, joined to the four corners at the
  cross-arm's top and the four at its bottom; it is SubDyn's interface joint.
- Supports: the four base corners, fixed.
- Sizes: the legs and the cross-arm's chords are ``leg_angle`` (leg width and
  thickness, m), every other member ``brace_angle``. An angle :math:`b \times t`
  without its root fillet has the area :math:`t(2b - t)`, the principal second
  moments of its two legs as rectangles meeting at the heel, the torsion constant
  :math:`(2b - t)t^3/3` and the least radius of gyration :math:`r_v`. The major
  principal axis is the angle's axis of symmetry, and each angle is turned
  (``MSpin``) so that it lies as on a tower: a leg's, or a cross-arm chord's,
  towards the corner it stands on, between its two faces' normals, and a face
  member's at 45 degrees to its face, one leg flat on the face and the other
  standing out of it, the flat leg towards the face's line of symmetry (or down,
  or out along the cross-arm, for a member centred on it), so that mirror images
  in the tower's planes of symmetry are turned alike (``angle_principal_axes``,
  the default). The faces are the
  true planes: the shaft's lean in with its taper, and the cross-arm's bottom
  and sides slope to its tip. The hangers and the stations' inner diagonals keep
  ``MSpin`` 0, as does every member with ``angle_principal_axes = false``, whose
  major axes then lie as SubDyn's default orientation puts them (for a tapering
  leg, horizontal and across its corner, so the weak axis points to the corner).
- Design data: the legs and chords are legs (bolted in both faces); every other
  member is bracing, or redundant for the struts between crossed diagonals,
  bolted by one leg with a normal framing eccentricity at both ends and no
  rotational restraint. The net area is the gross area (``NetArea`` 1): the
  tension checks take no bolt holes, which a member file can give (a typical
  bolted angle has 0.85).

One frame is built per tower type and cross-arm height (to the millimetre), for
the first tower of that height, and shared by the others, which share its first
natural frequency too. The frame's stiffness and mass are dense matrices held on
every rank, so a frame (generated or read) may have at most 6000 free degrees of
freedom, about 120 panels with crossed bracing (48 free degrees of freedom a
panel; single bracing, 24 a panel, stays below it up to the 200 panels
``frame_panels`` allows); a larger one is refused. It is written to
``<diagnostics_dir>/frame_<tower>.dat`` as a SubDyn input file, which SubDyn's
driver reads to the same stiffness, and its design data to
``frame_<tower>_members.dat`` in the member file layout of the next section.

Member checks
-------------

This section describes the strength check of every member
(``ERF_MemberChecks.H``), following ASCE 10-15 (Design of Latticed Steel
Transmission Structures) for angle members. A frame from ``frame_file`` is
checked when the type also gives ``member_file``; a generated frame always is.

- Forces: a member's elements' end forces in their local axes, at present: those
  under the frame's weight plus :math:`K_e(\mathbf u + a_1 \dot{\mathbf u})` of its
  motion (the elastic and the stiffness-proportional damping forces), or the
  static forces under the tower's present loads before the first step. The
  member's tension :math:`T` and compression :math:`C` are the largest along it.
- Slenderness: :math:`L/r` with :math:`L` the member's length joint to joint and
  :math:`r` the angle's :math:`r_v` (or, for a member that is not an angle, the
  least of its section's). Legs take :math:`KL/r = L/r`. Bracing and redundant
  members take, for :math:`L/r \le 120`, :math:`L/r` (concentric at both ends,
  curve 1), :math:`30 + 0.75 L/r` (eccentric at one end, curve 2) or
  :math:`60 + 0.5 L/r` (eccentric at both, curve 3); beyond 120, :math:`L/r`
  (unrestrained, curve 4), :math:`28.6 + 0.762 L/r` (partly restrained at one
  end, curve 5) or :math:`46.2 + 0.615 L/r` (at both, curve 6). The limits are
  :math:`L/r \le 150` for legs and :math:`KL/r \le 200` for bracing and 250 for
  redundant members; a member over its limit is flagged.
- Local buckling: an angle leg's :math:`w/t = (b - t)/t` against
  :math:`(w/t)_1 = 80\sqrt{E/F_y}/\sqrt{29000}` and
  :math:`(w/t)_2 = 144\sqrt{E/F_y}/\sqrt{29000}` (ASCE 10-15's
  :math:`80\psi/\sqrt{F_y}` and :math:`144\psi/\sqrt{F_y}`, written for any units and
  any :math:`E`, equal to them at :math:`E` = 29000 ksi): the yield stress stands up to
  :math:`(w/t)_1`, then :math:`F_{cr} = [1.677 - 0.677 (w/t)/(w/t)_1] F_y` up to
  :math:`(w/t)_2`, then :math:`F_{cr} = 0.0332 \pi^2 E/(w/t)^2`. :math:`w/t > 25` is flagged.
- Compression: :math:`F_a = [1 - (KL/r)^2/(2 C_c^2)] F_{cr}` up to
  :math:`C_c = \pi\sqrt{2E/F_{cr}}`, then :math:`\pi^2 E/(KL/r)^2`, on the gross area.
- Tension: :math:`F_y` on the net area (``NetArea`` times the gross area), or
  :math:`0.9 F_y` for an angle bolted by one leg.
- Utilisation: :math:`\max(T/(F_t A_n), C/(F_a A))`, 1 at the design strength. The
  strengths are nominal and the loads those computed, without load factors.
  Bending of the members, which the frame carries at its rigid joints, is not
  checked: ASCE 10-15's effective slenderness accounts for the end
  eccentricities of bolted angles.

A member file has one row per frame member, blank lines and lines starting with
``#`` or ``!`` skipped::

    # MemberID  Role  Fy(Pa)  b(m)  t(m)  NetArea(-)  Bolted  Ends  Restraint
    1  leg      3.45e8  0.15  0.012  0.85  both  concentric  none
    61 bracing  3.45e8  0.09  0.007  0.85  one   both        none

``Role`` is ``leg``, ``bracing`` or ``redundant``; ``Bolted`` is ``one`` or
``both`` (the angle's legs bolted at its ends); ``Ends`` is ``concentric``,
``one`` or ``both`` (the ends with a normal framing eccentricity, curves 1 to 3);
``Restraint`` is ``none``, ``one`` or ``both`` (the ends partly restrained
against rotation, curves 4 to 6). ``b = t = 0`` marks a member that is not an
angle. An angle's area must be within 10 % of its frame section's, so that the
strength and the stiffness describe the same member.

Steel at temperature
--------------------

This section describes how the steel's temperature softens and weakens the
frame, by EN 1993-1-2 (Eurocode 3, structural fire design) Table 3.1, linear
between its rows:

==================  ==  ===  ===  ===  ====  ====  ====  ====  ======  =====  ======  ====
:math:`\theta` (C)  20  200  300  400  500   600   700   800   900     1000   1100    1200
:math:`k_y`         1   1    1    1    0.78  0.47  0.23  0.11  0.06    0.04   0.02    0
:math:`k_E`         1   0.9  0.8  0.7  0.6   0.31  0.13  0.09  0.0675  0.045  0.0225  0
==================  ==  ===  ===  ===  ====  ====  ====  ====  ======  =====  ======  ====

(both 1 at 100 C and below). A member at :math:`\theta` has :math:`k_E E` and
:math:`k_E G` in its stiffness, and :math:`k_y F_y` and :math:`k_E E` in its checks; its
mass is unchanged, and thermal expansion is not modelled. Temperatures must stay
below 1200 C, where no stiffness is left.

``erf.conductors.<type>.steel_temperature`` sets every member of the type's
frames at one temperature for the run. A frame tower can also be heated member by
member while it moves (``FrameTower::set_temperature()``, the entry point for a
fire model): the frame is rebuilt and re-factored, the tower carries on from the
same position and velocity, and its static equilibrium under the weight, its
support reactions under it, its first natural frequency and its Newmark
factorisation follow the heated frame's stiffness; the Rayleigh coefficients stay those of
the frame at its first temperature.

Outputs
-------

This section lists the files a frame tower's driver writes to its diagnostics
directory from the frame's results (the conductor line model, with its keys
``erf.conductors.diagnostics_dir``, ``node_output_int`` and ``stats_start``;
the frame itself writes nothing).

- ``towers.dat``: for a tower whose members are checked, ``<tower>_utilisation``
  and ``<tower>_member``, the largest utilisation and the member it is in.
- ``tower_<tower>_stats.csv``: ``max_utilisation`` among the quantities.
- ``tower_<tower>_members_stats.csv``: per member, ``m<id>_axial`` (its tension, or
  minus its compression, whichever is larger; N) and ``m<id>_utilisation``, with
  their mean, rms, minimum and maximum since ``stats_start``.
- ``tower_<tower>_members.csv``: per member its role, joints, temperature, length,
  :math:`r`, :math:`L/r`, :math:`KL/r` and its limit, :math:`w/t`, :math:`F_y` and :math:`E` at
  its temperature, :math:`F_a`, its tensile and compressive strengths and the two
  flags.
- ``tower_<tower>_frame.dat``, every ``node_output_int`` steps: the frame nodes'
  displacements (frame axes, from where the frame stands under its weight at its
  first temperature) and the members' utilisations; its header lists the nodes'
  positions and the members' ids.

Reading SubDyn input files
--------------------------

This section lists the part of a SubDyn input file (OpenFAST 5.0 layout) the
reader takes, and what it refuses.

=========================================== ===============================================================
SubDyn section                              Read as
=========================================== ===============================================================
FEA and Craig-Bampton parameters            ``FEMMod`` 1 or 3 (2 and 4 refused), ``NDiv`` :math:`\ge 1`
Structure joints                            id and position; ``JointType`` must be 1 (cantilever)
Base reaction joints                        six flags (1 fixed, 0 free) and an optional SSI file, relative
                                            to the SubDyn file's directory
Interface joints                            one id: the cross-arm's centre, where the lines are tied
Members                                     ``MType`` 1/1c, 1r or 4 with the same property set at both
                                            ends and ``MSpin`` (deg); cables (2), rigid links (3) and
                                            springs (5) are refused
Circular, rectangular, arbitrary properties every row
Cable, rigid link, spring, cosine matrices  skipped
Joint additional concentrated masses        mass, inertia and the centre's offset
=========================================== ===============================================================

``write_subdyn()`` writes a frame back in the same layout, every value to 17
significant digits (and a support spring to an SSI file beside it), so that this
reader and SubDyn read the same frame; it writes no temperatures.

Every check names the file, and the line or the item: an unknown joint, a
property set missing from its member type's table, a zero-length member, a
joint on no member, a non-positive modulus or section value, a support listed
twice, an SSI file on a support with no free degree of freedom.

Verification
------------

This section lists what the unit tests ``DirectionCosines``, ``BeamElement``,
``Frame``, ``BeamMass``, ``FrameModes``, ``FrameDynamics``, ``RigidBodyMass``,
``FrameSubDyn``, ``FrameTower``, ``RigidLink``, ``LatticeFrame``,
``MemberChecks`` and ``HeatedFrame`` (``Tests/Unit/MovingBodies``) compare with.

- Beam theory, exact for these elements at the nodes, with one and five elements
  per member: a cantilever's tip deflection :math:`PL^3/(3EI)` (plus
  :math:`PL/(\kappa G A)` with shear deformation) in both bending planes, its tip
  rotation, :math:`PL/(EA)` and :math:`TL/(G J_t)`; the sag of a cantilever under its
  own weight, :math:`wL^4/(8EI)` (plus :math:`wL^2/(2\kappa G A)`); a vertical member's
  sink under its weight and a mass, :math:`\rho g L^2/(2E) + m g L/(EA)`.
- A spring base: the tip moves by the beam's bending plus the base's translation
  and tilt from the spring's own coupled solve.
- A fixed portal frame's sway, :math:`P h^2 (2a + 3b)/(12 a (a + 6b))` with
  :math:`a = E I_c/h` and :math:`b = E I_b/W` (slope-deflection).
- A triangulated truss of slender members: the pin-jointed member forces and
  reactions.
- Equilibrium: on a lattice tower with a spring base, the reactions balance
  arbitrary joint loads and gravity in force and in moment.
- SubDyn stiffness: three SubDyn input files (``Tests/test_files/FrameSubDynTower``,
  cases A and B of a 25 m lattice tower and case T of a 30 m one) read by ERF's
  reader give the stiffness condensed to the tower's interface joint, a 6 x 6
  matrix, that SubDyn writes as ``KBBt`` (OpenFAST 5.0.0, double precision), to
  the 7 digits SubDyn prints. Cases A and T have Euler-Bernoulli arbitrary
  sections; case B Timoshenko circular, rectangular and spun arbitrary sections,
  two elements per member and a base joint on a coupled spring.
- Mass: an element moves :math:`\rho A L` in every translation and :math:`\rho J_0 L` in a
  twist; a concentrated mass's 6 x 6 matrix gives the kinetic energy of the
  offset rigid body for any motion.
- Modes: a cantilever's bending (two modes in each plane), axial and torsional
  frequencies from beam theory; mass orthonormality and the eigenproblem's
  residual; a stiff post on a spring base with masses ringing at the base's own
  coupled frequencies; the Sturm count of modes below a shift.
- SubDyn modes: the 20 lowest natural frequencies of the two cases above and of a
  third with concentrated masses at the cross-arm tips (offset centres and
  products of inertia) and a mass on the spring base equal those SubDyn writes
  as ``Full_frequencies``, to its 7 digits; case B's cluster of eight modes within
  0.25 % of each other included.
- Coupling: a link keeps a force and its moment and does no work; a tower of
  case T (``Tests/test_files/FrameSubDynTower/caseT``, a 30 m lattice tower with a 12 m cross-arm) turned 30 degrees rings at SubDyn's first frequency
  (4.219267 Hz, 6800.110 kg), settles under steady drag and a line's pull where
  the frame's static solution of the linked loads puts it, and its footings
  carry the applied loads, their moment and the weight, the downwind legs the
  more compressed; a frame whose cross-arm is far from the lines' attachment is
  refused.
- Generated towers: a generated 30 m tower with a peak (case G,
  ``Tests/test_files/FrameSubDynTower/caseG``, written by ``write_subdyn()`` with
  ``angle_principal_axes = false``) has SubDyn's ``KBBt`` at the cross-arm's
  centre and its 10 lowest natural frequencies, to SubDyn's 7 digits (SubDyn's
  ``MSpin`` sense itself is checked on case B's turned diagonals); every
  generated angle lies on its principal axes as described, on the shaft and on
  the cross-arms; its geometry (corners, taper, crossings,
  tips) and member counts are those described; under a lateral load at the
  cross-arm, the members cut by any horizontal plane balance the load above it,
  the legs carrying more than three quarters of the overturning moment; a written file
  reads back to the same frame.
- Member checks: an angle's properties against the two rectangles by hand and a
  tabulated L100x100x10; the six slenderness curves and their meeting at 120; the
  local-buckling and column curves continuous at their limits; EN 1993-1-2's
  table; tensile and compressive strengths, flags and heated members; a column's
  axial force; the member file's round trip and its refusals. A frame tower's
  checks equal those of the frame's static solution, at rest before its first
  step and settled under steady loads.
- Steel at temperature: a frame at a uniform 600 C moves :math:`1/k_E` as far under
  the same loads and gravity, with the same member forces; a tower heated while
  settled keeps its place, rings at :math:`\sqrt{k_E}` times its frequency and
  settles where the heated frame's static solution puts it; a saved state
  carries the members' temperatures.
- Time response: a cantilever released in its first mode follows
  :math:`\cos(2 n \arctan(\omega h/2))` and keeps its energy to :math:`10^{-11}`; with Rayleigh
  damping the mode follows Newmark's recursion for the damped oscillator of its
  frequency and ratio; a damped tower under a sudden load and gravity settles at
  the static solution; a saved and restored state continues bit for bit.
