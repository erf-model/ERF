.. _sec:IBSEB:

Surface energy balance on immersed-boundary faces
=================================================

Buildings represented by immersed forcing (:cpp:`erf.buildings_type =
ImmersedForcing`) are resolved: every cell is either fluid or solid, and the
faces between them are the walls and roofs. The immersed-boundary surface
energy balance gives each such face its own skin temperature from

.. math::

   C \frac{dT_s}{dt} = SW_{abs} + LW_{net} - H - LE - G ,

with the absorbed shortwave, the net longwave, the sensible and latent heat
to the air, and the conduction into the wall. It is built in phases; the
sections below grow with them, and ``Source/ImmersedBoundarySEB/
IBSEB_DEVELOPMENT.md`` carries the plan and the phase log.

Face storage
------------

A wall face lies between a fluid cell, where the blanking is below one
half, and a solid cell, where it is at or above one half. Each face is
stored once, on the rank that owns the fluid cell, as a compact list rather
than a face-centred field: the fluid cell, the face direction, the side the
solid is on, the building id, the material id, the area, the skin
temperature, the slab temperatures, the view fractions and the current
fluxes, all in device arrays, contiguous per local box. Every input the
balance needs is in the fluid cell next to the face, so the per-step update
is one kernel over the list with no communication; only the reports reduce
across ranks.

Building ids come from labelling the solid columns of the blanking, four
connected in the horizontal and numbered in scan order; a face carries the
id of its solid neighbour's column. The per-building report,
:cpp:`erf.ibseb.csv_file`, lists every :cpp:`erf.ibseb.csv_int` steps the
number of faces, the area and the area-weighted mean skin temperature of
each building, and a ``[IBSEB]`` line prints the face counts per
direction, the number of buildings, the total area and the skin
temperature range.

For output the list is scattered into cell-centred fields: ``ibseb_nfaces``
and ``ibseb_tskin`` in the plotfile, and ``IBSEBState`` in the checkpoint,
which holds the skin and slab temperatures of up to six faces per fluid
cell. On restart the list is rebuilt from the blanking and refilled from
that field, so a restart does not depend on the number of ranks.
``Exec/CanonicalTests/SEB/Phase1_Storage`` checks the face counts against the
mask, the rank independence and the checkpoint round trip.

Shortwave and shadow
--------------------

The downwelling radiation reaches the balance through a provider,
:cpp:`erf.ibseb.radiation`. The built-in ``prescribed`` provider gives the
direct-normal irradiance, the diffuse irradiance on a horizontal surface and
the sun vector either as fixed inputs (:cpp:`erf.ibseb.sun_mode = fixed`,
for analytic tests) or from the site and time (``solar``): declination,
equation of time, hour angle, zenith and azimuth by the Spencer series, the
direct beam by the Bird transmission :math:`S_0 E_0 \tau^{1/\cos z}` and the
diffuse light as a fraction :math:`k_d` of what the beam lost on the way
down. The sun vector points from a surface toward the sun, with the azimuth
clockwise from north, so east is :math:`+x` and north :math:`+y`.

On a face with outward normal :math:`\mathbf{n}` the direct beam is
:math:`I_{dn}\, \max(0, \mathbf{n}\cdot\mathbf{s})` unless the face is in
shadow. Shadow is decided by a ray cast from the face centre toward the sun
against the height of every column of the level, a two-dimensional walk over
the columns the ray crosses: buildings stand on the ground and the ray only
rises, so the ray is blocked wherever its height on entering a column is
below that column's top. The column tops are a small array replicated on
every rank, so the test needs no communication and costs a few operations
per column crossed. The diffuse light on a face is the sky view fraction
times the horizontal diffuse plus the ground view fraction times the ground
albedo times the total horizontal irradiance, with the fractions from the
hemisphere sampling below. The absorbed shortwave is one minus the face
albedo times the sum.

View fractions and longwave
---------------------------

Once, at initialisation, every face samples a cosine-weighted hemisphere
around its outward normal: :cpp:`erf.ibseb.view_n_az` azimuths by
:cpp:`erf.ibseb.view_n_el` elevations, stratified with
:math:`\theta = \arcsin\sqrt{u}` from the normal so that every ray carries
the same weight and the counts are view factors. Each ray goes through the
same column walk as the shadow test, made direction-aware: a rising ray is
blocked by a column whose top is above its entry height and otherwise
reaches the sky; a falling ray is blocked by a column it descends into and
otherwise reaches the ground. The fractions of rays ending on the sky, the
ground and a building, :math:`f_{sky}`, :math:`f_{ground}` and
:math:`f_{bldg}`, sum to one. A roof sees no ground; a wall on open ground
sees half sky and half ground; a wall facing another building sees it.

The longwave arriving at a face is

.. math::

   LW_{in} = f_{sky}\, LW_{sky} + f_{ground}\, \varepsilon_g \sigma T_g^4
           + f_{bldg}\, \sigma T_s^4 ,

the sky term either the input :cpp:`erf.ibseb.lw_down` or, with
:cpp:`erf.ibseb.lw_mode = gray`, :math:`\varepsilon_{sky} \sigma T_{air}^4`
with the air temperature of the face's fluid cell; the ground term at the
input ground temperature; and the building term at the face's own skin
temperature, the isothermal-surroundings approximation under which a face
and the walls it sees exchange no net longwave. There are no face-to-face
view factors and no radiosity; the fractions are stored so a radiosity pass
can be added later without touching the balance. The net longwave,
:math:`\varepsilon (LW_{in} - \sigma T_s^4)`, is positive into the face.

``Exec/CanonicalTests/SEB/Phase3_Longwave`` checks the fractions of every
face against an independent hemisphere sampling, their closure and the
analytic values on the clean planes, and the longwave formulas on every
face, on one and four ranks.

Sensible and latent heat
------------------------

The sensible heat leaves a face through a wall function between its skin
temperature and the adjacent fluid cell, the same form the immersed forcing
uses for momentum. The tangential wind of the fluid cell, at half a cell
from the wall, gives the friction velocity with the roughness
:cpp:`erf.ibseb.z0_wall`, :math:`u_* = \kappa U_t / \ln(\delta/z_0)`; the
difference between the skin temperature, expressed as a potential
temperature with the cell's Exner function, and the cell's potential
temperature gives the temperature scale with :cpp:`erf.ibseb.z0h_wall`; and

.. math::

   H = \rho c_p u_* \theta_* = \rho c_p u_* \kappa
       \frac{\theta_s - \theta_{air}}{\ln(\delta/z_{0h})} ,

positive out of the face. Vertical faces are neutral, since buoyancy runs
along them; on roofs :cpp:`erf.ibseb.stability_correction` applies the
surface layer's similarity functions with a few fixed-point passes on the
Obukhov length. The latent flux is zero in this version and its argument
is in place for a wet-surface option.

The flux enters the atmosphere as a source: every face deposits
:math:`H A / (c_p V \Pi)` into the rho-theta equation of its fluid cell,
with :math:`A` the face area, :math:`V` the cell volume and :math:`\Pi`
the Exner function, so the cell warms by :math:`H A / (c_p V)` per second
in temperature. It is added after the sources are rebuilt at every slow
stage and never overwrites another term. :cpp:`erf.ibseb.couple_heat =
false` diagnoses the flux without applying it. Because the balance now
owns the temperature condition at the buildings, the immersed forcing's
surface-temperature inputs must not be set with it.

``Exec/CanonicalTests/SEB/Phase4_Sensible`` holds a cube's faces at 320 K in
an 8 m/s wind at 300 K and checks the wall function on every face against
the formulas, the rank independence, the heat budget of the closed domain
against the summed face flux, and a mass-inflow, pressure-outflow variant
whose wake is warmer than the inflow.

Conduction into the wall and materials
--------------------------------------

Every face carries a slab of :cpp:`erf.ibseb.n_slab_layers` uniform layers
between its skin and the building interior at :cpp:`erf.ibseb.T_interior`.
The conduction is solved implicitly with the Thomas algorithm, the form of
the SLUCM branch's slab solver with the skin temperature as the top
boundary instead of a flux; the layer centres sit half a layer below the
skin and half a layer above the interior, so both boundary fluxes use
:math:`2k/\Delta z`. The conduction into the slab through the skin,

.. math::

   G = \frac{2k}{\Delta z}\,(T_s - T_0) ,

is positive into the wall. The scheme is unconditionally stable, so thick
walls with a few layers and thin walls with many are both fine.

The conductivity, heat capacity and thickness of the slab, and the albedo
and emissivity of the face, come from a material library when
:cpp:`erf.ibseb.material_file` is given, a CSV in the SLUCM schema so one
file serves both models; :cpp:`erf.ibseb.material_default` applies to every
building and :cpp:`erf.ibseb.material_by_building` gives one id per
building. Without a file the uniform inputs apply to every face.

``Exec/CanonicalTests/SEB/Phase5_Ground`` checks a thick finely layered slab
against the semi-infinite erfc solution, a thin light slab against the
steady linear profile, the materials per building, and the slab through a
checkpoint restart.

The prognostic balance
----------------------

With :cpp:`erf.ibseb.prognostic` (the default once the balance is on) the
skin of every face is massless and its temperature at the end of each
atmospheric step closes the balance

.. math::

   SW_{abs} + \epsilon\,Q_{ext} + \epsilon\,\big[LW_{ext} - (1 - f_b)\,\sigma T_s^4\big]
   - C_H\,(T_s/\Pi - \theta_a) - LE - G(T_s) = 0 ,

the form of the SLUCM branch's facet solver with three changes. The
external incident flux :math:`Q_{ext}` is the hook through which a fire's
radiation will enter; it is absorbed with the longwave emissivity since
such sources are thermal, and :cpp:`erf.ibseb.Q_ext_uniform` sets it on
every face for tests. The wall term of the incoming longwave,
:math:`f_b \sigma T_s^4` under the isothermal-surroundings approximation,
depends on the unknown and folds into the emission as the factor
:math:`(1 - f_b)`. And the conduction is not the lagged
:math:`2k(T_s - T_0^n)/\Delta z` but the flux the implicit slab step will
actually take, :math:`G = a T_s - b`, whose two coefficients come from two
trial slab steps (the step is linear in :math:`T_s`); the balance and the
slab therefore agree to rounding instead of by one step's lag, and the
slab energy changes by exactly :math:`\Delta t\,(G - G_{bottom})`.

The sensible coefficient :math:`C_H = \rho c_p \kappa u_* / \ln(\delta/z_{0h})`
is frozen at the wind of the step, with the Exner function :math:`\Pi` of
the fluid cell turning the skin temperature into a potential temperature;
the latent flux is held at its stored value (zero so far). Newton's method
with the analytic Jacobian, every term of which is negative, converges in
two or three iterations from the previous step's temperature; a step cap
(:cpp:`erf.ibseb.newton_max_step_K`) and the bounds
:cpp:`erf.ibseb.T_skin_min` / :cpp:`erf.ibseb.T_skin_max` guard it. The
bounds are inputs because a face under a fire legitimately exceeds the
380 K of the urban canopy model; a face held at a bound leaves a non-zero
residual, which every face stores and the summary reports as
``resid_max``. After the solve the slab is advanced with the new skin
temperature and the longwave, sensible and conduction fluxes are rewritten
at it, so the heat the air receives during the step, the dumps and the
reports all describe the same closed balance. The atmosphere is seen at
the start of the step and the skin is implicit within it, the coupling of
a land-surface model.

Setting :cpp:`erf.ibseb.prognostic = false` keeps the skin at its initial
or restart value and diagnoses the terms around it, which is how the
regtests of the earlier phases check each term on its own.

``Exec/CanonicalTests/SEB/Phase6_Prognostic`` runs a cube under a fixed sun
with every step dumped and checks the residual on every face, the
consistency of every stored flux with the skin temperature, the slab
energy per step, the closure over the run (radiation in equals heat
convected, stored and conducted out to within the summed residual), and
an independent Python model with its own Newton and a dense slab solve
driven by the dumped forcing, which reproduces the skin temperature of
every face to 1e-9 K. It also exercises the external flux with the bound
raised, a checkpoint restart, and the sun rising over the cube at Boulder
on the solstice, where the east wall warms before the roof.

``Exec/CanonicalTests/SEB/Phase2_Shortwave`` puts a short box 40 m east of a
tall one and checks the shadow flag of every face against an independent
ray cast, the incidence on every orientation, the height to which the tall
box shadows the short one's west wall against
:math:`H - d \tan(\text{elevation})`, the agreement of one and four ranks,
and the solar mode against the solstice-noon zenith at Boulder. The
embedded-boundary reader steps each building edge over one cell, so the
boxes have a full-height core with a half-height rim, which the test
reads from the face dump rather than assuming.
