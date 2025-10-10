
 .. role:: cpp(code)
    :language: c++

.. _sec::surface_layer

Surface Layer Boundaries
------------------------
The surface layer provides an abstraction layer for users to directly specify
diffusive fluxes at a boundary via a multitude of methods. More specifically, the surface layer condition
applies an impenetrable condition for the boundary normal velocity but higher order extrapolation for
all the other variables and then allows a user to specify a method for calculating the diffusive fluxes

::

   erf.surface_layer.flux_type    = STRING    #flux types (donelan, moeng, custom)

The ``donelan`` flux type employs bulk drag coefficients to compute the diffusive stresses while the ``moeng`` type
employs Moeng's formulation for Monin-Obukhov similarity theory (MOST) and ``custom`` allows the user to directly
specify the fluxes through ``ustar; tstar; qstar``. Currently, the MOST pathway is the primary flux type employed
in ERF simulations and will be the focus in subsequent sections.


MOST Theory
~~~~~~~~~~~~~~~~~~~
Monin-Obukhov similarity theory (MOST) is used to describe the atmospheric surface layer (ASL), the lowest part of the atmospheric boundary layer.  The implementation of MOST in ERF follows that in `AMR-Wind <https://github.com/Exawind/amr-wind/>`_, which is based on the surface layer profiles presented in
`P. van der Laan, et al., Wind Energy, 2017 <https://onlinelibrary.wiley.com/doi/10.1002/we.2017>`_ and
`D. Etling, "Modeling the vertical ABL structure", 1999 <https://www.worldscientific.com/doi/abs/10.1142/9789814447164_0003>`_.
MOST theory assumes that the ASL is in a steady state and horizontally homogeneous, and kinematic fluxes due to turbulent transport (:math:`\overline{u^{'}w^{'}}`, :math:`\overline{v^{'}w^{'}}`, and :math:`\overline{\theta^{'}w^{'}}`) are constant with height.
:math:`\Phi_m` and :math:`\Phi_h` are the nondimensional wind shear and temperature gradient, respectively, which are assumed to follow universal similarity laws based on dimensional arguments.
With these assumptions, the MOST theory can be written as:

.. math::

  \overline{u^{'}} \overline{w^{'}} = {\rm const} = -u^{2}_{\star},

  \overline{w^{'}} \overline{\theta^{'}} = {\rm const} = -u_{\star}\theta_{\star},

  \overline{w^{'}} \overline{q_{v}^{'}} = {\rm const} = -u_{\star}q_{v,\star},

  \Phi_{m}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{u}(z)}{\partial z},

  \Phi_{h}(\zeta) = \frac{\kappa z}{u_{\star}} \frac{\partial \overline{\theta}(z)}{\partial z}

where the nondimensional gradients are expressed in terms of the MOST stability parameter, :math:`\zeta = \frac{z}{L} = -\frac{\kappa z}{u_{\star}^{3}} \frac{g}{\overline{\theta}} \overline{w^{'}\theta^{'}}`, which serves as a surface layer scaling parameter.
Here, :math:`L` is the Monin-Obukhov length,
:math:`u_{\star}` is the friction velocity (defined for :math:`u` aligned with the wind direction),
:math:`\theta_{\star}` is the surface layer temperature scale,
:math:`\overline{\theta}` is the reference virtual potential temperature for the ASL,
and :math:`\kappa` is the von Karman constant (taken to be :math:`0.41`).

Integration of the MOST assumption equations give the classical MOST profiles of mean velocity, potential temperature, and water vapor

.. math::

  \overline{u}(z)    &= \frac{u_{\star}}{\kappa} \left[ \mathrm{ln} \left(\frac{z}{z_0}\right) - \Psi_m(\zeta)\right],

  \overline{\theta}(z) - \theta_0 &= \frac{\theta_{\star}}{\kappa} \left[ \mathrm{ln}\left(\frac{z}{z_0}\right) - \Psi_{h}(\zeta) \right]

  \overline{q_{v}}(z) - q_{v,0} &= \frac{q_{v,\star}}{\kappa} \left[ \mathrm{ln}\left(\frac{z}{z_0}\right) - \Psi_{h}(\zeta) \right]


where :math:`\theta_0` and :math:`q_{v,0}` are the surface values and :math:`z_0` is a characteristic roughness height. The integrated similarity functions,

.. math::

  \Psi_{m}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\Phi_{m}(\zeta)]\frac{d\zeta}{\zeta},

  \Psi_{h}(\zeta) &= \int_{0} ^{\frac{z}{L}} [1-\Phi_{h}(\zeta)]\frac{d\zeta}{\zeta}

are calculated analytically from empirical gradient functions :math:`\Phi_m` and :math:`\Phi_h`, which are
defined piecewise for stable and unstable values of the stability parameter.

Unstable: :math:`(\zeta < 0)`

.. math::

  \Phi_{m} &= (1-\gamma_{1}\zeta)^{-\frac{1}{4}}, \quad
  \Psi_{m}    = \mathrm{ln}\left[\frac{1}{8}(1+\Phi_{m}^{-2})(1+\Phi_{m}^{-1})^{2}\right]-2\arctan(\Phi_{m}^{-1})+\frac{\pi}{2},

  \Phi_{h} &= \sigma_{\theta}(1-\gamma_{2}\zeta)^{-\frac{1}{2}}, \quad
  \Psi_{h}    = (1+\sigma_{\theta}) \mathrm{ln} \left[\frac{1}{2}(1+\Phi_{h}^{-1}) \right]+(1-\sigma_{\theta}) {\mathrm{ln}} \left[\frac{1}{2}(-1+\Phi_{h}^{-1})\right]

Stable: :math:`(\zeta > 0)`

.. math::
  \Phi_{m} &= 1+\beta \zeta, \quad \Psi_{m}=-\beta \zeta,

  \Phi_{h} &= \sigma_{\theta}+\beta \zeta, \quad \Psi_{h}=(1-\sigma_{\theta})\mathrm{ln}(\zeta)-\beta \zeta,

where the constants take the values proposed in `Dyer, Boundary Layer Meteorology, 1974
<https://link.springer.com/article/10.1007/BF00240838>`_:

.. math::
  \sigma_{\theta}=1, \quad \beta = 5, \quad \gamma_{1}=16, \quad \gamma_{2}=16

Inverting the equations above, the MOST stability parameter,

.. math::
  \zeta=\frac{z}{L} = -\kappa z \frac{g}{\bar{\theta}} \frac{\theta_{\star}}{u^{2}_{\star}}

is determined by the friction velocity

.. math::
  u_{\star} = \kappa \overline{u}/[\mathrm{ln}(z/z_0)-\Psi_{m}({z}/{L})]

and the characteristic surface layer temperature

.. math::
  \theta_{\star} = \kappa (\overline{\theta}-\theta_0)/[\mathrm{ln}(z / z_0)-\Psi_{h}(z/L)]


MOST Implementation
~~~~~~~~~~~~~~~~~~~

As noted in :ref:`sec:surface_layer`, the boundary conditions for velocity, temperature, and water vapor do not change
with the flux type. Therefore, the MOST implementation in ERF is a specific method for computing the diffusive fluxes, which are directly written into the stress tensor/vector. The MOST pathway is structured to allow either the surface temperature (:math:`\theta_0`) or surface temperature flux (:math:`\overline{w^{'}\theta^{'}}`) to be enforced. To compute the MOST flux, the following algorithm is applied:

#. Horizontal (planar) averages :math:`\bar{u}`, :math:`\bar{v}` and :math:`\overline{\theta}` are computed at a reference height :math:`z_{ref}` assumed to be within the surface layer.

#. Initially, neutral conditions (:math:`L=\infty, \zeta=0`) are assumed and used to compute a provisional :math:`u_{\star}` using the equation given above. If :math:`\theta_0` is specified, the above equation for :math:`\theta_{\star}` is applied and then the surface flux is computed :math:`\overline{w^{'}\theta^{'}} = -u_{\star} \theta_{\star}`. If :math:`\overline{w^{'}\theta^{'}}` is specified, :math:`\theta_{\star}` is computed as :math:`-\overline{w^{'}\theta^{'}}/u_{\star}` and the previous equation is inverted to compute :math:`\theta_0`.

#. The stability parameter :math:`\zeta` is recomputed using the equation given above based on the provisional values of :math:`u_{\star}` and :math:`\theta_{\star}`.

#. The previous two steps are repeated iteratively, sequentially updating the values of :math:`u_{\star}` and :math:`\zeta`, until the change in the value of :math:`u_{\star}` on each iteration falls below a specified tolerance.

#. Once the MOST iterations have converged, and the planar average surface flux values are known, the approach from `Moeng, Journal of the Atmospheric Sciences, 1984 <https://journals.ametsoc.org/view/journals/atsc/41/13/1520-0469_1984_041_2052_alesmf_2_0_co_2.xml>`_ is applied to consistently compute local surface-normal stress/flux values (e.g., :math:`\tau_{xz} = - \rho \overline{u^{'}w^{'}}`):

   .. math::

     \left. \frac{\tau_{xz}}{\rho} \right|_0 &= u_{\star}^{2} \frac{(u - \bar{u})|\mathbf{\bar{u}}| +  \bar{u}\sqrt{u^2 + v^2} }{|\mathbf{\bar{u}}|^2},

     \left. \frac{\tau_{yz}}{\rho}  \right|_0 &= u_{\star}^{2}  \frac{(v - \bar{v})|\mathbf{\bar{u}}| +  \bar{v}\sqrt{u^2 + v^2} }{|\mathbf{\bar{u}}|^2},

     \left.  \frac{\tau_{\theta z}}{\rho} \right|_0  &= \theta_\star u_{\star} \frac{|\mathbf{\bar{u}}| ({\theta} - \overline{\theta}) +
                                                \sqrt{u^2+v^2}  (\overline{\theta} - \theta_0) }{ |\mathbf{\bar{u}}| (\overline{\theta} -\theta_0) } =
                                                u_{\star} \kappa  \frac{|\mathbf{\bar{u}}| ({\theta} - \overline{\theta})  +
                                                \sqrt{u^2+v^2} (\overline{\theta} - \theta_0) }{ |\mathbf{\bar{u}}| [  \mathrm{ln}(z_{ref} / z_0)-\Psi_{h}(z_{ref}/L)] }

     \left.  \frac{\tau_{q_{v} z}}{\rho} \right|_0  &= q_{v,\star} u_{\star} \frac{|\mathbf{\bar{u}}| ({q_{v}} - \overline{q_{v}}) +
                                                \sqrt{u^2+v^2}  (\overline{q_{v}} - q_{v,0}) }{ |\mathbf{\bar{u}}| (\overline{q_{v}} - q_{v,0}) } =
                                                u_{\star} \kappa  \frac{|\mathbf{\bar{u}}| ({q_{v}} - \overline{q_{v}})  +
                                                \sqrt{u^2+v^2} (\overline{q_{v}} - q_{v,0}) }{ |\mathbf{\bar{u}}| [  \mathrm{ln}(z_{ref} / z_0)-\Psi_{h}(z_{ref}/L)] }

   where :math:`\bar{\psi}` are plane averaged values (at :math:`z_{ref}`) and
   :math:`|\mathbf{\bar{u}}|` is the plane averaged magnitude of horizontal velocity (plane averaged wind speed).
   We note a slight variation in the denominator of the velocity terms from the form of the
   equations presented in Moeng. This difference is due to how the stress components are computed
   --- i.e., the stress componen in Moeng's work, :math:`\tau_{xz,0}`, is given by :math:`\tau_{xz,0} = u_{\star}^{2} \bar{u}/|\mathbf{\bar{u}}|`,
   where :math:`\bar{u}/|\mathbf{\bar{u}}|` is the unit vector component applied to the total stress :math:`u_{\star}^{2}`.

   Finally, it must be noted that using terrain-fitted coorindates will modify the surface normal and tangent vectors.
   Consequently, the MOST implementation with terrain-fitted coorindates will formally require local vector rotations.
   Stress rotations with MOST are a work in progress but may be activated with ``erf.use_rotate_surface_flux = true``.
   Therefore, running with terrain (``erf.terrain_type = StaticFittedMesh``) and with MOST (``surface_flux.flux_type = "moeng"``) is not recommended.

MOST Inputs
~~~~~~~~~~~~~~~~~~~
To evaluate the fluxes with MOST, the surface rougness parameter :math:`z_{0}` must be specified. This quantity may be considered a constant or may be parameterized through the friction velocity :math:`u_{\star}`. ERF supports four methods for parameterizing the surface roughness: ``constant``, ``charnock``, ``modified_charnock``, and ``wave_coupled``. The latter three methods parameterize :math:`z_{0} = f(u_{\star})` and are described in `Jimenez & Dudhia, American Meteorological Society, 2018 <https://doi.org/10.1175/JAMC-D-17-0137.1>`_ and `Warner et. al, Ocean Modelling, 2010 <https://doi.org/10.1016/j.ocemod.2010.07.010>`_. The rougness calculation method may be specified with

::

   erf.most.roughness_type    = STRING    #Z_0 type (constant, charnock, modified_charnock, wave_couples)

If the ``charnock`` method is employed, the :math:`a` constant may be specified with ``erf.most.charnock_constant`` (defaults to 0.0185). If the ``modified_charnock`` method is employed, the depth :math:`d` may be specified with ``erf.most.modified_charnock_depth`` (defaults to 30 m). If the ``wave_coupled`` method is employed, the user must provide wave height and mean wavelength data.

While the MOST methods relevant to air-sea interfaces (``charnock``, ``modified_charnock``, and ``wave_coupled``) dynamically compute :math:`z_{0}`, the ``constant`` case can also allow for spatially varying :math:`z_{0}` if an inputs file is employed. More specifically, one may specify

::

   erf.most.roughness_file_name    = STRING    #Name of file that contains (x,y,z_0)

in the inputs file and ERF will populate the 2D :math:`z_{0}` array with values contained in the text file.

When computing an average :math:`\overline{\phi}` for the MOST boundary, where :math:`\phi` denotes a generic variable, ERF supports a variety of approaches. Specifically, ``planar averages`` and ``local region averages`` may be computed with or without ``time averaging``. With each averaging methodology, the query point :math:`z` may be determined from the following procedures: specified vertical distance :math:`z_{ref}` from the bottom surface, specified :math:`k_{index}`, or (when employing terrain-fitted coordinates) specified normal vector length :math:`z_{ref}`. The available inputs to the MOST boundary and their associated data types are

::

   erf.most.average_policy    = INT    #POLICY FOR AVERAGING
   erf.most.use_normal_vector = BOOL   #USE NORMAL VECTOR W/ TERRAIN?
   erf.most.use_interpolation = BOOL   #INTERPOLATE QUERY POINT W/ TERRAIN?
   erf.most.time_average      = BOOL   #USE TIME AVERAGING?
   erf.most.z0                = FLOAT  #SURFACE ROUGHNESS [m]
   erf.most.zref              = FLOAT  #QUERY DISTANCE (HEIGHT OR NORM LENGTH) [m]
   erf.most.surf_temp         = FLOAT  #SPECIFIED SURFACE TEMP [K]
   erf.most.surf_temp_flux    = FLOAT  #SPECIFIED SURFACE TEMP FLUX [K-m/s]
   erf.most.surf_heating_rate = FLOAT  #SPECIFIED RATE OF SURFACE TEMP CHANGE [K/h]
   erf.most.surf_moist        = FLOAT  #SPECIFIED SURFACE MOISTURE [-]
   erf.most.surf_moist_flux   = FLOAT  #SPECIFIED SURFACE MOISTURE FLUX [m/s]
   erf.most.k_arr_in          = INT    #SPECIFIED K INDEX ARRAY (MAXLEV)
   erf.most.radius            = INT    #SPECIFIED REGION RADIUS [grid cells]
   erf.most.time_window       = FLOAT  #WINDOW FOR TIME AVG [s]

We now consider two concrete examples. To employ an instantaneous ``planar average`` at a specified vertical height above the bottom surface, one would specify:

::

   erf.most.average_policy    = 0
   erf.most.use_normal_vector = false
   erf.most.time_average      = false
   erf.most.z0                = 0.1
   erf.most.zref              = 1.0
   erf.most.surf_temp_flux    = 0.0
   erf.most.surf_moist_flux   = 0.0

By contrast, ``local region averaging`` would be employed in conjunction with ``time averaging`` for the following inputs:

::

   erf.most.average_policy    = 1
   erf.most.use_normal_vector = true
   erf.most.use_interpolation = true
   erf.most.time_average      = true
   erf.most.z0                = 0.1
   erf.most.zref              = 1.0
   erf.most.surf_temp_flux    = 0.0
   erf.most.surf_moist_flux   = 0.0
   erf.most.radius            = 1
   erf.most.time_window       = 10.0

In the above case, ``use_normal_vector`` utilizes the a local surface-normal vector with length :math:`z_{ref}` to construct the positions of the query points. Each query point, and surrounding points that are within ``erf.most.radius`` from the query point, are interpolated to and averaged; for a radius of 1, 27 points are averaged. The ``time average`` is completed by way of an exponential filter function whose peak coincides with the current time step and tail extends backwards in time

.. math::

   \frac{1}{\tau} \int_{-\infty}^{0} \exp{\left(t/\tau\right)} \, f(t) \; \rm{d}t.

Due to the form of the above integral, it is advantageous to consider :math:`\tau` as a multiple of the simulation time step :math:`\Delta t`, which is specified by ``erf.most.time_window``. As ``erf.most.time_window`` is reduced to 0, the exponential filter function tends to a Dirac delta function (prior averages are irrelevant). Increasing ``erf.most.time_window`` extends the tail of the exponential and more heavily weights prior averages.

Low-speed corrections
~~~~~~~~~~~~~~~~~~~~~
The following options are available:

::

   erf.most.include_wstar       = true
   erf.most.include_subgrid_vel = true

These correspond to a mean surface velocity of

.. math::

   |\bar{\mathbf{u}}| = \sqrt{u^2 + v^2 + (\beta w_*)^2 + V_{sg}^2}

where :math:`w_*` is the (Deardorff) convective velocity scale and
:math:`\beta=1.2` (Beljaars 1995, QJRMS). The subgrid velocity scale
:math:`V_{sg}` handles weak large-scale flow that is underresolved (Mahrt & Sun
1995, MWR). This is parameterized as

.. math::

   V_{sg} = 0.32 \left(\frac{\Delta x}{5000} - 1 \right)^{0.33}

which vanishes for grid spacings of :math:`\Delta x < 5` km.
