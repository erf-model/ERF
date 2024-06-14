Wind farm models
==================

Introduction
-------------

ERF supports models for wind farm parametrization in which the effects of wind turbines are represented by imposing a momentum sink on the mean flow and/or turbulent kinetic energy (TKE). Currently only the Fitch model (`Fitch et al. 2012`_) and Explicit Wake Parametrization (EWP) model (`Volker et al. 2015`_) are supported.

.. _Fitch model:

Fitch model
------------

The Fitch model for wind farms introduced in `Fitch et al. 2012`_  models the effect of wind farms (See Fig. `1`) as source terms in the governing equations for the horizontal components of momentum (i.e., :math:`x` and :math:`y` momentum) and the turbulent kinetic energy (TKE). The wind turbine is discretized only in the vertical (ie. `z` direction). At a given cell :math:`(i,j,k)`, the source terms in the governing equations are

.. _`Fitch et al. 2012`: https://journals.ametsoc.org/view/journals/mwre/140/9/mwr-d-11-00352.1.xml

.. math::

    \frac{\partial u_{ijk}}{\partial t} &= \frac{u_{ijk}}{|V|_{ijk}}\frac{\partial |V|_{ijk}}{\partial t} \\
    \frac{\partial v_{ijk}}{\partial t} &= \frac{v_{ijk}}{|V|_{ijk}}\frac{\partial |V|_{ijk}}{\partial t} \\
    \frac{\partial \text{TKE}_{ijk}}{\partial t} &=  \frac{0.5N_{ij}C_{TKE}(|V|_{ijk})|V|_{ijk}^3A_{ijk}}{\Delta x \Delta y (z_{k+1}-z_k)}

where

.. math::

    \frac{\partial |V|_{ijk}}{\partial t} = \frac{0.5N_{ij}C_T(|V|_{ijk})|V|_{ijk}^2A_{ijk}}{\Delta x \Delta y (z_{k+1}-z_k)}

where `u` and `v` are horizontal components of velocity, `|V|` is the velocity magnitude, :math:`N_{ij}` is the number of turbines in cell :math:`(i,j)` (see Fig. :numref:`fig:WindFarm`), :math:`C_T` is the coefficient of thrust of the turbines, :math:`C_{TKE}` is the fraction of energy converted to turbulent kinetic energy -- both of which are functions of the velocity magnitude, and :math:`A_{ijk}` is the area intersected by the swept area of the turbine between levels :math:`z=z_k` and :math:`z= z_{k+1}`, and is explained in the next section.


Intersected area :math:`A_{ijk}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Consider :math:`A_k^{k+1}` -- the area intersected by the swept area of the wind turbine between :math:`z=z_k` and :math:`z = z_{k+1}`. We have (see Figs. :numref:`fig:WindTurbine_Fitch` and :numref:`fig:Fitch_Aijk` below)

.. math::

    A_k = \frac{\pi R^2}{2} - A_{ks}

where :math:`A_{ks}` is the area of the segment of the circle as shown in Fig. :numref:`fig:Fitch_Aijk`. We have from geometry, :math:`d_k = \min(|z_k - h|,R)` is the perpendicular distance of the center of the turbine to :math:`z = z_k`, where :math:`h` is the height of the center of the turbine from the ground. The area of the segment is

.. math::

    A_{ks} = R^2\theta_k - d_k\sqrt{R^2 - d_k^2}

where :math:`\theta_k = \cos^{-1}\left(\frac{d_k}{R}\right)`.

Hence, we have the intersected area :math:`A_{ijk}\equiv A_k^{k+1}` as

.. math::

    A_k^{k+1} =
    \begin{cases}
        |A_k - A_{k+1}| & \text{if } (z_k - h)(z_{k+1}-h) > 0 \\
        |A_k + A_{k+1}| & \text{if } (z_k - h)(z_{k+1}-h) \le 0 \\
    \end{cases}

An example of the Fitch model is in ``Exec/Fitch``

.. _fig:WindFarm:

.. figure:: ../figures/WindFarm_Fitch.png
   :width: 300
   :align: center

   Horizontal view of the wind farm in the Fitch model -- shows a wind farm in cell `(i,j)` with 5 wind turbines. The turbines are discretized only in the vertical direction.

.. _fig:WindTurbine_Fitch:

.. figure:: ../figures/WindTurbine_Fitch.png
   :width: 400
   :align: center

   The vertical discretization of the wind turbine in the Fitch model.

.. _fig:Fitch_Aijk:

.. figure:: ../figures/FitchModel_A_ijk.png
   :width: 400
   :align: center

   The area terminology in the Fitch model. The circle represents the area swept by the turbine blades.


.. _explicit-wake-parametrization-ewp-model:

Explicit Wake Parametrization (EWP) model
-----------------------------------------

The Explicit Wake Parametrization (EWP) model [`Volker et al. 2015`_] is very similar to the Fitch model, and models the effect of wind farms as source terms in the governing equations for the horizontal components of momentum (i.e., :math:`x` and :math:`y` momentum) and the turbulent kinetic energy (TKE). At a given cell :math:`(i,j,k)`, the source terms in the governing equations are:

.. math::
    \frac{\partial u_{ijk}}{\partial t} = -\sqrt{\frac{\pi}{8}}\frac{N_{ij}c_tr_0^2\overline{u}_0^2}{\Delta x \Delta y \sigma_e}
    \exp\left\{-\frac{1}{2}\left(\frac{z-h}{\sigma_e}\right)^2\right\}\cos(\phi(k))

.. math::
    \frac{\partial v_{ijk}}{\partial t} = -\sqrt{\frac{\pi}{8}}\frac{N_{ij}c_tr_0^2\overline{u}_0^2}{\Delta x \Delta y \sigma_e}
    \exp\left\{-\frac{1}{2}\left(\frac{z-h}{\sigma_e}\right)^2\right\}\sin(\phi(k))

.. math::
    \frac{\partial \text{TKE}_{ijk}}{\partial t} = -N_{ij}\rho A_rc_t\langle \overline{u}_{i,h}\overline{u'^2}_{i,h}\rangle

with

.. math::
    \sigma_e = \frac{\overline{u}_0}{3KL}\left[\left(\frac{2KL}{\overline{u}_0} + \sigma_0^2\right)^{\frac{3}{2}} - \sigma_0^3\right]

where :math:`N_{ij}` is the number of turbines in cell :math:`(i,j)`, :math:`c_t` is the thrust coefficient, :math:`r_0` is the rotor radius, :math:`\overline{u}_0` is the mean advection velocity at hub height, :math:`h` is the hub height, :math:`\sigma_0 \approx 1.7 r_0` [`Volker et al. 2017`_] is a length scale that accounts for near-wake expansion, :math:`L` is the downstream distance that the wake travels within the cell approximated as a fraction of the cell size, :math:`K` is the turbulence eddy diffusivity (:math:`m^2/s`), :math:`\Delta x` and :math:`\Delta y` are the mesh sizes in the horizontal directions, and :math:`\phi(k)` is the wind direction with respect to the x-axis. :math:`\overline{u}_{i,h}` and :math:`\overline{u'}_{i,h}` are the mean and the fluctuating values of the velocity components (subscript :math:`i` is the direction index) at the hub height :math:`h`, :math:`A_r = \pi r^2` is the swept area of the rotor and :math:`\rho` is the density of air.

The EWP model does not have a concept of intersected area by the turbine rotor like the Fitch model (see :ref:`Fitch model`). The exponential factor in the source terms for the velocities models the effect of the rotor for the momentum sinks (see Fig. :numref:`fig:WindTurbine_EWP`), and the turbulent kinetic energy source term only depends on the density, hub-height mean velocities and fluctuations, and the total swept area of the rotor :math:`A_r`.

.. _fig:WindTurbine_EWP:

.. figure:: ../figures/WindTurbine_EWP.png
   :width: 400
   :align: center

   In the EWP model, the exponential factor in the source terms for the velocities models the effect of the rotor for the momentum sinks unlike the Fitch model which computes the
   intersected area (see Fig. :numref:`fig:WindTurbine_Fitch`).

.. _`Volker et al. 2015`: https://gmd.copernicus.org/articles/8/3715/2015/
.. _`Volker et al. 2017`: https://iopscience.iop.org/article/10.1088/1748-9326/aa5d86


.. _Inputs:

Inputs for wind farm parametrization models
------------------------------------------------------------

Fitch, EWP
~~~~~~~~~~~

The following are the inputs required for simulations with Fitch, EWP models.

.. code-block:: cpp
	
    // The parametrization model - Fitch, EWP
    erf.windfarm_type = "Fitch"

    // How are the turbine locations specified? - using latitude-longitude
    // format or x-y format? lat_lon or x_y
    erf.windfarm_loc_type = "lat_lon"

    // If using lat_lon, then the latitude and longitude of
    // the lower bottom corner of the domain has to be specified
    // The following means 35 deg N, 100 deg W (note the negative sign)
    erf.latitude_lo      =   35.0 
    erf.longitude_lo     = -100.0 

    // Table having the wind turbine locations
    erf.windfarm_loc_table = "windturbines_1WT.txt"

    // Table having the specifications of the wind turbines. All turbines are assumed to 
    // have the same specifications
    erf.windfarm_spec_table = "wind-turbine_1WT.tbl"

| ``erf.windfarm_type`` has to be one of the supported models - ``Fitch``, ``EWP``.
| ``erf.windfarm_loc_type`` is a variable to specify how the wind turbine locations in the wind farm is specified. If using the latitude and longitude of the turbine location, this has to be ``lat_lon`` or if using x and y co-ordinates to specify the turbine locations, this input is ``x_y``. ``erf.latitude_lo`` and ``erf.longitude_lo`` specifies the latitude and longitude of the lower bottom corner of the domain box.  ie. if the domain box is specified as

.. code-block:: cpp

    geometry.prob_lo     = -25000.   0.  -10000
    geometry.prob_hi     =  25000. 10000. 10000.

| then ``(erf.latitude_lo, erf.longitude_lo)`` corresponds to ``(-25000, 0, -10000)``.
| The ``erf.windfarm_loc_table`` specifies the locations of the wind turbines in the wind farm. 
| For the latitude-longitude format, an example is as below. Each line specifies the latitude and longitude of the wind turbine location. The third entry simply has to be always 1 (WRF requires the third entry to be always 1, so maintaining same format here). The first entry means that the turbine is located at ``35.7828 deg N, 99.0168 deg W`` (note the negative sign in the entry corresponding to West).

.. code-block:: cpp

    35.7828828829 -99.0168 1
    35.8219219219 -99.0168 1
    35.860960961 -99.0168 1
    35.9 -99.0168 1
    35.939039039 -99.0168 1
    35.9780780781 -99.0168 1
    36.0171171171 -99.0168 1
    35.7828828829 -98.9705333333 1

For the x-y format, an example is as below. Each line specifies the x and y co-ordinates of the wind turbine location in metres

.. code-block:: cpp

	89264.99080053 91233.3333309002
	89259.1966417755 95566.6666710007
	89254.1277665419 99900.0000000001
	89249.7842982733 104233.333329
	89246.1663427532 108566.6666691
	89243.2739881117 112899.9999981
	93458.6633652711 86900.0000019001
	93450.4377452458 91233.3333309002
	93442.9032518779 95566.6666710007

| The ``erf.windfarm_spec_table`` gives the specifications of the wind turbines in the wind farm. All wind turbines are assumed to have the same specifications. Here is a sample specifications table.

.. code-block:: cpp

    4
    119.0 178.0 0.130 2.0
    9   0.805    50.0
    10   0.805    50.0
    11   0.805    50.0
    12   0.805    50.0

The first line is the number of pairs of entries for the power curve and thrust coefficient (there are 4 entries in this table which are in the last four lines of the table).
The second line gives the height in meters of the turbine hub, the diameter in
meters of the rotor, the standing thrust coefficient, and the nominal power of the turbine in MW.
The remaining lines (four in this case) contain the three values of: wind speed (m/s), thrust coefficient, and power production in kW.
