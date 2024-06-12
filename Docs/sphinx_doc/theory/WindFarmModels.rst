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

where `u` and `v` are horizontal components of velocity, `|V|` is the velocity magnitude, :math:`N_{ij}` is the number of turbines in cell :math:`(i,j)`, :math:`C_T` is the coefficient of thrust of the turbines, :math:`C_{TKE}` is the fraction of energy converted to turbulent kinetic energy -- both of which are functions of the velocity magnitude, and :math:`A_{ijk}` is the area intersected by the swept area of the turbine between levels :math:`z=z_k` and :math:`z= z_{k+1}`, and is explained in the next section.

Intersected area :math:`A_{ijk}`
_________________________________

Consider :math:`A_k^{k+1}` -- the area intersected by the swept area of the wind turbine between :math:`z=z_k` and :math:`z = z_{k+1}`. We have (see Figs. `2` and `3` below)

.. math::

    A_k = \frac{\pi R^2}{2} - A_{ks}

where :math:`A_{ks}` is the area of the segment of the circle as shown in Fig. `3`. We have from geometry, :math:`d_k = \min(|z_k - z_c|,R)` is the perpendicular distance of the center of the turbine to :math:`z = z_k`, where :math:`z_c` is the height of the center of the turbine from the ground. The area of the segment is

.. math::

    A_{ks} = R^2\theta_k - d_k\sqrt{R^2 - d_k^2}

where :math:`\theta_k = \cos^{-1}\left(\frac{d_k}{R}\right)`.

Hence, we have the intersected area :math:`A_{ijk}\equiv A_k^{k+1}` as

.. math::

    A_k^{k+1} =
    \begin{cases}
        |A_k - A_{k+1}| & \text{if } (z_k - z_c)(z_{k+1}-z_c) > 0 \\
        |A_k + A_{k+1}| & \text{if } (z_k - z_c)(z_{k+1}-z_c) \le 0 \\
    \end{cases}

An example of the Fitch model is in ``Exec/Fitch``

.. 1:

.. figure:: ../figures/WindFarm_Fitch.png
   :width: 300
   :align: center

   Horizontal view of the wind farm in the Fitch model -- shows a wind farm in cell `(i,j)` with 5 wind turbines. The turbines are discretized only in the vertical direction.

.. 2:

.. figure:: ../figures/WindTurbine_Fitch.png
   :width: 300
   :align: center

   The vertical discretization of the wind turbine in the Fitch model.

.. 3:

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

The EWP model does not have a concept of intersected area by the turbine rotor like the Fitch model (see :ref:`Fitch model`). The exponential factor in the source terms for the velocities models the effect of the rotor for the momentum source/sinks, and the turbulent kinetic energy source term only depends on the density, hub-height mean velocities and fluctuations, and the total swept area of the rotor :math:`A_r`.    


.. _`Volker et al. 2015`: https://gmd.copernicus.org/articles/8/3715/2015/
.. _`Volker et al. 2017`: https://iopscience.iop.org/article/10.1088/1748-9326/aa5d86
