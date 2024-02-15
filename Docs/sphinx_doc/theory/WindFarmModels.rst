Wind farm models
==================

Introduction
-------------

ERF supports models for wind farm parametrization in which the effects of wind turbines are represented by imposing a momentum sink on the mean flow and/or turbulent kinetic energy (TKE). Currently only the Fitch model (`Fitch et al. 2012`_) is supported.

Fitch model
------------

The Fitch model for wind farms introduced in `Fitch et al. 2012`_  models the effect of wind farms (See Fig. `1`) as source terms in the governing equations for the horizontal components of momentum (i.e., :math:`x` and :math:`y` momentum) and the turbulent kinetic energy (TKE). The wind turbine is discretized only in the vertical (ie. `z` direction). At a given cell :math:`(i,j,k)`, the source terms in the governing equations are

.. _`Fitch et al. 2012`: https://journals.ametsoc.org/view/journals/mwre/140/9/mwr-d-11-00352.1.xml

.. math::

    \frac{\partial u_{ijk}}{\partial t} &= \frac{u_{ijk}}{|V|_{ijk}}\frac{\partial |V|_{ijk}}{\partial t} \\
    \frac{\partial v_{ijk}}{\partial t} &= \frac{v_{ijk}}{|V|_{ijk}}\frac{\partial |V|_{ijk}}{\partial t} \\
    \frac{\partial \text{TKE}_{ijk}}{\partial t} &=  \frac{0.5N^{ij}C_{TKE}(|V|_{ijk})|V|_{ijk}^3A_{ijk}}{z_{k+1}-z_k}

where

.. math::

    \frac{\partial |V|_{ijk}}{\partial t} = \frac{0.5N^{ij}C_T(|V|_{ijk})|V|_{ijk}^2A_{ijk}}{z_{k+1}-z_k}

where `u` and `v` are horizontal components of velocity, `|V|` is the velocity magnitude, :math:`N^{ij}` is the number of turbines in cell :math:`(i,j)`, :math:`C_T` is the coefficient of thrust of the turbines, :math:`C_{TKE}` is the fraction of energy converted to turbulent kinetic energy -- both of which are functions of the velocity magnitude, and :math:`A_{ijk}` is the area intersected by the swept area of the turbine between levels :math:`z=z_k` and :math:`z= z_{k+1}`, and is explained in the next section.

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
