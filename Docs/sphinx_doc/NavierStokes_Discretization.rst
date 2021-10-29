#########################################
Discretization of Navier-Stokes Equations
#########################################
Last update: 2021-10-29

NOTE: For the sake of simplicity, the discretized equations mention time levels :math:`n` and :math:`n+1`. They should be treated as initial and final states of each RK3 stage.

Staggered Grids
===============
The staggered grids indicating where different variables are located.

XY Plane
--------
.. image:: figures/grid_discretization/stagger_XY.PNG
  :width: 400

YZ Plane
--------
.. image:: figures/grid_discretization/stagger_YZ.PNG
  :width: 400

Mass Conservation
=================

Difference Equation
-------------------

.. math::

   \begin{align}
   \rho_{i,j,k}^{n + 1} = \rho_{i,j,k}^{n} - \Delta t & \left\{         \frac{1}{\Delta x} \left\lbrack \left( \rho u \right)_{i + 1,j,k}^{n} - \left( \rho u \right)_{i.j,k}^{n} \right\rbrack \right. \\
                                                      & \hspace{-5pt} + \frac{1}{\Delta y} \left\lbrack \left( \rho v \right)_{i,j + 1,k}^{n} - \left( \rho v \right)_{i.j,k}^{n} \right\rbrack \\
                                                      & \left. +        \frac{1}{\Delta z} \left\lbrack \left( \rho w \right)_{i,j,k + 1}^{n} - \left( \rho w \right)_{i,j,k}^{n} \right\rbrack \right\}
   \end{align}


Contributions from different directions
---------------------------------------
.. image:: figures/grid_discretization/continuity_x.PNG
  :width: 400
.. image:: figures/grid_discretization/continuity_y.PNG
  :width: 400
.. image:: figures/grid_discretization/continuity_z.PNG
  :width: 400

Advection Contribution to DNS/LES
=================================

Momentum Conservation – U Momentum
----------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho u \right)_{i,j,k}^{n + 1} & = \left( \rho u \right)_{i,j,k}^{n} - \\
      \Delta t &  \left\{ \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \rho u \right)_{i + 1,j,k}^{n} + \left( \rho u \right)_{i,j,k}^{n}         \right)u_{i + \frac{1}{2},j,k}^{n} - \left( \left( \rho u \right)_{i,j,k}^{n} + \left( \rho u \right)_{i - 1,j,k}^{n} \right)u_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
               &        + \frac{1}{2\Delta y}  \left\lbrack \left( \left( \rho v \right)_{i,j + 1,k}^{n} + \left( \rho v \right)_{i - 1,j + 1,k}^{n} \right)u_{i,j + \frac{1}{2},k}^{n} - \left( \left( \rho v \right)_{i,j,k}^{n} + \left( \rho v \right)_{i - 1,j,k}^{n} \right)u_{i,j - \frac{1}{2},k}^{n} \right\rbrack          \\
               & + \left. \frac{1}{2\Delta z}  \left\lbrack \left( \left( \rho w \right)_{i,j,k + 1}^{n} + \left( \rho w \right)_{i - 1,j,k + 1}^{n} \right)u_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \rho w \right)_{i,j,k}^{n} + \left( \rho w \right)_{i - 1,j,k}^{n} \right)u_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
               & - \frac{\Delta t}{\Delta x}\left\lbrack p_{i,\ j,\ k}^{n} - p_{i - 1,\ j,\ k}^{n} \right\rbrack \\
   \end{align}

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/x_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/x_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/x_mom_advec_z.PNG
  :width: 400

Momentum Conservation – V Momentum
----------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho v \right)_{i,j,k}^{n + 1} & = \left( \rho v \right)_{i,j,k}^{n} - \\
    \Delta t & \left\{    \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \rho u \right)_{i + 1,j,k}^{n} + \left( \rho u \right)_{i + 1,j - 1,k}^{n} \right)v_{i + \frac{1}{2},j,k}^{n} - \left( \left( \rho u \right)_{i,j,k}^{n} + \left( \rho u \right)_{i,j - 1,k}^{n} \right)v_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
             & +          \frac{1}{2\Delta y}  \left\lbrack \left( \left( \rho v \right)_{i,j + 1,k}^{n} + \left( \rho v \right)_{i,j,k}^{n}         \right)v_{i,j + \frac{1}{2},k}^{n} - \left( \left( \rho v \right)_{i,j,k}^{n} + \left( \rho v \right)_{i,j - 1,k}^{n} \right)v_{i,j - \frac{1}{2},k}^{n} \right\rbrack          \\
             & + \left. \ \frac{1}{2\Delta z}  \left\lbrack \left( \left( \rho w \right)_{i,j,k + 1}^{n} + \left( \rho w \right)_{i,j - 1,k + 1}^{n} \right)v_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \rho w \right)_{i,j,k}^{n} + \left( \rho w \right)_{i,j - 1,k}^{n} \right)v_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
             & - \frac{\Delta t}{\Delta y}\left\lbrack p_{i,j,\ k}^{n} - p_{i,\ j - 1,\ k}^{n} \right\rbrack \\
   \end{align}

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/y_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/y_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/y_mom_advec_z.PNG
  :width: 400

Momentum Conservation – W Momentum
----------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho w \right)_{i,j,k}^{n + 1} & = \left( \rho w \right)_{i,j,k}^{n} - \\
   \Delta t & \left\{    \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \rho u \right)_{i + 1,j,k}^{n} + \left( \rho u \right)_{i + 1,j,k - 1}^{n} \right)w_{i + \frac{1}{2},j,k}^{n} - \left( \left( \rho u \right)_{i,j,k}^{n} + \left( \rho u \right)_{i,j,k - 1}^{n} \right)w_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
            & +          \frac{1}{2\Delta y}  \left\lbrack \left( \left( \rho v \right)_{i,j + 1,k}^{n} + \left( \rho v \right)_{i,j + 1,k - 1}^{n} \right)w_{i,j + \frac{1}{2},k}^{n} - \left( \left( \rho v \right)_{i,j,k}^{n} + \left( \rho v \right)_{i,j,k - 1}^{n} \right)w_{i,j - \frac{1}{2},k}^{n} \right\rbrack          \\
            & + \left. \ \frac{1}{2\Delta z}  \left\lbrack \left( \left( \rho w \right)_{i,j,k + 1}^{n} + \left( \rho w \right)_{i,j,k}^{n}         \right)w_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \rho w \right)_{i,j,k}^{n} + \left( \rho w \right)_{i,j,k - 1}^{n} \right)w_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
            & - \frac{\Delta t}{\Delta z}\left\lbrack p_{i,\ j,\ k}^{n} - p_{i,\ j,\ \ k - 1}^{n} \right\rbrack\  + \ \Delta t \left\lbrack \rho_{i,j,k - \ \frac{1}{2}}^{n} \right\rbrack g \\
   \end{align}

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/z_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/z_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/z_mom_advec_z.PNG
  :width: 400


Energy Conservation – Potential Temperature Advection
-----------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho \theta \right)_{i,j,k}^{n + 1}  = \left( \rho \theta \right)_{i,j,k}^{n} -
    \Delta t &   \left\{ \frac{1}{\Delta x}\ \left\lbrack \left( \rho u \right)_{i + 1,j,k}^{n} \theta_{i + \frac{1}{2},j,k}^{n} - \left( \rho u \right)_{i,j,k}^{n} \theta_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
             & +         \frac{1}{\Delta y}  \left\lbrack \left( \rho v \right)_{i,j + 1,k}^{n} \theta_{i,j + \frac{1}{2},k}^{n} - \left( \rho v \right)_{i,j,k}^{n} \theta_{i,j - \frac{1}{2},k}^{n} \right\rbrack          \\
             & + \left.  \frac{1}{\Delta z}  \left\lbrack \left( \rho w \right)_{i,j,k + 1}^{n} \theta_{i,j,k + \frac{1}{2}}^{n} - \left( \rho w \right)_{i,j,k}^{n} \theta_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
   \end{align}

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/temp_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/temp_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/temp_advec_z.PNG
  :width: 400

Scalar Conservation – Scalar Advection
--------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho C \right)_{i,j,k}^{n + 1} = \left( \rho C \right)_{i,j,k}^{n} -
   \Delta t & \left\{  \frac{1}{\Delta x}\ \left\lbrack \left( \rho u \right)_{i + 1,j,k}^{n} C_{i + \frac{1}{2},j,k}^{n} - \left( \rho u \right)_{i,j,k}^{n} C_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
            & +        \frac{1}{\Delta y}  \left\lbrack \left( \rho v \right)_{i,j + 1,k}^{n} C_{i,j + \frac{1}{2},k}^{n} - \left( \rho v \right)_{i,j,k}^{n} C_{i,j - \frac{1}{2},k}^{n} \right\rbrack          \\
            & + \left. \frac{1}{\Delta z}  \left\lbrack \left( \rho w \right)_{i,j,k + 1}^{n} C_{i,j,k + \frac{1}{2}}^{n} - \left( \rho w \right)_{i,j,k}^{n} C_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
   \end{align}


Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/scalar_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/scalar_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/scalar_advec_z.PNG
  :width: 400

Diagnostic Variables
--------------------

.. math::

  p_{i, j, k}^n = (\rho_{i, j, k}^n R_d \theta_{i, j, k}^n / p_0^{R_d / c_p} )^\gamma

.. math::

  T_{i, j, k}^n =  \frac{p_{i, j, k}^n}{  \rho_{i, j, k}^n R_d}

Here :math:`\rho_{i, j, k}^n, T_{i, j, k}^n, \theta_{i, j, k}^n`, and :math:`p_{i, j, k}^n` are the density, temperature, potential temperature and pressure, respectively;
these variables are all defined at cell centers of cell indexed by :math:`(i, j, k)` and at time level :math:`n`.

:math:`R_d` and :math:`c_p` are the gas constant and specific heat capacity for dry air respectively,
and :math:`\gamma = c_p / (c_p - R_d)` .  :math:`p_0` is a reference value for pressure.


Differencing of Different Orders
================================

:math:`[\rho, u, v, w, \theta]`,  :math:`m = i, j, k`,  and :math:`U_d = [u, v, w]` for :math:`[x, y, z]` directions respectively.

.. math::

   \begin{array}{lll}
   \left. q_{m + \frac{1}{2}} \right|^{2nd} = \frac{1}{2}\left( q_{m + 1} + q_{m} \right)   & & \\
   \left. q_{m + \frac{1}{2}} \right|^{4th} = \frac{7}{12}\left( q_{m + 1} + q_{m} \right)  & \hspace{-5pt} - \frac{1}{12}\left( q_{m + 2} + q_{m - 1} \right)                                                        & \\
   \left. q_{m + \frac{1}{2}} \right|^{6th} = \frac{37}{60}\left( q_{m + 1} + q_{m} \right) & \hspace{-5pt} - \frac{2}{15}\left( q_{m + 2} + q_{m - 1} \right)                                                        & \hspace{-5pt}+ \frac{1}{60}\left( q_{m + 3} + q_{m - 2} \right)          \\
                                                                                            & & \\
   \left. q_{m + \frac{1}{2}} \right|^{3rd} = \left. q_{m + \frac{1}{2}} \right|^{4th}      & \hspace{-5pt} + \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{12}\left\lbrack \left( q_{m + 2} + q_{m - 1} \right) \right. & \hspace{-5pt}- 3\left. \left( q_{m + 1} + q_{m} \right) \right\rbrack    \\
                                                                                            & & \\
   \left. q_{m + \frac{1}{2}} \right|^{5th} = \left. q_{m + \frac{1}{2}} \right|^{6th}      & \hspace{-5pt} - \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{60}\left\lbrack \left( q_{m + 3} + q_{m - 2} \right) \right. & \hspace{-5pt}- 5\left( q_{m + 2} + q_{m - 1} \right)                     + 10\left. \left( q_{m + 1} + q_{m} \right) \right\rbrack \\
   \end{array}


.. math::

   \begin{array}{lll}
   \left. q_{m - \frac{1}{2}} \right|^{2nd} = \frac{1}{2}\left( q_{m} + q_{m - 1} \right)   & & \\
   \left. q_{m - \frac{1}{2}} \right|^{4th} = \frac{7}{12}\left( q_{m} + q_{m - 1} \right)  & \hspace{-5pt} - \frac{1}{12}\left( q_{m + 1} + q_{m - 2} \right)                                                          & \\
   \left. q_{m - \frac{1}{2}} \right|^{6th} = \frac{37}{60}\left( q_{m} + q_{m - 1} \right) & \hspace{-5pt} - \frac{2}{15}\left( q_{m + 1} + q_{m - 2} \right)                                                          & \hspace{-5pt} + \frac{1}{60}\left( q_{m + 2} + q_{m - 3} \right) \\
    & & \\
   \left. q_{m - \frac{1}{2}} \right|^{3rd} = \left. q_{m - \frac{1}{2}} \right|^{4th}      & \hspace{-5pt} + \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{12}\left\lbrack \left( q_{m + 1} + q_{m - 2} \right) \right.\  & \hspace{-5pt} - 3\left. \ \left( q_{m} + q_{m - 1} \right) \right\rbrack \\
    & & \\
   \left. q_{m - \frac{1}{2}} \right|^{5th} = \left. q_{m - \frac{1}{2}} \right|^{6th}      & \hspace{-5pt} - \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{60}\left\lbrack \left( q_{m + 2} + q_{m - 1} \right) \right.\  & \hspace{-5pt} - 5\left( q_{m + 1} + q_{m - 2} \right) + 10\left. \left( q_{m} + q_{m - 1} \right) \right\rbrack
   \end{array}


Ref: Skamarock, W. C., Klemp, J. B., Dudhia, J., Gill, D. O., Liu, Z., Berner, J., ... Huang, X. -yu. (2019). A Description of the Advanced Research WRF Model Version 4 (No. NCAR/TN-556+STR). doi:10.5065/1dfh-6p97
`doi:10.5065/1dfh-6p97 <http://dx.doi.org/10.5065/1dfh-6p97>`_

Momentum, Thermal, and Scalar Diffusion Contribution to DNS
===========================================================

Strain Rate Tensor
------------------
The schematic below shows the definition of strain-rate components 
.. image:: figures/grid_discretization/StrainRate.PNG
  :width: 400

.. math::

   \begin{array}{ll}
   S_{11,i + \frac{1}{2}} = & \frac{1}{\Delta x}\left( u_{i + 1,j,k} - u_{i,j,k} \right) \\
   S_{11,i - \frac{1}{2}} = & \frac{1}{\Delta x}\left( u_{i,j,k} - u_{i - 1,j,k} \right) \\
   S_{12,j + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j + 1,k} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( v_{i,j + 1,k} - v_{i - 1,j + 1,k} \right) \right\rbrack \\
   S_{12,j - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   S_{13,k + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k + 1} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( w_{i,j,k + 1} - w_{i - 1,j,k + 1} \right) \right\rbrack \\
   S_{13,k - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   \end{array}


.. math::

   \begin{array}{ll}
   S_{21,i + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i + 1,j,k} - u_{i + 1,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i + 1,j,k} - v_{i,j,k} \right) \right\rbrack \\
   S_{21,i - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   S_{22,j + \frac{1}{2}} = & \frac{1}{\Delta y}\left( v_{i,j + 1,k} - v_{i,j,k} \right) \\
   S_{22,j - \frac{1}{2}} = & \frac{1}{\Delta y}\left( v_{i,j,k} - v_{i,j - 1,k} \right) \\
   S_{23,k + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k + 1} - v_{i,j,k} \right) + \frac{1}{\Delta y}\left( w_{i,j,k + 1} - w_{i,j - 1,k + 1} \right) \right\rbrack \\
   S_{23,k - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   \end{array


.. math::

   \begin{array}{ll}
   S_{31,i + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i + 1,j,k} - u_{i + 1,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i + 1,j,k} - w_{i,j,k} \right) \right\rbrack \\
   S_{31,i - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k}     - u_{i,j,k - 1}     \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   S_{32,j + \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j + 1,k} - v_{i,j + 1,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j + 1,k} - w_{i,j,k} \right) \right\rbrack \\
   S_{32,j - \frac{1}{2}} = & \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k}     - v_{i,j,k - 1}     \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   S_{33,k + \frac{1}{2}} = & \frac{1}{\Delta z}\left( w_{i,j,k + 1} - w_{i,j,k} \right) \\
   S_{33,k - \frac{1}{2}} = & \frac{1}{\Delta z}\left( w_{i,j,k} - w_{i,j,k - 1} \right) \\
   \end{array}

Momentum Conservation – U Momentum viscous stress divergence
------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho u \right)_{i,j,k}^{n + 1} = \left( \rho u \right)_{i,j,k}^{n} + & \\
     \Delta t &  \left\{ \frac{1}{\Delta x}\ \left\lbrack \tau_{11,i + \frac{1}{2}} - \tau_{11,i - \frac{1}{2}} \right\rbrack \right.\  \\
              &        + \frac{1}{\Delta y}\ \left\lbrack \tau_{12,j + \frac{1}{2}} - \tau_{12,j - \frac{1}{2}} \right\rbrack           \\
              & + \left. \frac{1}{\Delta z}\ \left\lbrack \tau_{13,k + \frac{1}{2}} - \tau_{13,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{align}

Note that LES equation has a similar format except how :math:`\tau_{11,i + \frac{1}{2}}`, :math:`\tau_{11,i - \frac{1}{2}}`, :math:`\tau_{12,j + \frac{1}{2}}`, :math:`\tau_{12,j - \frac{1}{2}}`, :math:`\tau_{13,k + \frac{1}{2}}`, and :math:`\tau_{13,k - \frac{1}{2}}`are defined.

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/x_mom_diff_a.PNG
  :width: 400
.. image:: figures/grid_discretization/x_mom_diff_b.PNG
  :width: 400

Momentum Conservation – V Momentum viscous stress divergence
------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho v \right)_{i,j,k}^{n + 1} = \left( \rho v \right)_{i,j,k}^{n} + & \\
     \Delta t & \left\{  \frac{1}{\Delta x} \left\lbrack \tau_{21,i + \frac{1}{2}} - \tau_{21,i - \frac{1}{2}} \right\rbrack \right.  \\
              &        + \frac{1}{\Delta y} \left\lbrack \tau_{22,j + \frac{1}{2}} - \tau_{22,j - \frac{1}{2}} \right\rbrack          \\
              & + \left. \frac{1}{\Delta z} \left\lbrack \tau_{23,k + \frac{1}{2}} - \tau_{23,k - \frac{1}{2}} \right\rbrack \right\}
   \end{align}

Note that LES equation has a similar format except how :math:`\tau_{21,i + \frac{1}{2}}`, :math:`\tau_{21,i - \frac{1}{2}}`, :math:`\tau_{22,j + \frac{1}{2}}`, :math:`\tau_{22,j - \frac{1}{2}}`, :math:`\tau_{23,k + \frac{1}{2}}`, and :math:`\tau_{23,k - \frac{1}{2}}`are defined.

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/y_mom_diff_a.PNG
  :width: 400
.. image:: figures/grid_discretization/y_mom_diff_b.PNG
  :width: 400

Momentum Conservation – W Momentum viscous stress divergence
------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho w \right)_{i,j,k}^{n + 1} = \left( \rho w \right)_{i,j,k}^{n} + & \\
    \Delta t &  \left\{ \frac{1}{\Delta x} \left\lbrack \tau_{31,i + \frac{1}{2}} - \tau_{31,i - \frac{1}{2}} \right\rbrack \right.  \\
             &        + \frac{1}{\Delta y} \left\lbrack \tau_{32,j + \frac{1}{2}} - \tau_{32,j - \frac{1}{2}} \right\rbrack            \\
             & + \left. \frac{1}{\Delta z} \left\lbrack \tau_{33,k + \frac{1}{2}} - \tau_{33,k - \frac{1}{2}} \right\rbrack \right\}
   \end{align}

Note that LES equation has a similar format except how :math:`\tau_{31,i + \frac{1}{2}}`, :math:`\tau_{31,i - \frac{1}{2}}`, :math:`\tau_{32,j + \frac{1}{2}}`, :math:`\tau_{32,j - \frac{1}{2}}`, :math:`\tau_{33,k + \frac{1}{2}}`, and :math:`\tau_{33,k - \frac{1}{2}}`are defined.

Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/z_mom_diff_a.PNG
  :width: 400
.. image:: figures/grid_discretization/z_mom_diff_b.PNG
  :width: 400

Energy Conservation – Potential Temperature Diffusion
-----------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{matrix}
   \left( \rho \theta \right)_{i,j,k}^{n + 1} & = & \left( \rho \theta \right)_{i,j,k}^{n} & + & \Delta t\rho_{i,j,k}\alpha_{T} & \left\{ \frac{1}{{\Delta x}^{2}}\ \left\lbrack \theta_{i + 1,j,k}^{n} - \ {2\theta}_{i,j,k}^{n} + \ \theta_{i - 1,j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{{\Delta y}^{2}}\left\lbrack \theta_{i,j + 1,k}^{n} - \ 2\theta_{i,j,k}^{n} + \ \theta_{i,j - 1,k}^{n} \right\rbrack \\
    & & & & & + \left. \frac{1}{{\Delta z}^{2}}\left\lbrack \theta_{i,j,k + 1}^{n} - \ {2\theta}_{i,j,k}^{n} + \ \theta_{i,j,k - 1}^{n} \right\rbrack \right\}
   \end{matrix}


Scalar Conservation – Scalar Diffusion
--------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{matrix}
   \left( \rho C \right)_{i,j,k}^{n + 1} & = & \left( \rho C \right)_{i,j,k}^{n} & + & \Delta t\rho_{i,j,k}\alpha_{S} & \left\{ \frac{1}{{\Delta x}^{2}}\ \left\lbrack C_{i + 1,j,k}^{n} - \ {2C}_{i,j,k}^{n} + \ C_{i - 1,j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{{\Delta y}^{2}}\left\lbrack C_{i,j + 1,k}^{n} - \ 2C_{i,j,k}^{n} + \ C_{i,j - 1,k}^{n} \right\rbrack \\
    & & & & & + \left. \frac{1}{{\Delta z}^{2}}\left\lbrack C_{i,j,k + 1}^{n} - \ {2C}_{i,j,k}^{n} + \ C_{i,j,k - 1}^{n} \right\rbrack \right\}
   \end{matrix}

Momentum, Thermal, and Scalar Diffusion Contribution to LES
===========================================================

Strain Rate and Eddy Viscosity
------------------------------

The goal is to compute eddy viscosity at the *cell centers* and interpolated them to the edges. Refer again to the strain rate tensor schematic.

.. image:: figures/grid_discretization/StrainRate.PNG
  :width: 400

.. math::

   \begin{array}{ll}
   S_{11} = & S_{11i + \frac{1}{2}} \\
   S_{22} = & S_{22j + \frac{1}{2}} \\
   S_{33} = & S_{33k + \frac{1}{2}}
   \end{array}


.. math::

   \begin{matrix}
   S_{12} = & \frac{1}{4}\left\lbrack S_{12i,j - \frac{1}{2}} + S_{12i,j + \frac{1}{2}} + S_{12i + 1,j - \frac{1}{2}} + S_{12i + 1,j + \frac{1}{2}} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix} \\
   S_{21} = & \frac{1}{4}\left\lbrack S_{21i - \frac{1}{2},j} + S_{21i + \frac{1}{2},j} + S_{21i - \frac{1}{2},j + 1} + S_{21i + \frac{1}{2},j + 1} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix} \\
   S_{13} = & \frac{1}{4}\left\lbrack S_{13i,k - \frac{1}{2}} + S_{13i,k + \frac{1}{2}} + S_{13i + 1,k - \frac{1}{2}} + S_{13i + 1,k + \frac{1}{2}} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix} \\
   S_{31} = & \frac{1}{4}\left\lbrack S_{31i - \frac{1}{2},k} + S_{31i + \frac{1}{2},k} + S_{31i - \frac{1}{2},k + 1} + S_{31i + \frac{1}{2},k + 1} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix} \\
   S_{23} = & \frac{1}{4}\left\lbrack S_{23j,k - \frac{1}{2}} + S_{23j,k + \frac{1}{2}} + S_{23j + 1,k - \frac{1}{2}} + S_{23j + 1,k + \frac{1}{2}} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix} \\
   S_{32} = & \frac{1}{4}\left\lbrack S_{32j - \frac{1}{2},k} + S_{32j + \frac{1}{2},k} + S_{32j - \frac{1}{2},k + 1} + S_{32j + \frac{1}{2},k + 1} \right\rbrack = \begin{smallmatrix} \text{Average of the 4 edges} \\ \text{surrouding the cell}\end{smallmatrix}
   \end{matrix}

Note that:

.. math::

   \begin{array}{ll}
   S_{12} = & S_{21} \\
   S_{13} = & S_{31} \\
   S_{23} = & S_{32}
   \end{array}

:math:`K_{i,j,k} = {2\left( C_{S} \right)^{2}\rho_{i,j,k}\left( {2S}_{mn}S_{mn} \right)}^{\frac{1}{2}}`,
where

.. math::

   S_{mn}S_{mn} = S_{11}^{2} + S_{22}^{2} + S_{33}^{2} + S_{12}^{2} + S_{13}^{2} + S_{23}^{2} + S_{21}^{2} + S_{31}^{2} + S_{32}^{2} \\

Owing to symmetry we need to compute 6 of the 9 tensor components.

.. image:: figures/grid_discretization/EddyViscosity.PNG
  :width: 400

The interpolated values of eddy-viscosity at the edges are the average
of the values at the centers of the 4 cells the edge is part of.

.. math::

   \begin{array}{c}
   K_{i + \frac{1}{2},j - \frac{1}{2},k} = \frac{1}{4}\left\lbrack K_{i,j - 1,k} + K_{i,j,k} + K_{i + 1,j - 1,k} + K_{i + 1,j,k} \right\rbrack \\
   K_{i + \frac{1}{2},j + \frac{1}{2},k} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j + 1,k} + K_{i + 1,j,k} + K_{i + 1,j + 1,k} \right\rbrack \\
   K_{i + \frac{1}{2},j,k - \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k - 1} + K_{i + 1,j,k} + K_{i + 1,j,k - 1} \right\rbrack \\
   K_{i + \frac{1}{2},j,k + \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k + 1} + K_{i,j,k} + K_{i + 1,j,k + 1} + K_{i + 1,j,k} \right\rbrack \\
   K_{i,j + \frac{1}{2},k - \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k - 1} + K_{i,j + 1,k} + K_{i,j + 1,k - 1} \right\rbrack \\
   K_{i,j + \frac{1}{2},k + \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k + 1} + K_{i,j + 1,k} + K_{i,j + 1,k + 1} \right\rbrack
   \end{array}

Momentum Conservation – U Momentum - subfilter stress divergence
----------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho u \right)_{i,j,k}^{n + 1} = \left( \rho u \right)_{i,j,k}^{n} + & \\
     \Delta t &  \left\{ \frac{1}{\Delta x}\ \left\lbrack \tau_{11,i + \frac{1}{2}} - \tau_{11,i - \frac{1}{2}} \right\rbrack \right.\  \\
              &        + \frac{1}{\Delta y}\ \left\lbrack \tau_{12,j + \frac{1}{2}} - \tau_{12,j - \frac{1}{2}} \right\rbrack           \\
              & + \left. \frac{1}{\Delta z}\ \left\lbrack \tau_{13,k + \frac{1}{2}} - \tau_{13,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{align}

.. math::

   \begin{array}{l}
   \tau_{11,i + \frac{1}{2}} = K_{i,j,k}                            \ S_{11,i + \frac{1}{2}} = K_{i,j,k}                                                    \frac{1}{\Delta x}\left( u_{i + 1,j,k} - u_{i,j,k} \right) \\
   \tau_{11,i - \frac{1}{2}} = K_{i - 1,j,k}                        \ S_{11,i - \frac{1}{2}} = K_{i - 1,j,k}                                                \frac{1}{\Delta x}\left( u_{i,j,k} - u_{i - 1,j,k} \right) \\
   \tau_{12,j + \frac{1}{2}} = K_{i - \frac{1}{2},j + \frac{1}{2},k}\ S_{12,j + \frac{1}{2}} = K_{i - \frac{1}{2},j + \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j + 1,k} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( v_{i,j + 1,k} - v_{i - 1,j + 1,k} \right) \right\rbrack \\
   \tau_{12,j - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\ S_{12,j - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   \tau_{13,k + \frac{1}{2}} = K_{i - \frac{1}{2},j,k + \frac{1}{2}}\ S_{13,k + \frac{1}{2}} = K_{i - \frac{1}{2},j,k + \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k + 1} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( w_{i,j,k + 1} - w_{i - 1,j,k + 1} \right) \right\rbrack \\
   \tau_{13,k - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\ S_{13,k - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack
   \end{array}

Momentum Conservation – V Momentum - subfilter stress divergence
----------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho v \right)_{i,j,k}^{n + 1} = \left( \rho v \right)_{i,j,k}^{n} + & \\
     \Delta t & \left\{  \frac{1}{\Delta x} \left\lbrack \tau_{21,i + \frac{1}{2}} - \tau_{21,i - \frac{1}{2}} \right\rbrack \right.  \\
              &        + \frac{1}{\Delta y} \left\lbrack \tau_{22,j + \frac{1}{2}} - \tau_{22,j - \frac{1}{2}} \right\rbrack          \\
              & + \left. \frac{1}{\Delta z} \left\lbrack \tau_{23,k + \frac{1}{2}} - \tau_{23,k - \frac{1}{2}} \right\rbrack \right\}
   \end{align}

.. math::

   \begin{array}{ll}
   \tau_{21,i + \frac{1}{2}} = K_{i + \frac{1}{2},j - \frac{1}{2},k}\ S_{21,i + \frac{1}{2}} & \hspace{-5pt} = K_{i + \frac{1}{2},j - \frac{1}{2},k} \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i + 1,j,k} - u_{i + 1,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i + 1,j,k} - v_{i,j,k} \right) \right\rbrack \\
   \tau_{21,i - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\ S_{21,i - \frac{1}{2}} & \hspace{-5pt} = K_{i - \frac{1}{2},j - \frac{1}{2},k} \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} -     u_{i,j - 1,k} \right)     + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   \tau_{22,j + \frac{1}{2}} = K_{i,j,k}                            \ S_{22,j + \frac{1}{2}} & \hspace{-5pt} = K_{i,j,k}                                                     \frac{1}{\Delta y}\left( v_{i,j + 1,k} - v_{i,j,k} \right) \\
   \tau_{22,j - \frac{1}{2}} = K_{i,j - 1,k}                        \ S_{22,j - \frac{1}{2}} & \hspace{-5pt} = K_{i,j - 1,k}                                                 \frac{1}{\Delta y}\left( v_{i,j,k} -     v_{i,j - 1,k} \right) \\
   \tau_{23,k + \frac{1}{2}} = K_{i,j - \frac{1}{2},k + \frac{1}{2}}\ S_{23,k + \frac{1}{2}} & \hspace{-5pt} = K_{i,j - \frac{1}{2},k + \frac{1}{2}} \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k + 1} - v_{i,j,k} \right)     + \frac{1}{\Delta y}\left( w_{i,j,k + 1} - w_{i,j - 1,k + 1} \right) \right\rbrack \\
   \tau_{23,k - \frac{1}{2}} = K_{i,j - \frac{1}{2}k - \frac{1}{2}} \ S_{23,k - \frac{1}{2}} & \hspace{-5pt} = K_{i,j - \frac{1}{2},k - \frac{1}{2}} \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} -     v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k}         \right) \right\rbrack \\
   \end{array}

Momentum Conservation – W Momentum - subfilter stress divergence
----------------------------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

.. math::

   \begin{align}
   \left( \rho w \right)_{i,j,k}^{n + 1} = \left( \rho w \right)_{i,j,k}^{n} + & \\
    \Delta t &  \left\{ \frac{1}{\Delta x} \left\lbrack \tau_{31,i + \frac{1}{2}} - \tau_{31,i - \frac{1}{2}} \right\rbrack \right.  \\
             &        + \frac{1}{\Delta y} \left\lbrack \tau_{32,j + \frac{1}{2}} - \tau_{32,j - \frac{1}{2}} \right\rbrack            \\
             & + \left. \frac{1}{\Delta z} \left\lbrack \tau_{33,k + \frac{1}{2}} - \tau_{33,k - \frac{1}{2}} \right\rbrack \right\}
   \end{align}

.. math::

   \begin{array}{ll}
   \tau_{31,i + \frac{1}{2}} = K_{i + \frac{1}{2},j,k - \frac{1}{2}}\ S_{31,i + \frac{1}{2}} & \hspace{-5pt} = K_{i + \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i + 1,j,k} - u_{i + 1,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i + 1,j,k} - w_{i,j,k} \right) \right\rbrack \\
   \tau_{31,i - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\ S_{31,i - \frac{1}{2}} & \hspace{-5pt} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   \tau_{32,j + \frac{1}{2}} = K_{i,j + \frac{1}{2},k - \frac{1}{2}}\ S_{32,j + \frac{1}{2}} & \hspace{-5pt} = K_{i,j + \frac{1}{2},k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j + 1,k} - v_{i,j + 1,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j + 1,k} - w_{i,j,k} \right) \right\rbrack \\
   \tau_{32,j - \frac{1}{2}} = K_{i,j - \frac{1}{2},k - \frac{1}{2}}\ S_{32,j - \frac{1}{2}} & \hspace{-5pt} = K_{i,j - \frac{1}{2},k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   \tau_{33,k + \frac{1}{2}} = K_{i,j,k}\ S_{33,k + \frac{1}{2}}                             & \hspace{-5pt} = K_{i,j,k}\frac{1}{\Delta z}\left( w_{i,j,k + 1} - w_{i,j,k} \right) \\
   \tau_{33,k - \frac{1}{2}} = K_{i,j, k - 1}\ S_{33,k - \frac{1}{2}}                        & \hspace{-5pt} = K_{i,j,k - 1}\frac{1}{\Delta z}\left( w_{i,j,k} - w_{i,j,k - 1} \right)
   \end{array}

Energy Conservation- Subgrid heat flux
--------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~
This section is yet to be implemented in the code.

.. math::

   \begin{align}
   \left( \rho \theta \right)_{i,j,k}^{n + 1}  =  \left( \rho \theta \right)_{i,j,k}^{n}  - & \\
     \Delta t & \left\{   \frac{1}{\Delta x} \left\lbrack {\vartheta_{1}}_{i + \frac{1}{2},j,k}^{n} - {\vartheta_{1}}_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.  \\
              &         + \frac{1}{\Delta y} \left\lbrack {\vartheta_{2}}_{i,j + \frac{1}{2},k}^{n} - {\vartheta_{2}}_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
              & + \left.  \frac{1}{\Delta z} \left\lbrack {\vartheta_{3}}_{i,j,k + \frac{1}{2}}^{n} - {\vartheta_{3}}_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\}
   \end{align}

.. math::

   \begin{array}{lcl}
   \vartheta_{1i,j,k}               & = & K_{i,j,k}\frac{1}{2\Delta x}\ \left\lbrack \theta_{i + 1,j,k}^{n} - \theta_{i - 1,j,k}^{n} \right\rbrack \\
   \vartheta_{2i,j,k}               & = & K_{i,j,k}\frac{1}{2\Delta y}\ \left\lbrack \theta_{i,j + 1,k}^{n} - \theta_{i,j - 1,k}^{n} \right\rbrack \\
   \vartheta_{3i,jk}                & = & K_{i,j,k}\frac{1}{2\Delta z}\ \left\lbrack \theta_{i,j,k + 1}^{n} - \theta_{i,j,k - 1}^{n} \right\rbrack \\
   \vartheta_{1i + \frac{1}{2},j,k} & = & \frac{1}{2}\left( \vartheta_{1i + 1,j,k} + \vartheta_{1i,j,k} \right) \\
   \vartheta_{1i - \frac{1}{2},j,k} & = & \frac{1}{2}\left( \vartheta_{1i,j,k}     + \vartheta_{1i - 1,j,k} \right) \\
   \vartheta_{1i,j + \frac{1}{2},k} & = & \frac{1}{2}\left( \vartheta_{1i,j + 1,k} + \vartheta_{1i,j,k} \right) \\
   \vartheta_{1i,j - \frac{1}{2},k} & = & \frac{1}{2}\left( \vartheta_{1i,j,k}     + \vartheta_{1i1,j - 1,k} \right) \\
   \vartheta_{1i,j,k + \frac{1}{2}} & = & \frac{1}{2}\left( \vartheta_{1i,j,k +}   + \vartheta_{1i,j,k} \right) \\
   \vartheta_{1i,j,k - \frac{1}{2}} & = & \frac{1}{2}\left( \vartheta_{1i,j,k}     + \vartheta_{1i1,j,k - 1} \right)
   \end{array}

Scalar Conservation- Subgrid scalar flux
----------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~

Prognostic Equation for Subgrid Kinetic Energy
----------------------------------------------

Difference Equation
~~~~~~~~~~~~~~~~~~~
This section is yet to be implemented in the code.

.. math::

   \begin{align}
   \left( \rho e \right)_{i,j,k}^{n + 1} = \left( \rho e \right)_{i,j,k}^{n} - & \\
     \Delta t & \left\{  \frac{1}{\Delta x}\left\lbrack \left( \rho u \right)_{i + 1,j,k}^{n}e_{i + \frac{1}{2},j,k}^{n} - \left( \rho u \right)_{i,j,k}^{n}e_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right. \\
              & +        \frac{1}{\Delta y}\left\lbrack \left( \rho v \right)_{i,j + 1,k}^{n}e_{i,j + \frac{1}{2},k}^{n} - \left( \rho v \right)_{i,j,k}^{n}e_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
              & + \left. \frac{1}{\Delta z}\left\lbrack \left( \rho w \right)_{i,j,k + 1}^{n}e_{i,j,k + \frac{1}{2}}^{n} - \left( \rho w \right)_{i,j,k}^{n}e_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
              & + \frac{g}{\Theta}\vartheta_{3} - \tau_{mn}\frac{\partial u_{m}}{\partial x_{n}} - \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} - \epsilon
   \end{align}

.. math::

   \begin{array}{l}
   \vartheta_{i} = K_{H}\frac{\partial\theta}{\partial x_{i}} \\
    K_{H} = \left( 1 + 2\frac{\mathcal{l}}{\Delta s} \right)K_{M} \\
    K_{M} = 0.1\mathcal{l}e^{\frac{1}{2}} = 0.1\mathcal{l}e_{i,j,k}^{\frac{1}{2}} \\
    K_{{M}_{i,j,k}} = 0.1\mathcal{l}e_{i,j,k}^{\frac{1}{2}}
    \end{array}

:math:`\mathcal{l} = \Delta s = \sqrt[3]{\Delta x \Delta y \Delta z}`,

Convective case:

.. math::

   \begin{align}
   \mathcal{l} & = 0.76\ e^{\frac{1}{2}}\left( \frac{g}{\Theta}\frac{\partial\theta}{\partial z} \right) \\
               & = 0.76e_{i,j,k}^{\frac{1}{2}}\left\lbrack \frac{g}{\Theta}\frac{1}{2\Delta z}\left( \theta_{i,j,k + 1}^{n} - \theta_{i,j,k - 1}^{n} \right) \right\rbrack
   \end{align}

.. math::

   \begin{align}
   \vartheta_{1} = {K_{H}}_{i,j,k}\frac{1}{2\Delta x} \left\lbrack \theta_{i + 1,j,k}^{n} - \theta_{i - 1,j,k}^{n} \right\rbrack \\
   \vartheta_{2} = {K_{H}}_{i,j,k}\frac{1}{2\Delta y} \left\lbrack \theta_{i,j + 1,k}^{n} - \theta_{i,j - 1,k}^{n} \right\rbrack \\
   \vartheta_{3} = {K_{H}}_{i,j,k}\frac{1}{2\Delta z} \left\lbrack \theta_{i,j,k + 1}^{n} - \theta_{i,j,k - 1}^{n} \right\rbrack
   \end{align}

.. math::

   \begin{align}
   \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} & = \\
          K_{i,j,k} & \left\{  \frac{1}{2\Delta x}\left\lbrack e_{i + 1,j,k}^{n} - e_{i - 1,j,k}^{n} \right\rbrack \right. \\
                    &        + \frac{1}{2\Delta y}\left\lbrack e_{i,j + 1,k}^{n} - e_{i,j - 1,k}^{n} \right\rbrack  \\
                    & + \left. \frac{1}{2\Delta z}\left\lbrack e_{i,j,k + 1}^{n} - e_{i,j,k - 1}^{n} \right\rbrack \right\}
   \end{align}

.. math::

   \epsilon = C_{\epsilon}\rho_{i,j,k}\frac{\left( e_{i,j,k} \right)^{\frac{3}{2}}}{\mathcal{l}}

.. math::

   C_{\epsilon} = 0.19 + 0.51\frac{\mathcal{l}}{\Delta s}

.. math::

   \begin{align}
   \tau_{mn}\frac{\partial u_{m}}{\partial x_{n}} & = KS_{mn}\frac{\partial u_{m}}{\partial x_{n}} \\
                                                  & = KS_{mn}S_{mn} \\
                                                  & = K(S_{11}^{2} + S_{22}^{2} + S_{33}^{2} + S_{12}^{2} + S_{13}^{2} + S_{23}^{2} + S_{21}^{2} + S_{31}^{2} + S_{32}^{2})
   \end{align}




Contributions from different directions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. image:: figures/grid_discretization/TKE_x.PNG
  :width: 400
.. image:: figures/grid_discretization/TKE_y.PNG
  :width: 400
.. image:: figures/grid_discretization/TKE_z.PNG
  :width: 400
