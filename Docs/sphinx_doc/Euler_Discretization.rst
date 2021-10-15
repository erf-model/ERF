.. highlight:: rst

###################################################
Finite Difference Discretization of Euler Equations
###################################################

Last update: 2021-03-04

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

    \begin{align*}
        \rho_{i, j, k}^{n+1} = \rho_{i, j, k}^{n}
        - \Delta t
        \{& \frac{1}{\Delta x} [ {(\rho u)}_{i+1, j, k}^{n} - {(\rho u)}_{i, j, k}^{n}]\\
        + & \frac{1}{\Delta y} [ {(\rho v)}_{i, j+1, k}^{n} - {(\rho v)}_{i, j, k}^{n}]\\
        + & \frac{1}{\Delta z} [ {(\rho w)}_{i, j, k+1}^{n} - {(\rho w)}_{i, j, k}^{n}] \}
    \end{align*}

Divergence Compoments
---------------------
.. image:: figures/grid_discretization/continuity_x.PNG
  :width: 400
.. image:: figures/grid_discretization/continuity_y.PNG
  :width: 400
.. image:: figures/grid_discretization/continuity_z.PNG
  :width: 400

X-Momentum Conservation
=======================

Difference Equation
-------------------

.. math::

    \begin{align*}
        (\rho u)_{i, j, k}^{n+1} = (\rho u)_{i, j, k}^{n}
        - \Delta t
        \{& \frac{1}{2 \Delta x} [( {(\rho u)}_{i+1, j, k}^{n} + {(\rho u)}_{i, j, k}^{n}) u_{i+\frac{1}{2},j,k}^n
                                 -( {(\rho u)}_{i, j, k}^{n} + {(\rho u)}_{i-1, j, k}^{n}) u_{i-\frac{1}{2},j,k}^n]\\
        + & \frac{1}{2 \Delta y} [( {(\rho v)}_{i, j+1, k}^{n} + {(\rho v)}_{i-1, j+1, k}^{n}) u_{i, j+\frac{1}{2},k}^n
                                 -( {(\rho v)}_{i, j, k}^{n} + {(\rho v)}_{i-1, j, k}^{n}) u_{i, j-\frac{1}{2},k}^n]\\
        + & \frac{1}{2 \Delta z} [( {(\rho w)}_{i, j, k+1}^{n} + {(\rho w)}_{i-1, j, k+1}^{n}) u_{i, j, k+\frac{1}{2}}^n
                                 -( {(\rho w)}_{i, j, k}^{n} + {(\rho w)}_{i-1, j, k}^{n}) u_{i, j, k-\frac{1}{2}}^n] \}\\
        - & \frac{\Delta t}{\Delta x}[p_{i, j, k}^{n} - p_{i-1, j, k}^{n}]    
    \end{align*}

Divergence Compoments
---------------------
.. image:: figures/grid_discretization/x_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/x_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/x_mom_advec_z.PNG
  :width: 400

Y-Momentum Conservation
=======================

Difference Equation
-------------------

.. math::

    \begin{align*}
    (\rho v)_{i, j, k}^{n+1} = (\rho v)_{i, j, k}^{n}
    - \Delta t
        \{& \frac{1}{2 \Delta x} [( {(\rho u)}_{i+1, j, k}^{n} + {(\rho u)}_{i+1, j-1, k}^{n}) v_{i+\frac{1}{2},j,k}^n
                                 -( {(\rho u)}_{i, j, k}^{n} + {(\rho u)}_{i, j-1, k}^{n}) v_{i-\frac{1}{2},j,k}^n]\\
        + & \frac{1}{2 \Delta y} [( {(\rho v)}_{i, j+1, k}^{n} + {(\rho v)}_{i, j, k}^{n}) v_{i, j+\frac{1}{2},k}^n
                              -( {(\rho v)}_{i, j, k}^{n} + {(\rho v)}_{i, j-1, k}^{n}) v_{i, j-\frac{1}{2},k}^n]\\
        + & \frac{1}{2 \Delta z} [( {(\rho w)}_{i, j, k+1}^{n} + {(\rho w)}_{i, j-1, k+1}^{n}) v_{i, j, k+\frac{1}{2}}^n
                              -( {(\rho w)}_{i, j, k}^{n} + {(\rho w)}_{i, j-1, k}^{n}) v_{i, j, k-\frac{1}{2}}^n] \}\\
        - & \frac{\Delta t}{\Delta y}[p_{i, j, k}^{n} - p_{i, j-1, k}^{n}]
    \end{align*}

Divergence Compoments
---------------------
.. image:: figures/grid_discretization/y_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/y_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/y_mom_advec_z.PNG
  :width: 400

Z-Momentum Conservation
=======================

Difference Equation
-------------------
  
.. math::

    \begin{align*}
    (\rho w)_{i, j, k}^{n+1} = (\rho w)_{i, j, k}^{n}
    - \Delta t
        \{& \frac{1}{2 \Delta x} [( {(\rho u)}_{i+1, j, k}^{n} + {(\rho u)}_{i+1, j, k-1}^{n}) w_{i+\frac{1}{2},j,k}^n
                                 -( {(\rho u)}_{i, j, k}^{n} + {(\rho u)}_{i, j, k-1}^{n}) w_{i-\frac{1}{2},j,k}^n]\\
        + & \frac{1}{2 \Delta y} [( {(\rho v)}_{i, j+1, k}^{n} + {(\rho v)}_{i, j+1, k-1}^{n}) w_{i, j+\frac{1}{2},k}^n
                                 -( {(\rho v)}_{i, j, k}^{n} + {(\rho v)}_{i, j, k-1}^{n}) w_{i, j-\frac{1}{2},k}^n]\\
        + & \frac{1}{2 \Delta z} [( {(\rho w)}_{i, j, k+1}^{n} + {(\rho w)}_{i, j, k}^{n}) w_{i, j, k+\frac{1}{2}}^n
                                 -( {(\rho w)}_{i, j, k}^{n} + {(\rho w)}_{i, j, k-1}^{n}) w_{i, j, k-\frac{1}{2}}^n] \}\\
        - & \frac{\Delta t}{\Delta z}[p_{i, j, k}^{n} - p_{i, j, k-1}^{n}] + \Delta t g \rho_{i, j, k-\frac{1}{2}}^n
    \end{align*}


Divergence Compoments
---------------------
.. image:: figures/grid_discretization/z_mom_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/z_mom_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/z_mom_advec_z.PNG
  :width: 400


Energy Conservation
===================

Difference Equation
-------------------
   
.. math::

    \begin{align*}
        (\rho \theta)_{i, j, k}^{n+1} = (\rho \theta)_{i, j, k}^{n}
        - \Delta t
        \{& \frac{1}{\Delta x} [{(\rho u)}_{i+1, j, k}^{n} \theta_{i+\frac{1}{2},j,k}^n
                                 -{(\rho u)}_{i, j, k}^{n}  \theta_{i-\frac{1}{2},j,k}^n]\\
        + & \frac{1}{\Delta y} [{(\rho v)}_{i, j+1, k}^{n} \theta_{i, j+\frac{1}{2},k}^n
                                 -{(\rho v)}_{i, j, k}^{n} \theta_{i, j-\frac{1}{2},k}^n]\\
        + & \frac{1}{\Delta z} [{(\rho w)}_{i, j, k+1}^{n} \theta_{i, j, k+\frac{1}{2}}^n
                                 -{(\rho w)}_{i, j, k}^{n} \theta_{i, j, k-\frac{1}{2}}^n] \}
    \end{align*}


Divergence Compoments
---------------------
.. image:: figures/grid_discretization/temp_advec_x.PNG
  :width: 400
.. image:: figures/grid_discretization/temp_advec_y.PNG
  :width: 400
.. image:: figures/grid_discretization/temp_advec_z.PNG
  :width: 400

Diagnostic Variables
====================
 
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

:math:`[\rho, u, v, w, \theta]`,  :math:`m = i, j, k`,  and :math:`U_d = [u, v, w]` for :math:`[x, y, z]` differences respectively.

.. math::

    \begin{align*}
        q_{m+\frac{1}{2}}|^{2^{nd}} &= \frac{1}{2} (q_{m+1} + q_{m})\\
        q_{m+\frac{1}{2}}|^{4^{th}} &=   \frac{7}{12} (q_{m+1} + q_{m})
                                      - \frac{1}{12} (q_{m+2} + q_{m-1})\\
        q_{m+\frac{1}{2}}|^{6^{th}} &=   \frac{37}{60} (q_{m+1} + q_{m})
                                      - \frac{2}{15} (q_{m+2} + q_{m-1})
                                      + \frac{1}{60} (q_{m+1} + q_{m-2})\\
        q_{m+\frac{1}{2}}|^{3^{rd}} &=   q_{m+\frac{1}{2}}|^{4^{th}}
                                       + \frac{U_d}{|U_d|} \frac{1}{12} [ 
                                               (q_{m+2} + q_{m-1})
                                              - 3(q_{m+1} + q_{m})]\\
        q_{m+\frac{1}{2}}|^{5^{th}} &=   q_{m+\frac{1}{2}}|^{6^{th}}
                                       - \frac{U_d}{|U_d|} \frac{1}{60} [
                                                (q_{m+3} + q_{m-2}) 
                                          -  5(q_{m+2} + q_{m-1})
                                          + 10(q_{m+1} + q_{m})]\\
    \end{align*}

Ref: Skamarock, W. C., Klemp, J. B., Dudhia, J., Gill, D. O., Liu, Z., Berner, J., ... Huang, X. -yu. (2019). A Description of the Advanced Research WRF Model Version 4 (No. NCAR/TN-556+STR). doi:10.5065/1dfh-6p97
`doi:10.5065/1dfh-6p97 <http://dx.doi.org/10.5065/1dfh-6p97>`_
