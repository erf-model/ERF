Finite difference discretization of Euler/Navier-Stokes equations for the ERF model
===================================================================================

+------------------------------------------------------------------------------------------------+------------+
| |image0|                                                                                       | |image1|   |
+================================================================================================+============+
| Figure 1. X-Y (left) and Y-Z (right) staggered grids indicating where variables are located.   |
+------------------------------------------------------------------------------------------------+------------+

Mass Conservation
=================

.. math::

   \begin{matrix}
   \rho_{i,j,k}^{n + 1} & = & \rho_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{\Delta x} \right.\  & \left\lbrack \left( \text{ρu} \right)_{i + 1,j,k}^{n} \right.\  & - & \left. \ \left( \text{ρu} \right)_{i.j,k}^{n} \right\rbrack \\
    & & & & & + \frac{1}{\Delta y} & \left\lbrack \left( \text{ρv} \right)_{i,j + 1,k}^{n} \right.\  & - & \left. \ \left( \text{ρv} \right)_{i.j,k}^{n} \right\rbrack \\
    & & & & & + \frac{1}{\Delta z} & \left\lbrack \left( \text{ρw} \right)_{i,j,k + 1}^{n} \right.\  & - & \left. \ \left. \ \left( \text{ρw} \right)_{i,j,k}^{n} \right\rbrack \right\} \\
   \end{matrix}

+---------------------------------------------------------------------------------------------------+------------+------------+
| |image2|                                                                                          | |image3|   | |image4|   |
+===================================================================================================+============+============+
| Figure 2. Divergence components: x direction (left), y direction (center), z direction (right).   |
+---------------------------------------------------------------------------------------------------+------------+------------+

Advection Contribution to DNS/LES
=================================

Momentum Conservation – U Momentum
----------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρu} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρu} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \text{ρu} \right)_{i + 1,j,k}^{n} + \left( \text{ρu} \right)_{i,j,k}^{n} \right)u_{i + \frac{1}{2},j,k}^{n} - \left( \left( \text{ρu} \right)_{i,j,k}^{n} + \left( \text{ρu} \right)_{i - 1,j,k}^{n} \right)u_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{2\Delta y}\left\lbrack \left( \left( \text{ρv} \right)_{i,j + 1,k}^{n} + \left( \text{ρv} \right)_{i - 1,j + 1,k}^{n} \right)u_{i,j + \frac{1}{2},k}^{n} - \left( \left( \text{ρv} \right)_{i,j,k}^{n} + \left( \text{ρv} \right)_{i - 1,j,k}^{n} \right)u_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & \  + \ \left. \ \frac{1}{2\Delta z}\left\lbrack \left( \left( \text{ρw} \right)_{i,j,k + 1}^{n} + \left( \text{ρw} \right)_{i - 1,j,k + 1}^{n} \right)u_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \text{ρw} \right)_{i,j,k}^{n} + \left( \text{ρw} \right)_{i - 1,j,k}^{n} \right)u_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
    & & & & & - \frac{\Delta t}{\Delta x}\left\lbrack p_{i,\ j,\ k}^{n} - p_{i - 1,\ j,\ k}^{n} \right\rbrack \\
   \end{matrix}

+----------------------------------------------------------------------------------------+------------+------------+
| |image5|                                                                               | |image6|   | |image7|   |
+========================================================================================+============+============+
| Figure 3. U momentum: x advection (left), y advection (center), z advection (right).   |
+----------------------------------------------------------------------------------------+------------+------------+

Momentum Conservation – V Momentum
----------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρv} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρv} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \text{ρu} \right)_{i + 1,j,k}^{n} + \left( \text{ρu} \right)_{i + 1,j - 1,k}^{n} \right)v_{i + \frac{1}{2},j,k}^{n} - \left( \left( \text{ρu} \right)_{i,j,k}^{n} + \left( \text{ρu} \right)_{i,j - 1,k}^{n} \right)v_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{2\Delta y}\left\lbrack \left( \left( \text{ρv} \right)_{i,j + 1,k}^{n} + \left( \text{ρv} \right)_{i,j,k}^{n} \right)v_{i,j + \frac{1}{2},k}^{n} - \left( \left( \text{ρv} \right)_{i,j,k}^{n} + \left( \text{ρv} \right)_{i,j - 1,k}^{n} \right)v_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{2\Delta z}\left\lbrack \left( \left( \text{ρw} \right)_{i,j,k + 1}^{n} + \left( \text{ρw} \right)_{i,j - 1,k + 1}^{n} \right)v_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \text{ρw} \right)_{i,j,k}^{n} + \left( \text{ρw} \right)_{i,j - 1,k}^{n} \right)v_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
    & & & & & - \frac{\Delta t}{\Delta y}\left\lbrack p_{i,j,\ k}^{n} - p_{i,\ j - 1,\ k}^{n} \right\rbrack \\
   \end{matrix}

+----------------------------------------------------------------------------------------+------------+-------------+
| |image8|                                                                               | |image9|   | |image10|   |
+========================================================================================+============+=============+
| Figure 4. V momentum: x advection (left), y advection (center), z advection (right).   |
+----------------------------------------------------------------------------------------+------------+-------------+

Momentum Conservation – W Momentum
----------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρw} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρw} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{2\Delta x}\ \left\lbrack \left( \left( \text{ρu} \right)_{i + 1,j,k}^{n} + \left( \text{ρu} \right)_{i + 1,j,k - 1}^{n} \right)w_{i + \frac{1}{2},j,k}^{n} - \left( \left( \text{ρu} \right)_{i,j,k}^{n} + \left( \text{ρu} \right)_{i,j,k - 1}^{n} \right)w_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{2\Delta y}\left\lbrack \left( \left( \text{ρv} \right)_{i,j + 1,k}^{n} + \left( \text{ρv} \right)_{i,j + 1,k - 1}^{n} \right)w_{i,j + \frac{1}{2},k}^{n} - \left( \left( \text{ρv} \right)_{i,j,k}^{n} + \left( \text{ρv} \right)_{i,j,k - 1}^{n} \right)w_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{2\Delta z}\left\lbrack \left( \left( \text{ρw} \right)_{i,j,k + 1}^{n} + \left( \text{ρw} \right)_{i,j,k}^{n} \right)w_{i,j,k + \frac{1}{2}}^{n} - \left( \left( \text{ρw} \right)_{i,j,k}^{n} + \left( \text{ρw} \right)_{i,j,k - 1}^{n} \right)w_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
    & & & & & - \frac{\Delta t}{\Delta z}\left\lbrack p_{i,\ j,\ k}^{n} - p_{i,\ j,\ \ k - 1}^{n} \right\rbrack\  + \ \Delta\text{t }\left\lbrack \rho_{i,j,k - \ \frac{1}{2}}^{n} \right\rbrack\text{g} \\
   \end{matrix}

+----------------------------------------------------------------------------------------+-------------+-------------+
| |image11|                                                                              | |image12|   | |image13|   |
+========================================================================================+=============+=============+
| Figure 5. W momentum: x advection (left), y advection (center), z advection (right).   |
+----------------------------------------------------------------------------------------+-------------+-------------+

Energy Conservation – Potential Temperature Advection 
------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρθ} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρθ} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{\Delta x}\ \left\lbrack \left( \text{ρu} \right)_{i + 1,j,k}^{n}\text{ θ}_{i + \frac{1}{2},j,k}^{n} - \left( \text{ρu} \right)_{i,j,k}^{n}\text{ θ}_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{\Delta y}\left\lbrack \left( \text{ρv} \right)_{i,j + 1,k}^{n}\text{ θ}_{i,j + \frac{1}{2},k}^{n} - \left( \text{ρv} \right)_{i,j,k\ }^{n}\theta_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{\Delta z}\left\lbrack \left( \text{ρw} \right)_{i,j,k + 1\ }^{n}\theta_{i,j,k + \frac{1}{2}}^{n} - \left( \text{ρw} \right)_{i,j,k}^{n}\text{ θ}_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
   \end{matrix}

+------------------------------------------------------------------------------------------------------------+-------------+-------------+
| |image14|                                                                                                  | |image15|   | |image16|   |
+============================================================================================================+=============+=============+
| Figure 6. Potential temperature equation: x advection (left), y advection (center), z advection (right).   |
+------------------------------------------------------------------------------------------------------------+-------------+-------------+

Scalar Conservation – Scalar Advection 
---------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρS} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρS} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{\Delta x}\ \left\lbrack \left( \text{ρu} \right)_{i + 1,j,k}^{n}\text{ S}_{i + \frac{1}{2},j,k}^{n} - \left( \text{ρu} \right)_{i,j,k}^{n}\text{ S}_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{\Delta y}\left\lbrack \left( \text{ρv} \right)_{i,j + 1,k}^{n}\text{ S}_{i,j + \frac{1}{2},k}^{n} - \left( \text{ρv} \right)_{i,j,k\ }^{n}S_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{\Delta z}\left\lbrack \left( \text{ρw} \right)_{i,j,k + 1\ }^{n}S_{i,j,k + \frac{1}{2}}^{n} - \left( \text{ρw} \right)_{i,j,k}^{n}\text{ S}_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
   \end{matrix}

+----------------------------------------------------------------------------------------------------------+-------------+-------------+
| |image17|                                                                                                | |image18|   | |image19|   |
+==========================================================================================================+=============+=============+
| Figure 7. Scalar conservation equation: x advection (left), y advection (center), z advection (right).   |
+----------------------------------------------------------------------------------------------------------+-------------+-------------+

Pressure Diagnostic
-------------------

    This is included here but pressure can be diagnosed from other
    variables solved for irrespective of whether the effect of advection
    or diffusion or both or none is considered.

.. math:: \operatorname{}{\text{ }p_{i,\ j,\ k}^{n}} = \rho_{i,\ j,\ k}^{n}R_{d}\theta_{i,\ j,\ k}^{n}\left( \frac{p_{i,\ j,\ k}^{n}}{p_{0}} \right)^{\frac{R_{d}}{c_{p}}}

.. math:: \operatorname{}{\text{ }p_{i,\ j,\ k}^{n}} = \left\lbrack \rho_{i,\ j,\ k}^{n}R_{d}\theta_{i,\ j,\ k}^{n}\left( \frac{1}{p_{0}} \right)^{\frac{R_{d}}{c_{p}}} \right\rbrack^{\gamma}

.. math:: q = \left\lbrack \rho,\ u,v,w,\theta \right\rbrack\ \ \ \ \ \ \ \ \ \ \ \ \ \ \ m = i,j,k\ \ \ \ \ \ \ \ \ \ \ \ U_{d} = \left\lbrack u,v,w \right\rbrack\mathrm{\text{\ for\ }}\left\lbrack x,y,z \right\rbrack\ \mathrm{\text{differences}},\ \mathrm{\text{respecively}}\text{\ \ \ }

.. math::

   \begin{matrix}
   \left. \ q_{m + \frac{1}{2}} \right|^{2nd} & = & \frac{1}{2}\left( q_{m + 1} + q_{m} \right) & & & & & & & \\
   \left. \ q_{m + \frac{1}{2}} \right|^{4th} & = & \frac{7}{12}\left( q_{m + 1} + q_{m} \right) & - & \frac{1}{12}\left( q_{m + 2} + q_{m - 1} \right) & & & & & \\
   \left. \ q_{m + \frac{1}{2}} \right|^{6th} & = & \frac{37}{60}\left( q_{m + 1} + q_{m} \right) & - & \frac{2}{15}\left( q_{m + 2} + q_{m - 1} \right) & + & \frac{1}{60}\left( q_{m + 3} + q_{m - 2} \right) & & & \\
    & & & & & & & & & \\
   \left. \ q_{m + \frac{1}{2}} \right|^{3rd} & = & \left. \ q_{m + \frac{1}{2}} \right|^{4th} & + & \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{12}\left\lbrack \left( q_{m + 2} + q_{m - 1} \right) \right.\  & - & 3\left. \ \left( q_{m + 1} + q_{m} \right) \right\rbrack & & & \\
    & & & & & & & & & \\
   \left. \ q_{m + \frac{1}{2}} \right|^{5th} & = & \left. \ q_{m + \frac{1}{2}} \right|^{6th} & - & \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{60}\left\lbrack \left( q_{m + 3} + q_{m - 2} \right) \right.\  & - & 5\left( q_{m + 2} + q_{m - 1} \right) & + & 10\left. \ \left( q_{m + 1} + q_{m} \right) \right\rbrack & \\
   \end{matrix}

.. math::

   \begin{matrix}
   \left. \ q_{m - \frac{1}{2}} \right|^{2nd} & = & \frac{1}{2}\left( q_{m} + q_{m - 1} \right) & & & & & & & \\
   \left. \ q_{m - \frac{1}{2}} \right|^{4th} & = & \frac{7}{12}\left( q_{m} + q_{m - 1} \right) & - & \frac{1}{12}\left( q_{m + 1} + q_{m - 2} \right) & & & & & \\
   \left. \ q_{m - \frac{1}{2}} \right|^{6th} & = & \frac{37}{60}\left( q_{m} + q_{m - 1} \right) & - & \frac{2}{15}\left( q_{m + 1} + q_{m - 2} \right) & + & \frac{1}{60}\left( q_{m + 2} + q_{m - 3} \right) & & & \\
    & & & & & & & & & \\
   \left. \ q_{m - \frac{1}{2}} \right|^{3rd} & = & \left. \ q_{m - \frac{1}{2}} \right|^{4th} & + & \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{12}\left\lbrack \left( q_{m + 1} + q_{m - 2} \right) \right.\  & - & 3\left. \ \left( q_{m} + q_{m - 1} \right) \right\rbrack & & & \\
    & & & & & & & & & \\
   \left. \ q_{m - \frac{1}{2}} \right|^{5th} & = & \left. \ q_{m - \frac{1}{2}} \right|^{6th} & - & \frac{U_{d}}{\left| U_{d} \right|}\frac{1}{60}\left\lbrack \left( q_{m + 2} + q_{m - 1} \right) \right.\  & - & 5\left( q_{m + 1} + q_{m - 2} \right) & + & 10\left. \ \left( q_{m} + q_{m - 1} \right) \right\rbrack & \\
   \end{matrix}

Momentum, Thermal, and Scalar Diffusion Contribution to DNS
===========================================================

Strain Rate Tensor
------------------

+-------------------------------------------+
| |image20|                                 |
+===========================================+
| Figure 8. Strain rate tensor schematic.   |
+-------------------------------------------+

Momentum Conservation – U Momentum viscous stress divergence
------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρu} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρu} \right)_{i,j,k}^{n} & + & \Delta t\ 2\rho_{i,j,k}\nu & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack S_{11,i + \frac{1}{2}} - S_{11,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack S_{12,j + \frac{1}{2}} - S_{12,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack S_{13,k + \frac{1}{2}} - S_{13,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   S_{11,i + \frac{1}{2}} = \frac{1}{\Delta x}\left( u_{i + 1,j,k} - u_{i,j,k} \right) \\
   S_{11,i - \frac{1}{2}} = \frac{1}{\Delta x}\left( u_{i,j,k} - u_{i - 1,j,k} \right) \\
   S_{12,j + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j + 1,k} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( v_{i,j + 1,k} - v_{i - 1,j + 1,k} \right) \right\rbrack \\
   S_{12,j - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   S_{13,k + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k + 1} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( w_{i,j,k + 1} - w_{i - 1,j,k + 1} \right) \right\rbrack \\
   S_{13,k - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   \end{matrix}

+----------------------------------------------------+-------------+
| |image21|                                          | |image22|   |
+====================================================+=============+
| Figure 9. Viscous stress divergence – U momentum   |
+----------------------------------------------------+-------------+

Momentum Conservation – V Momentum viscous stress divergence
------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρv} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρv} \right)_{i,j,k}^{n} & + & \Delta t\ 2\rho_{i,j,k}\nu & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack S_{21,i + \frac{1}{2}} - S_{21,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack S_{22,j + \frac{1}{2}} - S_{22,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack S_{23,k + \frac{1}{2}} - S_{23,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   S_{21,i + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i + 1,j,k} - u_{i + 1,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i + 1,j,k} - v_{i,j,k} \right) \right\rbrack \\
   S_{21,i - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   S_{22,j + \frac{1}{2}} = \frac{1}{\Delta y}\left( v_{i,j + 1,k} - v_{i,j,k} \right) \\
   S_{22,j - \frac{1}{2}} = \frac{1}{\Delta y}\left( v_{i,j,k} - v_{i,j - 1,k} \right) \\
   S_{23,k + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k + 1} - v_{i,j,k} \right) + \frac{1}{\Delta y}\left( w_{i,j,k + 1} - w_{i,j - 1,k + 1} \right) \right\rbrack \\
   S_{23,k - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   \end{matrix}

+-----------------------------------------------------+-------------+
| |image23|                                           | |image24|   |
+=====================================================+=============+
| Figure 10. Viscous stress divergence – V momentum   |
+-----------------------------------------------------+-------------+

Momentum Conservation – W Momentum viscous stress divergence
------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρw} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρw} \right)_{i,j,k}^{n} & + & \Delta t\ 2\rho_{i,j,k}\nu & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack S_{31,i + \frac{1}{2}} - S_{31,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack S_{32,j + \frac{1}{2}} - S_{32,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack S_{33,k + \frac{1}{2}} - S_{33,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   S_{31,i + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i + 1,j,k} - u_{i + 1,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i + 1,j,k} - w_{i,j,k} \right) \right\rbrack \\
   S_{31,i - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   S_{32,j + \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j + 1,k} - v_{i,j + 1,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j + 1,k} - w_{i,j,k} \right) \right\rbrack \\
   S_{32,j - \frac{1}{2}} = \frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   S_{33,k + \frac{1}{2}} = \frac{1}{\Delta z}\left( w_{i,j,k + 1} - w_{i,j,k} \right) \\
   S_{33,k - \frac{1}{2}} = \frac{1}{\Delta z}\left( w_{i,j,k} - w_{i,j,k - 1} \right) \\
   \end{matrix}

+-----------------------------------------------------+-------------+
| |image25|                                           | |image26|   |
+=====================================================+=============+
| Figure 11. Viscous stress divergence – W momentum   |
+-----------------------------------------------------+-------------+

Energy Conservation – Potential Temperature Diffusion 
------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρθ} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρθ} \right)_{i,j,k}^{n} & + & \Delta t\rho_{i,j,k}\alpha_{T} & \left\{ \frac{1}{{\Delta x}^{2}}\ \left\lbrack \theta_{i + 1,j,k}^{n} - \ {2\theta}_{i,j,k}^{n} + \ \theta_{i - 1,j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{{\Delta y}^{2}}\left\lbrack \theta_{i,j + 1,k}^{n} - \ 2\theta_{i,j,k}^{n} + \ \theta_{i,j - 1,k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{{\Delta z}^{2}}\left\lbrack \theta_{i,j,k + 1}^{n} - \ {2\theta}_{i,j,k}^{n} + \ \theta_{i,j,k - 1}^{n} \right\rbrack \right\} \\
   \end{matrix}

Scalar Conservation – Scalar Diffusion 
---------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρC} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρC} \right)_{i,j,k}^{n} & + & \Delta t\rho_{i,j,k}\alpha_{S} & \left\{ \frac{1}{{\Delta x}^{2}}\ \left\lbrack C_{i + 1,j,k}^{n} - \ {2C}_{i,j,k}^{n} + \ C_{i - 1,j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{{\Delta y}^{2}}\left\lbrack C_{i,j + 1,k}^{n} - \ 2C_{i,j,k}^{n} + \ C_{i,j - 1,k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{{\Delta z}^{2}}\left\lbrack C_{i,j,k + 1}^{n} - \ {2C}_{i,j,k}^{n} + \ C_{i,j,k - 1}^{n} \right\rbrack \right\} \\
   \end{matrix}

Momentum, Thermal, and Scalar Diffusion Contribution to LES
===========================================================

Strain Rate and Eddy Viscosity
------------------------------

    The goal is to compute eddy viscosity at the *cell centers* and
    interpolated them to the edges. Refer again to the strain rate
    tensor schematic.

    |image27|

.. math::

   \begin{matrix}
   S_{11} = S_{11i + \frac{1}{2}} \\
   S_{22} = S_{22j + \frac{1}{2}} \\
   \begin{matrix}
   S_{33} = S_{33k + \frac{1}{2}} \\
   S_{12} = \frac{1}{4}\left\lbrack S_{12i,j - \frac{1}{2}} + S_{12i,j + \frac{1}{2}} + S_{12i + 1,j - \frac{1}{2}} + S_{12i + 1,j + \frac{1}{2}} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   \begin{matrix}
   S_{21} = \frac{1}{4}\left\lbrack S_{21i - \frac{1}{2},j} + S_{21i + \frac{1}{2},j} + S_{21i - \frac{1}{2},j + 1} + S_{21i + \frac{1}{2},j + 1} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   \begin{matrix}
   S_{13} = \frac{1}{4}\left\lbrack S_{13i,k - \frac{1}{2}} + S_{13i,k + \frac{1}{2}} + S_{13i + 1,k - \frac{1}{2}} + S_{13i + 1,k + \frac{1}{2}} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   S_{31} = \frac{1}{4}\left\lbrack S_{31i - \frac{1}{2},k} + S_{31i + \frac{1}{2},k} + S_{31i - \frac{1}{2},k + 1} + S_{31i + \frac{1}{2},k + 1} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   \begin{matrix}
   S_{23} = \frac{1}{4}\left\lbrack S_{23j,k - \frac{1}{2}} + S_{23j,k + \frac{1}{2}} + S_{23j + 1,k - \frac{1}{2}} + S_{23j + 1,k + \frac{1}{2}} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   S_{32} = \frac{1}{4}\left\lbrack S_{32j - \frac{1}{2},k} + S_{32j + \frac{1}{2},k} + S_{32j - \frac{1}{2},k + 1} + S_{32j + \frac{1}{2},k + 1} \right\rbrack = Average\ of\ the\ 4\ edges\ surrouding\ the\ cell \\
   \end{matrix}

    
    Note that:

.. math:: S_{12} = S_{21}

.. math:: S_{13} = S_{31}

.. math:: S_{23} = S_{32}

:math:`K_{i,j,k} = {- 2\left( C_{S} \right)^{2}\rho_{i,j,k}\left( {2S}_{\text{mn}}S_{\text{mn}} \right)}^{\frac{1}{2}}`,
where

.. math::

   \begin{matrix}
   S_{\text{mn}}S_{\text{mn}} = S_{11}^{2} + S_{22}^{2} + S_{33}^{2} + S_{12}^{2} + S_{13}^{2} + S_{23}^{2} + S_{21}^{2} + S_{31}^{2} + S_{32}^{2} \\
   \end{matrix}

Owing to symmetry we need to compute 6 of the 9 tensor components.

+------------------------------+
| |image28|                    |
+==============================+
| Figure 12. Eddy viscosity.   |
+------------------------------+

The interpolated values of eddy-viscosity at the edges are the average
of the values at the centers of the 4 cells the edge is part of.

.. math::

   \begin{matrix}
   \begin{matrix}
   K_{i + \frac{1}{2},j - \frac{1}{2},k} = \frac{1}{4}\left\lbrack K_{i,j - 1,k} + K_{i,j,k} + K_{i + 1,j - 1,k} + K_{i + 1,j,k} \right\rbrack \\
   \begin{matrix}
   K_{i + \frac{1}{2},j + \frac{1}{2},k} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j + 1,k} + K_{i + 1,j,k} + K_{i + 1,j + 1,k} \right\rbrack \\
   \begin{matrix}
   K_{i + \frac{1}{2},j,k - \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k - 1} + K_{i + 1,j,k} + K_{i + 1,j,k - 1} \right\rbrack \\
   \begin{matrix}
   K_{i + \frac{1}{2},j,k + \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k + 1} + K_{i,j,k} + K_{i + 1,j,k + 1} + K_{i + 1,j,k} \right\rbrack \\
   \begin{matrix}
   K_{i,j + \frac{1}{2},k - \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k - 1} + K_{i,j + 1,k} + K_{i,j + 1,k - 1} \right\rbrack \\
   \begin{matrix}
   K_{i,j + \frac{1}{2},k + \frac{1}{2}} = \frac{1}{4}\left\lbrack K_{i,j,k} + K_{i,j,k + 1} + K_{i,j + 1,k} + K_{i,j + 1,k + 1} \right\rbrack \\
   \end{matrix} \\
   \end{matrix} \\
   \end{matrix} \\
   \end{matrix} \\
   \end{matrix} \\
   \end{matrix} \\
   \end{matrix}

Momentum Conservation – U Momentum - subfilter stress divergence
----------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρu} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρu} \right)_{i,j,k}^{n} & - & \Delta\text{t } & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack \tau_{11,i + \frac{1}{2}} - \tau_{11,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack \tau_{12,j + \frac{1}{2}} - \tau_{12,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack \tau_{13,k + \frac{1}{2}} - \tau_{13,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   \tau_{11,i + \frac{1}{2}}{= K_{i,j,k}\text{ S}}_{11,i + \frac{1}{2}} = K_{i,j,k}\frac{1}{\Delta x}\left( u_{i + 1,j,k} - u_{i,j,k} \right) \\
   \tau_{11,i - \frac{1}{2}} = {K_{i - 1,j,k}\text{ S}}_{11,i - \frac{1}{2}} = K_{i - 1,j,k}\frac{1}{\Delta x}\left( u_{i,j,k} - u_{i - 1,j,k} \right) \\
   \tau_{12,j + \frac{1}{2}} = {K_{i - \frac{1}{2},j + \frac{1}{2},k}\text{ S}}_{12,j + \frac{1}{2}} = K_{i - \frac{1}{2},j + \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j + 1,k} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( v_{i,j + 1,k} - v_{i - 1,j + 1,k} \right) \right\rbrack \\
   \tau_{12,j - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\ S_{12,j - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   \tau_{13,k + \frac{1}{2}} = {K_{i - \frac{1}{2},j,k + \frac{1}{2}}\text{ S}}_{13,k + \frac{1}{2}} = K_{i - \frac{1}{2},j,k + \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k + 1} - u_{i,j,k} \right) + \frac{1}{\Delta x}\left( w_{i,j,k + 1} - w_{i - 1,j,k + 1} \right) \right\rbrack \\
   {\tau_{13,k - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\text{ S}}_{13,k - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   \end{matrix}

Momentum Conservation – V Momentum - subfilter stress divergence
----------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρv} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρv} \right)_{i,j,k}^{n} & - & \Delta\text{t } & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack \tau_{21,i + \frac{1}{2}} - \tau_{21,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack \tau_{22,j + \frac{1}{2}} - \tau_{22,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack \tau_{23,k + \frac{1}{2}} - \tau_{23,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   {{\tau_{21,i + \frac{1}{2}} = K}_{i + \frac{1}{2},j - \frac{1}{2},k}\text{ S}}_{21,i + \frac{1}{2}} = K_{i + \frac{1}{2},j - \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i + 1,j,k} - u_{i + 1,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i + 1,j,k} - v_{i,j,k} \right) \right\rbrack \\
   {{\tau_{21,i - \frac{1}{2}} = K}_{i - \frac{1}{2},j - \frac{1}{2},k}\text{ S}}_{21,i - \frac{1}{2}} = K_{i - \frac{1}{2},j - \frac{1}{2},k}\frac{1}{2}\left\lbrack \frac{1}{\Delta y}\left( u_{i,j,k} - u_{i,j - 1,k} \right) + \frac{1}{\Delta x}\left( v_{i,j,k} - v_{i - 1,j,k} \right) \right\rbrack \\
   {\tau_{22,j + \frac{1}{2}} = K}_{i,j,k}\ S_{22,j + \frac{1}{2}} = K_{i,j,k}\frac{1}{\Delta y}\left( v_{i,j + 1,k} - v_{i,j,k} \right) \\
   \tau_{22,j - \frac{1}{2}} = K_{i,j - 1,k}\ S_{22,j - \frac{1}{2}} = K_{i,j - 1,k}\frac{1}{\Delta y}\left( v_{i,j,k} - v_{i,j - 1,k} \right) \\
   \tau_{23,k + \frac{1}{2}} = K_{i,j - \frac{1}{2},k + \frac{1}{2}}\ S_{23,k + \frac{1}{2}} = K_{i,j - \frac{1}{2},k + \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k + 1} - v_{i,j,k} \right) + \frac{1}{\Delta y}\left( w_{i,j,k + 1} - w_{i,j - 1,k + 1} \right) \right\rbrack \\
   \tau_{23,k - \frac{1}{2}} = K_{i,j - \frac{1}{2}k - \frac{1}{2}}\text{ S}_{23,k - \frac{1}{2}} = K_{i,j - \frac{1}{2},k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   \end{matrix}

Momentum Conservation – W Momentum - subfilter stress divergence
----------------------------------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρw} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρw} \right)_{i,j,k}^{n} & - & \Delta\text{t } & \left. \ \left\{ \frac{1}{\Delta x}\ \left\lbrack \tau_{31,i + \frac{1}{2}} - \tau_{31,i - \frac{1}{2}} \right\rbrack \right.\  + \frac{1}{\Delta y}\ \left\lbrack \tau_{32,j + \frac{1}{2}} - \tau_{32,j - \frac{1}{2}} \right\rbrack + \frac{1}{\Delta z}\ \left\lbrack \tau_{33,k + \frac{1}{2}} - \tau_{33,k - \frac{1}{2}} \right\rbrack \right\} \\
   \end{matrix}

.. math::

   \begin{matrix}
   {K_{i + \frac{1}{2},j,k - \frac{1}{2}}\text{ S}}_{31,i + \frac{1}{2}} = K_{i + \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i + 1,j,k} - u_{i + 1,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i + 1,j,k} - w_{i,j,k} \right) \right\rbrack \\
   K_{i - \frac{1}{2},j,k - \frac{1}{2}}\ S_{31,i - \frac{1}{2}} = K_{i - \frac{1}{2},j,k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( u_{i,j,k} - u_{i,j,k - 1} \right) + \frac{1}{\Delta x}\left( w_{i,j,k} - w_{i - 1,j,k} \right) \right\rbrack \\
   K_{i,j + \frac{1}{2},k - \frac{1}{2}}\ S_{32,j + \frac{1}{2}} = K_{i,j + \frac{1}{2},k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j + 1,k} - v_{i,j + 1,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j + 1,k} - w_{i,j,k} \right) \right\rbrack \\
   K_{i,j - \frac{1}{2},k - \frac{1}{2}}S_{32,j - \frac{1}{2}} = K_{i,j - \frac{1}{2},k - \frac{1}{2}}\frac{1}{2}\left\lbrack \frac{1}{\Delta z}\left( v_{i,j,k} - v_{i,j,k - 1} \right) + \frac{1}{\Delta y}\left( w_{i,j,k} - w_{i,j - 1,k} \right) \right\rbrack \\
   K_{i,j,k}S_{33,k + \frac{1}{2}} = K_{i,j,k}\frac{1}{\Delta z}\left( w_{i,j,k + 1} - w_{i,j,k} \right) \\
   K_{i,jk - 1}S_{33,k - \frac{1}{2}} = K_{i,j,k - 1}\frac{1}{\Delta z}\left( w_{i,j,k} - w_{i,j,k - 1} \right) \\
   \end{matrix}

Energy Conservation- Subgrid heat flux
--------------------------------------

.. math::

   \begin{matrix}
   \left( \text{ρθ} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρθ} \right)_{i,j,k}^{n} & - & \Delta t & \left\{ \frac{1}{\Delta x}\ \left\lbrack {\vartheta_{1}}_{i + \frac{1}{2},j,k}^{n} - {\ \vartheta_{1}}_{i - \frac{1}{2},j,k}^{n} \right\rbrack \right.\  \\
    & & & & & + \frac{1}{\Delta y}\left\lbrack {\ \vartheta_{2}}_{i,j + \frac{1}{2},k}^{n} - {\vartheta_{2}}_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \left. \ \frac{1}{\Delta z}\left\lbrack {\vartheta_{3}}_{i,j,k + \frac{1}{2}}^{n} - {\ \vartheta_{3}}_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \right\} \\
   \end{matrix}

.. math:: \vartheta_{1i,j,k} = K_{i,j,k}\frac{1}{2\Delta x}\ \left\lbrack \text{ θ}_{i + 1,j,k}^{n} - \text{ θ}_{i - 1,j,k}^{n} \right\rbrack

.. math:: \vartheta_{2i,j,k} = K_{i,j,k}\frac{1}{2\Delta y}\ \left\lbrack \text{ θ}_{i,j + 1,k}^{n} - \text{ θ}_{i,j - 1,k}^{n} \right\rbrack

.. math:: \vartheta_{3i,jk} = K_{i,j,k}\frac{1}{2\Delta z}\ \left\lbrack \text{ θ}_{i,j,k + 1}^{n} - \text{ θ}_{i,j,k - 1}^{n} \right\rbrack

.. math:: \vartheta_{1i + \frac{1}{2},j,k} = \frac{1}{2}\left( \vartheta_{1i + 1,j,k} + \vartheta_{1i,j,k} \right)

.. math:: \vartheta_{1i - \frac{1}{2},j,k} = \frac{1}{2}\left( \vartheta_{1i,j,k} + \vartheta_{1i - 1,j,k} \right)

.. math:: \vartheta_{1i,j + \frac{1}{2},k} = \frac{1}{2}\left( \vartheta_{1i,j + 1,k} + \vartheta_{1i,j,k} \right)

.. math:: \vartheta_{1i,j - \frac{1}{2},k} = \frac{1}{2}\left( \vartheta_{1i,j,k} + \vartheta_{1i1,j - 1,k} \right)

.. math:: \vartheta_{1i,j,k + \frac{1}{2}} = \frac{1}{2}\left( \vartheta_{1i,j,k +} + \vartheta_{1i,j,k} \right)

.. math:: \vartheta_{1i,j,k - \frac{1}{2}} = \frac{1}{2}\left( \vartheta_{1i,j,k} + \vartheta_{1i1,j,k - 1} \right)

i.  .. rubric:: Scalar Conservation- Subgrid scalar flux
       :name: scalar-conservation--subgrid-scalar-flux

ii. .. rubric:: Prognostic Equation for Subgrid Kinetic Energy
       :name: prognostic-equation-for-subgrid-kinetic-energy

.. math::

   \begin{matrix}
   \left( \text{ρe} \right)_{i,j,k}^{n + 1} & = & \left( \text{ρe} \right)_{i,j,k}^{n} & - & \Delta t & \frac{1}{\Delta x}\ \left\lbrack \left( \text{ρu} \right)_{i + 1,j,k}^{n}\text{ e}_{i + \frac{1}{2},j,k}^{n} - \left( \text{ρu} \right)_{i,j,k}^{n}\text{ e}_{i - \frac{1}{2},j,k}^{n} \right\rbrack \\
    & & & & & + \frac{1}{\Delta y}\left\lbrack \left( \text{ρv} \right)_{i,j + 1,k}^{n}\text{ e}_{i,j + \frac{1}{2},k}^{n} - \left( \text{ρv} \right)_{i,j,k\ }^{n}e_{i,j - \frac{1}{2},k}^{n} \right\rbrack \\
    & & & & & + \frac{1}{\Delta z}\left\lbrack \left( \text{ρw} \right)_{i,j,k + 1\ }^{n}e_{i,j,k + \frac{1}{2}}^{n} - \left( \text{ρw} \right)_{i,j,k}^{n}\text{ e}_{i,j,k - \frac{1}{2}}^{n} \right\rbrack \\
    & & & & & + \frac{g}{\Theta}\vartheta_{3} - \tau_{\text{mn}}\frac{\partial u_{m}}{\partial x_{n}} - \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} - \epsilon \\
   \end{matrix}

.. math:: \vartheta_{i} = K_{H}\frac{\partial\theta}{\partial x_{i}}

.. math:: K_{H} = \left( 1 + 2\frac{\mathcal{l}}{\Delta s} \right)K_{M}

.. math:: K_{M} = 0.1\mathcal{l}e^{\frac{1}{2}} = 0.1\mathcal{l}e_{i,j,k}^{\frac{1}{2}}

.. math:: K_{Mi,j,k} = 0.1\mathcal{l}e_{i,j,k}^{\frac{1}{2}}

:math:`\mathcal{l} = \Delta s = \sqrt[3]{\Delta\text{x }\Delta\text{y }\Delta\text{z }}`,
convective case

.. math:: \mathcal{l} = 0.76\ e^{\frac{1}{2}}\left( \frac{g}{\Theta}\frac{\partial\theta}{\partial z} \right) = 0.76e_{i,j,k}^{\frac{1}{2}}\left\lbrack \frac{g}{\Theta}\frac{1}{2\Delta z}\left( \text{ θ}_{i,j,k + 1}^{n} - \text{ θ}_{i,j,k - 1}^{n} \right) \right\rbrack

.. math:: \vartheta_{1} = {K_{H}}_{i,j,k}\frac{1}{2\Delta x}\ \left\lbrack \text{ θ}_{i + 1,j,k}^{n} - \text{ θ}_{i - 1,j,k}^{n} \right\rbrack

.. math:: \vartheta_{2} = {K_{H}}_{i,j,k}\frac{1}{2\Delta y}\ \left\lbrack \text{ θ}_{i,j + 1,k}^{n} - \text{ θ}_{i,j - 1,k}^{n} \right\rbrack

.. math:: \vartheta_{3} = {K_{H}}_{i,j,k}\frac{1}{2\Delta z}\ \left\lbrack \text{ θ}_{i,j,k + 1}^{n} - \text{ θ}_{i,j,k - 1}^{n} \right\rbrack

.. math:: \frac{\partial\left\langle \left( u_{n}^{'}\rho e + u_{n}^{'}p^{'} \right) \right\rangle}{\partial x_{n}} = K_{i,j,k}\left\{ \frac{1}{2\Delta x}\ \left\lbrack \text{ e}_{i + 1,j,k}^{n} - \text{ e}_{i - 1j,k}^{n} \right\rbrack + \right.\ \frac{1}{2\Delta y}\left\lbrack \text{ e}_{i,j + 1,k}^{n} - e_{i,j - 1,k}^{n} \right\rbrack + \left. \ \frac{1}{2\Delta z}\left\lbrack e_{i,j,k + 1}^{n} - \text{ e}_{i,j,k - 1}^{n} \right\rbrack \right\}

.. math:: \epsilon = C_{\epsilon}\rho_{i,j,k}\frac{\left( e_{i,j,k} \right)^{\frac{3}{2}}}{\mathcal{l}}

.. math:: C_{\epsilon} = 0.19 + 0.51\frac{\mathcal{l}}{\Delta s}

.. math:: \tau_{\text{mn}}\frac{\partial u_{m}}{\partial x_{n}} = KS_{\text{mn}}\frac{\partial u_{m}}{\partial x_{n}} = KS_{\text{mn}}S_{\text{mn}} = K(S_{11}^{2} + S_{22}^{2} + S_{33}^{2} + S_{12}^{2} + S_{13}^{2} + S_{23}^{2} + S_{21}^{2} + S_{31}^{2} + S_{32}^{2})

+-------------+-------------+-------------+
| |image29|   | |image30|   | |image31|   |
+=============+=============+=============+
+-------------+-------------+-------------+

Figure 13. Subgrid kinetic energy.

.. |image0| image:: media/image1.png
   :width: 3.88814in
   :height: 2.18708in
.. |image1| image:: media/image2.png
   :width: 3.84709in
   :height: 2.16399in
.. |image2| image:: media/image3.PNG
   :width: 3.79259in
   :height: 2.13333in
.. |image3| image:: media/image4.PNG
   :width: 3.81927in
   :height: 2.14834in
.. |image4| image:: media/image5.PNG
   :width: 3.81889in
   :height: 2.14813in
.. |image5| image:: media/image6.png
   :width: 3.29000in
   :height: 2.18000in
.. |image6| image:: media/image7.png
   :width: 3.29000in
   :height: 2.18000in
.. |image7| image:: media/image8.png
   :width: 3.29000in
   :height: 2.18000in
.. |image8| image:: media/image9.png
   :width: 3.29000in
   :height: 2.18000in
.. |image9| image:: media/image10.png
   :width: 3.29000in
   :height: 2.18000in
.. |image10| image:: media/image11.png
   :width: 3.29000in
   :height: 2.18000in
.. |image11| image:: media/image12.png
   :width: 3.29000in
   :height: 2.18000in
.. |image12| image:: media/image13.png
   :width: 3.29000in
   :height: 2.18000in
.. |image13| image:: media/image14.png
   :width: 3.29000in
   :height: 2.18000in
.. |image14| image:: media/image15.png
   :width: 3.87552in
   :height: 2.17998in
.. |image15| image:: media/image16.png
   :width: 3.85514in
   :height: 2.16851in
.. |image16| image:: media/image17.png
   :width: 3.75371in
   :height: 2.11146in
.. |image17| image:: media/image18.png
   :width: 3.87551in
   :height: 2.17998in
.. |image18| image:: media/image19.png
   :width: 3.85514in
   :height: 2.16851in
.. |image19| image:: media/image20.png
   :width: 3.75371in
   :height: 2.11146in
.. |image20| image:: media/image21.png
   :width: 3.85000in
   :height: 1.97000in
.. |image21| image:: media/image22.png
   :width: 3.46000in
   :height: 2.18000in
.. |image22| image:: media/image23.png
   :width: 3.37000in
   :height: 2.19000in
.. |image23| image:: media/image24.png
   :width: 3.46000in
   :height: 2.18000in
.. |image24| image:: media/image25.png
   :width: 3.37000in
   :height: 2.19000in
.. |image25| image:: media/image26.png
   :width: 3.37000in
   :height: 2.19000in
.. |image26| image:: media/image27.png
   :width: 3.37000in
   :height: 2.19000in
.. |image27| image:: media/image21.png
   :width: 3.85000in
   :height: 1.97000in
.. |image28| image:: media/image28.png
   :width: 3.34000in
   :height: 1.94000in
.. |image29| image:: media/image29.png
   :width: 3.08000in
   :height: 2.04000in
.. |image30| image:: media/image30.png
   :width: 3.08000in
   :height: 2.04000in
.. |image31| image:: media/image31.png
   :width: 3.08000in
   :height: 2.05000in
