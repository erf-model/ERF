
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran


.. _Equations:

Map Factors (Dry)
=============================

ERF supports the use of isotropic projections only, e.g. Lambert conformal, polar stereographic and Mercator.
Like WRF, ERF implements the projections using map factors, which are defined as the ratio of the distance
in computational space to the corresponding distance on the earth's surface.

With map factors, the dry governing equations have the following form

.. math::
  \frac{\partial \rho}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u} m^{-1})
                                      - m   \frac{\partial (\rho w m^{-1})}{\partial z}

  \frac{\partial (\rho u})}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} u m^{-1})
                                           - m   \frac{\partial (\rho w            u m^{-1})}{\partial z}
                                           - m \frac{\partial p^\prime}{\partial x}
                                           + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho v})}{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} v m^{-1})
                                           - m   \frac{\partial (\rho w            v m^{-1})}{\partial z}
                                           - m   \frac{\partial p^\prime}{\partial y}
                                           + \nabla \cdot \tau + \mathbf{F},

  \frac{\partial (\rho w) }{\partial t} &= - m^2 \nabla_H \cdot (\rho \mathbf{u_H} w m^{-1})
                                           - m   \frac{\partial (\rho w            w m^{-1})}{\partial z}
                                           -     \frac{\partial p^\prime}{\partial z}
                                           + \rho^\prime \mathbf{g}
                                           + \nabla \cdot \tau + F_w,

  \frac{\partial (\rho \theta)}{\partial t} &= - m^2 \nabla \cdot (\rho \mathbf{u_H} \theta m^{-1})
                                               - m   \frac{\partial (\rho w \theta m^{-1})}{\partial z}
                                               + \nabla \cdot ( \rho \alpha_{T}\ \nabla \theta) + F_{\rho \theta},

  \frac{\partial (\rho C)}{\partial t} &= - m^2 \nabla \cdot (\rho \mathbf{u_H} C m^{-1})
                                          - m   \frac{\partial (\rho w C m^{-1})}{\partial z}
                                          + \nabla \cdot (\rho \alpha_{C}\ \nabla C)

where
:math:`u_H` is the velocity in the horizontal (lateral) only,
:math:`\nabla_H` is the horizontal (lateral) gradient,
:math:`m` is the map factor at the appropriate spatial location and :math:`m^{-1} = 1 / m`.
