
 .. role:: cpp(code)
    :language: c++

 .. role:: f(code)
    :language: fortran

.. _Buoyancy:

Buoyancy
=========

ERF has three options for how to define the buoyancy force.  Even in the absence of moisture these
expressions are not equivalent.

Type 1
------

One version of the buoyancy force is expressed simply as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g}

which follows most directly from the Navier-Stokes equations with no additional approximations.

Type 2
------

One version of the buoyancy force is given in Marat F. Khairoutdinov and David A. Randall's paper (J. Atm Sciences, 607, 1983):

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx -\rho_0 \mathbf{g} ( \frac{T^\prime}{\bar{T}}
                 + 0.61 q_v^\prime - q_c - q_i - q_p - \frac{p^\prime}{\bar{p}} )

This can be derived by starting with

.. math::
   p = \rho (R_d q_d + R_v q_v) T = \rho R_d T (q_d + \frac{R_v}{R_d} q_v) =
        \rho R_d T [1 + (\frac{R_v}{R_d} − 1) q_v − q_c − q_i - q_p ]

where we have replaced :math:`q_d` by :math:`1 - q_v - q_c - q_i - q_p`.

Then, assuming the perturbations of :math:`p^\prime`, :math:`T^\prime`, and :math:`\rho^\prime`
are small compared with the total pressure, temperature, and density, respectively,
and :math:`\rho = \rho_d + \rho_v + \rho_c + \rho_i + \rho_p`
we can write

.. math::
   \frac{p^\prime}{p} = \frac{\rho^\prime}{\rho} + \frac{T^\prime}{T} + \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+(\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p)}

which allows us to write

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx \rho \left( \frac{p^\prime}{p} - \frac{T^\prime}{T} -
         \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+(\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p)} \right)

We can re-write

.. math::
     \frac{(\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime}{1+ ( (\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p) ) } )
     \approx
     (\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime) (1 - ( (\frac{R_v}{R_d}-1)q_v - q_c - q_i - q_p) )

     \approx
     (\frac{R_v}{R_d}-1) q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime)

if we retain only first-order terms in the mixing ratios.

Thus the buoyancy term can be finally written as

.. math::
     \mathbf{B} = \rho^\prime \mathbf{g} \approx \rho \left( \frac{p^\prime}{p} - \frac{T^\prime}{T} -
         ( \frac{R_v}{R_d}-1) q_v^\prime + q_c^\prime + q_i^\prime + q_p^\prime \right)

If we assume that :math:`q_c = q_c^\prime`, :math:`q_i = q_i^\prime`, and :math:`q_p = q_p^\prime`
because :math:`\bar{q_i} = \bar{q_c} = \bar{q_p} = 0`,  and we replace :math:`\frac{R_v}{R_d}` by 1.61,
then this is identical to the expression at the top of this section.

Type 3
------

The third option for buoyancy force is

.. math::
   \mathbf{B} = -\rho_0 \mathbf{g} ( 0.61 q_v^\prime - q_c^\prime - q_i^\prime - q_p^\prime
                  + \frac{T^\prime}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) )


To implement the buoyancy force term, we assume :math:`T_v = T (1 + (\frac{R_v}{R_d} − 1 ) q_v − q_c − q_i - q_p) \approx T`, then we have

.. math::
    p = \rho (R_d q_d + R_v q_v) T = \rho R_d T [1 + (\frac{R_v}{R_d} − 1) q_v − q_c − q_i - q_p ] = \rho R_d T_v


so the perturbation of :math:`\rho` can be written as

.. math::
   \frac{p^\prime}{p} = \frac{\rho^\prime}{\rho} + \frac{T_v^\prime}{T_v}


then :math:`\frac{\rho^\prime}{\rho}` is

.. math::
   \frac{\rho^\prime}{\rho} = \frac{p^\prime}{p} - \frac{T_v^\prime}{T_v}

if we ignore the term :math:`\frac{p^\prime}{p}`, the equation implemented can be written as

.. math::
   \frac{T_v^\prime}{T_v} \approx \frac{\bar{T} [ (\frac{R_v}{R_d}-1) (q_v-\bar{q_v}) - (q_c + q_i + q_p - \bar{q_c} - \bar{q_i} - \bar{q_p})] +
                           (T - \bar{T})[1+(\frac{R_v}{R_d}-1) \bar{q_v} - \bar{q_c} - \bar{q_i} - \bar{q_p} ]}{\bar{T}}


after reorganizing the terms, we get

.. math::
   \mathbf{B} = \mathbf{g} \rho^\prime = -\rho \mathbf{g} [ 0.61 (q_v - \bar{q_v}) - (q_c - \bar{q_c} + q_i - \bar{q_i} + q_p - \bar{q_p})
                  + \frac{T - \bar{T}}{\bar{T}} (1.0 + 0.61 \bar{q_v} - \bar{q_i} - \bar{q_c} - \bar{q_p}) ]

where the overbar represents a horizontal average of the current state.

