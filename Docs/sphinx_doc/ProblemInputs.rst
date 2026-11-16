.. _sec:ProblemInputs:

Problem-Specific Inputs
=======================

Inputs with the ``prob.`` prefix are **not** global ERF options.  Each one is
read by the problem setup that the executable was built with, so the same name
means different things in different problems -- ``prob.T_0`` is an initial
temperature in one setup and a reference temperature in another, and
``prob.x_c`` is a bubble centre in one and a scalar-blob centre in another.  A
``prob.`` input that a given problem does not read is silently ignored.

For that reason this page is an index rather than a reference: it lists, for
each problem setup in ``Source/Prob``, the ``prob.`` inputs that setup reads.
The meaning, acceptable values and defaults are in the header named beside it,
where each input is read with its default on the same line.

Problems in ``Exec`` may define further ``prob.`` inputs of their own, in that
problem's ``prob.cpp``; those are not listed here.

.. note::

   This index is derived from the ``ParmParse`` calls in ``Source/Prob``.  A
   setup that reads an input through a helper shared with another setup will
   show that input under both.

Initial state perturbation
--------------------------

+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Problem setup           | Defined in                                       | Inputs read (prefix ``prob.``)                                                     |
+=========================+==================================================+====================================================================================+
| ABL,                    | ``ERF_InitCustomPert_ABL.H``                     | ``A_0``, ``KE_0``, ``KE_decay_height``, ``KE_decay_order``, ``T_0``,               |
| constant_rhotheta_src   |                                                  | ``T_0_Pert_Mag``, ``bubble_dT``, ``bubble_x_c``, ``bubble_x_r``, ``bubble_y_c``,   |
|                         |                                                  | ``bubble_y_r``, ``bubble_z_c``, ``bubble_z_r``, ``pert_deltaT``, ``pert_deltaU``,  |
|                         |                                                  | ``pert_deltaV``, ``pert_periods_T``, ``pert_ref_height``, ``pert_rhotheta``,       |
|                         |                                                  | ``rhoKE_0``, ``rho_0`` (the ``bubble_*`` inputs: see :ref:`sec:ABLWarmBubble`)     |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| AnelasticWallDiffusion  | ``ERF_InitCustomPert_AnelasticWallDiffusion.H``  | ``axis``, ``theta_hi``, ``theta_lo``                                               |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Bomex                   | ``ERF_InitCustomPert_Bomex.H``                   | ``A_0``, ``KE_0``, ``T_0``, ``T_0_Pert_Mag``, ``custom_TKE``, ``pert_ref_height``, |
|                         |                                                  | ``qv_0_Pert_Mag``, ``rho_0``                                                       |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Bubble                  | ``ERF_InitCustomPert_Bubble.H``                  | ``T_0``, ``T_pert``, ``T_pert_is_airtemp``, ``do_moist_bubble``, ``eq_pot_temp``,  |
|                         |                                                  | ``perturb_rho``, ``qt_init``, ``theta_pert``, ``use_empirical_psat``, ``x_c``,     |
|                         |                                                  | ``x_r``, ``y_c``, ``y_r``, ``z_c``, ``z_r``                                        |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| DataAssimilation_ISV    | ``ERF_InitCustomPert_DataAssimilation_ISV.H``    | ``M_inf``, ``R``, ``T_inf``, ``alpha``, ``beta``, ``gamma``, ``p_inf``, ``sigma``, |
|                         |                                                  | ``xc``, ``yc``                                                                     |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| DensityCurrent          | ``ERF_InitCustomPert_DensityCurrent.H``          | ``T_pert``, ``x_c``, ``x_r``, ``z_c``, ``z_r``                                     |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| EBPoiseuille            | ``ERF_InitCustomPert_EBPoiseuille.H``            | ``A_0``, ``B_0``, ``prob_type``, ``rad_0``, ``rho_0``, ``xc_frac``, ``xradius``,   |
|                         |                                                  | ``yc_frac``, ``zc_frac``, ``zradius``                                              |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| FlowInABox              | ``ERF_InitCustomPert_FlowInABox.H``              | ``T_0_Pert_Mag``                                                                   |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| GATE                    | ``ERF_InitCustomPert_GATE.H``                    | ``A_0``, ``KE_0``, ``T_0``, ``T_0_Pert_Mag``, ``custom_TKE``, ``pert_deltaQV``,    |
|                         |                                                  | ``pert_deltaT``, ``pert_periods_QV``, ``pert_periods_T``, ``pert_ref_height``,     |
|                         |                                                  | ``qv_0_Pert_Mag``, ``rho_0``                                                       |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| IsentropicVortex        | ``ERF_InitCustomPert_IsentropicVortex.H``        | ``M_inf``, ``R``, ``T_inf``, ``alpha``, ``beta``, ``gamma``, ``p_inf``, ``sigma``, |
|                         |                                                  | ``xc``, ``yc``                                                                     |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| MovingTerrain           | ``ERF_InitCustomPert_MovingTerrain.H``           | ``Ampl``, ``wavelength``                                                           |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| MultiSpeciesBubble      | ``ERF_InitCustomPert_MultiSpeciesBubble.H``      | ``T_0``, ``T_pert``, ``T_pert_is_airtemp``, ``do_moist_bubble``, ``eq_pot_temp``,  |
|                         |                                                  | ``perturb_rho``, ``qt_init``, ``theta_pert``, ``use_empirical_psat``, ``x_c``,     |
|                         |                                                  | ``x_r``, ``y_c``, ``y_r``, ``z_c``, ``z_r``                                        |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| ParticleTests           | ``ERF_InitCustomPert_ParticleTests.H``           | ``T_pert``, ``x_c``, ``x_r``, ``z_c``, ``z_r``                                     |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| RICO                    | ``ERF_InitCustomPert_RICO.H``                    | ``A_0``, ``KE_0``, ``T_0_Pert_Mag``, ``custom_TKE``, ``pert_deltaQV``,             |
|                         |                                                  | ``pert_deltaT``, ``pert_periods_QV``, ``pert_periods_T``, ``pert_ref_height``,     |
|                         |                                                  | ``qv_0_Pert_Mag``, ``wbar_sub_max``                                                |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3D          | ``ERF_InitCustomPert_SDMCongestus3D.H``          | ``A_0``, ``KE_0``, ``Q_0_Pert_Mag``, ``T_0``, ``T_0_Pert_Mag``, ``pert_deltaQV``,  |
|                         |                                                  | ``pert_deltaT``, ``pert_periods_QV``, ``pert_periods_T``, ``pert_ref_height``,     |
|                         |                                                  | ``rho_0``                                                                          |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3DCold      | ``ERF_InitCustomPert_SDMCongestus3DCold.H``      | ``A_0``, ``KE_0``, ``Q_0_Pert_Mag``, ``T_0``, ``T_0_Pert_Mag``, ``pert_deltaQV``,  |
|                         |                                                  | ``pert_deltaT``, ``pert_periods_QV``, ``pert_periods_T``, ``pert_ref_height``,     |
|                         |                                                  | ``rho_0``                                                                          |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| ScalarAdvDiff           | ``ERF_InitCustomPert_ScalarAdvDiff.H``           | ``A_0``, ``B_0``, ``prob_type``, ``rad_0``, ``rho_0``, ``theta_pert_amp``,         |
|                         |                                                  | ``xc_frac``, ``yc_frac``, ``zc_frac``                                              |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SquallLine              | ``ERF_InitCustomPert_SquallLine.H``              | ``T_tr``, ``eq_pot_temp``, ``height``, ``qt_init``, ``theta_0``, ``theta_c``,      |
|                         |                                                  | ``theta_tr``, ``use_empirical_psat``, ``x_c``, ``x_r``, ``z_c``, ``z_r``, ``z_tr`` |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SuperCell               | ``ERF_InitCustomPert_SuperCell.H``               | ``T_tr``, ``eq_pot_temp``, ``height``, ``qt_init``, ``theta_0``, ``theta_c``,      |
|                         |                                                  | ``theta_tr``, ``use_empirical_psat``, ``x_c``, ``x_r``, ``y_c``, ``y_r``, ``z_c``, |
|                         |                                                  | ``z_r``, ``z_tr``                                                                  |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| TaylorGreenVortex       | ``ERF_InitCustomPert_TaylorGreenVortex.H``       | ``M_0``, ``T_0``, ``V_0``, ``rho_0``                                               |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| TurbulentInflow         | ``ERF_InitCustomPert_TurbulentInflow.H``         | ``A_0``, ``KE_0``, ``KE_decay_height``, ``KE_decay_order``, ``T_0``,               |
|                         |                                                  | ``T_0_Pert_Mag``, ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``,   |
|                         |                                                  | ``W_0_Pert_Mag``, ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``,            |
|                         |                                                  | ``pert_periods_V``, ``pert_ref_height``, ``pert_rhotheta``, ``rho_0``              |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| WPS, Metgrid            | ``ERF_InitCustomPert_KE.H``                      | ``KE_0``, ``KE_decay_height``, ``KE_decay_order``                                  |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+

.. _sec:ABLWarmBubble:

ABL warm bubble
~~~~~~~~~~~~~~~

The ABL setup (and ``constant_rhotheta_src``, which shares its perturbation file, keys and
start-up checks) can start the run with one warm bubble, made as WRF's idealized fire
cases make it (``delt_perturbation``, ``xrad_perturbation``, ``yrad_perturbation``,
``zrad_perturbation`` and ``hght_perturbation`` in WRF's ``module_initialize_fire.F``).
The potential temperature of each cell is raised by

.. math::

   \theta' = \Delta\theta \cos^2\left(\frac{\pi r}{2}\right), \qquad
   r = \sqrt{\left(\frac{x - x_c}{x_r}\right)^2 + \left(\frac{y - y_c}{y_r}\right)^2
            + \left(\frac{z - z_c}{z_r}\right)^2} \le 1,

with :math:`(x, y)` the cell centre and :math:`z` the cell-centre height, terrain included
(WRF uses the geopotential height of the mass point, also terrain included). With EB terrain
the cells are Cartesian and :math:`z` is the cell's own height, so the bubble sits at a fixed
altitude whatever the ground below it; an ``eb2.geometry = plane`` boundary shifts the heights
so that the plane lies at zero. Cells with :math:`r > 1` are not changed, and the bubble is not
wrapped across periodic boundaries (WRF's is not either): a bubble that crosses a periodic
boundary is cut there.

- **Compressible:** the bubble is made at constant pressure, as WRF makes it: :math:`\rho\theta`
  and the mixing ratio :math:`q_v` are kept and the density drops to
  :math:`\rho\theta / (\theta + \theta')`, so the bubble starts buoyant and with no pressure
  pulse. Other moisture species are not rescaled. WRF also re-integrates the geopotential of
  the column above the bubble; on ERF's fixed heights that is a change of about 1 Pa and is
  left out.
- **Fixed-in-time density (anelastic):** the density is held, so the bubble is carried by
  :math:`\rho\theta` alone (:math:`\rho\theta \mathrel{+}= \rho\,\theta'`).

The bubble is added to the random (``T_0_Pert_Mag``) and sinusoidal (``pert_deltaT``)
perturbations when those are set as well. WRF places the centre on the mass point
:math:`i_c = (\mathrm{ide}-\mathrm{ids})/2` (rounded down), whose centre is
:math:`x_c = (i_c - 1/2)\,\Delta x` from the domain's west edge: for an even number of cells,
half a cell short of the domain centre (and likewise in y). Use that value to match WRF exactly.

The following table lists the inputs. A start-up line reports the bubble at each level whose
grids hold the domain's lowest corner (always level 0; a refined patch away from that corner
gets the bubble without the line).

+---------------------+-------------------------------------------------------------+--------------------+---------+
| Input               | Meaning                                                     | Acceptable values  | Default |
+=====================+=============================================================+====================+=========+
| ``prob.bubble_dT``  | potential-temperature excess at the bubble centre [K]       | finite Real        | 0 (off) |
|                     | (WRF ``delt_perturbation``); 0 adds no bubble               |                    |         |
+---------------------+-------------------------------------------------------------+--------------------+---------+
| ``prob.bubble_x_c`` | bubble centre x, y and height z [m] (WRF centres x and y on | finite Real        | 0       |
| ``prob.bubble_y_c`` | the domain; z is ``hght_perturbation``)                     |                    |         |
| ``prob.bubble_z_c`` |                                                             |                    |         |
+---------------------+-------------------------------------------------------------+--------------------+---------+
| ``prob.bubble_x_r`` | bubble radii in x, y and z [m] (WRF ``xrad_perturbation``,  | Real > 0 when      | 0       |
| ``prob.bubble_y_r`` | ``yrad_perturbation``, ``zrad_perturbation``)               | ``bubble_dT`` is   |         |
| ``prob.bubble_z_r`` |                                                             | not 0              |         |
+---------------------+-------------------------------------------------------------+--------------------+---------+

A non-finite value or a radius that is not positive stops the run when ``bubble_dT`` is not 0.
A centre or radius given without ``bubble_dT`` stops the run; with ``bubble_dT = 0`` it gives a
warning that they are not used.

Initial velocity perturbation
-----------------------------

+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| Problem setup         | Defined in                                         | Inputs read (prefix ``prob.``)                                                     |
+=======================+====================================================+====================================================================================+
| ABL                   | ``ERF_InitCustomPertVels_ABL.H``                   | ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``, ``W_0_Pert_Mag``,   |
|                       |                                                    | ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``, ``pert_periods_V``,          |
|                       |                                                    | ``pert_ref_height``                                                                |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| Bomex                 | ``ERF_InitCustomPertVels_Bomex.H``                 | ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``, ``W_0_Pert_Mag``,   |
|                       |                                                    | ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``, ``pert_periods_V``,          |
|                       |                                                    | ``pert_ref_height``, ``qv_0_Pert_Mag``                                             |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| ConstantU             | ``ERF_InitCustomPertVels_ConstantU.H``             | ``U_0``, ``V_0``, ``W_0``                                                          |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| CouettePoiseuille     | ``ERF_InitCustomPertVels_CouettePoiseuille.H``     | ``U_0``, ``V_0``, ``W_0``, ``pert_delta_u``, ``pert_delta_v``,                     |
|                       |                                                    | ``pert_periods_u``, ``pert_periods_v``, ``prob_type``                              |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| DataAssimilation_ISV  | ``ERF_InitCustomPertVels_DataAssimilation_ISV.H``  | ``M_inf``, ``R``, ``T_inf``, ``alpha``, ``beta``, ``gamma``, ``sigma``, ``xc``,    |
|                       |                                                    | ``yc``                                                                             |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| EkmanSpiral           | ``ERF_InitCustomPertVels_EkmanSpiral.H``           | ``rho_0``                                                                          |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| IsentropicVortex      | ``ERF_InitCustomPertVels_IsentropicVortex.H``      | ``M_inf``, ``R``, ``T_inf``, ``alpha``, ``beta``, ``gamma``, ``sigma``, ``xc``,    |
|                       |                                                    | ``yc``                                                                             |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| MovingTerrain         | ``ERF_InitCustomPertVels_MovingTerrain.H``         | ``Ampl``, ``wavelength``                                                           |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| ParticleTests         | ``ERF_InitCustomPertVels_ParticleTests.H``         | ``U_0``                                                                            |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| RICO                  | ``ERF_InitCustomPertVels_RICO.H``                  | ``T_0_Pert_Mag``, ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``,   |
|                       |                                                    | ``W_0_Pert_Mag``, ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``,            |
|                       |                                                    | ``pert_periods_V``, ``pert_ref_height``, ``qv_0_Pert_Mag``                         |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3D        | ``ERF_InitCustomPertVels_SDMCongestus3D.H``        | ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``, ``W_0_Pert_Mag``,   |
|                       |                                                    | ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``, ``pert_periods_V``,          |
|                       |                                                    | ``pert_ref_height``                                                                |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| ScalarAdvDiff         | ``ERF_InitCustomPertVels_ScalarAdvDiff.H``         | ``U_0``, ``V_0``, ``W_0``, ``prob_type``, ``uRef``, ``z0``, ``zRef``               |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| TaylorGreenVortex     | ``ERF_InitCustomPertVels_TaylorGreenVortex.H``     | ``V_0``                                                                            |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| Terrain3DHemisphere   | ``ERF_InitCustomPertVels_Terrain3DHemisphere.H``   | ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0_Pert_Mag``,            |
|                       |                                                    | ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``, ``pert_periods_V``,          |
|                       |                                                    | ``pert_ref_height``                                                                |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| TurbulentInflow       | ``ERF_InitCustomPertVels_TurbulentInflow.H``       | ``A_0``, ``KE_0``, ``KE_decay_height``, ``KE_decay_order``, ``T_0``,               |
|                       |                                                    | ``T_0_Pert_Mag``, ``U_0``, ``U_0_Pert_Mag``, ``V_0``, ``V_0_Pert_Mag``, ``W_0``,   |
|                       |                                                    | ``W_0_Pert_Mag``, ``pert_deltaU``, ``pert_deltaV``, ``pert_periods_U``,            |
|                       |                                                    | ``pert_periods_V``, ``pert_ref_height``, ``pert_rhotheta``, ``rho_0``              |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| UserDefined           | ``ERF_InitCustomPertVels_UserDefined.H``           | ``U_0``, ``V_0``, ``W_0``                                                          |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+
| WitchOfAgnesi         | ``ERF_InitCustomPertVels_WitchOfAgnesi.H``         | ``U_0``, ``V_0``                                                                   |
+-----------------------+----------------------------------------------------+------------------------------------------------------------------------------------+

Heating source
--------------

+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| Problem setup       | Defined in                                          | Inputs read (prefix ``prob.``)                                                     |
+=====================+=====================================================+====================================================================================+
| Bomex               | ``ERF_UpdateRhoThetaSources_Bomex.H``               | ``advection_heating_rate``, ``source_cutoff``, ``source_cutoff_transition``        |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| Constant            | ``ERF_UpdateRhoThetaSources_Constant.H``            | ``advection_heating_rate``                                                         |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| RICO                | ``ERF_UpdateRhoThetaSources_RICO.H``                | ``advection_heating_rate``                                                         |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3D      | ``ERF_UpdateRhoThetaSources_SDMCongestus3D.H``      | ``advection_heating_rate``, ``advection_heating_rate_base``                        |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3DCold  | ``ERF_UpdateRhoThetaSources_SDMCongestus3DCold.H``  | ``advection_heating_rate``, ``source_cutoff``                                      |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+
| SineMassFlux        | ``ERF_UpdateRhoThetaSources_SineMassFlux.H``        | ``advection_heating_rate``, ``source_cutoff``, ``source_cutoff_transition``        |
+---------------------+-----------------------------------------------------+------------------------------------------------------------------------------------+

Moisture source
---------------

+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Problem setup       | Defined in                                       | Inputs read (prefix ``prob.``)                                                     |
+=====================+==================================================+====================================================================================+
| Bomex               | ``ERF_UpdateRhoQtSources_Bomex.H``               | ``advection_moisture_rate``, ``moisture_source_cutoff``,                           |
|                     |                                                  | ``moisture_source_cutoff_transition``                                              |
+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| RICO                | ``ERF_UpdateRhoQtSources_RICO.H``                | ``advection_moisture_rate``, ``moisture_source_cutoff``                            |
+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3D      | ``ERF_UpdateRhoQtSources_SDMCongestus3D.H``      | ``advection_moisture_rate``, ``advection_moisture_rate_base``                      |
+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SDMCongestus3DCold  | ``ERF_UpdateRhoQtSources_SDMCongestus3DCold.H``  | ``advection_moisture_rate``, ``source_cutoff``                                     |
+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| SineMassFlux        | ``ERF_UpdateRhoQtSources_SineMassFlux.H``        | ``advection_moisture_rate``, ``moisture_source_cutoff``,                           |
|                     |                                                  | ``moisture_source_cutoff_transition``                                              |
+---------------------+--------------------------------------------------+------------------------------------------------------------------------------------+

Subsidence
----------

+----------+------------------------------------+------------------------------------------------------------------------------------+
| Problem  | Defined in                         | Inputs read (prefix ``prob.``)                                                     |
| setup    |                                    |                                                                                    |
+==========+====================================+====================================================================================+
| Bomex    | ``ERF_UpdateWSubsidence_Bomex.H``  | ``wbar_cutoff_max``, ``wbar_cutoff_min``, ``wbar_sub_max``                         |
+----------+------------------------------------+------------------------------------------------------------------------------------+
| RICO     | ``ERF_UpdateWSubsidence_RICO.H``   | ``wbar_cutoff_max``                                                                |
+----------+------------------------------------+------------------------------------------------------------------------------------+

Other
-----

+----------------------+--------------------------------+------------------------------------------------------------------------------------+
| Problem setup        | Defined in                     | Inputs read (prefix ``prob.``)                                                     |
+======================+================================+====================================================================================+
| CloudChamber         | ``ERF_CloudChamber.H``         | ``initial_relative_humidity``, ``initial_temperature_bottom``,                     |
|                      |                                | ``initial_temperature_top``, ``p_inf``, ``perturbation_mode``,                     |
|                      |                                | ``physical_temperature``, ``qv_bottom``, ``qv_top``, ``relative_humidity``,        |
|                      |                                | ``rh``, ``roughness``, ``temperature_perturbation_amplitude``,                     |
|                      |                                | ``thermodynamic_initialization``, ``theta_bottom``,                                |
|                      |                                | ``theta_perturbation_amplitude``, ``theta_top``                                    |
+----------------------+--------------------------------+------------------------------------------------------------------------------------+
| InitDensityHSE       | ``ERF_InitDensityHSE.H``       | ``T_from_theta_in_moist_init``, ``T_tr``, ``dtheta_dz``, ``eq_pot_temp``,          |
|                      |                                | ``height``, ``qt_init``, ``theta_0``, ``theta_tr``, ``use_empirical_psat``,        |
|                      |                                | ``z_tr``                                                                           |
+----------------------+--------------------------------+------------------------------------------------------------------------------------+
| InitRayleighDamping  | ``ERF_InitRayleighDamping.H``  | ``rayleigh_T_0``, ``rayleigh_U_0``, ``rayleigh_V_0``, ``rayleigh_W_0``             |
+----------------------+--------------------------------+------------------------------------------------------------------------------------+
