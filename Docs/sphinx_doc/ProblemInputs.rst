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
| ABL                     | ``ERF_InitCustomPert_ABL.H``                     | ``A_0``, ``KE_0``, ``KE_decay_height``, ``KE_decay_order``, ``T_0``,               |
|                         |                                                  | ``T_0_Pert_Mag``, ``pert_deltaU``, ``pert_deltaV``, ``pert_ref_height``,           |
|                         |                                                  | ``pert_rhotheta``, ``rhoKE_0``, ``rho_0``                                          |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| AnelasticWallDiffusion  | ``ERF_InitCustomPert_AnelasticWallDiffusion.H``  | ``axis``, ``theta_hi``, ``theta_lo``                                               |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Bomex                   | ``ERF_InitCustomPert_Bomex.H``                   | ``A_0``, ``KE_0``, ``T_0``, ``T_0_Pert_Mag``, ``custom_TKE``, ``pert_ref_height``, |
|                         |                                                  | ``qv_0_Pert_Mag``, ``rho_0``                                                       |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| Bubble                  | ``ERF_InitCustomPert_Bubble.H``                  | ``T_0``, ``T_pert``, ``T_pert_is_airtemp``, ``do_moist_bubble``, ``eq_pot_temp``,  |
|                         |                                                  | ``perturb_rho``, ``qt_init``, ``theta_pert``, ``use_empircal_psat``, ``x_c``,      |
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
| KE                      | ``ERF_InitCustomPert_KE.H``                      | ``KE_0``, ``KE_decay_height``, ``KE_decay_order``                                  |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| MovingTerrain           | ``ERF_InitCustomPert_MovingTerrain.H``           | ``Ampl``, ``wavelength``                                                           |
+-------------------------+--------------------------------------------------+------------------------------------------------------------------------------------+
| MultiSpeciesBubble      | ``ERF_InitCustomPert_MultiSpeciesBubble.H``      | ``T_0``, ``T_pert``, ``T_pert_is_airtemp``, ``do_moist_bubble``, ``eq_pot_temp``,  |
|                         |                                                  | ``perturb_rho``, ``qt_init``, ``theta_pert``, ``use_empircal_psat``, ``x_c``,      |
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
| ScalarAdvDiff           | ``ERF_InitCustomPert_ScalarAdvDiff.H``           | ``A_0``, ``B_0``, ``prob_type``, ``rad_0``, ``rho_0``, ``xc_frac``, ``yc_frac``,   |
|                         |                                                  | ``zc_frac``                                                                        |
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
| CouettePoiseuille     | ``ERF_InitCustomPertVels_CouettePoiseuille.H``     | ``U_0``, ``V_0``, ``W_0``, ``pert_delta_u``, ``pert_hi``, ``pert_lo``,             |
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
