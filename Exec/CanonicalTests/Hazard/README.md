# Hazard

## Overview
Canonical hazard tests covering coupled fire+dust scenarios, terrain effects, smoke transport, health diagnostics, and haboob-style multi-hazard setups.

## Subcases
| Subdirectory | Description | Character |
|---|---|---|
| `DustGaussianHill` | This hazard case uses smooth Gaussian-hill terrain to test terrain-modified winds, plume lifting, and dust concentration response over idealized topography. | analytical |
| `DustGaussianPit` | This case uses a smooth Gaussian depression to examine sheltered flow, recirculation tendencies, or reduced ventilation effects on dust transport. | analytical |
| `DustHealthIntegration` | This hazard case integrates several dust-health diagnostics, such as silica, visibility, or occupational exposure, within a single workflow. | empirical / regression |
| `DustSTEL` | This case validates rolling or short-window occupational exposure metrics analogous to a short-term exposure limit (STEL). | empirical / regression |
| `DustSilica` | This case focuses on silica-related hazard metrics derived from the dust field, supporting occupational or health-risk assessments. | empirical / regression |
| `DustTerrainFlat` | This case provides a terrain-free hazard reference for comparing the effects of hills, pits, and slopes on dust transport and exposure. | empirical / regression |
| `DustTerrainSlopeEffect` | This case examines how a simple slope modifies near-surface winds and dust response compared with the flat-terrain baseline. | empirical / regression |
| `DustVisibility` | This hazard case tests the conversion from suspended dust concentrations to visibility-oriented diagnostics relevant to environmental and safety analysis. | empirical / regression |
| `FireDustBaseline` | This case provides a simple fire+dust hazard reference against which stronger interaction, lofting, or terrain effects can be compared. | empirical / regression |
| `FireDustInteraction1` | This case introduces a basic fire-induced dust interaction mechanism, such as plume lofting or source enhancement, beyond the baseline hazard setup. | empirical / regression |
| `FireDustInteraction2` | This case adds a different or stronger fire-dust coupling mechanism so its effects can be separated from the first interaction scenario. | empirical / regression |
| `FireDustInteraction3` | This case extends the interaction suite with an additional coupling pathway or stronger forcing, broadening the hazard regression matrix. | empirical / regression |
| `FireDustInteractions12` | This case verifies that two interaction mechanisms can be enabled together and produce a coherent coupled hazard response. | empirical / regression |
| `FireDustInteractions123` | This case verifies that two interaction mechanisms can be enabled together and produce a coherent coupled hazard response. | empirical / regression |
| `FireDustLoftingScaling` | This case focuses on the hazard parameters that translate fire-plume buoyancy into dust lofting strength or depth. | empirical / regression |
| `FireDustMassConservation` | This hazard case verifies that fire-driven dust emission and transport remain mass-consistent in a coupled simulation, especially when dry deposition is disabled so airborne mass should not decrease spuriously. | analytical |
| `FireDustTerrainCoupled` | This case combines terrain-modified winds with fire-dust coupling to verify the integrated hazard response in non-flat topography. | empirical / regression |
| `FireDustWindStrength` | This case examines how stronger or weaker background winds alter the combined fire and dust hazard response. | empirical / regression |
| `FireSmokeBaseline` | This case provides a reference smoke-enabled hazard configuration without extra dust interaction complexity. | empirical / regression |
| `FireSmokeDustCoupled` | This case validates simultaneous operation of thermal fire forcing, smoke-tracer transport, and dust evolution in one hazard scenario. | empirical / regression |
| `FireSmokeWindTransport` | This case highlights advection and dispersion of smoke under imposed winds within the hazard framework. | empirical / regression |
| `GaussianTerrain` | This directory collects related canonical ERF test cases under the GaussianTerrain theme. It provides an organizational overview for the child cases listed below. | analytical |
| `HaboobFireFlat` | This case couples an idealized cold-pool / dust-storm initialization with fire on flat ground, verifying the basic haboob-fire workflow. | empirical / regression |
| `HaboobFireHill` | This case adds hill topography to the haboob-fire setup to examine how terrain reshapes the cold-pool flow and resulting fire response. | empirical / regression |
| `HaboobFirePit` | This case uses bowl-shaped terrain to test confined or sheltered haboob-fire interactions and resulting hazard evolution. | empirical / regression |
| `Visualization` | Documentation and supporting assets. | empirical / regression |

## Notes
- See each subdirectory README for case-specific purpose, expected results, and key parameters.
- See `inputs_hazard_master_reference` for module-level reference settings.
