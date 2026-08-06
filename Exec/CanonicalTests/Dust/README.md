# Dust

## Overview
Canonical dust tests covering threshold/emission physics, transport, deposition, diagnostics, schedules, and optional chemistry/suppression workflows.

## Subcases
| Subdirectory | Description | Character |
|---|---|---|
| `DustAtmCoupling` | This case exercises dust evolution under atmospheric forcing with coupling enabled, focusing on whether the atmospheric state drives dust transport consistently. | empirical / regression |
| `DustAtmReturn` | This case documents the return-coupling pathway in which dust can feed back or be reintroduced according to the configured return logic. | empirical / regression |
| `DustBlast` | This case validates time-dependent dust release from a prescribed blast schedule, checking source activation timing and integration with the rest of the dust model. | empirical / regression |
| `DustCriticalMaterials` | This case provides a configuration for dust scenarios that track specialized material properties or composition-sensitive behavior relevant to critical-materials applications. | empirical / regression |
| `DustDeposition` | This case exercises deposition physics that remove dust from the atmosphere and transfer mass toward the surface or deposition budget. | empirical / regression |
| `DustEmission` | This case checks the vertical dust-emission flux calculation once friction velocity exceeds threshold, including soil-property dependence and bin-wise emission behavior. | empirical / regression |
| `DustGrid` | This case verifies that a dust grid with its own refinement ratio is configured correctly relative to the atmospheric mesh. | empirical / regression |
| `DustIntegration` | This case combines several dust features such as schedules, diagnostics, and transport in a single end-to-end workflow representative of application use. | empirical / regression |
| `DustMRFDiffusion` | This case checks that the MRF boundary-layer parameterization mixes dust as expected under turbulent atmospheric conditions. | empirical / regression |
| `DustMSHA` | This case exercises mining-style occupational dust exposure calculations and summary diagnostics relative to MSHA thresholds. | empirical / regression |
| `DustMultiSite` | This case verifies that multiple independent dust-emission sites can be configured and represented within one simulation. | empirical / regression |
| `DustNAAQS` | This case documents ambient-air-quality diagnostics for particulate matter, including comparison against NAAQS-style thresholds. | empirical / regression |
| `DustOutput` | This case focuses on dust plotfile, statistics, and diagnostic output controls to confirm the module writes the intended fields. | empirical / regression |
| `DustPHREEQCFeedback` | This case tests ingestion of PHREEQC-derived geochemical fields and application of those fields to dust behavior or diagnostics. | empirical / regression |
| `DustParticles` | This case verifies particle-size-bin configuration, including transport or diagnostics across multiple dust classes. | empirical / regression |
| `DustRoadSchedule` | This case verifies CSV-driven activation of dust sources associated with road activity or haul-road operation. | empirical / regression |
| `DustScaffold` | This case is a compact baseline for dust-module setup, useful for smoke tests and initial configuration validation before adding specialized features. | empirical / regression |
| `DustSettling` | This case tests gravitational settling of airborne dust following injection, with expected motion and tendencies based on Stokes settling theory for the configured particle sizes. | analytical |
| `DustSuppression` | This case checks that suppression factors or maps reduce emission strength as configured, representing mitigation measures such as water spray or suppressants. | empirical / regression |
| `DustSurfaceReader` | This case validates ingestion of gridded surface-property files such as soil type and silt fraction that influence threshold and emission calculations. | empirical / regression |
| `DustThreshold` | This case verifies the threshold-friction-velocity calculation used to activate dust emission and confirms that emission begins only above threshold. | analytical |
| `DustWindExtract` | This case checks that near-surface winds and related friction-velocity inputs are diagnosed and passed correctly into dust threshold and emission routines. | empirical / regression |
| `PhreeqcReader` | This case focuses on reading and parsing PHREEQC-style external chemistry data needed for optional dust geochemistry coupling. | empirical / regression |

## Notes
- See each subdirectory README for case-specific purpose, expected results, and key parameters.
- See `inputs_dust_master_reference` for module-level reference settings.
