# Radiation Regression and Verification Cases

This directory collects canonical regression, verification, and benchmark cases for ERF's radiation capabilities. The cases are organized to validate both isolated radiative physics and coupled pathways such as thermodynamic forcing, surface-property wiring, boundary-layer interaction, and simplified surface energy balance workflows.

Use these cases as the primary index for manual verification, regression maintenance, and documentation cross-reference when extending the TwoStream implementation in `Source/Radiation/`.

## Test Cases

| Test Case | Description | Prerequisites |
|---|---|---|
| [`SW_ClearSky_Analytical/`](SW_ClearSky_Analytical/) | Clear-sky shortwave Beer-Lambert verification against an analytical attenuation profile. | Radiation diagnostics enabled; single-column-style analytic comparison. |
| [`LW_Isothermal/`](LW_Isothermal/) | Longwave gray-gas consistency test in an isothermal atmosphere with zero net heating. | Radiation diagnostics enabled; isothermal test mode. |
| [`SW_Cloud_Layer/`](SW_Cloud_Layer/) | Shortwave absorption test with a prescribed cloud optical-depth layer. | Cloud optical-depth parameters configured in the inputs file. |
| [`SW_Scattering_Cloud/`](SW_Scattering_Cloud/) | Shortwave scattering validation using the two-stream scattering approximation. | Scattering enabled with cloud optical properties. |
| [`TwoStream_RhoTheta_Coupling/`](TwoStream_RhoTheta_Coupling/) | Verifies deposition of radiative heating into `qheating_rates` and coupling into `RhoTheta`. | TwoStream enabled; thermodynamic source-term coupling active. |
| [`TwoStream_TimeIntegration/`](TwoStream_TimeIntegration/) | Confirms diagnostic cadence and time-integration consistency for repeated TwoStream calls. | TwoStream enabled; multi-step run with diagnostic call-site tracking. |
| [`TwoStream_Benchmark_Suite/`](TwoStream_Benchmark_Suite/) | Aggregates benchmark metrics across representative shortwave, longwave, cloudy, and coupled test cases. | Existing case outputs or manually executed benchmark inputs; Python validation scripts. |
| [`TwoStream_NonuniformDZ/`](TwoStream_NonuniformDZ/) | Validates heating-rate behavior when physical vertical spacing varies with height. | Nonuniform vertical geometry configured in the case inputs. |
| [`TwoStream_SurfaceHeterogeneity/`](TwoStream_SurfaceHeterogeneity/) | Verifies column-varying albedo, emissivity, and surface temperature with robust scalar fallback behavior. | Surface-property fields or scalar fallback parameters configured. |
| [`TwoStream_DynamicTau_MoistCloud/`](TwoStream_DynamicTau_MoistCloud/) | Exercises moisture- and cloud-aware dynamic optical-depth diagnosis. | Moisture, density, and cloud-related radiation inputs enabled. |
| [`TwoStream_PBL_MRF_YSU_Coupling/`](TwoStream_PBL_MRF_YSU_Coupling/) | Validates YSUNew-focused radiative tendency coupling and limiter/smoother safeguards. | TwoStream enabled; YSUNew selected in the test inputs. |
| [`TwoStream_ProgCloudFraction/`](TwoStream_ProgCloudFraction/) | Checks RH/`qc`-based diagnosed cloud fraction with bounded, backward-compatible behavior. | Prognostic cloud-fraction option enabled in `RadChoice`. |
| [`TwoStream_Aerosol_Turbidity/`](TwoStream_Aerosol_Turbidity/) | Verifies prescribed aerosol or turbidity optical-depth contributions and safe fallbacks. | Aerosol/turbidity parameters configured in the inputs file. |
| [`TwoStream_DiurnalSolarGeometry/`](TwoStream_DiurnalSolarGeometry/) | Exercises time-varying solar zenith geometry with fixed-angle fallback retained. | Dynamic solar-geometry controls configured when feature-on case is desired. |
| [`TwoStream_SEB_MultiFabInfra/`](TwoStream_SEB_MultiFabInfra/) | Validates allocation and population of simplified surface energy balance support fields. | SEB infrastructure enabled; surface-property defaults or LSM passthrough available. |
| [`TwoStream_SEB_Diagnostic/`](TwoStream_SEB_Diagnostic/) | Computes and reports diagnostic surface energy balance residuals without prognostic feedback. | SEB infrastructure enabled; diagnostic mode active for feature-on case. |
| [`TwoStream_SEB_Prognostic/`](TwoStream_SEB_Prognostic/) | Advances prognostic surface temperature and moisture with limiter/clamp safeguards. | SEB infrastructure and prognostic mode enabled in the input variant. |

## Shared Resources

- `sounding_us_standard_atm` — Reference atmospheric sounding used by the canonical radiation cases when an analytic or idealized initialization is sufficient.

## Running Tests

Most cases follow the same pattern:

```bash
cd Exec/CanonicalTests/Radiation/<CaseName>
mpirun -np 1 erf.ex <inputs-file>
python3 <check-script>.py
```

Case-specific input files, checker names, and expected artifacts are documented in the corresponding test-case subdirectory when a dedicated README is present.

## Related Documentation

- `RAD_DEVELOPMENT.md` — design history, capability notes, and validation mapping for the radiation implementation.
- `Source/DataStructs/ERF_RadStruct.H` — runtime control structure and parameter definitions for radiation options.

## References

- **Beer-Lambert Law:** Beer, A., 1852. Bestimmung der Absorption des rothen Lichts in farbigen Flüssigkeiten. *Annalen der Physik und Chemie*, 86(78), 78-88.
- **Two-Stream Model:** Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam, 1989. Rapid calculation of radiative heating rates and photodissociation rates in inhomogeneous multiple scattering atmospheres. *J. Geophys. Res.*, 94, 16387–16405. https://doi.org/10.1029/JD094iD13p16387
- **Stefan-Boltzmann Law:** Kirchhoff, G., 1860. Ueber den Zusammenhang zwischen den Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung. *Monatsberichte der Akademie der Wissenschaften zu Berlin*, 783-787.
