# Phase 1 Two-Stream Radiation Module RegTests

This directory contains regression tests (RegTests) for Phase 1 of the two-stream
atmospheric radiation model in ERF.

## Overview

Phase 1 implements a minimal, clear-sky radiation solver with two test cases:

1. **SW_ClearSky_Analytical** — Validates shortwave (solar) radiation using the
   Beer-Lambert direct-beam model against an exact analytical solution.

2. **LW_Isothermal** — Validates longwave (thermal) radiation using a gray-gas
   two-stream model in a special isothermal configuration where net heating
   is zero (perfect self-consistency test).

## Shared Resources

### `sounding_us_standard_atm`

A 1976 U.S. Standard Atmosphere temperature/wind sounding used to initialize
both test cases. Format:
```
z(m)    T(K)      u(m/s) v(m/s) q(kg/kg)
0       288.15    5.0    0.0    0.010
1000    284.9     5.0    0.0    0.008
...
```

## Test Case 1: SW_ClearSky_Analytical

**Physics:** Beer-Lambert direct-beam solar radiation

**Configuration:** `SW_ClearSky_Analytical/inputs`

**Key Parameters:**
- Domain: 1000 × 1000 m horizontal, 10 km vertical (20 layers)
- Simulation: 1 timestep, 0.1 s duration
- Solar constant: S₀ = 1361 W/m²
- Solar zenith angle: 60°
- Optical depth per layer: τ = 0.05

**Analytical Solution:**

The downwelling direct-beam flux at height z is:

```
F_dir(z) = S₀ · cos(zenith) · exp(-τ_cumulative(z) / cos(zenith))
```

where τ_cumulative = τ · k (k = number of layers above z).

**Validation:**

The Python script `check_flux_accuracy.py` verifies:
1. TOA flux matches S₀ · cos(zenith) = 1361 × cos(60°) ≈ 680.5 W/m²
2. Surface flux matches analytical attenuation through 20 layers
3. All fluxes are non-negative
4. Relative accuracy within 5%

**Running:**

```bash
cd Exec/CanonicalTests/Radiation/SW_ClearSky_Analytical
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Test Case 2: LW_Isothermal

**Physics:** Gray-gas two-stream longwave radiation (isothermal mode)

**Configuration:** `LW_Isothermal/inputs`

**Key Parameters:**
- Domain: 1000 × 1000 m horizontal, 10 km vertical (20 layers)
- Simulation: 1 timestep, 0.1 s duration
- Isothermal temperature: T_iso = 288.15 K
- Longwave optical depth per layer: τ_lw = 1.0
- Isothermal test mode: ENABLED (forces all T = T_iso everywhere)

**Analytical Solution:**

In isothermal mode, all layers have the same temperature, so:

- Upwelling flux: F_up(z) = σ · T_iso⁴ everywhere (constant)
- Downwelling flux: F_down(z) = σ · T_iso⁴ everywhere (constant)
- Net flux: F_net = F_up - F_down = 0 everywhere
- Heating rate: dT/dt = 0 everywhere

where σ = 5.670374419 × 10⁻⁸ W/(m²·K⁴) is the Stefan-Boltzmann constant.

For T_iso = 288.15 K:
```
F = 5.670e-8 × (288.15)⁴ ≈ 384.4 W/m²
```

**Validation:**

The Python script `check_flux_accuracy.py` verifies:
1. F_up matches σ · T_iso⁴
2. F_down matches σ · T_iso⁴
3. F_up and F_down are equal (within machine precision)
4. Maximum heating rate is exactly zero (within round-off)
5. Energy conservation is satisfied

This is a stringent consistency test: if the two-stream solver is
correct, isothermal conditions must remain isothermal.

**Running:**

```bash
cd Exec/CanonicalTests/Radiation/LW_Isothermal
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Diagnostic Output

Both tests generate a CSV file with radiation diagnostics:
- `radiation_sw_diag.dat` (SW test)
- `radiation_lw_diag.dat` (LW test)

Columns: `step, time, SW_surface, SW_TOA, F_up_surface, F_down_toa, heating_rate_max`

The ERF solver also prints lines starting with `RADIATION_DIAG:` that the
check scripts parse.

## References

- **Beer-Lambert Law:** Beer, A., 1852. Bestimmung der Absorption des rothen
  Lichts in farbigen Flüssigkeiten. *Annalen der Physik und Chemie*, 86(78), 78-88.

- **Two-Stream Model:** Toon, O. B., C. P. McKay, T. P. Ackerman, and K. Santhanam,
  1989. Rapid calculation of radiative heating rates and photodissociation rates
  in inhomogeneous multiple scattering atmospheres. *J. Geophys. Res.*, 94,
  16387–16405. https://doi.org/10.1029/JD094iD13p16387

- **Stefan-Boltzmann Law:** Kirchhoff, G., 1860. Ueber den Zusammenhang zwischen
  den Emissionsvermögen und den Absorptionsvermögen der Körper für Wärmestrahlung.
  *Monatsberichte der Akademie der Wissenschaften zu Berlin*, 783-787.

## Related Documentation

- **Phase 1 Specification:** See problem_statement.md in the repository root
- **Radiation Module:** `Source/Radiation/` directory
- **Control Parameters:** `Source/DataStructs/ERF_RadStruct.H`
