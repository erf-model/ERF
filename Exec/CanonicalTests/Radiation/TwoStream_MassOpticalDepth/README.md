# Mass-Based Optical Depth Test

## Objective

Validate the `tau_model = mass` option: the shortwave and longwave optical depth of every layer
follows its mass path, so the column optical depth is a property of the atmosphere rather than of
the vertical grid.

## Test Design

The same moist column (qv = 8 g/kg, dry-adiabatic-ish theta profile, no clouds) is run on 32 and 64
layers over 1024 m with

```
erf.radiation.tau_model = "mass"
erf.radiation.sw_kabs_dry   = 4.0e-6   # m^2/kg
erf.radiation.sw_kscat_dry  = 3.0e-6   # m^2/kg, Rayleigh (omega = 1, g = 0)
erf.radiation.sw_kabs_vapor = 4.0e-3   # m^2/kg
erf.radiation.lw_kabs_dry   = 1.0e-4   # m^2/kg
erf.radiation.lw_kabs_vapor = 0.1      # m^2/kg
```

With the per-layer model the 64-layer column would carry twice the optical depth of the 32-layer one;
with the mass model the two columns have the same optical depth and the same fluxes.

## Files

- `inputs_coarse`, `inputs_fine` — 32- and 64-layer runs
- `input_sounding_mass` — moist sounding
- `check_mass_optical_depth.py` — compares the two diagnostics files

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_MassOpticalDepth
mpirun -np 1 erf.ex inputs_coarse
mpirun -np 1 erf.ex inputs_fine
python3 check_mass_optical_depth.py
```

## Validation Criteria

1. `SW_surface`, `SW_up_TOA`, `LW_up_TOA` agree within 1% between the grids
2. `LW_net_surface` agrees within 2 W/m²
3. `SW_up_TOA` is positive (Rayleigh scattering plus surface reflection)
