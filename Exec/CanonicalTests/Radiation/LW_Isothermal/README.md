# Longwave Isothermal Column Test

## Objective

Validate the gray-gas longwave solver on a physically isothermal column: an atmosphere at a
uniform absolute temperature over a black surface at the same temperature.

## Test Design

### Configuration

- **Domain**: 3000 m × 3000 m horizontal, 1024 m vertical (64 layers)
- **Time**: one slow step (0.25 s)
- **Sounding**: `theta(z) = T0 exp(g z / (c_p T0))` with `T0 = 300 K`, so `T = theta * pi = T0`
  at every level (the sounding is piecewise linear between five levels, accurate to ~0.05 K)
- **Surface**: `surface_temp_k = 300`, `surface_emissivity_lw = 1`
- **Longwave optical depth per layer**: `tau_lw = 1.0` (the 64-layer column is opaque)
- **Shortwave**: disabled

### Key Physics

Every layer emits `sigma T0^4`. In an opaque isothermal column the upward flux is `sigma T0^4`
on every interface and the downward flux approaches `sigma T0^4` at the surface, so

```
LW_up_TOA      = sigma T0^4              ≈ 459.3 W/m²
LW_net_surface = sigma T0^4 exp(-tau_col) ≈ 0 W/m²
```

The column still cools to space (the downward flux vanishes at the top), so the heating rate is
non-zero and strongest in the top layers. This case previously used an `isothermal_test` override
that forced `F_up = F_down` and zero heating; it now exercises the real solver.

## Files

- `inputs` — control file
- `input_sounding_lw_isothermal` — isothermal-temperature sounding
- `check_flux_accuracy.py` — validation script
- `radiation_lw_diag.dat` — reference diagnostics

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/LW_Isothermal
mpirun -np 1 erf.ex inputs
python3 check_flux_accuracy.py
```

## Validation Criteria

1. `LW_up_TOA` within 0.5% of `sigma T0^4`
2. `|LW_net_surface|` below 0.5% of `sigma T0^4`
3. `heating_rate_max` finite and non-zero
