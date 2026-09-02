# Flame_Temperature

## Purpose
Covers the two flame-temperature formulations not exercised elsewhere. The
`byram_radiant` default is used throughout the other fire cases, so these two
complete the set.

## Physics / Model Features Exercised
- Flame temperature diagnostic from fireline intensity
- McAlpine heat-balance and Nelson emissivity formulations
- Flame tilt diagnostic (Nelson case)

## Expected Results
Flame temperature should be finite and above ambient wherever the fire is
burning, and equal to ambient in unburned cells. Peak temperature should track
fireline intensity. The three methods differ in magnitude but should agree in
sign and rough range; a large spread between them points at the intensity input
rather than at the method.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.flame_temp_method` | `"mcalpine_heat"` / `"nelson_emissivity"` | Selects the formulation. |
| `erf.fire.compute_flame_tilt` | `true` (Nelson case) | Also computes the flame tilt angle. |

## References
- McAlpine & Wakimoto 1991, The acceleration of fire from point source to equilibrium spread.
- Nelson 2003, Power of the fire - a thermodynamic analysis.
