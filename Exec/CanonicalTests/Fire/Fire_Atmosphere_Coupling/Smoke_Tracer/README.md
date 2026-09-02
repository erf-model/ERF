# Smoke_Tracer

## Purpose
Exercises passive smoke-tracer emission driven by fire sensible heat release.

## Physics / Model Features Exercised
- Smoke emission proportional to fire heat flux
- Transport of the smoke scalar by the dycore
- Lagged fire-atmosphere coupling as the heat source

## Expected Results
Smoke should appear only where the fire is burning and be advected downwind.
Column-integrated smoke grows while fuel is being consumed and levels off once
the fuel behind the front is exhausted. Smoke mass should never be negative.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.smoke_enable` | `true` | Enables the smoke tracer component. |
| `erf.fire.smoke_emission_factor` | `0.015` | Smoke mass emitted per unit fuel burned [kg/kg]. |
| `erf.fire.smoke_heat_of_comb` | `1.86e7` | Heat of combustion used to convert heat release to fuel mass [J/kg]. |
| `erf.fire.coupling_type` | `"lagged"` | Fire heat is injected into the atmosphere. |

## Notes
- Enabling smoke adds a conserved component, so plotfile component counts differ
  from the other coupling cases.

## References
- Urbanski 2014, Wildland fire emissions, carbon, and climate: Emission factors.
