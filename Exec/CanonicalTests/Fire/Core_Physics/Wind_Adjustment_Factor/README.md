# Wind_Adjustment_Factor

## Purpose
Compares the two wind adjustment factor formulations that reduce the reference
wind to midflame height before it enters the spread calculation.

## Physics / Model Features Exercised
- Wind adjustment factor applied to the effective fire wind
- Andrews (2012) unsheltered and BehavePlus formulations

## Expected Results
Both cases should give a midflame wind below the reference wind, and a spread
rate below an unadjusted run. The two formulations differ for the same fuel bed
depth, so the pair isolates the effect of the formula itself.

An unrecognised `waf_formula` prints a warning and silently falls back to a WAF
of 0.4, so a result matching neither reference value means the string was not
recognised rather than that the formula behaved oddly.

## Key Parameters
| Parameter | Value | Description |
|-----------|-------|-------------|
| `erf.fire.use_waf` | `true` | Enables the wind adjustment factor. |
| `erf.fire.waf_formula` | `"andrews"` / `"behaviorplus"` | Selects the formulation. |
| `erf.fire.fuel_model_id` | `1` | Fuel bed depth from this model drives the WAF. |

## References
- Andrews 2012, Modeling wind adjustment factor and midflame wind speed for Rothermel's surface fire spread model.
