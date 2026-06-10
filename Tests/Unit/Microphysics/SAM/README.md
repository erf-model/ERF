# SAM Derivation Notes

This note records the symbolic identities that back the SAM source-term and
component-flux tests in this directory.

The expressions below were executed with SymPy during this task against the
algebra used by the helper surface in `Source/Microphysics/SAM/ERF_SAMUtils.H`.
They are included here as executed reductions, not as hypothetical derivations.

## Accretion formulas

Executed SymPy output:

```text
dprc = accrrc*dt*qcc*qpr**powr1
dpsc = accrsc*dt*qcc*qps**pows1
dpgc = accrgc*dt*qcc*qpg**powg1
dpsi = accrsi*dt*qii*qps**pows1
dpgi = accrgi*dt*qii*qpg**powg1
```

## Accretion elasticities

Executed SymPy output:

```text
dprc: donor_ratio=lambda, collector_ratio=lambda**powr1
dpsc: donor_ratio=lambda, collector_ratio=lambda**pows1
dpgc: donor_ratio=lambda, collector_ratio=lambda**powg1
dpsi: donor_ratio=lambda, collector_ratio=lambda**pows1
dpgi: donor_ratio=lambda, collector_ratio=lambda**powg1
```

These are the scaling identities behind the donor-linear and collector power-law
contracts.

## Source partition total-water residual

Partition formulas:

```text
dqpr = (dqca + dqia)*omp + dprc
dqps = (dqca + dqia)*(1 - omp)*(1 - omg) + dpsc + dpsi
dqpg = (dqca + dqia)*(1 - omp)*omg + dpgc + dpgi
```

Executed SymPy reduction:

```text
partition_total_water_residual = 0
```

This is the identity

```text
(dqpr + dqps + dqpg)
- (dqca + dqia + dprc + dpsc + dpgc + dpsi + dpgi)
= 0
```

and, with `dqc = dqca + dprc + dpsc + dpgc` and
`dqi = dqia + dpsi + dpgi`, equivalently

```text
(dqpr + dqps + dqpg) - (dqc + dqi) = 0
```

## Latent-proxy residual for source partitioning

The source-step latent proxy is

```text
tabs + fac_cond*qv - fac_fus*(qci + qps + qpg)
```

For the corrected temperature update

```text
tabs' = tabs + fac_fus*(dqca*(1 - omp) - dqia*omp + dpsc + dpgc)
```

Executed SymPy reduction:

```text
corrected = 0
```

For the older update that omitted `dpsc + dpgc`

```text
tabs_old' = tabs + fac_fus*(dqca*(1 - omp) - dqia*omp)
```

Executed SymPy reduction:

```text
missing_dpsc_dpgc = -dpgc*fac_fus - dpsc*fac_fus
```

This is the missing fusion-heating residual that the riming test protects.

## Evaporation conservation and derivative

Executed SymPy output:

```text
rain_rate = evapr1*sqrt(qpr) + evapr2*qpr**powr2
snow_rate = evaps1*sqrt(qps) + evaps2*qps**pows2
graupel_rate = evapg1*sqrt(qpg) + evapg2*qpg**powg2
evap_total_water_residual = 0
evap_latent_proxy_residual = 0
generic_evap_derivative = a/(2*sqrt(q)) + b*p*q**p/q
```

The derivative line is the simplified form of

```text
d/dq[a*sqrt(q) + b*q**p]
```

which is equivalent to

```text
0.5*a/sqrt(q) + p*b*q**(p - 1)
```

The latent-proxy residual reduces to zero after substituting

```text
fac_sub = fac_cond + fac_fus
```

## Limiter proportionality

For the active branch,

```text
scalec = qcl_lim / (dqc + eps)
scalei = qci_lim / (dqi + eps)
```

Executed SymPy output:

```text
active_liquid_ratio_residual = 0
active_ice_ratio_residual = 0
active_liquid_total_residual = eps*qcl_lim/(dqc + eps)
active_ice_total_residual = eps*qci_lim/(dqi + eps)
```

So every liquid sink shares `scalec`, every ice sink shares `scalei`, pairwise
ratios are preserved, and the gap to the available cloud mass is the expected
epsilon-scale residual.

For the inactive branch,

```text
scalec_inactive = dqc / (dqc + eps)
```

Executed SymPy output:

```text
inactive_liquid_identity_residual = -dprc*eps/(dqc + eps)
```

The same epsilon-scale identity residual applies component-wise to the other
liquid and ice sinks.

## Component-flux derivative

Executed SymPy output:

```text
F = sqrt(rho0)*v*(q*rho)**(c + 1)/sqrt(rho)
dF/dq = q**(c + 1)*rho**(c + 1)*sqrt(rho0)*v*(c + 1)/(q*sqrt(rho))
dF/dq residual = 0
dF/drho = q**(c + 1)*rho**(c + 1)*sqrt(rho0)*v*(2*c + 1)/(2*rho**(3/2))
dF/drho residual = 0
```

The zero residuals certify the compact forms used in the tests:

```text
dF/dq   = (1 + c) * F / q
dF/drho = (c + 1/2) * F / rho
```

## Sedimentation telescoping budget

Executed SymPy reduction:

```text
((f1 - f0) + (f2 - f1) + (f3 - f2)) - (f3 - f0) = 0
```

This is the column-budget identity behind the component-wise `PrecipFall`
budgets. After multiplying each cell contribution by `rho*detJ` and using the
cell update's `dJinv`, the same telescoping structure closes the detJ-weighted
budget.

## Mapping To Tests

| Derivation or contract | Tests |
| --- | --- |
| Exact accretion formulas | `SAMScalar.AccretionRatesMatchIndependentFormulas` |
| Accretion elasticities | `SAMScalar.AccretionRatesHaveExpectedElasticities` |
| Accretion selector thresholds | `SAMScalar.AccretionSelectorBoundaries`, `SAMScalar.AccretionActivationThresholds` |
| Source partition total-water residual | `SAMScalar.PartitionedPrecipConservesSources`, `SAMPhysicalProperties.PrecipSources_AllSpeciesNonzeroExpectedDiagnostics` |
| Corrected latent-proxy residual and missing `dpsc + dpgc` term | `SAMPhysicalProperties.PrecipSources_ConserveLatentProxyOrExposeRimingHeatingBug`, `SAMPhysicalProperties.PrecipSources_AllSpeciesNonzeroExpectedDiagnostics` |
| Evaporation exact formulas | `SAMScalar.PrecipEvaporationRatesMatchIndependentFormulas` |
| Evaporation derivative | `SAMScalar.PrecipEvaporationRateDerivativesMatchAnalyticFormulas` |
| Cloud-sink limiter proportionality and epsilon residual | `SAMScalar.CloudSinkLimiterPreservesProportionsWhenActive`, `SAMScalar.CloudSinkLimiterIdentityWhenInactive`, `SAMScalar.CloudSinkLimiterCapsSinks` |
| Exact all-source active diagnostics | `SAMPhysicalProperties.PrecipSources_AllSpeciesNonzeroExpectedDiagnostics`, `SAMPhysicalProperties.PrecipSources_AllSpeciesNonzeroAllSourceTermsActive` |
| Total-water conservation across representative source cases | `SAMPhysicalProperties.PrecipSources_ConserveTotalWaterAcrossRepresentativeCases` |
| Public constant-pressure latent-proxy closure | `SAMPhysicalProperties.PublicPrecipConstantPressureLatentProxySingleCell`, `SAMParallel.PrecipConstantPressureLatentProxyBudgetDecompositionInvariant`, `SAMParallel.PrecipDetJWeightedConstantPressureLatentProxyBudgetDecompositionInvariant` |
| Component-flux derivative | `SAMScalar.PrecipComponentFluxDerivativeInPositiveDomain` |
| Component-flux tiny-positive activation | `SAMScalar.PrecipComponentFluxTinyPositiveContract`, `SAMScalar.PrecipComponentFlux_ZeroWhenComponentIsZero` |
| Component-wise sedimentation telescoping budget | `SAMPhysicalProperties.PrecipFall_RainSnowGraupelComponentBudgetsClose`, `SAMPhysicalProperties.PrecipFall_DetJWeightedComponentBudgetsClose`, `SAMParallel.PrecipFallComponentBudgetsIndependentOfDecomposition` |

## Scope Notes

This README is intentionally limited to the algebra exercised by the tests in
this directory. It does not attempt to restate the full SAM microphysics model.