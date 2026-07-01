# MRF PBL Model – WRF-Style Configuration Tests

This directory contains test cases for the MRF (Medium Range Forecast) PBL model as
implemented following WRF's `module_bl_mrf.F`, covering neutral, unstable, and stable
atmospheric boundary layer conditions.

WRF reference:
https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F


## Test Cases


### 1. Neutral Boundary Layer (`inputs_neutral`)
**Reference:** Sorbjan (1989): Structure of the Atmospheric Boundary Layer.

Neutral conditions with zero surface heat flux and uniform potential temperature profile.
Shear-driven turbulence dominates. Tests MRF in its simplest stable configuration.

```bash
./erf inputs_neutral
```

**Key characteristics:**
- Zero surface heat flux (adiabatic)
- Uniform potential temperature (no buoyancy)
- Wind shear as only turbulence source
- PBL height determined by u*/f scaling
- Countergradient corrections disabled (not applicable)

---

### 2. Unstable Boundary Layer (`inputs_unstable`)
Reference: BOMEX (Barbados Oceanographic and Meteorological Experiment)
Siebesma et al. (2003): A large eddy simulation intercomparison study of shallow cumulus 
convection. J. Atmos. Sci., 60, 1201-1219.

Strong surface heating with buoyancy-driven convection. Tests MRF's countergradient
flux (HGAMT/HGAMQ) corrections for PBL height prediction.

```bash
./erf inputs_unstable
```

Key characteristics:
- Strong surface heating (5 K/h)
- Superadiabatic temperature near surface
- Rapid PBL growth via convection
- Countergradient corrections important (HGAMT positive, VPERT limiting critical)
- Moisture diffusivity reflects latent heat effects
- Shallow cumulus cloud formation at PBL top

---

### 3. Stable Boundary Layer (`inputs_stable`)
Reference: GABLS1 (GEWEX Atmospheric Boundary Layer Study 1)
Beare et al. (2006): An intercomparison of large-eddy simulations of the stable boundary layer.
Boundary-Layer Meteorol., 118, 247-272.

Strong surface cooling with stable stratification. Tests MRF in weak mixing regime where
Richardson number effects dominate.

```bash
./erf inputs_stable
```

Key characteristics:
- Strong surface cooling (-0.25 K/h)
- Statically stable stratification (dθ/dz = 0.01 K/m)
- Very shallow mixed layer (h ~ 50-100 m)
- Inertial oscillations from Coriolis force
- Countergradient corrections suppressed (downgradient fluxes in stable)
- Free atmosphere mixing via YSU scheme important
- Richardson number dependent diffusivity essential

---

### 4. Cloud-Topped Boundary Layer (`inputs_cloud_topped`)
Reference: BOMEX (Barbados Oceanographic and Meteorological Experiment)
Siebesma et al. (2003): A large eddy simulation intercomparison study of shallow cumulus 
convection. J. Atmos. Sci., 60, 1201-1219.

Realistic warm-season conditions with shallow cumulus clouds forming at PBL top. Tests MRF's 
ability to represent cloud-topped boundary layers with realistic entrainment and cloud-top 
evaporative cooling effects.

```bash
./erf inputs_cloud_topped
```

Key characteristics:
- Moderate surface heating (8 K/h) and moisture supply
- Moist well-mixed layer with shallow cumulus formation
- Cloud formation at PBL top via vertical convection
- Sharp entrainment zone between PBL and free troposphere
- Strong capping inversion limiting vertical extent
- Cloud-top evaporative cooling and entrainment effects
- Countergradient corrections important (strong HGAMT/HGAMQ in this case)
- Moisture effects maximize virtual potential temperature (latent heat dominates)
- Cloud-aware stability adjustments improve mixing representation
- Radiative effects can be enabled for realistic cloud forcing

---

### 5. Baseline MRF (`inputs_baseline`)
Standard MRF without any optional features. Serves as a reference for neutral/stable
boundary layers.

```bash
./erf inputs_baseline
```

---

### 6. Full WRF-Style (`inputs_full_enhancements`)
Enables both WRF-aligned features:
- HGAMT + HGAMQ countergradient correction for PBL height (heat and moisture via q_*)
- Moisture turbulent diffusivity with Prq ≈ Prt

```bash
./erf inputs_full_enhancements
```

Key parameters:
```
pbl.enable_mrf_countergradient = true   # HGAMT/HGAMQ for PBL height finding
pbl.mrf_moistvars               = true   # moisture diffusivity (Prq ~ Prt)
```

---

### 7. Diurnal Cycle Transition (`inputs_diurnal`)
Simulates a dynamic 24-to-48-hour diurnal cycle driven by time-varying surface conditions. This test validates the smooth transition of the Monin-Obukhov length $L$ and the corresponding activation/deactivation of countergradient fluxes ($\text{HGAMT}/\text{HGAMQ}$) without causing spurious numerical oscillations.

```bash
./erf inputs_diurnal
```

---

### 8. Marine Boundary Layer (`inputs_marine`)
Represents marine environments characterized by high latent heat flux (evaporation) but weak sensible heat flux. It specifically validates the land/water mask discrimination safeguard which zeroes out the moisture countergradient term $HGAMQ$ over water bodies while keeping moisture buoyancy corrections active.

```bash
./erf inputs_marine
```

---

### 9. Complex Terrain (`inputs_terrain_mrf`)
A simulation over a simple 2D bell-shaped (Witch of Agnesi) hill utilizing terrain-fitted coordinates (`StaticFittedMesh`). It ensures that vertical coordinate metrics and cell-centered operations function stably under grid distortion without introducing unphysical localized shear or artificial mixing coefficients.

```bash
./erf inputs_terrain_mrf
```

---

### 10. Strong Wind / High Shear (`inputs_high_shear`)
A shear-dominated neutral boundary layer representing extreme storm conditions (geostrophic wind of 35 m/s). It evaluates the performance and safety bounds of the stability-corrected Prandtl number ($0.5 \le Pr_t \le 4.0$) and division-by-zero guards under massive wind shear.

```bash
./erf inputs_high_shear
```

---

### 11. Arctic Polar Stable (`inputs_arctic`)
An extremely stable boundary layer under intense polar radiative cooling (-1 K/h). This test evaluates the Richardson-number-dependent stable free-atmosphere mixing routines and the minimum conservative diffusivity limit ($K_{min} = 0.1\text{ m}^2/\text{s}$) under conditions of completely collapsed boundary layer turbulence.

```bash
./erf inputs_arctic
```

---

## Advanced/Edge Case Tests (Enhanced Coverage)

### 12. Weak Convection Transition (`inputs_weak_convection_transition`)

A 24-hour simulation starting from neutral conditions transitioning to weak convection as surface heating gradually increases (0 → 2 K/h). This test validates smooth activation of countergradient corrections without spurious oscillations or negative HGAMT/HGAMQ values.

Key Physics Tested:
- Smooth activation of HGAMT from 0 to positive values (no oscillations)
- Monotonic PBL height growth without abrupt jumps
- Stable Richardson number evolution during stability transitions
- Realistic Monin-Obukhov length calculations across neutral-to-unstable transition

```bash
./erf inputs_weak_convection_transition
```

Expected Behavior:
- PBL height grows monotonically from ~100 m to ~300-500 m over 24 hours
- HGAMT increases smoothly with heating rate (no negative values)
- No numerical oscillations in diffusivity coefficients

---

### 13. Very Low Wind Speed / Calm Conditions (`inputs_calm_conditions`)

Tropical calm-air morning scenario with geostrophic wind < 0.5 m/s and convective heating (5 K/h). Tests wind shear safety threshold (≥ 1.0e-8) and proper fallback to minimum diffusivity bounds (Kmin = 0.1 m²/s).

Key Physics Tested:
- Numerical stability with wind shear → 0 (tests division-by-zero guards)
- Proper minimum diffusivity bounds (Kmin) when shear is weak
- No NaN/Inf in K_m, K_h, K_q coefficients
- Richardson number calculations remain stable near zero wind

```bash
./erf inputs_calm_conditions
```

Expected Behavior:
- No NaN or Infinity values in output
- Diffusivity coefficients bounded by [Kmin, Kmax] = [0.1, 300] m²/s
- Model runs to completion without crashing
- Shear-independent convective mixing occurs

---

### 14. Saturated Boundary Layer (`inputs_saturated_layer`)

Fog or marine stratus layer with high relative humidity (85% → 99% near surface) and 0.1 g/kg cloud water. Validates the saturation-aware HGAMQ limiter which smoothly ramps down HGAMQ to zero as RH > 95% to prevent moisture pumping and grid-point storms.

Key Physics Tested:
- HGAMQ smoothly ramps down as RH approaches 95-100%
- No moisture pumping instabilities in saturated conditions
- Proper interaction between latent heat effects and countergradient mixing
- Realistic PBL height in highly moist conditions

```bash
./erf inputs_saturated_layer
```

Expected Behavior:
- HGAMQ → 0 as surface RH > 95% (smooth ramp, not abrupt)
- Stable moisture profiles with no grid-point oscillations
- Cloud-top mixing remains stable despite high RH
- PBL height remains physically reasonable despite saturation

---

### 15. Rapid Surface Cooling / Inversion Collapse (`inputs_rapid_cooling`)

Rapid surface temperature drop simulating inversion formation (initial +5 K/h heating → rapid -3 K/h cooling → -0.1 K/h stable recovery). Tests numerical stability and smooth evolution of VPERT and HOL during sign changes in heat flux.

Key Physics Tested:
- Smooth VPERT evolution from positive to zero as HGAMT deactivates
- Stable HOL calculations when switching from unstable to stable regimes
- No oscillations in Ri_g during rapid Monin-Obukhov sign changes
- Physical PBL height collapse without numerical shocks

```bash
./erf inputs_rapid_cooling
```

Expected Behavior:
- Model remains numerically stable throughout 6-hour transition
- PBL height decreases smoothly from ~500 m to ~50 m
- HGAMT → 0 smoothly as heat flux becomes negative
- No spurious diffusivity oscillations during transition

---

### 16. Extreme Heat Flux (`inputs_extreme_heating`)

Extreme surface heating (15-20 K/h) representing intense heat wave or desert thermals. Validates that HGAMT properly limits to GAMCRT = 3 K maximum and VPERT limits to [0, GAMCRT] to prevent unphysical PBL overshoot.

Key Physics Tested:
- HGAMT limited to GAMCRT = 3 K (physical upper bound)
- VPERT bounded to [0, 3 K] range
- PBL height grows realistically despite extreme forcing (doesn't overshoot domain)
- Model remains stable under extreme atmospheric conditions

```bash
./erf inputs_extreme_heating
```

Expected Behavior:
- HGAMT saturates at 3 K (no values > 3 K)
- VPERT saturates at 3 K despite extreme heat flux
- PBL height growth slows as it approaches PBL height cap (0.9 × z_max)
- Conservative diffusivity bounds (Kmax = 300 m²/s) prevent excessive mixing

---


### 17. Fine Temporal Resolution (`inputs_fine_dt_stable`)

Very short time step (0.01 s vs standard 0.05-0.1 s) simulation over 3 hours in stable conditions. Tests accumulation of temporal discretization errors and validates conservation laws at high CFL stringency.

Key Physics Tested:
- No accumulation of temporal truncation errors over 3 hours
- Numerical stability maintained with very small dt
- Energy conservation validated at fine temporal resolution
- Model performance characteristics under extreme CFL constraints

```bash
./erf inputs_fine_dt_stable
```

Expected Behavior:
- Model runs to completion without divergence
- Energy and momentum remain well-conserved
- Diffusivity coefficients remain stable and bounded
- Results consistent with standard dt simulations (when scaled appropriately)

---

