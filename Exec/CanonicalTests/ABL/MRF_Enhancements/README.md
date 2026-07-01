# MRF PBL Model – WRF-Style Configuration Tests

This directory contains test cases for the MRF (Medium Range Forecast) PBL model as
implemented following WRF's `module_bl_mrf.F`, covering neutral, unstable, and stable
atmospheric boundary layer conditions.

**WRF reference:**
https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F

## Implementation Notes

Recent updates to ERF's MRF implementation align key physics with WRF's `module_bl_mrf.F`:

| Feature | Change | Reference |
|---------|--------|-----------|
| **Stability function φ_m** | Businger-Dyer: `(1-16·HOL)^(-1/4)` | WRF L857 |
| **SFCFLG stable-side gating** | Disable countergradient when stable (obuk > 0) | WRF L808, 867-884 |
| **WSCALE convective velocity** | Bounds: `u*/5 ≤ wstar ≤ 16·u*` | WRF L863-865 |
| **Terrain-aware coordinates** | Cell-centered heights with `Compute_Zrel_AtCellCenter()` | ✅ Verified |
| **YSU free-atmosphere mixing** | Richardson-number-dependent above PBL | ✅ Already correct |

The implementation features:
- **HGAMT / HGAMQ** (WRF lines 872–879): Countergradient corrections for PBL height finding
- **Moisture diffusivity** with Prq ≈ Prt (WRF lines 968–986) when `mrf_moistvars = true`
- **No explicit entrainment** (consistent with WRF MRF formulation)

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
**Reference:** BOMEX (Barbados Oceanographic and Meteorological Experiment)
Siebesma et al. (2003): A large eddy simulation intercomparison study of shallow cumulus 
convection. J. Atmos. Sci., 60, 1201-1219.

Strong surface heating with buoyancy-driven convection. Tests MRF's critical countergradient
flux (HGAMT/HGAMQ) corrections for accurate PBL height prediction.

```bash
./erf inputs_unstable
```

**Key characteristics:**
- Strong surface heating (5 K/h)
- Superadiabatic temperature near surface
- Rapid PBL growth via convection
- Countergradient corrections ESSENTIAL (HGAMT positive, VPERT limiting critical)
- Moisture diffusivity reflects latent heat effects
- Shallow cumulus cloud formation at PBL top

---

### 3. Stable Boundary Layer (`inputs_stable`)
**Reference:** GABLS1 (GEWEX Atmospheric Boundary Layer Study 1)
Beare et al. (2006): An intercomparison of large-eddy simulations of the stable boundary layer.
Boundary-Layer Meteorol., 118, 247-272.

Strong surface cooling with stable stratification. Tests MRF in weak mixing regime where
Richardson number effects dominate.

```bash
./erf inputs_stable
```

**Key characteristics:**
- Strong surface cooling (-0.25 K/h)
- Statically stable stratification (dθ/dz = 0.01 K/m)
- Very shallow mixed layer (h ~ 50-100 m)
- Inertial oscillations from Coriolis force
- Countergradient corrections suppressed (downgradient fluxes in stable)
- Free atmosphere mixing via YSU scheme critical
- Richardson number dependent diffusivity essential

---

### 4. Cloud-Topped Boundary Layer (`inputs_cloud_topped`)
**Reference:** BOMEX (Barbados Oceanographic and Meteorological Experiment)
Siebesma et al. (2003): A large eddy simulation intercomparison study of shallow cumulus 
convection. J. Atmos. Sci., 60, 1201-1219.

Realistic warm-season conditions with shallow cumulus clouds forming at PBL top. Tests MRF's 
ability to represent cloud-topped boundary layers with realistic entrainment and cloud-top 
evaporative cooling effects.

```bash
./erf inputs_cloud_topped
```

**Key characteristics:**
- Moderate surface heating (8 K/h) and moisture supply
- Moist well-mixed layer with shallow cumulus formation
- Cloud formation at PBL top via vertical convection
- Sharp entrainment zone between PBL and free troposphere
- Strong capping inversion limiting vertical extent
- Cloud-top evaporative cooling and entrainment effects
- Countergradient corrections ESSENTIAL (strongest HGAMT/HGAMQ in this case)
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
**HIGH PRIORITY: Tests smooth HGAMT/HGAMQ activation**

A 24-hour simulation starting from neutral conditions transitioning to weak convection as surface heating gradually increases (0 → 2 K/h). This test validates smooth activation of countergradient corrections without spurious oscillations or negative HGAMT/HGAMQ values.

**Key Physics Tested:**
- Smooth activation of HGAMT from 0 to positive values (no oscillations)
- Monotonic PBL height growth without abrupt jumps
- Stable Richardson number evolution during stability transitions
- Realistic Monin-Obukhov length calculations across neutral-to-unstable transition

```bash
./erf inputs_weak_convection_transition
```

**Expected Behavior:**
- PBL height grows monotonically from ~100 m to ~300-500 m over 24 hours
- HGAMT increases smoothly with heating rate (no negative values)
- No numerical oscillations in diffusivity coefficients

---

### 13. Very Low Wind Speed / Calm Conditions (`inputs_calm_conditions`)
**HIGH PRIORITY: Tests numerical stability in division-by-zero prone calculations**

Tropical calm-air morning scenario with geostrophic wind < 0.5 m/s and convective heating (5 K/h). Tests wind shear safety threshold (≥ 1.0e-8) and proper fallback to minimum diffusivity bounds (Kmin = 0.1 m²/s).

**Key Physics Tested:**
- Numerical stability with wind shear → 0 (tests division-by-zero guards)
- Proper minimum diffusivity bounds (Kmin) when shear is weak
- No NaN/Inf in K_m, K_h, K_q coefficients
- Richardson number calculations remain stable near zero wind

```bash
./erf inputs_calm_conditions
```

**Expected Behavior:**
- No NaN or Infinity values in output
- Diffusivity coefficients bounded by [Kmin, Kmax] = [0.1, 300] m²/s
- Model runs to completion without crashing
- Shear-independent convective mixing occurs

---

### 14. Saturated Boundary Layer (`inputs_saturated_layer`)
**HIGH PRIORITY: Validates saturation-aware HGAMQ limiter (ERF UNIQUE ENHANCEMENT)**

Fog or marine stratus layer with high relative humidity (85% → 99% near surface) and 0.1 g/kg cloud water. Validates the ERF-unique saturation-aware HGAMQ limiter which smoothly ramps down HGAMQ to zero as RH > 95% to prevent moisture pumping and grid-point storms.

**Key Physics Tested:**
- HGAMQ smoothly ramps down as RH approaches 95-100% (ERF unique safeguard!)
- No moisture pumping instabilities in saturated conditions
- Proper interaction between latent heat effects and countergradient mixing
- Realistic PBL height in highly moist conditions

```bash
./erf inputs_saturated_layer
```

**Expected Behavior:**
- HGAMQ → 0 as surface RH > 95% (smooth ramp, not abrupt)
- Stable moisture profiles with no grid-point oscillations
- Cloud-top mixing remains stable despite high RH
- PBL height remains physically reasonable despite saturation

---

### 15. Rapid Surface Cooling / Inversion Collapse (`inputs_rapid_cooling`)
**MEDIUM PRIORITY: Tests robustness during rapid stability transitions**

Rapid surface temperature drop simulating inversion formation (initial +5 K/h heating → rapid -3 K/h cooling → -0.1 K/h stable recovery). Tests numerical stability and smooth evolution of VPERT and HOL during sign changes in heat flux.

**Key Physics Tested:**
- Smooth VPERT evolution from positive to zero as HGAMT deactivates
- Stable HOL calculations when switching from unstable to stable regimes
- No oscillations in Ri_g during rapid Monin-Obukhov sign changes
- Physical PBL height collapse without numerical shocks

```bash
./erf inputs_rapid_cooling
```

**Expected Behavior:**
- Model remains numerically stable throughout 6-hour transition
- PBL height decreases smoothly from ~500 m to ~50 m
- HGAMT → 0 smoothly as heat flux becomes negative
- No spurious diffusivity oscillations during transition

---

### 16. Extreme Heat Flux (`inputs_extreme_heating`)
**MEDIUM PRIORITY: Tests upper bounds on HGAMT (GAMCRT = 3 K limit)**

Extreme surface heating (15-20 K/h) representing intense heat wave or desert thermals. Validates that HGAMT properly limits to GAMCRT = 3 K maximum and VPERT limits to [0, GAMCRT] to prevent unphysical PBL overshoot.

**Key Physics Tested:**
- HGAMT limited to GAMCRT = 3 K (physical upper bound)
- VPERT bounded to [0, 3 K] range
- PBL height grows realistically despite extreme forcing (doesn't overshoot domain)
- Model remains stable under extreme atmospheric conditions

```bash
./erf inputs_extreme_heating
```

**Expected Behavior:**
- HGAMT saturates at 3 K (no values > 3 K)
- VPERT saturates at 3 K despite extreme heat flux
- PBL height growth slows as it approaches PBL height cap (0.9 × z_max)
- Conservative diffusivity bounds (Kmax = 300 m²/s) prevent excessive mixing

---

### 17. Complex Multi-Scale Terrain (`inputs_complex_terrain_mrf`)
**MEDIUM PRIORITY: Extends terrain validation to realistic multi-scale topography**

Extends the simple 2D Witch-of-Agnesi hill test to more realistic multi-scale terrain with multiple peaks/valleys and varying slope angles. Validates that MRF diffusivity and vertical derivatives remain stable and consistent across terrain transitions.

**Key Physics Tested:**
- Consistent K_m, K_h, K_q across terrain transitions (no artifacts at slopes)
- Stable h_zeta and vertical derivative metrics over complex terrain
- No spurious shear generation from terrain coordinate distortion
- Physical PBL height evolution independent of local topography

```bash
./erf inputs_complex_terrain_mrf
```

**Expected Behavior:**
- Diffusivity profiles smooth across terrain transitions
- No localized shear artifacts near steep slopes
- Richardson number calculations remain stable over complex terrain
- PBL height diagnosis unaffected by underlying terrain complexity

---

### 18. Fine Temporal Resolution (`inputs_fine_dt_stable`)
**LOW PRIORITY: Tests numerical stability under stringent CFL conditions**

Very short time step (0.01 s vs standard 0.05-0.1 s) simulation over 3 hours in stable conditions. Tests accumulation of temporal discretization errors and validates conservation laws at high CFL stringency.

**Key Physics Tested:**
- No accumulation of temporal truncation errors over 3 hours
- Numerical stability maintained with very small dt
- Energy conservation validated at fine temporal resolution
- Model performance characteristics under extreme CFL constraints

```bash
./erf inputs_fine_dt_stable
```

**Expected Behavior:**
- Model runs to completion without divergence
- Energy and momentum remain well-conserved
- Diffusivity coefficients remain stable and bounded
- Results consistent with standard dt simulations (when scaled appropriately)

---

## Model Safeguards and Hardening

The ERF MRF implementation includes critical model safeguards to address numerical and physical vulnerabilities beyond standard WRF MRF configurations:

1. **Predictor Monin-Obukhov Length / $HOL$ Bounding Safeguard:**
   Limits the Monin-Obukhov ratio $HOL = \text{sf} \times h / L$ in the initial predictor calculation of `phiM` to prevent floating-point overflow and Division-by-Zero under highly convective, low-wind conditions when $L \to 0$.
2. **Absolute Bounds on Diagnosed PBL Height ($h$):**
   Clamps the diagnosed PBL height to a maximum fraction of the domain height ($0.9 \times z_{max}$) and bounds it from below by the first vertical grid cell center (minimum 10 m). This prevents unphysical runaway PBL heights and division-by-zero in cells close to the surface.
3. **Saturation-Aware Moisture Countergradient ($HGAMQ$) Limiter:**
   Smoothly ramps down the nonlocal moisture countergradient contribution to zero as local relative humidity exceeds 95%. This prevents unphysical moisture pumping and runaway grid-point storm instabilities in saturated conditions.
4. **Upper Bound on Gradient Richardson Number ($Ri_g$) above PBL:**
   Enforces both upper (100.0) and lower (-100.0) bounds on the diagnosed Richardson number in the free atmosphere to prevent extreme scale heights from inducing numerical shocks in high-stability conditions.

## Available MRF Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `pbl.pbl_mrf_Ribcr` | 0.5 | Critical bulk Richardson number (WRF BRCR) |
| `pbl.pbl_mrf_const_b` | 7.8 | Surface layer factor (WRF CFAC) |
| `pbl.pbl_mrf_sf` | 0.1 | Surface layer fraction (WRF SFAC) |
| `pbl.enable_mrf_countergradient` | false | Apply HGAMT/HGAMQ to PBL height finding |
| `pbl.mrf_moistvars` | false | Enable moisture diffusivity (Prq ≈ Prt) |

## Test Case Summary: Stability Regimes

| Case | Stability | Surface Heat Flux | dθ/dz | Turbulence | Key Physics |
|------|-----------|-------------------|-------|-----------|-------------|
| Neutral | Neutral | 0 K/h | 0 K/m | Shear-driven | u*/f scaling, no buoyancy |
| Unstable | Unstable | +5 K/h | -0.02 K/m | Convective | HGAMT/HGAMQ critical, VPERT limits |
| Cloud-Topped | Very Unstable | +8 K/h | -0.02 K/m | Convective with clouds | Moisture effects, entrainment, HGAMT/HGAMQ strongest |
| Stable | Stable | -0.25 K/h | +0.01 K/m | Suppressed | Richardson number, inertial oscillations |

## How to Use These Test Cases

### Validation of MRF Implementation
```bash
# Test neutral PBL (baseline physics)
./erf inputs_neutral

# Test unstable PBL (countergradient corrections)
./erf inputs_unstable

# Test cloud-topped PBL (realistic warm-season conditions)
./erf inputs_cloud_topped

# Test stable PBL (Richardson number effects)
./erf inputs_stable
```

### Regression Testing
These cases should be used to validate that MRF changes don't break:
1. **inputs_neutral**: Shear-driven mixing and u*/f PBL height scaling
2. **inputs_unstable**: Convective PBL growth with HGAMT/HGAMQ corrections
3. **inputs_cloud_topped**: Cloud formation, entrainment, and moisture effects on PBL height
4. **inputs_stable**: Weak mixing and inertial oscillations

### Performance Analysis
Compare against reference metrics:
- **Neutral**: PBL height h ≈ 0.16 * u* / f (Garratt 1994)
- **Unstable**: Rapid growth to 500-1000 m within 3-6 hours
- **Cloud-Topped**: PBL height 800-1200 m with cloud top at 950-1050 m, cloud fraction 0.3-0.6
- **Stable**: Shallow layer h ≈ 50-100 m with inertial oscillations

## References

### MRF Core Implementation
- Hong, S.-Y. and H.-L. Pan (1996): Nonlocal boundary layer vertical diffusion in a
  medium-range forecast model. *Mon. Wea. Rev.*, 124, 2322–2339.
  https://doi.org/10.1175/1520-0493(1996)124<2322:NBLVDI>2.0.CO;2

### YSU Free Atmosphere Scheme (used above PBL in MRF)
- Hong, S.-Y., Y. Noh, and J. Dudhia (2006): A new vertical diffusion package with 
  an explicit treatment of entrainment processes. *Mon. Wea. Rev.*, 134, 2318–2341.
  https://doi.org/10.1175/MWR3250.1

### Test Case References

#### Neutral Boundary Layer
- Sorbjan, Z. (1989): Structure of the Atmospheric Boundary Layer.
  Kluwer Academic Publishers.
- Garratt, J. R. (1994): The Atmospheric Boundary Layer.
  *Cambridge University Press*, 2nd ed.

#### Unstable Boundary Layer
- Siebesma, A. P., et al. (2003): A large eddy simulation intercomparison study 
  of shallow cumulus convection. *J. Atmos. Sci.*, 60, 1201-1219.
  https://doi.org/10.1175/1520-0469(2003)60<1201:ALESIUS>2.0.CO;2

#### Stable Boundary Layer
- Beare, R. J., et al. (2006): An intercomparison of large-eddy simulations of the 
  stable boundary layer. *Boundary-Layer Meteorol.*, 118, 247-272.
  https://doi.org/10.1007/s10546-005-9009-5

