# MRF PBL Model – WRF-Style Configuration Tests

This directory contains test cases for the MRF (Medium Range Forecast) PBL model as
implemented following WRF's `module_bl_mrf.F`, covering neutral, unstable, and stable
atmospheric boundary layer conditions.

**WRF reference:**
https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F

The implementation uses:
- **HGAMT / HGAMQ** (WRF lines 872–879) only for PBL height finding — not stored in the
  diffusion coefficients.
- **Moisture diffusivity** with Prq ≈ Prt (WRF lines 968–986) when `mrf_moistvars = true`.
- No entrainment terms (WRF MRF does not include explicit entrainment fluxes).

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

