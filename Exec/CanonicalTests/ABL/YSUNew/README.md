# YSUNew PBL Model – Test Cases

This directory contains test cases for the YSUNew (Yang-Seng-Ueno New) PBL model as
implemented in ERF, covering neutral, unstable, stable, and moist atmospheric boundary layer conditions.

## Overview

The YSUNew scheme is an enhanced version of the YSU scheme with improved treatment of:
- Countergradient fluxes (HGAMT/HGAMQ) for unstable conditions
- Entrainment effects at PBL top
- Moisture effects on turbulent mixing
- Stability-dependent mixing parameters

## Test Cases

### 1. Stable Boundary Layer (`inputs_stable`)

**Physical scenario:** Stable conditions with strong surface cooling and statically stable stratification.

**Key parameters:**
- `erf.pbl_type = YSUNew`
- `erf.pbl_ysu_land_Ribcr = 0.25`
- `erf.enable_ysu_countergradient = false` (not applicable in stable conditions)
- `erf.enable_ysu_entrainment = false` (weak entrainment in stable layer)
- `erf.most.surf_temp = 265.0 K` (cold surface)
- `erf.most.surf_heating_rate = -0.25 K/h` (cooling)
- Domain: 4000×4000×2000 m, 4×4×64 cells
- Time steps: 100 (total ~5 seconds simulation time)

**Expected behavior:**
- PBLH (turbulence.Turb_lengthscale) should grow to 50–300 m after 100 steps
- Weak turbulent kinetic energy throughout domain
- Inertial oscillations visible in wind profiles
- K_turb > 0 inside PBL

**How to run:**
```bash
./erf inputs_stable
```

**Expected output files:**
- `plt_stable/` (plotfiles)
- `chk_stable*` (checkpoint files)
- `stable_hist.dat`, `stable_profiles.dat`, `stable_covar.dat` (diagnostics)

---

### 2. Unstable/Convective Boundary Layer (`inputs_unstable`)

**Physical scenario:** Unstable conditions with strong surface heating, buoyancy-driven convection, and rapid PBL growth.

**Key parameters:**
- `erf.pbl_type = YSUNew`
- `erf.pbl_ysu_land_Ribcr = 0.25`
- `erf.enable_ysu_countergradient = true` (critical for accurate PBL height)
- `erf.enable_ysu_entrainment = true` (important for convective PBL top)
- `erf.pbl_ysunew_phifac = 8.0` (stability function parameter)
- `erf.most.surf_temp = 308.0 K` (warm surface)
- `erf.most.surf_heating_rate = 0.1 K*m/s` (strong heating, positive = unstable)
- `erf.ysu_moistvars = true` (moisture effects)
- Domain: 4000×4000×2000 m, 4×4×64 cells
- Time steps: 100

**Expected behavior:**
- PBLH should grow to 500–2000 m after 100 steps (depending on resolution/time stepping)
- Strong convective updrafts visible in vertical velocity field
- Countergradient heat flux active in lower half of PBL
- Virtual potential temperature well-mixed inside PBL
- K_turb significantly elevated compared to stable case

**How to run:**
```bash
./erf inputs_unstable
```

**Expected output files:**
- `plt_unstable/` (plotfiles)
- `chk_unstable*` (checkpoint files)
- `unstable_hist.dat`, `unstable_profiles.dat` (diagnostics)

---

### 3. Neutral Boundary Layer (`inputs_neutral`)

**Physical scenario:** Neutral conditions with zero surface heat flux and uniform potential temperature. Shear-driven turbulence dominates.

**Key parameters:**
- `erf.pbl_type = YSUNew`
- `erf.enable_ysu_countergradient = false` (not relevant for neutral conditions)
- `erf.enable_ysu_entrainment = false`
- `erf.most.surf_temp = 300.0 K`
- `erf.most.surf_heating_rate = 0.0 K/h` (no surface heat flux)
- `erf.ysu_moistvars = false`
- Domain: 4000×4000×2000 m, 4×4×64 cells
- Time steps: 100

**Expected behavior:**
- PBLH should be 100–500 m after 100 steps
- Moderate turbulent kinetic energy with shear-driven structure
- Wind profile logarithmic in surface layer
- Well-defined mixed layer with inversion at top
- K_turb determined primarily by wind shear

**How to run:**
```bash
./erf inputs_neutral
```

**Expected output files:**
- `plt_neutral/` (plotfiles)
- `chk_neutral*` (checkpoint files)
- `neutral_hist.dat`, `neutral_profiles.dat` (diagnostics)

---

### 4. Moist Unstable Boundary Layer (`inputs_moist`)

**Physical scenario:** Moist unstable conditions with high lower-tropospheric moisture content. Latent heating effects amplify buoyancy-driven convection.

**Key parameters:**
- `erf.pbl_type = YSUNew`
- `erf.enable_ysu_countergradient = true`
- `erf.enable_ysu_entrainment = true`
- `erf.pbl_ysunew_phifac = 8.0`
- `erf.ysu_moistvars = true` (moisture effects on diffusivity)
- `erf.use_moisture = true`
- `erf.moisture_model = Kessler` (cloud/precipitation physics)
- Initial sounding includes moisture: q_v = 0.01 kg/kg at surface, tapering to 0.002 at 2000 m
- Domain: 4000×4000×2000 m, 4×4×64 cells
- Time steps: 100

**Expected behavior:**
- PBLH should reach 500–2000 m (similar to unstable case, but with enhanced latent heating)
- Potential for shallow cumulus cloud formation near PBL top
- Virtual potential temperature (θ_v) well-mixed, typically 1–2 K higher than dry case
- K_turb elevated due to both sensible and latent heating
- Moisture flux significant in lower half of PBL

**How to run:**
```bash
./erf inputs_moist
```

**Expected output files:**
- `plt_moist/` (plotfiles)
- `chk_moist*` (checkpoint files)
- `moist_hist.dat`, `moist_profiles.dat` (diagnostics)

---

## Physics Validation

### Expected PBLH Ranges (after 100 steps at dt=0.05 s):

| Scenario       | PBLH Range | Notes                                 |
|----------------|-----------|---------------------------------------|
| Stable         | 50–300 m  | Weak mixing, Richardson number limited|
| Unstable       | 500–2000 m| Strong convection, rapid growth      |
| Neutral        | 100–500 m | Shear-driven, well-defined inversion |
| Moist Unstable | 500–2000 m| Latent heating enhances convection   |

### Turbulent Kinetic Energy (K_turb):

- **K_turb > 0** everywhere inside the PBL for all test cases
- **Stable case:** K_turb ~ 0.001–0.01 m²/s² (weak mixing)
- **Neutral case:** K_turb ~ 0.01–0.1 m²/s² (moderate shear-driven mixing)
- **Unstable cases:** K_turb ~ 0.1–1.0 m²/s² (strong buoyancy-driven mixing)

### Eddy Diffusivity (K_eddy):

- **K_eddy > 0** everywhere inside the PBL
- Must not produce NaN or Inf values in output
- Should vary smoothly with height and time
- Typical ranges:
  - Stable: 0.1–1 m²/s
  - Neutral: 1–10 m²/s
  - Unstable: 10–100 m²/s

---

## Diagnostic Output

Each test case produces:

1. **Plotfiles** (`plt_*/`): Full 3D fields at diagnostic intervals
   - Key fields: `density`, `temp`, `theta`, velocity components, `EddyDiff_Mom_v`

2. **Checkpoint files** (`chk_*/`): Full solution state for restart capability

3. **History/Profile files** (`.dat`): Planar-averaged quantities
   - Profiles of θ, q_v, wind components, TKE
   - Domain-averaged sums of mass, energy, momentum

---

## Running Tests Individually

For example, to run just the unstable case and examine diagnostics:

```bash
./erf inputs_unstable

# View PBLH evolution
tail -100 unstable_profiles.dat | column -t

# Check for NaN values in plotfiles (requires amrex_plot_tool or similar)
# amrex_plot_tool plt_unstable00015 | grep -i nan
```

---

## Automated Testing

A shell script `check_results.sh` is provided to validate all four test cases:

```bash
./check_results.sh
```

This script:
- Runs each test case to max_step = 100
- Checks for NaN values in output
- Validates PBLH is physically reasonable (10–5000 m range)
- Reports PASS/FAIL for each case

---

## Reference Literature

**YSU Scheme:**
- Hong, S.-Y., Y. Noh, and J. Dudhia, 2006: A new vertical diffusion package with an explicit treatment of entrainment processes. *Mon. Wea. Rev.*, 134, 2318–2341.

**Stable Boundary Layer:**
- Beare, R. J., et al., 2006: An intercomparison of large-eddy simulations of the stable boundary layer. *Boundary-Layer Meteorol.*, 118, 247–272.

**Convective Boundary Layer:**
- Siebesma, A. P., et al., 2003: A large eddy simulation intercomparison study of shallow cumulus convection. *J. Atmos. Sci.*, 60, 1201–1219.

---

## Troubleshooting

**Issue:** Model aborts with NaN in first few steps
- **Solution:** Check input sounding file format (must have pressure in first line, then (z, θ, dθ/dz, u, v) for each level)

**Issue:** PBLH unreasonably high (> 3000 m in stable case)
- **Solution:** Verify `pbl_ysu_land_Ribcr` and `enable_ysu_countergradient` settings

**Issue:** Plotfiles missing variables
- **Solution:** Check `erf.plot_vars_1` parameter includes desired diagnostics

---
