# MRF PBL Model – WRF-Style Configuration Tests

This directory contains test cases for the MRF (Medium Range Forecast) PBL model as
implemented following WRF's `module_bl_mrf.F`.

**WRF reference:**
https://github.com/wrf-model/WRF/blob/master/phys/module_bl_mrf.F

The implementation uses:
- **HGAMT / HGAMQ** (WRF lines 872–879) only for PBL height finding — not stored in the
  diffusion coefficients.
- **Moisture diffusivity** with Prq ≈ Prt (WRF lines 968–986) when `mrf_moistvars = true`.
- No entrainment terms (WRF MRF does not include explicit entrainment fluxes).

## Test Cases

### 1. Baseline MRF (`inputs_baseline`)
Standard MRF without any optional features. Serves as a reference for neutral/stable
boundary layers.

```bash
./erf inputs_baseline
```

### 2. Full WRF-Style (`inputs_full_enhancements`)
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

## References

- Hong, S.-Y. and H.-L. Pan (1996): Nonlocal boundary layer vertical diffusion in a
  medium-range forecast model. *Mon. Wea. Rev.*, 124, 2322–2339.
- Hong, S.-Y. et al. (2006): A new vertical diffusion package with an explicit treatment
  of entrainment processes. *Mon. Wea. Rev.*, 134, 2318–2341.
