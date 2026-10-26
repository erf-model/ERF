# Prognostic Cloud Fraction Test

## Purpose

Run the TwoStream radiation with `erf.radiation.cloud_fraction_prog_enable = true`,
so the cloud fraction of the cloud band is diagnosed from relative humidity (RH)
and cloud water (qc) instead of the fixed `erf.radiation.cloud_fraction`.

The case is a smoke test of that code path. Its sounding is dry (qv = 0) and the
deck sets no moisture model, so RH = 0 and qc = 0 in every cell, the diagnosed
cloud fraction is 0, and the cloud band adds no optical depth. It does not test
the RH value itself; the unit tests in
`Tests/Unit/Radiation/ERF_GTestTwoStreamColumn.cpp` do
(`RelativeHumidityOfAHalfSaturatedParcel`,
`RelativeHumidityStaysBoundedWhenSaturationPressureExceedsP`,
`CloudFractionLiquidWaterTermIsAThreshold`).

## Configuration

- **Domain**: 3000 m × 3000 m × 2710 m, 8 × 8 × 64 cells, stretched in z
  (`erf.initial_dz = 4`, ratio 1.06), periodic in x and y
- **Flow**: ABL with a 15 m/s geostrophic wind, MRF PBL, surface layer at the bottom
- **Time**: `stop_time = 5` s with `erf.fixed_dt = 0.5` s
- **Radiation**: `cloud_layer` optical-depth profile, cloud band 500-2000 m,
  `cloud_tau_per_layer = 0.5`, fixed sun (cos zenith 0.5)
- **Cloud fraction**: `rh_min = 0.7`, `rh_max = 1.0`, `qc_scale = 1e-3` kg/kg

## Key Physics

The cloud fraction of each layer is (see `Docs/sphinx_doc/theory/Radiation.rst`,
"Cloud Fraction Diagnosis"):

```
cf = min(1, clamp((RH - rh_min) / (rh_max - rh_min), 0, 1) + min(1, qc / qc_scale))
```

RH = qv / qsat, with qsat from the Magnus saturation vapour pressure over liquid
water. Inside the cloud band the layer optical depth gains `cf * cloud_tau_per_layer`.

## Files

- `inputs` — control file with the prognostic cloud fraction enabled
- `input_sounding` — dry sounding (qv = 0)
- `check_progcf.py` — reads the radiation diagnostics and the plotfiles

## How to run

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_ProgCloudFraction
mpirun -np 1 erf.ex inputs
python3 check_progcf.py
```

## Expected Results

`check_progcf.py` prints `RESULT: PASS` when the radiation diagnostics parse,
their values are finite, and every plotfile has a Header.

## CTest

None. The case is not registered in CTest.
