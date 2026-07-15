# SineMassFlux

A column with prescribed Eulerian-transported water vapour that condenses and forms cloud
and rain, following the Shipway and Hill (2012) kinematic warm-rain driver.

## References

- Shipway, B.J. and Hill, A.A. (2012), *Diagnosis of systematic differences between
  multiple parametrizations of warm rain microphysics using a kinematic framework*,
  Q.J.R. Meteorol. Soc., **138**, 2196-2211. <https://doi.org/10.1002/qj.1913>
- PySDM (<https://github.com/open-atmos/PySDM>):
  `examples/PySDM_examples/Shipway_and_Hill_2012`

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics, NH42SO4 aerosol, run as a
  one-dimensional (vertical) column; the reference configuration.
- `inputs_Kessler` - same column with the bulk Kessler scheme.
- `inputs_SAM` - convective variant with prescribed surface fluxes and the SAM bulk
  scheme on a 32 x 32 x 100 grid.

## Sounding file

- `input_sounding` - theta/qv/(u, v) profile to 3260 m, used by all of the above. The
  `inputs_SDM` domain top is 2975 m; the profile extends above it and is interpolated
  onto the grid, matching the `SDM_SineMassFlux` regression test.
