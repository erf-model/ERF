# SineMassFlux

A column with prescribed Eulerian-transported water vapor that condenses and forms cloud
and rain, following the Shipway and Hill (2012) kinematic warm-rain driver.

## References

- Beydoun, H., Ghosh, D., de Jong, E. K., McGuffin, D., Kendrick, C., Lee, J., Gardner, D. J.,
  and Lundquist, K. A. (2026), *Super Droplet Cloud Microphysics in the Energy Research and
  Forecasting Model*, ESS Open Archive (preprint).
  <https://doi.org/10.22541/essoar.176858434.43447149/v1>
- Shipway, B.J. and Hill, A.A. (2012), *Diagnosis of systematic differences between
  multiple parametrizations of warm rain microphysics using a kinematic framework*,
  Q.J.R. Meteorol. Soc., **138**, 2196-2211. <https://doi.org/10.1002/qj.1913>
- PySDM (<https://github.com/open-atmos/PySDM>):
  `examples/PySDM_examples/Shipway_and_Hill_2012`

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics, NH42SO4 aerosol, run as a
  one-dimensional (vertical) column: the case of Figure 1 of Beydoun et al. (2026), rebuilt
  from the paper because its deck is not available, and set up as PySDM's Shipway and Hill
  (2012) example (theta held at the sounding, geometric collision kernel with efficiency 1,
  128 super-droplets per cell, 5 s time step). The vertical grid (dz = 59.5 m) is that of the
  plotted run; the text of the paper gives 50 m and 256 super-droplets per cell. The deck gives
  more rain water and surface precipitation than the ERF run of the paper.
- `inputs_Kessler` - same column with the bulk Kessler scheme.
- `inputs_SAM` - convective variant with prescribed surface fluxes and the SAM bulk
  scheme on a 32 x 32 x 100 grid.

## Sounding file

- `input_sounding_SDM` - theta/qv/(u, v) profile of `inputs_SDM`, inferred from Figure 1 of
  Beydoun et al. (2026): 1015 hPa, 295.825 K and 15 g/kg at the ground, supersaturated at
  about 600-900 m at the start, as in that figure.
- `input_sounding` - theta/qv/(u, v) profile to 3260 m, used by `inputs_Kessler`, `inputs_SAM`
  and the `SDM_SineMassFlux` regression test.
