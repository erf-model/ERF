# RICO: Rain In Cumulus over the Ocean

Idealized large-eddy simulation of precipitating shallow trade-wind cumulus, based on
the GCSS/RICO model intercomparison case. The domain is 12.8 x 12.8 x 4 km on a
128 x 128 x 100 grid, periodic in the horizontal. The flow is driven by prescribed
large-scale forcing (geostrophic wind, subsidence, and advective/radiative tendencies)
together with bulk-aerodynamic surface fluxes (the `rico` surface-layer type). The
default run is 24 hours.

## Reference

- van Zanten, M.C. et al. (2011), *Controls on precipitation and cloudiness in
  simulations of trade-wind cumulus as observed in a composite case based on RICO*,
  J. Adv. Model. Earth Syst., **3**, M06001.
  <https://doi.org/10.1029/2011MS000056>

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics with NH4HSO4 aerosol and two
  log-normal initialization modes; the reference configuration.
- `inputs_Kessler` - same case with the bulk Kessler warm-rain scheme.
- `inputs_Kessler_regtest` - short, deterministic-perturbation variant of
  `inputs_Kessler` for a quick regression check.
- `inputs_SAM` - same forcing with the SAM bulk scheme on a smaller
  6.4 x 6.4 x 4 km / 64 x 64 x 100 grid and prescribed surface fluxes.
- `input_sounding` - initial theta, qv, and (u, v) profiles.
