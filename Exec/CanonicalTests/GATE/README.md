# GATE

Idealized tropical maritime convection on a 12.8 x 12.8 x 4 km domain (128 x 128 x 100),
periodic in the horizontal. The large-scale forcing (geostrophic wind, subsidence, and
advective tendencies) and the bulk-aerodynamic surface fluxes are configured as in the
RICO case (see `../RICO`); the case is named for the GARP Atlantic Tropical Experiment.

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics with NH4HSO4 aerosol and two
  log-normal initialization modes; the reference configuration.
- `inputs_Kessler` - same case with the bulk Kessler warm-rain scheme.
- `inputs_Kessler_regtest` - short variant of `inputs_Kessler` for a quick regression check.
- `inputs_SAM` - same case with the SAM bulk scheme.
- `inputs_Morrison` - same case with the Morrison two-moment scheme.
- `input_sounding` - initial theta, qv, and (u, v) profiles.
