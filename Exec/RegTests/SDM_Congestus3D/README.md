# SDM_Congestus3D

Idealized three-dimensional congestus (towering cumulus) convection driven by prescribed
large-scale advective tendencies of heat and moisture. The domain is 6.4 x 6.4 x 10 km on
a 128 x 128 x 200 grid, periodic in the horizontal, capped by a slip wall. The surface
layer supplies momentum drag only (zero sensible and latent heat flux), so convection is
forced entirely by the spatially varying `rhotheta` and moisture source terms; their
magnitudes follow the values annotated after Trapp in the input files.

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics with NH4HSO4 aerosol and two
  log-normal initialization modes; the reference configuration.
- `inputs_Kessler` - same dynamics and forcing with the bulk Kessler scheme.
- `inputs_Kessler_regtest` - short, deterministic-perturbation variant of
  `inputs_Kessler` for a quick regression check.
- `input_sounding` - initial theta and qv profiles (read because
  `erf.init_type = input_sounding`).
