# DYCOMS-II RF01

Idealized large-eddy simulation of nocturnal marine stratocumulus, based on the DYCOMS-II
RF01 (research flight 1) intercomparison case. The domain is 1.6 x 1.6 x 1.6 km on a
32 x 32 x 320 grid (50 m horizontal, 5 m vertical), periodic in the horizontal, with
prescribed surface fluxes, geostrophic and subsidence forcing, simple longwave radiation,
and super-droplet microphysics.

## Reference

- Stevens, B. et al. (2005), *Evaluation of Large-Eddy Simulations via Observations of
  Nocturnal Marine Stratocumulus*, Mon. Wea. Rev., **133**, 1443-1462.
  <https://doi.org/10.1175/MWR2930.1>

## Input files

- `inputs_DYCOMS2RF01` - super-droplet (Lagrangian) microphysics with NH4HSO4 aerosol and
  simple longwave radiation; the reference configuration.
- `input_sounding` - initial theta and qv profiles (idealized initialization).
