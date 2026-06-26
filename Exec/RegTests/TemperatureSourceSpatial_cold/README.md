# TemperatureSourceSpatial_cold

Idealized convection driven by prescribed, spatially varying source terms of heat and
moisture (after Trapp), the same forcing approach as `SDM_Congestus3D`. The "cold"
configuration extends the column high enough to activate ice processes. Two reference
decks are provided: a warm case and a deeper case that exercises the ice microphysics.

## Input files

- `inputs_SDM` - warm super-droplet case on a 6.4 x 6.4 x 10 km domain (128 x 128 x 200).
- `inputs_SDM_ice` - deep super-droplet case with ice microphysics on a 6.4 x 0.2 x 15 km
  domain (128 x 4 x 300), tracking the six condensate species and the per-species
  (water/ice) super-droplet mass densities.
- `inputs_Kessler` - warm case with the bulk Kessler scheme.
- `inputs_Kessler_regtest` - short variant of `inputs_Kessler` for a quick regression check.
- `input_sounding` - initial theta and qv profiles (idealized initialization).
