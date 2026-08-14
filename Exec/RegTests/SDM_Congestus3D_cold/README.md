# SDM_Congestus3D_cold

Idealized convection driven by prescribed, spatially varying source terms of heat and
moisture (after Trapp), the same forcing approach as `SDM_Congestus3D`. The "cold"
configuration extends the column high enough to activate ice processes. Alongside the
convection decks, the same problem source is used (with `kinematic_mode` and a prescribed
sounding) to drive a set of single-process ice-microphysics verification columns.

## References

- Mitra, S.K., Vohl, O., Ahr, M. and Pruppacher, H.R. (1990), *A wind tunnel and theoretical
  study of the melting behavior of atmospheric ice particles*, J. Atmos. Sci., **47**,
  584-591.
- Niemand, M. et al. (2012), *A particle-surface-area-based parameterization of immersion
  freezing on desert dust particles*, J. Atmos. Sci., **69**, 3077-3092.
  <https://doi.org/10.1175/JAS-D-11-0249.1>
- Beard, K.V. and Grover, S.N. (1974), *Numerical collision efficiencies for small raindrops
  colliding with micron size particles*, J. Atmos. Sci., **31**, 543-550.

## Input files

- `inputs_SDM` - warm super-droplet case on a 6.4 x 6.4 x 10 km domain (128 x 128 x 200).
- `inputs_SDM_ice` - deep super-droplet case with ice microphysics on a 6.4 x 0.2 x 15 km
  domain (128 x 4 x 300), tracking the six condensate species and the per-species
  (water/ice) super-droplet mass densities.
- `inputs_Kessler` - warm case with the bulk Kessler scheme.
- `inputs_Kessler_regtest` - short variant of `inputs_Kessler` for a quick regression check.

## Ice-microphysics verification columns

Single-process tests in a prescribed (kinematic) 1D column, each with its own sounding:

- `inputs_SDM_ice_sublimation` - monodisperse ice injected at the top of an isothermal
  (-15 C), subsaturated-over-ice column falls and loses mass by vapour diffusion.
- `inputs_SDM_ice_melting` - low-density ice flakes injected at the 0 C top fall into warmer
  air and melt, with the mixed-phase terminal velocity blending from ice to water.
- `inputs_SDM_ice_freezing` - supercooled droplets carrying an insoluble INP, injected at the
  warm bottom and lifted by a prescribed updraft, freeze by immersion (Niemand INAS) as they
  cool.
- `inputs_SDM_ice_riming` - solid ice injected at the top of a -10 C water-saturated column
  falls through a supercooled cloud of 10 um droplets and grows by collecting them.

## Sounding files

- `input_sounding` - theta and qv profiles for the convection decks (idealized initialization).
- `input_sounding_sublimation` - isothermal -15 C, RH over ice 82%, 0-1000 m.
- `input_sounding_melting` - 0 C at the 800 m top warming downward to +4.8 C (lapse
  0.6 C/100 m), water-saturated.
- `input_sounding_freezing` - -8 C at the base to -40 C aloft (lapse 8 K/km), water-saturated,
  0-4000 m.
- `input_sounding_riming` - isothermal -10 C, water-saturated, 0-1000 m.
