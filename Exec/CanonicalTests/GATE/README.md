# GATE

Idealized tropical convection following the GATE Phase III cloud-resolving-model case of
Grabowski, Wu, and Moncrieff (1996), on a 12.8 x 12.8 x 4 km domain (128 x 128 x 100),
periodic in the horizontal. The large-scale forcing is prescribed height-dependent profiles
of the apparent heat source (a theta tendency), the apparent moisture sink (a qv tendency),
and the geostrophic wind, fit to the GATE composite data; the large-scale subsidence is set
to zero (w_bar = 0, as required with periodic lateral boundaries). Only the bulk-aerodynamic
surface fluxes use the RICO formulation (see `../RICO`). The case is named for the GARP
Atlantic Tropical Experiment (GATE, 1974).

## Reference

- Grabowski, W.W., Wu, X. and Moncrieff, M.W. (1996), *Cloud-resolving modeling of tropical cloud
  systems during Phase III of GATE. Part I: Two-dimensional experiments*, J. Atmos. Sci., **53**,
  3684-3709. `10.1175/1520-0469(1996)053<3684:CRMOTC>2.0.CO;2`

## Input files

- `inputs_SDM` - super-droplet (Lagrangian) microphysics with NH4HSO4 aerosol and two
  log-normal initialization modes; the reference configuration.
- `inputs_Kessler` - same case with the bulk Kessler warm-rain scheme.
- `inputs_Kessler_regtest` - short variant of `inputs_Kessler` for a quick regression check.
- `inputs_SAM` - same case with the SAM bulk scheme.
- `inputs_Morrison` - same case with the Morrison two-moment scheme.
- `input_sounding` - initial theta, qv, and (u, v) profiles.
