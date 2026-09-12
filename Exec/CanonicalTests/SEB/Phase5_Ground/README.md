# SEB/Phase5_Ground

Phase 5 of the immersed-boundary surface energy balance: conduction into
the walls and roofs through a slab per face, and the material library. The
skin temperature is still prescribed; the slab beneath it now evolves.

```
./run_ground.sh /path/to/erf_exec        # NP=4 by default
python3 plot_ground.py faces_thin        # G of every face against height
```

## The slab

Each face carries `n_slab_layers` uniform layers between its skin and the
building interior at `T_interior`. The solver is the SLUCM branch's implicit
tridiagonal conduction with one change: the skin temperature is the top
boundary instead of a flux, since that is what the balance holds, and the
conduction into the slab, G, is what comes out. Layer centres sit half a
layer below the skin and half a layer above the interior, so both boundary
fluxes use 2k/dz. Unconditionally stable; up to 256 layers.

## Materials

`erf.ibseb.material_file` names a CSV in the SLUCM schema (`mat_id, name,
albedo, emissivity, k_therm_W_per_mK, rho_cp_J_per_m3K, thickness_m,
description`); `material_default` applies to every building and
`material_by_building` gives one id per building in building order. The
albedo and emissivity now enter the shortwave and longwave per face. Without
a file the uniform `k_therm`, `rho_cp`, `thickness`, `albedo` and
`emissivity` inputs apply.

## The scenario

Two 40 m cubes on a 10 m grid in a light wind, faces held at 320 K over
slabs that start at 320 K with a 300 K interior; the flux to the air is
switched off so only the slab does anything and the runs are cheap.

## What is checked

1. **Thick slab, semi-infinite response** (200 mm concrete in 1 mm layers,
   50 s): the wave from the cold interior travels about 6 mm, so the flux
   through the skin stays zero, the top layer stays at the skin temperature,
   and the bottom layer follows the erfc solution for a boundary step.
2. **Thin slab, steady state** (20 mm, light, 100 s): G equals k dT / L on
   every face within 1 percent and the layer temperatures lie on the
   straight line between skin and interior.
3. **Materials by building**: building 1 carries the concrete and building 2
   the timber of the CSV, with all five properties.
4. **Restart**: the thin deck checkpointed at step 200 and restarted to 400
   reproduces the straight run's slab, conduction flux and geometry exactly.
   The atmosphere-derived columns (wind, density, air temperature, sensible
   flux) are only required to be close, because the immersed-forcing
   atmosphere of the development branch does not restart bit-for-bit: the
   wind at the faces differs by about one part in ten thousand after the
   restart. That is outside the balance and worth a separate look.

## Reference

```
== thick slab, semi-infinite response (4 ranks, 50 s)
  skin flux stays zero while the wave from the interior is far away: PASS (max |G| 0.00e+00 W/m2)
  bottom layer vs semi-infinite erfc solution at t = 50.0 s: expected 300.921 K, got 300.925-300.925 K -> PASS
  top layer still at the skin temperature: PASS
== thin slab, steady linear profile (4 ranks, 100 s)
  thin slab at steady state: G = k dT / L = 1000.0 W/m2, got 1000.0-1000.0 -> PASS
  linear profile: top layer 317.50 K, bottom 302.50 K -> PASS
== materials by building (4 ranks)
  building 1 carries material 1 (albedo 0.25, k 1.5 W/m/K, 0.3 m): PASS
  building 2 carries material 2 (albedo 0.4, k 0.15 W/m/K, 0.15 m): PASS
== slab through a checkpoint (4 ranks)
  slab, G and geometry after restart equal the straight run: PASS (max |diff| 0.0e+00)
  atmosphere-derived columns close (development IF restart is not bit-exact): PASS (max rel diff 2.5e-05)
ALL PASS

```
