#!/usr/bin/env python3
"""Write the terrain, sounding and profile files of the Terrain_Inflow_Profile cases.

    python3 gen_inputs.py

Every terrain puts the ground on the inflow face at 150 m or more, so a profile
applied by level or by height above the floor differs clearly from one applied
by height above the ground. Terrain files are in ERF's text format (nx, ny, the x and the y node
coordinates, then the heights with x outer and y inner) on 50 m nodes over
the 1600 m x 800 m domain:

  terrain_plateau.txt     flat ground raised to 150 m: a level inflow face off the floor
  terrain_incline.txt     ground rising downstream from 150 m at the inflow face to 390 m (slope 0.15)
  terrain_decline.txt     ground falling downstream from 290 m at the inflow face to 50 m
  terrain_crossridge.txt  a ridge across the inflow face: 150 m at y = 0 and 800 m, 230 m at y = 400 m
  terrain_flat.txt        flat ground on the floor, for the dirichlet_file parity case

Profiles, the log law of inflow_log_law.* in inputs_inflow (10 m/s at 10 m over
z0 = 0.1 m from 270 degrees, tke from u*^2/Cmu0^2 tapering to 700 u*):

  inflow_profile_crossridge.txt  "# z u v T tke" for xlo.inflow_profile = file
  inflow_levels_flat.txt         z u v w rows for xlo.dirichlet_file
  inflow_profile_flat.txt        "# z u v w", the same rows for xlo.inflow_profile = file
  input_sounding                 a uniform 10 m/s westerly at 300 K for the interior
"""
import math

LX, LY, DX = 1600.0, 800.0, 50.0
NX, NY = int(LX / DX) + 1, int(LY / DX) + 1
KAPPA, CMU0 = 0.41, 0.5562
SPEED, HREF, Z0, TKE_ZSCALE = 10.0, 10.0, 0.1, 700.0
USTAR = KAPPA * SPEED / math.log((HREF + Z0) / Z0)


def write_terrain(name, height):
    with open(name, "w") as f:
        f.write(f"{NX}\n{NY}\n")
        for i in range(NX):
            f.write(f"{i * DX:.2f}\n")
        for j in range(NY):
            f.write(f"{j * DX:.2f}\n")
        for i in range(NX):
            for j in range(NY):
                f.write(f"{height(i * DX, j * DX):.4f}\n")


write_terrain("terrain_plateau.txt", lambda x, y: 150.0)
write_terrain("terrain_incline.txt", lambda x, y: 150.0 + 0.15 * x)
write_terrain("terrain_decline.txt", lambda x, y: 290.0 - 0.15 * x)
write_terrain("terrain_crossridge.txt", lambda x, y: 190.0 - 40.0 * math.cos(2.0 * math.pi * y / LY))
write_terrain("terrain_flat.txt", lambda x, y: 0.0)


def speed(z):
    return USTAR / KAPPA * math.log((z + Z0) / Z0)


def tke(z):
    return USTAR ** 2 / CMU0 ** 2 * max((USTAR * TKE_ZSCALE - z) / (max(USTAR, 0.01) * TKE_ZSCALE), 0.01)


# the ground, then 59 heights spaced geometrically from 1 cm to 1000 m
zs = [0.0] + [0.01 * (1000.0 / 0.01) ** (m / 58.0) for m in range(59)]

with open("inflow_profile_crossridge.txt", "w") as f:
    f.write("# Heights above the local ground; the log law of inputs_inflow\n")
    f.write("# z u v T tke\n")
    for z in zs:
        f.write(f"{z:12.5f} {speed(z):10.5f} {0.0:8.4f} {300.0:7.2f} {tke(z):10.6f}\n")

with open("inflow_levels_flat.txt", "w") as f:
    for z in zs:
        f.write(f"{z:12.5f} {speed(z):10.5f} {0.0:8.4f} {0.0:6.2f}\n")

with open("inflow_profile_flat.txt", "w") as f:
    f.write("# z u v w\n")
    for z in zs:
        f.write(f"{z:12.5f} {speed(z):10.5f} {0.0:8.4f} {0.0:6.2f}\n")

with open("input_sounding", "w") as f:
    f.write("1000.0 300.0 0.0\n")
    f.write("   0.0 300.0 0.0 10.0 0.0\n")
    f.write("1000.0 300.0 0.0 10.0 0.0\n")
