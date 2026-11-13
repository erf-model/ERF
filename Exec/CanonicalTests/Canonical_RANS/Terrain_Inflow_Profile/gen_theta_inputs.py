#!/usr/bin/env python3
"""Write the sounding and profile files of the theta-above-ground case (inputs_theta).

    python3 gen_theta_inputs.py

Theta is 300 K up to 300 m above the ground, rises by 8 K to 400 m, and by
3 K/km above that, with a uniform 10 m/s westerly. Over the incline (ground
150 m at the inflow face to 390 m at the outflow) this inversion sits 150 to
390 m higher in absolute height at the outflow end than at the inflow face, so
a sounding read at physical heights and one read above the ground differ by up
to 8 K in the interior.

  input_sounding_inversion        z theta qv u v, heights read above the ground with
                                  erf.input_sounding_theta_above_ground = true
  inflow_profile_inversion.txt    "# z u v T", the same column for xlo.inflow_profile = file
"""

Z_INV, DTH_INV, DZ_INV, LAPSE = 300.0, 8.0, 100.0, 0.003
SPEED = 10.0


def theta(z):
    if z <= Z_INV:
        return 300.0
    if z <= Z_INV + DZ_INV:
        return 300.0 + DTH_INV * (z - Z_INV) / DZ_INV
    return 300.0 + DTH_INV + LAPSE * (z - Z_INV - DZ_INV)


zs = [0.0, 100.0, 200.0, 300.0, 325.0, 350.0, 375.0, 400.0, 600.0, 800.0, 1000.0, 1200.0]

with open("input_sounding_inversion", "w") as f:
    f.write("1000.0 300.0 0.0\n")
    for z in zs:
        f.write(f"{z:8.1f} {theta(z):8.3f} 0.0 {SPEED:5.1f} 0.0\n")

with open("inflow_profile_inversion.txt", "w") as f:
    f.write("# Heights above the local ground: the column of input_sounding_inversion\n")
    f.write("# z u v T\n")
    for z in zs:
        f.write(f"{z:8.1f} {SPEED:6.2f} {0.0:5.2f} {theta(z):8.3f}\n")
