#!/usr/bin/env python3
"""Restart from the checkpoints older versions of ERF wrote, and check how their SST and
TSK are brought to the absolute temperature [K] they are now kept as.

Without a surface-temperature contract marker, a WRFInput run stored SST_0 and TSK_0 as
potential temperature converted with the PSFC of its wrfinput, and a Metgrid run stored
absolute temperature; with the version 1 marker, both stored potential temperature, a
Metgrid run converting with the surface pressure of its met_em input.  Version 2 stores
absolute temperature.  Each case copies chk00000 (version 2), turns the copy into such a
checkpoint, restarts from it for one step, and compares SST_0 and TSK_0 of the
checkpoint the restart writes with what the migration must give:

  exact        the SST and TSK of chk00000, recovered with the PSFC of wrfinput_d01
  approximate  theta (p_restart/p0)^(R_d/c_p), with p_restart the surface pressure of
               the restarted state; not the SST of the input, which it misses by
               (p_restart/p_input)^(R_d/c_p)
  unchanged    the SST and TSK of chk00000, bit for bit
  abort        the restart stops with the messages given

The mixed cases start from mixed_chk00000 instead, written from wrfinput_mixed, whose west
half is land with the placeholder SST WRF writes there, zero.  The surface layer takes no
SST over land, nor anywhere at a level with TSK, which supplies every cell, so a restart
keeps such placeholders as they are: wherever the legacy checkpoint holds no positive
temperature, the restart must give back what it holds, bit for bit.  An invalid value
where the surface layer takes the field, or a surface pressure that cannot convert a
stored temperature, aborts it.  Some cases remove TSK_0 (as checkpoints written before
it was) or LMASK_0, or zero a cell, before restarting.

The nest cases start from nest_chk00000, written by nest_files.i, whose level 1 lists two
init files, nest_d02 and nest_d03, which starts later and had not been read.  Older
versions restarted such a checkpoint, but which file each level-1 grid took its SST and
TSK from is not recorded, so the pressure to convert them back with is not known, and the
restart aborts.  Restarted once with nest_d02 alone, the only file read, it recovers them
exactly and writes a version-2 checkpoint, which the deck as it was restarts from.

Usage: legacy_restarts.py ERF_EXECUTABLE RUN_DIR
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np


R_D = 287.0
C_P = 1004.5
P_0 = 1.0e5
GAMMA = 1.4
GRAV = 9.81
KAPPA = R_D / C_P
RHO_COMP = 0
RHO_THETA_COMP = 1
# The surface pressure of the met_em input of the Metgrid run that wrote the version 1
# checkpoint of the "approximate" case: no restart can read it again
P_METGRID_INPUT = 101325.0
SURFACE_FIELDS = ("SST_0", "TSK_0")
ASSERTION = "erf.legacy_surface_temperature_init_type"
# The initial checkpoint each case starts from, and the wrfinput it was written from
BASES = {"d01": ("chk00000", "wrfinput_d01"), "mixed": ("mixed_chk00000", "wrfinput_mixed")}
# Each level of nest_files.i: the file it took its SST and TSK from, and the cells (from,
# number) that file covers
NEST_LEVELS = ((0, "nest_d01", (0, 0), (16, 16)), (1, "nest_d02", (6, 6), (12, 12)))
# Cells (i, j) of wrfinput_mixed: land, where wrfinput_mixed_badpsfc has PSFC -1 Pa, and
# water
LAND_CELL = (3, 9)
WATER_CELL = (12, 5)
# The cells of wrfinput_mixed whose SST is a placeholder: its land
LAND_CELLS = 128


class Case:
    def __init__(self, name, marker, init_type, theta_pressure, options, expect, messages,
                 base="d01", remove=(), zero=()):
        self.name = name
        self.marker = marker                  # contract version to write, None to remove it
        self.init_type = init_type            # erf.init_type in job_info, None to remove it
        self.theta_pressure = theta_pressure  # "wrfinput", a pressure [Pa], or None
        self.options = options                # extra restart inputs
        self.expect = expect
        # must all appear in the restart output
        self.messages = (messages,) if isinstance(messages, str) else tuple(messages)
        self.base = base                      # the initial checkpoint it starts from (BASES)
        self.remove = remove                  # fields to remove from the legacy checkpoint
        self.zero = zero                      # (field, (i, j)) to set to zero in it


def kept(cells):
    """What the restart says of the SST placeholders it keeps at level 0"""
    return (f"converting the SST at level 0 to absolute temperature: kept {cells} cells "
            "where the surface layer does not take the SST as they are")


EXACT_MIXED = "EXACT recovery at level 0 with the PSFC of wrfinput_mixed"


CASES = [
    # A WRFInput checkpoint from before the marker, restarted as it was run
    Case("markerless_wrfinput", None, "WRFInput", "wrfinput", [], "exact",
         "EXACT recovery at level 0 with the PSFC of wrfinput_d01"),
    # The same, restarted with two level-0 init files: level 0 read its SST, TSK and PSFC
    # from the first only, and its PSFC must be read from that one, not left unread
    Case("two_level0_init_files", None, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_d01 wrfinput_d01b"], "exact",
         "EXACT recovery at level 0 with the PSFC of wrfinput_d01"),
    # A version 1 checkpoint written by a Metgrid run, restarted by a WRFInput run: the
    # PSFC of these wrfinput files is not the pressure it was converted with
    Case("v1_metgrid_on_wrfinput_restart", 1, "Metgrid", P_METGRID_INPUT, [], "approximate",
         "APPROXIMATE migration at level 0"),
    # A version 1 checkpoint written by a WRFInput run
    Case("v1_wrfinput", 1, "WRFInput", "wrfinput", [], "exact",
         "EXACT recovery at level 0 with the PSFC of wrfinput_d01"),
    # A markerless checkpoint whose job_info does not record the init type, with and
    # without the assertion
    Case("markerless_asserted", None, None, "wrfinput", [ASSERTION + "=WRFInput"], "exact",
         "erf.legacy_surface_temperature_init_type asserts that a WRFInput run converted"),
    Case("markerless_unasserted", None, None, "wrfinput", [], "abort",
         "Set erf.legacy_surface_temperature_init_type"),
    # A WRFInput checkpoint restarted by a run that reads no wrfinput
    Case("wrfinput_on_non_wrfinput_restart", None, "WRFInput", "wrfinput",
         ["erf.init_type=Uniform"], "abort", "restart with erf.init_type = WRFInput"),
    # Version 2: absolute temperature, restored as it is
    Case("v2", 2, "WRFInput", None, [], "unchanged", "Reading SST data"),

    # Part land, part water.  With TSK, which supplies every cell, the SST over land is a
    # placeholder no surface layer takes: kept as it is, while the SST over water and the
    # TSK are recovered
    Case("mixed_markerless", None, "WRFInput", "wrfinput", [], "exact",
         [EXACT_MIXED, kept(LAND_CELLS)], base="mixed"),
    Case("mixed_v1", 1, "WRFInput", "wrfinput", [], "exact",
         [EXACT_MIXED, kept(LAND_CELLS)], base="mixed"),
    # A zero SST over water is a placeholder too at a level with TSK
    Case("mixed_v1_water_sst_zero_with_tsk", 1, "WRFInput", "wrfinput", [], "exact",
         kept(LAND_CELLS + 1), base="mixed", zero=[("SST_0", WATER_CELL)]),
    # Without TSK, as checkpoints written before TSK_0 was, land takes the prescribed
    # temperature and keeps its SST placeholders, and water takes the SST ...
    Case("mixed_v1_without_tsk", 1, "WRFInput", "wrfinput", [], "exact",
         kept(LAND_CELLS), base="mixed", remove=["TSK_0"]),
    # ... so a zero SST over water aborts
    Case("mixed_v1_without_tsk_water_sst_zero", 1, "WRFInput", "wrfinput", [], "abort",
         [kept(LAND_CELLS),
          f"({WATER_CELL[0]},{WATER_CELL[1]}), a water cell where the surface layer takes "
          "the SST, has theta = 0,"],
         base="mixed", remove=["TSK_0"], zero=[("SST_0", WATER_CELL)]),
    # TSK supplies every cell, so a zero TSK over land aborts
    Case("mixed_v1_land_tsk_zero", 1, "WRFInput", "wrfinput", [], "abort",
         f"({LAND_CELL[0]},{LAND_CELL[1]}), a land cell where the surface layer takes the "
         "TSK, has theta = 0,", base="mixed", zero=[("TSK_0", LAND_CELL)]),
    # A surface pressure of -1 Pa leaves the SST placeholder of its cell as it is, but
    # cannot convert the TSK there, which aborts
    Case("mixed_v1_invalid_psfc", 1, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_mixed_badpsfc"], "abort",
         [kept(LAND_CELLS),
          f"({LAND_CELL[0]},{LAND_CELL[1]}), a land cell where the surface layer takes the "
          "TSK", "surface pressure = -1 Pa"], base="mixed"),
    # Without LMASK_0 the restart takes the land mask from erf.is_land, land by default,
    # as the surface layer will ...
    Case("mixed_v1_without_lmask", 1, "WRFInput", "wrfinput", [], "exact",
         kept(LAND_CELLS), base="mixed", remove=["LMASK_0"]),
    # ... and water with erf.is_land = 0, where, without TSK, every placeholder is taken
    Case("mixed_v1_without_lmask_or_tsk_on_water", 1, "WRFInput", "wrfinput",
         ["erf.is_land=0"], "abort",
         ["a water cell where the surface layer takes the SST, has theta = 0,",
          f"{LAND_CELLS} cells could not be converted"],
         base="mixed", remove=["LMASK_0", "TSK_0"]),
    # Version 2, placeholders and all, restored as it is
    Case("mixed_v2", 2, "WRFInput", None, [], "unchanged", "Reading SST data", base="mixed"),
]


def multifab_fabs(level_dir, name):
    """Each FAB of a VisMF MultiFab: data file, offset of its data, dtype, lo, hi, ncomp,
    and the number of ghost cells"""
    header = (level_dir / (name + "_H")).read_text()
    lines = header.splitlines()
    ncomp = int(lines[2])
    nghost = [int(v) for v in re.findall(r"-?\d+", lines[3])]
    nghost = (nghost * 3)[:3] if len(nghost) == 1 else nghost[:3]
    fabs = []
    for filename, offset in re.findall(r"FabOnDisk:\s+(\S+)\s+(\d+)", header):
        path = level_dir / filename
        with open(path, "rb") as fh:
            fh.seek(int(offset))
            fab_header = fh.readline().decode()
            start = fh.tell()
        # FAB ((8, (64 11 52 0 1 12 0 1023)),(8, (8 7 6 5 4 3 2 1)))((lo) (hi) (type)) n
        descriptor = re.search(r"\((\d+), \(([\d ]+)\)\)\)\(", fab_header)
        nbytes = int(descriptor.group(1))
        little = descriptor.group(2).split()[0] == str(nbytes)
        dtype = np.dtype(("<" if little else ">") + ("f8" if nbytes == 8 else "f4"))
        lo, hi = re.findall(r"\(([-\d,]+)\)", fab_header.split(")))")[1])[:2]
        fabs.append(dict(path=path, start=start, dtype=dtype, ncomp=ncomp, nghost=nghost,
                         lo=[int(v) for v in lo.split(",")],
                         hi=[int(v) for v in hi.split(",")]))
    return fabs


def fab_data(fab):
    """The data of a FAB as an array indexed [n, k, j, i]"""
    shape = (fab["ncomp"],) + tuple(h - l + 1 for l, h in zip(reversed(fab["lo"]),
                                                              reversed(fab["hi"])))
    return np.fromfile(fab["path"], dtype=fab["dtype"], count=int(np.prod(shape)),
                       offset=fab["start"]).reshape(shape)


def write_fab_data(fab, data):
    with open(fab["path"], "r+b") as fh:
        fh.seek(fab["start"])
        fh.write(np.ascontiguousarray(data, dtype=fab["dtype"]).tobytes())


def surface_field(checkpoint, name, shape, level=0, origin=(0, 0)):
    """The valid cells of a planar field of a checkpoint at a level, as an array [j, i] over
    the cells from origin (i, j)"""
    result = np.full(shape, np.nan)
    for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
        data = fab_data(fab)[0, 0]
        ng = fab["nghost"]
        lo = [fab["lo"][d] + ng[d] - origin[d] for d in range(2)]
        hi = [fab["hi"][d] - ng[d] - origin[d] for d in range(2)]
        if min(lo) < 0 or hi[0] >= shape[1] or hi[1] >= shape[0]:
            raise AssertionError(f"{checkpoint}/{name} has a grid outside the cells checked")
        result[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = \
            data[ng[1]:ng[1] + hi[1] - lo[1] + 1, ng[0]:ng[0] + hi[0] - lo[0] + 1]
    if np.isnan(result).any():
        raise AssertionError(f"{checkpoint}/{name} does not cover the cells checked")
    return result


def to_potential_temperature(checkpoint, pressure, level=0, origin=(0, 0)):
    """Store SST_0 and TSK_0 of a checkpoint as potential temperature, converted as older
    versions converted them: each cell with the surface pressure of the nearest cell of the
    file, pressure being an array [j, i] over the cells from origin (i, j), or a constant"""
    for name in SURFACE_FIELDS:
        for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
            data = fab_data(fab).astype(np.float64)
            if np.isscalar(pressure):
                p = pressure
            else:
                ny, nx = pressure.shape
                i = np.clip(np.arange(fab["lo"][0], fab["hi"][0] + 1) - origin[0], 0, nx - 1)
                j = np.clip(np.arange(fab["lo"][1], fab["hi"][1] + 1) - origin[1], 0, ny - 1)
                p = pressure[np.ix_(j, i)]
            data[0, 0] = data[0, 0] * (P_0 / p) ** KAPPA
            write_fab_data(fab, data)


def restart_surface_pressure(checkpoint, shape):
    """The surface pressure the restarted state gives, from its lowest cell, as the
    restart computes it: p(rho theta) + rho g (z_cell - z_surface), with no moisture"""
    ny, nx = shape
    rho = np.full(shape, np.nan)
    rho_theta = np.full(shape, np.nan)
    for fab in multifab_fabs(checkpoint / "Level_0", "Cell"):
        data = fab_data(fab)
        lo, hi = fab["lo"], fab["hi"]
        if lo[2] == 0:
            rho[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = data[RHO_COMP, 0]
            rho_theta[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = data[RHO_THETA_COMP, 0]
    z = np.full((2, ny + 1, nx + 1), np.nan)
    for fab in multifab_fabs(checkpoint / "Level_0", "Z_Phys_nd"):
        data = fab_data(fab)[0]
        lo, hi = fab["lo"], fab["hi"]
        for k in (0, 1):
            if not lo[2] <= k <= hi[2]:
                continue
            i0, i1 = max(lo[0], 0), min(hi[0], nx)
            j0, j1 = max(lo[1], 0), min(hi[1], ny)
            z[k, j0:j1 + 1, i0:i1 + 1] = \
                data[k - lo[2], j0 - lo[1]:j1 - lo[1] + 1, i0 - lo[0]:i1 - lo[0] + 1]
    if np.isnan(rho).any() or np.isnan(z).any():
        raise AssertionError(f"{checkpoint} does not cover the lowest cells")
    corners = z[:, :-1, :-1] + z[:, :-1, 1:] + z[:, 1:, :-1] + z[:, 1:, 1:]
    delta_z = 0.125 * (corners[1] - corners[0])
    return P_0 * (R_D * rho_theta / P_0) ** GAMMA + rho * GRAV * delta_z


def edit_job_info(checkpoint, init_type):
    """Record init_type as the erf.init_type of the job_info, or remove it if None"""
    path = checkpoint / "job_info"
    lines = path.read_text().splitlines(keepends=True)
    assignment = re.compile(r"^\s*erf\.init_type(\s|=)")
    found = any(assignment.match(line) for line in lines)
    if not found:
        raise AssertionError(f"{path} records no erf.init_type")
    edited = []
    for line in lines:
        if assignment.match(line):
            if init_type is not None:
                edited.append(f'erf.init_type = "{init_type}"\n')
        else:
            edited.append(line)
    path.write_text("".join(edited))


def zero_cell(checkpoint, name, cell):
    """Set a planar field of a checkpoint to zero at cell (i, j), ghost cells included"""
    i, j = cell
    found = False
    for fab in multifab_fabs(checkpoint / "Level_0", name):
        lo, hi = fab["lo"], fab["hi"]
        if lo[0] <= i <= hi[0] and lo[1] <= j <= hi[1]:
            data = fab_data(fab).copy()
            data[:, :, j - lo[1], i - lo[0]] = 0.0
            write_fab_data(fab, data)
            found = True
    if not found:
        raise AssertionError(f"{checkpoint}/{name} does not hold the cell {cell}")


def remove_field(checkpoint, name):
    """Remove a MultiFab from a checkpoint, as if it had not been written"""
    level_dir = checkpoint / "Level_0"
    paths = [level_dir / (name + "_H")] + sorted(level_dir.glob(name + "_D_*"))
    if not paths[0].exists():
        raise AssertionError(f"{checkpoint} holds no {name}")
    for path in paths:
        path.unlink()


def restart(exe, run_dir, inputs, checkpoint, name, options, env, max_step=1):
    """Restart from a checkpoint to max_step, writing name_chk<step>: the exit status and
    the output"""
    command = [str(exe), str(run_dir / inputs), f"erf.restart={checkpoint}",
               f"max_step={max_step}", f"erf.check_file={name}_chk", "erf.check_int=1",
               "amrex.call_addr2line=0"] + options
    log_path = run_dir / (name + ".log")
    with open(log_path, "w") as log:
        status = subprocess.run(command, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT,
                                env=env).returncode
    return status, log_path.read_text()


def run_case(case, exe, run_dir, bases, env):
    wrfinput = BASES[case.base][1]
    original, psfc = bases[case.base]
    legacy = run_dir / ("legacy_" + case.name)
    output = run_dir / (case.name + "_chk00001")
    log_path = run_dir / (case.name + ".log")
    for path in (legacy, output):
        shutil.rmtree(path, ignore_errors=True)
    shutil.copytree(original, legacy)

    marker = legacy / "surface_temperature_contract"
    if case.marker is None:
        marker.unlink()
    else:
        marker.write_text(f"{case.marker}\n")
    edit_job_info(legacy, case.init_type)
    if case.theta_pressure == "wrfinput":
        to_potential_temperature(legacy, psfc)
    elif case.theta_pressure is not None:
        to_potential_temperature(legacy, case.theta_pressure)
    for name in case.remove:
        remove_field(legacy, name)
    for name, cell in case.zero:
        zero_cell(legacy, name, cell)
    # The legacy fields as the restart reads them, None where removed
    stored = {name: None if name in case.remove else surface_field(legacy, name, psfc.shape)
              for name in SURFACE_FIELDS}

    options = list(case.options)
    if not any(option.startswith("erf.nc_init_file_0=") for option in options):
        options.append(f"erf.nc_init_file_0={wrfinput}")
    status, log_text = restart(exe, run_dir, "WRFInput_LegacySurfaceTemperature.i",
                               legacy.name, case.name, options, env)
    for message in case.messages:
        if message not in log_text:
            raise AssertionError(f"{case.name}: the restart output ({log_path}) does not "
                                 f"say \"{message}\"")
    if case.expect == "abort":
        if status == 0:
            raise AssertionError(f"{case.name}: the restart did not abort; see {log_path}")
        return "aborted as it must"
    if status != 0:
        raise AssertionError(f"{case.name}: the restart failed ({status}); see {log_path}")

    report = []
    for name in SURFACE_FIELDS:
        if stored[name] is None:
            if (output / "Level_0" / (name + "_H")).exists():
                raise AssertionError(f"{case.name}: the restart wrote the {name} removed")
            report.append(f"{name} absent")
            continue
        expected = surface_field(original, name, psfc.shape)
        actual = surface_field(output, name, psfc.shape)
        if case.expect == "unchanged":
            if not np.array_equal(actual, expected):
                raise AssertionError(f"{case.name}: {name} changed on restart")
            continue
        # Placeholders, which hold no positive temperature, come back as they are
        placeholder = ~(stored[name] > 0.0)
        if not np.array_equal(actual[placeholder], stored[name][placeholder]):
            raise AssertionError(f"{case.name}: {name} changed the placeholders it holds")
        actual = actual[~placeholder]
        expected = expected[~placeholder]
        if case.expect == "approximate":
            # Its effective surface temperature at the restart, not the SST of the input
            theta = expected * (P_0 / case.theta_pressure) ** KAPPA
            p_restart = restart_surface_pressure(original, psfc.shape)[~placeholder]
            bias = np.abs(actual - expected).min()
            if bias < 5.0:
                raise AssertionError(f"{case.name}: {name} is within {bias} K of the SST of "
                                     "the input, which this migration cannot recover")
            expected = theta * (p_restart / P_0) ** KAPPA
            report.append(f"{name} {bias:.1f} K or more from the input")
        error = np.abs(actual / expected - 1.0).max()
        if error > 1.0e-10:
            raise AssertionError(f"{case.name}: {name} is off by a relative {error:.3e}")
        report.append(f"{name} to {error:.1e}")
        if placeholder.any():
            report.append(f"{placeholder.sum()} {name} placeholders kept")
    return ", ".join(report) if report else "bit for bit"


def nest_fields(checkpoint):
    """SST_0 and TSK_0 of a checkpoint of nest_files.i at each level"""
    return {(level, name): surface_field(checkpoint, name, shape, level, origin)
            for level, _, origin, shape in NEST_LEVELS for name in SURFACE_FIELDS}


def run_nest_cases(exe, run_dir, env):
    """The nest cases (see above): the name and the result of each"""
    original = run_dir / "nest_chk00000"
    if (original / "surface_temperature_contract").read_text().split() != ["2"]:
        raise AssertionError("nest_chk00000 is not a version 2 checkpoint")
    expected = nest_fields(original)
    legacy = run_dir / "legacy_nest_v1"
    shutil.rmtree(legacy, ignore_errors=True)
    shutil.copytree(original, legacy)
    (legacy / "surface_temperature_contract").write_text("1\n")
    for level, wrfinput, origin, _ in NEST_LEVELS:
        with nc.Dataset(run_dir / wrfinput) as ds:
            psfc = np.asarray(ds["PSFC"][0], dtype=np.float64)
        to_potential_temperature(legacy, psfc, level, origin)
    results = []

    # The deck as it was
    name = "nest_v1_two_level1_init_files"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name, [], env)
    message = ("at level 1 cannot be recovered: the level lists 2 init files "
               "(erf.nc_init_file_1), and the checkpoint does not record which of them")
    if status == 0 or message not in log:
        raise AssertionError(f"{name}: the restart did not abort saying \"{message}\"; "
                             f"see {name}.log")
    results.append((name, "aborted as it must"))

    # Once with nest_d02 alone, the only file read
    name = "nest_v1_first_level1_init_file"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name,
                          ["erf.nc_init_file_1=nest_d02"], env)
    message = "EXACT recovery at level 1 with the PSFC of nest_d02"
    if status != 0 or message not in log:
        raise AssertionError(f"{name}: the restart failed ({status}) or did not say "
                             f"\"{message}\"; see {name}.log")
    migrated = run_dir / (name + "_chk00001")
    actual = nest_fields(migrated)
    error = max(np.abs(actual[key] / expected[key] - 1.0).max() for key in expected)
    if error > 1.0e-10:
        raise AssertionError(f"{name}: the SST or TSK is off by a relative {error:.3e}")
    results.append((name, f"SST_0 and TSK_0 at both levels to {error:.1e}"))

    # Then the deck as it was, from the version-2 checkpoint that wrote
    name = "nest_v2_two_level1_init_files"
    status, log = restart(exe, run_dir, "nest_files.i", migrated.name, name, [], env,
                          max_step=2)
    if status != 0:
        raise AssertionError(f"{name}: the restart failed ({status}); see {name}.log")
    continued = nest_fields(run_dir / (name + "_chk00002"))
    if any(not np.array_equal(continued[key], actual[key]) for key in actual):
        raise AssertionError(f"{name}: the SST or TSK changed on restart")
    results.append((name, "bit for bit"))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exe", type=Path)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    bases = {}
    for base, (checkpoint, wrfinput) in BASES.items():
        original = run_dir / checkpoint
        if (original / "surface_temperature_contract").read_text().split() != ["2"]:
            raise AssertionError(f"{checkpoint} is not a version 2 checkpoint")
        with nc.Dataset(run_dir / wrfinput) as ds:
            psfc = np.asarray(ds["PSFC"][0], dtype=np.float64)
        bases[base] = (original, psfc)

    env = dict(os.environ)
    if "CRAY_MPICH_DIR" in env:
        env["MPICH_GPU_SUPPORT_ENABLED"] = "0"

    failures = 0
    for case in CASES:
        try:
            print(f"{case.name}: {run_case(case, args.exe, run_dir, bases, env)}")
        except AssertionError as error:
            print(f"FAILED {error}")
            failures += 1
    try:
        for name, result in run_nest_cases(args.exe, run_dir, env):
            print(f"{name}: {result}")
    except AssertionError as error:
        print(f"FAILED {error}")
        failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
