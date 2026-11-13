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
  abort        the restart stops with the message given

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


class Case:
    def __init__(self, name, marker, init_type, theta_pressure, options, expect, message):
        self.name = name
        self.marker = marker                  # contract version to write, None to remove it
        self.init_type = init_type            # erf.init_type in job_info, None to remove it
        self.theta_pressure = theta_pressure  # "wrfinput", a pressure [Pa], or None
        self.options = options                # extra restart inputs
        self.expect = expect
        self.message = message                # must appear in the restart output


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


def surface_field(checkpoint, name, shape):
    """The valid cells of a planar field of a checkpoint, as an array [j, i]"""
    result = np.full(shape, np.nan)
    for fab in multifab_fabs(checkpoint / "Level_0", name):
        data = fab_data(fab)[0, 0]
        ng = fab["nghost"]
        lo = [fab["lo"][d] + ng[d] for d in range(2)]
        hi = [fab["hi"][d] - ng[d] for d in range(2)]
        result[lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = \
            data[ng[1]:ng[1] + hi[1] - lo[1] + 1, ng[0]:ng[0] + hi[0] - lo[0] + 1]
    if np.isnan(result).any():
        raise AssertionError(f"{checkpoint}/{name} does not cover the domain")
    return result


def to_potential_temperature(checkpoint, pressure):
    """Store SST_0 and TSK_0 of a checkpoint as potential temperature, converted as older
    versions converted them: each cell with the surface pressure of the nearest cell of the
    domain, pressure being an array [j, i] or a constant"""
    for name in SURFACE_FIELDS:
        for fab in multifab_fabs(checkpoint / "Level_0", name):
            data = fab_data(fab).astype(np.float64)
            if np.isscalar(pressure):
                p = pressure
            else:
                ny, nx = pressure.shape
                i = np.clip(np.arange(fab["lo"][0], fab["hi"][0] + 1), 0, nx - 1)
                j = np.clip(np.arange(fab["lo"][1], fab["hi"][1] + 1), 0, ny - 1)
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


def run_case(case, exe, run_dir, original, psfc, env):
    legacy = run_dir / ("legacy_" + case.name)
    output_prefix = case.name + "_chk"
    output = run_dir / (output_prefix + "00001")
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

    command = [str(exe), str(run_dir / "WRFInput_LegacySurfaceTemperature.i"),
               f"erf.restart={legacy.name}", "max_step=1", f"erf.check_file={output_prefix}",
               "erf.check_int=1", "amrex.call_addr2line=0"] + case.options
    with open(log_path, "w") as log:
        status = subprocess.run(command, cwd=run_dir, stdout=log, stderr=subprocess.STDOUT,
                                env=env).returncode
    log_text = log_path.read_text()
    if case.message not in log_text:
        raise AssertionError(f"{case.name}: the restart output ({log_path}) does not say "
                             f"\"{case.message}\"")
    if case.expect == "abort":
        if status == 0:
            raise AssertionError(f"{case.name}: the restart did not abort; see {log_path}")
        return "aborted as it must"
    if status != 0:
        raise AssertionError(f"{case.name}: the restart failed ({status}); see {log_path}")

    report = []
    for name in SURFACE_FIELDS:
        expected = surface_field(original, name, psfc.shape)
        actual = surface_field(output, name, psfc.shape)
        if case.expect == "unchanged":
            if not np.array_equal(actual, expected):
                raise AssertionError(f"{case.name}: {name} changed on restart")
            continue
        if case.expect == "approximate":
            # Its effective surface temperature at the restart, not the SST of the input
            theta = expected * (P_0 / case.theta_pressure) ** KAPPA
            p_restart = restart_surface_pressure(original, psfc.shape)
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
    return ", ".join(report) if report else "bit for bit"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exe", type=Path)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    original = run_dir / "chk00000"
    if (original / "surface_temperature_contract").read_text().split() != ["2"]:
        raise AssertionError("chk00000 is not a version 2 checkpoint")
    with nc.Dataset(run_dir / "wrfinput_d01") as ds:
        psfc = np.asarray(ds["PSFC"][0], dtype=np.float64)

    env = dict(os.environ)
    if "CRAY_MPICH_DIR" in env:
        env["MPICH_GPU_SUPPORT_ENABLED"] = "0"

    failures = 0
    for case in CASES:
        try:
            print(f"{case.name}: {run_case(case, args.exe, run_dir, original, psfc, env)}")
        except AssertionError as error:
            print(f"FAILED {error}")
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
