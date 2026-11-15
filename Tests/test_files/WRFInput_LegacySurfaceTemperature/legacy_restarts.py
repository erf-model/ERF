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

  input        the SST and TSK of chk00000, recovered with the available wrfinput PSFC
  approximate  theta (p_restart/p0)^(R_d/c_p), using diagnosed ground pressure or parent
               ground pressure on an aloft-only fine level
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

At a level without an SST the surface layer takes neither field, so with SST_0 removed a
zero TSK is a placeholder too, kept as it is, while every positive TSK is still converted;
with the SST, the same TSK is taken at every cell and aborts the restart.  The version 2
cases with the same edits check that the run takes the TSK just so: with the SST it
rejects such a TSK, and without it runs, restoring the TSK as it is.  These cases restart
with a surface layer below (SURFACE_LAYER), which the decks do not have, so that the run
takes the SST and TSK wherever it takes them at all.

The nest cases start from nest_chk00000, written by nest_files.i, whose level 1 lists two
init files, nest_d02 and nest_d03, which starts 10 s later. A run reads a later file only
once it starts, so until then level 1 holds what nest_d02 alone gave it. In this no-regrid
fixture, the selected first-file PSFC reproduces those values numerically; that does not
establish exact provenance for arbitrary restarted fine grids. The version-2 checkpoint it
writes restarts bit for bit. Once a later file
had started by the time of the checkpoint, the run that wrote it may have read that file
as well, or regridded the level for it, which is not recorded, and the restart aborts:
with nest_d03_now, nest_d03 starting at the start, listed second, and from the checkpoint
nest_later_chk00010, at t = 10 s, of the deck run on, while its nest_later_chk00005, at
t = 5 s, uses the first-file pressure policy.

A level 1 without an SST of its own is given level 0's on restart, so the surface layer
takes its TSK, and a zero there aborts the restart, of a version 1 checkpoint and the run
from a version 2 one alike; with no SST at any level, it is a placeholder, kept.  These
cases restart with a surface layer below too.

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
BASES = {
    "d01": ("chk00000", "wrfinput_d01", (0, 0), (16, 16)),
    "mixed": ("mixed_chk00000", "wrfinput_mixed", (0, 0), (16, 16)),
    "crop": ("crop_chk00000", "wrfinput_crop", (0, 0), (8, 8)),
}
# Each level of nest_files.i: the file it took its SST and TSK from, and the cells (from,
# number) that file covers
NEST_LEVELS = ((0, "nest_d01", (0, 0), (16, 16)), (1, "nest_d02", (6, 6), (12, 12)))
# Cells (i, j) of wrfinput_mixed: land, where wrfinput_mixed_badpsfc has PSFC -1 Pa, and
# water
LAND_CELL = (3, 9)
WATER_CELL = (12, 5)
# The cells of wrfinput_mixed whose SST is a placeholder: its land
LAND_CELLS = 128
# A zero TSK over land and over water
TSK_ZEROS = ("TSK_0", LAND_CELL), ("TSK_0", WATER_CELL)
# A cell (i, j) of level 1 of nest_files.i, which nest_d02 covers
NEST_CELL = (10, 13)
# Restart inputs that put a surface layer below, which takes its surface temperature from
# the SST and TSK (the decks have a slip wall there), with the closure it needs
SURFACE_LAYER = ["zlo.type=surface_layer", "erf.les_type=Smagorinsky", "erf.Cs=0.1"]


class Case:
    def __init__(self, name, marker, init_type, theta_pressure, options, expect, messages,
                 base="d01", remove=(), zero=(), input_file="WRFInput_LegacySurfaceTemperature.i"):
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
        self.input_file = input_file


def kept(cells, field="SST", level=0):
    """What the restart says of the placeholders of a field it keeps at a level"""
    return (f"converting the {field} at level {level} to absolute temperature: kept {cells} "
            f"cells where the surface layer does not take the {field} as they are")


INPUT_MIXED = "pressure recovery at level 0 with the PSFC of wrfinput_mixed"


CASES = [
    # A WRFInput checkpoint from before the marker, restarted as it was run
    Case("markerless_wrfinput", None, "WRFInput", "wrfinput", [], "input",
         "pressure recovery at level 0 with the PSFC of wrfinput_d01"),
    # The same, restarted with two level-0 init files: level 0 read its SST, TSK and PSFC
    # from the first only, and its PSFC must be read from that one, not left unread
    Case("two_level0_init_files", None, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_d01 wrfinput_d01b"], "input",
         "pressure recovery at level 0 with the PSFC of wrfinput_d01"),
    # A version 1 checkpoint written by a Metgrid run, restarted by a WRFInput run: the
    # PSFC of these wrfinput files is not the pressure it was converted with
    Case("v1_metgrid_on_wrfinput_restart", 1, "Metgrid", P_METGRID_INPUT, [], "approximate",
         "APPROXIMATE migration at level 0"),
    # A version 1 checkpoint written by a WRFInput run
    Case("v1_wrfinput", 1, "WRFInput", "wrfinput", [], "input",
         "pressure recovery at level 0 with the PSFC of wrfinput_d01"),
    # A small level-zero checkpoint restarted with a larger WRFInput file must recover
    # only the PSFC region mapped to its local grids, with spatially varying values.
    Case("cropped_markerless_wrfinput", None, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_crop"], "input",
         "pressure recovery at level 0 with the PSFC of wrfinput_crop", base="crop",
         input_file="WRFInput_LegacySurfaceTemperature_Cropped.i"),
    Case("cropped_v1_wrfinput", 1, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_crop"], "input",
         "pressure recovery at level 0 with the PSFC of wrfinput_crop", base="crop",
         input_file="WRFInput_LegacySurfaceTemperature_Cropped.i"),
    # The fallback is decided across all covered valid samples on this level. Uncovered
    # parent-seeded cells are excluded from the decision.
    Case("allzero_psfc_uses_p0", 1, "WRFInput", P_0,
         ["erf.nc_init_file_0=wrfinput_allzero"], "input",
         ["all zero on the covered level-0 domain; replacing its covered cells with p_0",
          "all covered valid PSFC samples read for the level were zero; covered cells were "
          "replaced by p_0"]),
    # One actual checkpoint valid box has zero PSFC while others have positive PSFC; the
    # level-wide fallback must not accept that zero region by itself.
    Case("zero_valid_box_in_mixed_psfc_does_not_use_p0", 1, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_zero_region"], "abort",
         ["surface pressure = 0 Pa", "cells cannot be converted"]),
    Case("nonfinite_psfc_does_not_use_p0", 1, "WRFInput", "wrfinput",
         ["erf.nc_init_file_0=wrfinput_nanpsfc"], "abort",
         ["surface pressure = nan Pa", "cells cannot be converted"]),
    # A markerless checkpoint whose job_info does not record the init type, with and
    # without the assertion
    Case("markerless_asserted", None, None, "wrfinput", [ASSERTION + "=WRFInput"], "input",
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
    Case("mixed_markerless", None, "WRFInput", "wrfinput", [], "input",
         [INPUT_MIXED, kept(LAND_CELLS)], base="mixed"),
    Case("mixed_v1", 1, "WRFInput", "wrfinput", [], "input",
         [INPUT_MIXED, kept(LAND_CELLS)], base="mixed"),
    # A zero SST over water is a placeholder too at a level with TSK
    Case("mixed_v1_water_sst_zero_with_tsk", 1, "WRFInput", "wrfinput", [], "input",
         kept(LAND_CELLS + 1), base="mixed", zero=[("SST_0", WATER_CELL)]),
    # Without TSK, as checkpoints written before TSK_0 was, land takes the prescribed
    # temperature and keeps its SST placeholders, and water takes the SST ...
    Case("mixed_v1_without_tsk", 1, "WRFInput", "wrfinput", [], "input",
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
    Case("mixed_v1_without_lmask", 1, "WRFInput", "wrfinput", [], "input",
         kept(LAND_CELLS), base="mixed", remove=["LMASK_0"]),
    # ... and water with erf.is_land = 0, where, without TSK, every placeholder is taken
    Case("mixed_v1_without_lmask_or_tsk_on_water", 1, "WRFInput", "wrfinput",
         ["erf.is_land=0"], "abort",
         ["a water cell where the surface layer takes the SST, has theta = 0,",
          f"{LAND_CELLS} cells (0 ghosts) could not be converted"],
         base="mixed", remove=["LMASK_0", "TSK_0"]),
    # Version 2, placeholders and all, restored as it is
    Case("mixed_v2", 2, "WRFInput", None, [], "unchanged", "Reading SST data", base="mixed"),

    # Without SST_0 the surface layer takes neither field: a zero TSK is a placeholder, over
    # land and water alike, kept as it is, while every positive TSK is converted ...
    Case("mixed_markerless_without_sst_tsk_zero", None, "WRFInput", "wrfinput",
         SURFACE_LAYER, "input", [INPUT_MIXED, kept(2, "TSK")], base="mixed",
         remove=["SST_0"], zero=TSK_ZEROS),
    Case("mixed_v1_without_sst_tsk_zero", 1, "WRFInput", "wrfinput", SURFACE_LAYER, "input",
         [INPUT_MIXED, kept(2, "TSK")], base="mixed", remove=["SST_0"], zero=TSK_ZEROS),
    Case("mixed_v1_without_sst", 1, "WRFInput", "wrfinput", SURFACE_LAYER, "input",
         INPUT_MIXED, base="mixed", remove=["SST_0"]),
    # ... while with the SST it takes the TSK at every cell, and the same zeros abort
    Case("mixed_markerless_tsk_zero", None, "WRFInput", "wrfinput", SURFACE_LAYER, "abort",
         ["converting the TSK at level 0 to absolute temperature: 2 cells (0 ghosts) could "
          "not be converted", "where the surface layer takes the TSK, has theta = 0,"],
         base="mixed", zero=TSK_ZEROS),
    # The run takes the TSK just so: with the SST it rejects such a TSK ...
    Case("mixed_v2_tsk_zero", 2, "WRFInput", None, SURFACE_LAYER, "abort",
         "SurfaceLayer surface humidity update: surface conversion failed at level 0",
         base="mixed", zero=TSK_ZEROS),
    # ... and without it runs, restoring the TSK as it is
    Case("mixed_v2_without_sst_tsk_zero", 2, "WRFInput", None, SURFACE_LAYER, "unchanged",
         "Reading TSK data", base="mixed", remove=["SST_0"], zero=TSK_ZEROS),
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


def zero_cell(checkpoint, name, cell, level=0):
    """Set a planar field of a checkpoint to zero at cell (i, j) of a level, ghost cells
    included"""
    i, j = cell
    found = False
    for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
        lo, hi = fab["lo"], fab["hi"]
        if lo[0] <= i <= hi[0] and lo[1] <= j <= hi[1]:
            data = fab_data(fab).copy()
            data[:, :, j - lo[1], i - lo[0]] = 0.0
            write_fab_data(fab, data)
            found = True
    if not found:
        raise AssertionError(f"{checkpoint}/{name} does not hold the cell {cell}")


def remove_field(checkpoint, name, level=0):
    """Remove a MultiFab of a level from a checkpoint, as if it had not been written"""
    level_dir = checkpoint / f"Level_{level}"
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


def check_restart(name, status, log, messages, abort=False):
    """That a restart aborted, or ran, saying each of messages"""
    for message in messages:
        if message not in log:
            raise AssertionError(f"{name}: the restart output ({name}.log) does not say "
                                 f"\"{message}\"")
    if abort and status == 0:
        raise AssertionError(f"{name}: the restart did not abort; see {name}.log")
    if not abort and status != 0:
        raise AssertionError(f"{name}: the restart failed ({status}); see {name}.log")


def check_version_2(name, checkpoint):
    """That a restart wrote its checkpoint as version 2, absolute temperature"""
    marker = checkpoint / "surface_temperature_contract"
    if not marker.exists() or marker.read_text().split() != ["2"]:
        raise AssertionError(f"{name}: {checkpoint.name} is not a version 2 checkpoint")


def edited_copy(checkpoint, name, remove=(), zero=()):
    """A copy of a checkpoint without the fields (level, field) of remove, and with the
    cells (level, field, (i, j)) of zero set to zero"""
    copy = checkpoint.parent / ("legacy_" + name)
    shutil.rmtree(copy, ignore_errors=True)
    shutil.copytree(checkpoint, copy)
    for level, field in remove:
        remove_field(copy, field, level)
    for level, field, cell in zero:
        zero_cell(copy, field, cell, level)
    return copy


def run_case(case, exe, run_dir, bases, env):
    wrfinput = BASES[case.base][1]
    original, psfc = bases[case.base]
    legacy = run_dir / ("legacy_" + case.name)
    output = run_dir / (case.name + "_chk00001")
    log_path = run_dir / (case.name + ".log")
    for path in (legacy, output):
        shutil.rmtree(path, ignore_errors=True)
    shutil.copytree(original, legacy)

    if case.base == "crop" and len(valid_boxes(original, "SST_0", 0)) < 4:
        raise AssertionError(f"{case.name}: the cropped fixture does not have four disjoint "
                             "local regions for the PSFC read")
    if case.name == "zero_valid_box_in_mixed_psfc_does_not_use_p0":
        with nc.Dataset(run_dir / "wrfinput_zero_region") as ds:
            fixture_psfc = np.asarray(ds["PSFC"][0], dtype=np.float64)
        assert_zero_and_positive_psfc_boxes(original, fixture_psfc)

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
    status, log_text = restart(exe, run_dir, case.input_file,
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
    check_version_2(case.name, output)

    report = []
    for name in SURFACE_FIELDS:
        if stored[name] is None:
            if (output / "Level_0" / (name + "_H")).exists():
                raise AssertionError(f"{case.name}: the restart wrote the {name} removed")
            report.append(f"{name} absent")
            continue
        actual = surface_field(output, name, psfc.shape)
        if case.expect == "unchanged":
            # As the restart read it
            if not np.array_equal(actual, stored[name], equal_nan=True):
                raise AssertionError(f"{case.name}: {name} changed on restart")
            report.append(f"{name} bit for bit")
            continue
        expected = surface_field(original, name, psfc.shape)
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


def nest_fields(checkpoint, names=SURFACE_FIELDS):
    """The fields names (SST_0 and TSK_0) of a checkpoint of nest_files.i at each level"""
    return {(level, name): surface_field(checkpoint, name, shape, level, origin)
            for level, _, origin, shape in NEST_LEVELS for name in names}


def check_no_sst(name, checkpoint):
    """That a checkpoint of nest_files.i holds no SST at any level"""
    for level, *_ in NEST_LEVELS:
        if (checkpoint / f"Level_{level}" / "SST_0_H").exists():
            raise AssertionError(f"{name}: the restart wrote an SST at level {level}")


def legacy_nest_copy(checkpoint, name):
    """A version 1 copy, legacy_<name>, of a checkpoint of nest_files.i, whose SST and TSK
    are potential temperature converted with the PSFC of the file of each level"""
    copy = checkpoint.parent / ("legacy_" + name)
    shutil.rmtree(copy, ignore_errors=True)
    shutil.copytree(checkpoint, copy)
    (copy / "surface_temperature_contract").write_text("1\n")
    for level, wrfinput, origin, _ in NEST_LEVELS:
        with nc.Dataset(checkpoint.parent / wrfinput) as ds:
            psfc = np.asarray(ds["PSFC"][0], dtype=np.float64)
        to_potential_temperature(copy, psfc, level, origin)
    return copy


def check_nest_recovered(name, checkpoint, expected):
    """That the known no-regrid fixture matches its source values numerically; this is not
    a general guarantee of exact fine-level recovery"""
    check_version_2(name, checkpoint)
    actual = nest_fields(checkpoint)
    error = max(np.abs(actual[key] / expected[key] - 1.0).max() for key in expected)
    if error > 1.0e-10:
        raise AssertionError(f"{name}: the SST or TSK is off by a relative {error:.3e}")
    return actual, error


def valid_boxes(checkpoint, name, level):
    """The valid boxes of a checkpoint MultiFab at a level"""
    boxes = set()
    for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
        ng = fab["nghost"]
        boxes.add(tuple(fab["lo"][d] + ng[d] for d in range(2)) +
                  tuple(fab["hi"][d] - ng[d] for d in range(2)))
    return boxes


def assert_zero_and_positive_psfc_boxes(checkpoint, psfc):
    """Assert that mixed-zero PSFC contains both zero and positive actual valid boxes."""
    boxes = valid_boxes(checkpoint, "SST_0", 0)
    zero_boxes = []
    positive_boxes = []
    for ilo, jlo, ihi, jhi in boxes:
        values = psfc[jlo:jhi + 1, ilo:ihi + 1]
        if values.size and np.all(values == 0.0):
            zero_boxes.append((ilo, jlo, ihi, jhi))
        if values.size and np.all(np.isfinite(values) & (values > 0.0)):
            positive_boxes.append((ilo, jlo, ihi, jhi))
    if not zero_boxes or not positive_boxes:
        raise AssertionError("mixed-zero PSFC fixture must contain an all-zero valid box "
                             "and a separate all-positive valid box; actual boxes are "
                             f"{sorted(boxes)}")


def valid_boxes_3d(checkpoint, name, level):
    """The valid 3D boxes of a checkpoint MultiFab at a level"""
    boxes = set()
    for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
        ng = fab["nghost"]
        boxes.add(tuple(fab["lo"][d] + ng[d] for d in range(3)) +
                  tuple(fab["hi"][d] - ng[d] for d in range(3)))
    return boxes


def valid_cells(checkpoint, name, level):
    """A planar checkpoint MultiFab as {(i, j): value} over its valid cells"""
    result = {}
    for fab in multifab_fabs(checkpoint / f"Level_{level}", name):
        data = fab_data(fab)[0, 0]
        ng = fab["nghost"]
        ilo, jlo = fab["lo"][0] + ng[0], fab["lo"][1] + ng[1]
        ihi, jhi = fab["hi"][0] - ng[0], fab["hi"][1] - ng[1]
        for j in range(jlo, jhi + 1):
            for i in range(ilo, ihi + 1):
                result[(i, j)] = data[j - fab["lo"][1], i - fab["lo"][0]]
    return result


def conserved_cells(checkpoint, level):
    """The valid checkpoint Cell state as {(i,j,k): component values}"""
    result = {}
    for fab in multifab_fabs(checkpoint / f"Level_{level}", "Cell"):
        data = fab_data(fab)
        ng = fab["nghost"]
        ilo, jlo, klo = (fab["lo"][d] + ng[d] for d in range(3))
        ihi, jhi, khi = (fab["hi"][d] - ng[d] for d in range(3))
        for k in range(klo, khi + 1):
            for j in range(jlo, jhi + 1):
                for i in range(ilo, ihi + 1):
                    result[(i, j, k)] = data[:, k - fab["lo"][2],
                                             j - fab["lo"][1], i - fab["lo"][0]]
    return result


def diagnosed_ground_pressure(checkpoint):
    """Flat-deck ground pressure diagnosed from level-0 rho, rho-theta and cell height"""
    state = conserved_cells(checkpoint, 0)
    return { (i, j): P_0 * (R_D * values[1] / P_0) ** GAMMA +
                     values[0] * GRAV * 50.0
             for (i, j, k), values in state.items() if k == 0 }


def elevated_surface_pressure_at_cell(rho_theta, rho):
    """Ground-face pressure proxy at the lowest cell of an aloft 100 m layer"""
    return P_0 * (R_D * rho_theta / P_0) ** GAMMA + rho * GRAV * 50.0


def set_aloft_checkpoint_pressure(checkpoint, target_cell_pressure):
    """Set only the lowest valid rho-theta plane of aloft level-1 boxes to a target pressure"""
    target_rho_theta = (P_0 / R_D) * (target_cell_pressure / P_0) ** (1.0 / GAMMA)
    for fab in multifab_fabs(checkpoint / "Level_1", "Cell"):
        ng = fab["nghost"]
        klo = fab["lo"][2] + ng[2]
        if klo <= 0:
            raise AssertionError("the Metgrid pressure fixture's level 1 reaches the ground")
        data = fab_data(fab).copy()
        for j in range(fab["lo"][1] + ng[1], fab["hi"][1] - ng[1] + 1):
            for i in range(fab["lo"][0] + ng[0], fab["hi"][0] - ng[0] + 1):
                data[1, klo - fab["lo"][2], j - fab["lo"][1],
                     i - fab["lo"][0]] = target_rho_theta
        write_fab_data(fab, data)
    return target_rho_theta


def set_valid_ground_component(checkpoint, level, comp, value):
    """Set one conserved component on every valid ground-plane cell of a level"""
    found = False
    for fab in multifab_fabs(checkpoint / f"Level_{level}", "Cell"):
        ng = fab["nghost"]
        klo = fab["lo"][2] + ng[2]
        if klo != 0:
            continue
        data = fab_data(fab).copy()
        ilo, jlo = fab["lo"][0] + ng[0], fab["lo"][1] + ng[1]
        ihi, jhi = fab["hi"][0] - ng[0], fab["hi"][1] - ng[1]
        data[comp, klo - fab["lo"][2], jlo - fab["lo"][1]:jhi - fab["lo"][1] + 1,
             ilo - fab["lo"][0]:ihi - fab["lo"][0] + 1] = value
        write_fab_data(fab, data)
        found = True
    if not found:
        raise AssertionError(f"{checkpoint}: level {level} has no valid ground cells")


def cell_cons_prolong(coarse, ratio=3):
    """AMReX's limited cell-consistent linear interpolation for a smooth 2D field"""
    ny, nx = coarse.shape
    padded = np.pad(coarse, 1, mode="edge")
    sx = np.zeros_like(coarse)
    sy = np.zeros_like(coarse)
    for j in range(ny):
        for i in range(nx):
            center = padded[j + 1, i + 1]
            left, right = padded[j + 1, i], padded[j + 1, i + 2]
            down, up = padded[j, i + 1], padded[j + 2, i + 1]
            dcx = 0.5 * (right - left)
            dfx, dbx = 2.0 * (right - center), 2.0 * (center - left)
            dcy = 0.5 * (up - down)
            dfy, dby = 2.0 * (up - center), 2.0 * (center - down)
            if dfx * dbx >= 0.0:
                sx[j, i] = np.copysign(min(abs(dfx), abs(dbx), abs(dcx)), dcx)
            if dfy * dby >= 0.0:
                sy[j, i] = np.copysign(min(abs(dfy), abs(dby), abs(dcy)), dcy)
            neighborhood = padded[j:j + 3, i:i + 3]
            dmax = (abs(sx[j, i]) + abs(sy[j, i])) * (ratio - 1) / (2.0 * ratio)
            if dmax > 0.0:
                alpha = min(1.0, (neighborhood.max() - center) / dmax,
                            (center - neighborhood.min()) / dmax)
                sx[j, i] *= alpha
                sy[j, i] *= alpha

    fine = np.empty((ny * ratio, nx * ratio), dtype=np.float64)
    for j in range(fine.shape[0]):
        cj, rj = divmod(j, ratio)
        yoff = (rj + 0.5) / ratio - 0.5
        for i in range(fine.shape[1]):
            ci, ri = divmod(i, ratio)
            xoff = (ri + 0.5) / ratio - 0.5
            fine[j, i] = coarse[cj, ci] + xoff * sx[cj, ci] + yoff * sy[cj, ci]
    return fine


def make_regridded_legacy_copy(checkpoint, name, parent_pressure):
    """Encode a post-RemakeLevel checkpoint as old theta, retaining its regridded layout"""
    copy = checkpoint.parent / ("legacy_" + name)
    shutil.rmtree(copy, ignore_errors=True)
    shutil.copytree(checkpoint, copy)
    (copy / "surface_temperature_contract").write_text("1\n")
    edit_job_info(copy, "WRFInput")
    parent_pressure_fine = cell_cons_prolong(parent_pressure)

    for field in SURFACE_FIELDS:
        parent_temperature = surface_field(checkpoint, field, parent_pressure.shape)
        parent_theta = parent_temperature * (P_0 / parent_pressure) ** KAPPA
        fine_theta = cell_cons_prolong(parent_theta)
        for fab in multifab_fabs(copy / "Level_1", field):
            data = fab_data(fab).copy()
            ii = np.arange(fab["lo"][0], fab["hi"][0] + 1)
            jj = np.arange(fab["lo"][1], fab["hi"][1] + 1)
            fi = np.clip(ii, 0, fine_theta.shape[1] - 1)
            fj = np.clip(jj, 0, fine_theta.shape[0] - 1)
            data[0, 0] = fine_theta[np.ix_(fj, fi)]
            write_fab_data(fab, data)

        # Independently read the fabricated checkpoint field back and compare each valid
        # fine index with a fresh parent-theta prolongation. This catches index scaling errors
        # in the fixture itself before the migration test can accidentally validate them.
        expected_parent_theta = surface_field(checkpoint, field, parent_pressure.shape) * \
            (P_0 / parent_pressure) ** KAPPA
        expected_fine_theta = cell_cons_prolong(expected_parent_theta)
        fabricated_theta = valid_cells(copy, field, 1)
        relative_errors = []
        for (i, j), value in fabricated_theta.items():
            if not (0 <= i < expected_fine_theta.shape[1] and
                    0 <= j < expected_fine_theta.shape[0]):
                raise AssertionError(f"{copy.name}: fine valid cell {(i, j)} is outside the "
                                     "parent-prolonged theta field")
            expected = expected_fine_theta[j, i]
            relative_errors.append(abs(value / expected - 1.0))
        if not relative_errors or max(relative_errors) > 2.0e-12:
            error = max(relative_errors) if relative_errors else float("inf")
            raise AssertionError(f"{copy.name}: fabricated {field} is not the independently "
                                 f"prolonged parent theta field (relative error {error:.3e})")

    return copy, parent_pressure_fine


def run_nest_regrid_pressure_case(exe, run_dir, env):
    """Exercise legacy migration after RemakeLevel changed a WRFInput fine footprint"""
    checkpoint = run_dir / "nest_regrid_chk00002"
    shutil.rmtree(checkpoint, ignore_errors=True)
    log_path = run_dir / "nest_regrid.log"
    status, _ = restart(
        exe, run_dir, "nest_files.i", "nest_chk00000", "nest_regrid",
        ["erf.nc_init_file_1=nest_d02_full", "erf.regrid_int=1",
         "erf.refinement_indicators=box1", "erf.box1.max_level=1",
         "erf.box1.start_time=0", "erf.box1.in_box_lo=0 0 0",
         "erf.box1.in_box_hi=8000 8000 800"], env, max_step=2)
    if status != 0:
        raise AssertionError(f"nest_regrid: run failed ({status}); see {log_path}")
    if not checkpoint.exists():
        raise AssertionError("nest_regrid: the post-regrid checkpoint was not written")
    if valid_boxes(run_dir / "nest_chk00000", "SST_0", 1) == \
            valid_boxes(checkpoint, "SST_0", 1):
        raise AssertionError("nest_regrid: level 1 did not change grids through RemakeLevel")

    with nc.Dataset(run_dir / "nest_d01") as ds:
        parent_pressure = np.asarray(ds["PSFC"][0], dtype=np.float64)
    with nc.Dataset(run_dir / "nest_d02") as ds:
        fine_pressure = np.asarray(ds["PSFC"][0], dtype=np.float64)
    legacy, parent_pressure_fine = make_regridded_legacy_copy(
        checkpoint, "nest_regrid_v1", parent_pressure)
    origin = NEST_LEVELS[1][2]

    name = "nest_regrid_v1_restart"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name,
                          ["erf.nc_init_file_1=nest_d02"], env, max_step=3)
    check_restart(name, status, log,
                  ["APPROXIMATE migration: the PSFC of nest_d02", "do not commute"])
    output = run_dir / (name + "_chk00003")
    check_version_2(name, output)

    differences = []
    for field in SURFACE_FIELDS:
        actual = valid_cells(output, field, 1)
        parent_temp = surface_field(checkpoint, field, parent_pressure.shape)
        parent_theta = parent_temp * (P_0 / parent_pressure) ** KAPPA
        fine_theta = cell_cons_prolong(parent_theta)
        fine_file_on_erf = np.full_like(parent_pressure_fine, np.nan)
        # The level-1 WRFInput file begins at the nonzero nest origin. The migration overlays
        # it there and uses recovered parent pressure where the regridded patch extends beyond.
        for j in range(fine_pressure.shape[0]):
            for i in range(fine_pressure.shape[1]):
                fine_file_on_erf[origin[1] + j, origin[0] + i] = fine_pressure[j, i]
        errors = []
        policy_delta = []
        stored_theta = valid_cells(legacy, field, 1)
        for (i, j), value in actual.items():
            if not (0 <= i < fine_theta.shape[1] and 0 <= j < fine_theta.shape[0]):
                raise AssertionError(f"{name}: fine valid cell {(i, j)} is outside the parent")
            p_file = fine_file_on_erf[j, i]
            p_parent = parent_pressure_fine[j, i]
            p_used = p_file if np.isfinite(p_file) else p_parent
            theta = stored_theta[(i, j)]
            expected = theta * (p_used / P_0) ** KAPPA
            parent_recovery = theta * (p_parent / P_0) ** KAPPA
            errors.append(abs(value / expected - 1.0))
            if np.isfinite(p_file):
                policy_delta.append(abs(expected - parent_recovery))
        error = max(errors)
        if error > 2.0e-10:
            raise AssertionError(f"{name}: {field} is off the selected pressure policy by "
                                 f"a relative {error:.3e}")
        if not policy_delta or max(policy_delta) < 0.5:
            raise AssertionError(f"{name}: the fine and parent PSFC produced no measurable "
                                 "regrid migration difference")
        differences.append(f"{field} policy error {error:.1e}, fine/parent pressure "
                           f"difference up to {max(policy_delta):.2f} K")
    return ", ".join(differences)


def run_nest_aloft_metgrid_case(exe, run_dir, env):
    """Regrid a WRFInput checkpoint aloft, migrate Metgrid-v1 arrays, then reach the ground"""
    checkpoint = run_dir / "nest_aloft_chk00002"
    shutil.rmtree(checkpoint, ignore_errors=True)
    log_path = run_dir / "nest_aloft.log"
    status, _ = restart(
        exe, run_dir, "nest_files.i", "nest_chk00000", "nest_aloft",
        ["erf.init_type=Uniform", "erf.check_file=nest_aloft_chk", "erf.check_int=2",
         "erf.regrid_int=1",
         "amr.refine_whole_domain_dir=0", "erf.refinement_indicators=box1",
         "erf.box1.max_level=1", "erf.box1.start_time=1",
         "erf.box1.in_box_lo=3000 3000 400", "erf.box1.in_box_hi=5000 5000 800"],
        env, max_step=2)
    if status != 0:
        raise AssertionError(f"nest_aloft: run failed ({status}); see {log_path}")
    if not checkpoint.exists():
        raise AssertionError("nest_aloft: the checkpoint was not written")
    aloft_boxes = valid_boxes_3d(checkpoint, "Cell", 1)
    if not aloft_boxes or min(box[2] for box in aloft_boxes) <= 0:
        raise AssertionError("nest_aloft: the regridded level-1 boxes reach the ground")

    legacy = run_dir / "legacy_nest_aloft_metgrid_v1"
    shutil.rmtree(legacy, ignore_errors=True)
    shutil.copytree(checkpoint, legacy)
    (legacy / "surface_temperature_contract").write_text("1\n")
    edit_job_info(legacy, "Metgrid")
    to_potential_temperature(legacy, 50000.0, level=0)
    to_potential_temperature(legacy, 50000.0, level=1)

    # Keep the synthetic elevated state physically warm enough for the one-step restart,
    # while remaining at least 10 kPa away from the true ground pressure.
    target_cell_pressure = 75000.0
    target_rho_theta = set_aloft_checkpoint_pressure(legacy, target_cell_pressure)
    legacy_state = conserved_cells(legacy, 1)
    first_aloft = next((key, values) for key, values in legacy_state.items()
                       if key[2] == min(box[2] for box in aloft_boxes))
    elevated_pressure = elevated_surface_pressure_at_cell(
        target_rho_theta, first_aloft[1][0])
    ground_pressure = diagnosed_ground_pressure(checkpoint)
    if not ground_pressure:
        raise AssertionError("nest_aloft: no level-0 ground pressure could be diagnosed")
    ground_values = np.asarray(list(ground_pressure.values()))
    if np.min(np.abs(ground_values - elevated_pressure)) < 10000.0:
        raise AssertionError("nest_aloft: elevated pressure is not substantially different "
                             "from ground pressure")

    # The parent has no surface arrays, but its ground pressure must still be prepared for
    # this aloft child. Also make the failure mode explicit if that parent pressure itself
    # cannot be diagnosed.
    parent_missing = run_dir / "legacy_nest_aloft_metgrid_parent_without_surface"
    shutil.rmtree(parent_missing, ignore_errors=True)
    shutil.copytree(legacy, parent_missing)
    for field in SURFACE_FIELDS:
        remove_field(parent_missing, field, level=0)
    parent_name = "nest_aloft_metgrid_parent_without_surface_restart"
    status, log = restart(exe, run_dir, "nest_files.i", parent_missing.name, parent_name,
                          [], env, max_step=3)
    check_restart(parent_name, status, log,
                  ["APPROXIMATE migration at level 1", "coarser level's ground pressure"])
    parent_output = run_dir / (parent_name + "_chk00003")
    check_version_2(parent_name, parent_output)
    parent_pressure = np.full((16, 16), np.nan, dtype=np.float64)
    for (i, j), value in ground_pressure.items():
        parent_pressure[j, i] = value
    parent_pressure_fine = cell_cons_prolong(parent_pressure)
    parent_errors = []
    for field in SURFACE_FIELDS:
        stored_theta = valid_cells(parent_missing, field, 1)
        actual = valid_cells(parent_output, field, 1)
        if actual.keys() != stored_theta.keys():
            raise AssertionError(f"{parent_name}: {field} changed layout")
        parent_errors.extend(
            abs(value / (stored_theta[(i, j)] * (parent_pressure_fine[j, i] / P_0) ** KAPPA)
                - 1.0)
            for (i, j), value in actual.items())
    parent_error = max(parent_errors) if parent_errors else float("inf")
    if parent_error > 2.0e-10:
        raise AssertionError(f"{parent_name}: child did not use prepared parent ground pressure "
                             f"(relative error {parent_error:.3e})")

    invalid_parent = run_dir / "legacy_nest_aloft_metgrid_parent_bad_pressure"
    shutil.rmtree(invalid_parent, ignore_errors=True)
    shutil.copytree(parent_missing, invalid_parent)
    set_valid_ground_component(invalid_parent, 0, RHO_COMP, 0.0)
    invalid_name = "nest_aloft_metgrid_parent_bad_pressure_restart"
    status, log = restart(exe, run_dir, "nest_files.i", invalid_parent.name, invalid_name,
                          [], env, max_step=3)
    check_restart(invalid_name, status, log,
                  ["legacy checkpoint surface-pressure recovery",
                   "density of the state cell is not usable"], abort=True)

    name = "nest_aloft_metgrid_v1_restart"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name, [], env,
                          max_step=3)
    check_restart(name, status, log,
                  ["APPROXIMATE migration at level 1", "coarser level's ground pressure"])
    migrated = run_dir / (name + "_chk00003")
    check_version_2(name, migrated)

    pressure_parent = np.full((16, 16), np.nan, dtype=np.float64)
    for (i, j), value in ground_pressure.items():
        pressure_parent[j, i] = value
    if np.isnan(pressure_parent).any():
        raise AssertionError("nest_aloft: level 0 does not cover its full ground plane")
    pressure_fine = cell_cons_prolong(pressure_parent)
    migration_errors = []
    elevated_differences = []
    for field in SURFACE_FIELDS:
        theta = valid_cells(legacy, field, 1)
        actual = valid_cells(migrated, field, 1)
        if actual.keys() != theta.keys():
            raise AssertionError(f"{name}: {field} changed layout before the regrid")
        for (i, j), value in actual.items():
            expected = theta[(i, j)] * (pressure_fine[j, i] / P_0) ** KAPPA
            wrong = theta[(i, j)] * (elevated_pressure / P_0) ** KAPPA
            migration_errors.append(abs(value / expected - 1.0))
            elevated_differences.append(abs(value - wrong))
    migration_error = max(migration_errors)
    if migration_error > 2.0e-10:
        raise AssertionError(f"{name}: migration did not use parent ground pressure "
                             f"(relative error {migration_error:.3e})")
    if max(elevated_differences) < 5.0:
        raise AssertionError(f"{name}: elevated stand-in pressure caused no measurable "
                             "temperature difference")

    # A later regrid expands level 1 down to the ground. RemakeLevel interpolates already
    # migrated absolute temperatures from the parent; those values must remain kelvin.
    next_name = "nest_aloft_metgrid_regridded"
    status, log = restart(
        exe, run_dir, "nest_files.i", migrated.name, next_name,
        ["erf.init_type=Uniform", "erf.regrid_int=1", "amr.refine_whole_domain_dir=2",
         "erf.refinement_indicators=box1", "erf.box1.max_level=1",
         "erf.box1.start_time=0", "erf.box1.in_box_lo=3000 3000 0",
         "erf.box1.in_box_hi=5000 5000 800"], env, max_step=4)
    if status != 0:
        raise AssertionError(f"{next_name}: restart failed ({status}); see {next_name}.log")
    regridded = run_dir / (next_name + "_chk00004")
    boxes = valid_boxes_3d(regridded, "Cell", 1)
    if not boxes or min(box[2] for box in boxes) != 0:
        raise AssertionError(f"{next_name}: level-1 regrid did not bring coverage to ground")
    for field in SURFACE_FIELDS:
        parent_temperature = surface_field(regridded, field, (16, 16), level=0)
        expected = cell_cons_prolong(parent_temperature)
        actual = valid_cells(regridded, field, 1)
        errors = [abs(value - expected[j, i]) for (i, j), value in actual.items()]
        if not errors or max(errors) > 2.0e-10:
            raise AssertionError(f"{next_name}: {field} is not kelvin after ground regrid")

    return (f"parent-pressure migration error {migration_error:.1e}; elevated-pressure bias "
            f"up to {max(elevated_differences):.1f} K; regridded ground coverage in kelvin")


def started_message(file, start, checkpoint_time):
    """What a restart says when an init file of level 1 after its first had started"""
    return ("at level 1 cannot be recovered: the level lists 2 init files "
            f"(erf.nc_init_file_1), and {file}, after the first, had started by the time "
            f"of the checkpoint (it starts at t = {start} s at level 0, and the checkpoint "
            f"is at t = {checkpoint_time} s): the run that wrote it may have read that file "
            "as well")


def run_nest_cases(exe, run_dir, env):
    """The nest cases (see above): the name and the result of each"""
    original = run_dir / "nest_chk00000"
    if (original / "surface_temperature_contract").read_text().split() != ["2"]:
        raise AssertionError("nest_chk00000 is not a version 2 checkpoint")
    expected = nest_fields(original)
    legacy = legacy_nest_copy(original, "nest_v1")
    results = []
    from_first = ("pressure recovery at level 1 with the PSFC of nest_d02, the first of its 2 "
                  "init files: none of the others had started by the time of the checkpoint")

    # The deck as it was: nest_d03 had not started, so level 1 had read nest_d02 alone
    name = "nest_v1_two_level1_init_files"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name, [], env)
    check_restart(name, status, log,
                  ["pressure recovery at level 0 with the PSFC of nest_d01", from_first,
                   "APPROXIMATE migration: the PSFC of nest_d02"])
    migrated = run_dir / (name + "_chk00001")
    actual, error = check_nest_recovered(name, migrated, expected)
    results.append((name, f"SST_0 and TSK_0 at both levels to {error:.1e}"))

    # The same as with nest_d02 alone
    name = "nest_v1_first_level1_init_file"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name,
                          ["erf.nc_init_file_1=nest_d02"], env)
    check_restart(name, status, log,
                  ["pressure recovery at level 1 with the PSFC of nest_d02\n"])
    alone, error = check_nest_recovered(name, run_dir / (name + "_chk00001"), expected)
    if any(not np.array_equal(alone[key], actual[key]) for key in actual):
        raise AssertionError(f"{name}: the SST or TSK differs from {migrated.name}'s")
    results.append((name, f"SST_0 and TSK_0 at both levels to {error:.1e}, the same as "
                          "with both files, bit for bit"))

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

    # A file listed after nest_d02 that starts at the start may have been read as well
    name = "nest_v1_level1_init_file_started"
    status, log = restart(exe, run_dir, "nest_files.i", legacy.name, name,
                          ["erf.nc_init_file_1=nest_d02 nest_d03_now"], env)
    check_restart(name, status, log, [started_message("nest_d03_now", 0, 0)], abort=True)
    results.append((name, "aborted as it must"))

    # The deck run on to t = 10 s, when nest_d03 starts, with no regridding to read it: its
    # checkpoint at t = 5 s selects the first file. This run has no regrid, so the
    # expected field remains the known no-regrid fixture value ...
    status, _ = restart(exe, run_dir, "nest_files.i", original.name, "nest_later",
                        ["erf.check_int=5"], env, max_step=10)
    if status != 0:
        raise AssertionError(f"nest_later: the run to t = 10 s failed ({status}); "
                             "see nest_later.log")
    name = "nest_v1_before_level1_init_file_starts"
    written = run_dir / "nest_later_chk00005"
    status, log = restart(exe, run_dir, "nest_files.i", legacy_nest_copy(written, name).name,
                          name, [], env, max_step=6)
    check_restart(name, status, log, [from_first])
    _, error = check_nest_recovered(name, run_dir / (name + "_chk00006"),
                                    nest_fields(written))
    results.append((name, f"SST_0 and TSK_0 at both levels to {error:.1e}"))

    # ... and the one at t = 10 s not, as the run may have read nest_d03 by then
    name = "nest_v1_once_level1_init_file_starts"
    written = run_dir / "nest_later_chk00010"
    status, log = restart(exe, run_dir, "nest_files.i", legacy_nest_copy(written, name).name,
                          name, [], env, max_step=11)
    check_restart(name, status, log, [started_message("nest_d03", 10, 10)], abort=True)
    results.append((name, "aborted as it must"))

    # Level 1 without an SST of its own is given level 0's on restart, so the surface layer
    # takes its TSK, and a zero there aborts the restart ...
    name = "nest_v1_level1_without_sst_tsk_zero"
    edited = edited_copy(legacy, name, remove=[(1, "SST_0")], zero=[(1, "TSK_0", NEST_CELL)])
    status, log = restart(exe, run_dir, "nest_files.i", edited.name, name,
                          ["erf.nc_init_file_1=nest_d02"] + SURFACE_LAYER, env)
    check_restart(name, status, log,
                  ["converting the TSK at level 1 to absolute temperature: 1 cells "
                   "(0 ghosts) could not be converted",
                   f"({NEST_CELL[0]},{NEST_CELL[1]}), a water cell where the surface layer "
                   "takes the TSK, has theta = 0,"], abort=True)
    results.append((name, "aborted as it must"))

    # ... as the run from a version 2 checkpoint does
    name = "nest_v2_level1_without_sst_tsk_zero"
    edited = edited_copy(original, name, remove=[(1, "SST_0")],
                         zero=[(1, "TSK_0", NEST_CELL)])
    status, log = restart(exe, run_dir, "nest_files.i", edited.name, name, SURFACE_LAYER,
                          env)
    check_restart(name, status, log,
                  ["SurfaceLayer surface humidity update: surface conversion failed at "
                   "level 1"], abort=True)
    results.append((name, "the run aborted as it must"))

    # With no SST at any level, the zero is a placeholder, kept as it is, and every positive
    # TSK is converted ...
    tsk = ("TSK_0",)
    name = "nest_v1_without_sst_tsk_zero"
    edited = edited_copy(legacy, name, remove=[(0, "SST_0"), (1, "SST_0")],
                         zero=[(1, "TSK_0", NEST_CELL)])
    stored = nest_fields(edited, tsk)
    status, log = restart(exe, run_dir, "nest_files.i", edited.name, name,
                          ["erf.nc_init_file_1=nest_d02"] + SURFACE_LAYER, env)
    check_restart(name, status, log,
                  ["pressure recovery at level 0 with the PSFC of nest_d01",
                   "pressure recovery at level 1 with the PSFC of nest_d02",
                   kept(1, "TSK", 1)])
    output = run_dir / (name + "_chk00001")
    check_version_2(name, output)
    check_no_sst(name, output)
    actual = nest_fields(output, tsk)
    placeholders = 0
    error = 0.0
    for key, value in actual.items():
        placeholder = ~(stored[key] > 0.0)
        if not np.array_equal(value[placeholder], stored[key][placeholder]):
            raise AssertionError(f"{name}: the TSK placeholder at level {key[0]} changed")
        placeholders += int(placeholder.sum())
        relative = value[~placeholder] / expected[key][~placeholder] - 1.0
        error = max(error, np.abs(relative).max())
    if placeholders != 1:
        raise AssertionError(f"{name}: {placeholders} TSK placeholders, not 1")
    if error > 1.0e-10:
        raise AssertionError(f"{name}: the TSK is off by a relative {error:.3e}")
    results.append((name, f"TSK_0 at both levels to {error:.1e}, 1 placeholder kept, "
                          "no SST"))

    # ... and the run from a version 2 checkpoint restores the TSK as it is
    name = "nest_v2_without_sst_tsk_zero"
    edited = edited_copy(original, name, remove=[(0, "SST_0"), (1, "SST_0")],
                         zero=[(1, "TSK_0", NEST_CELL)])
    stored = nest_fields(edited, tsk)
    status, log = restart(exe, run_dir, "nest_files.i", edited.name, name, SURFACE_LAYER,
                          env)
    check_restart(name, status, log, ["Reading TSK data"])
    output = run_dir / (name + "_chk00001")
    check_version_2(name, output)
    check_no_sst(name, output)
    actual = nest_fields(output, tsk)
    if any(not np.array_equal(actual[key], stored[key]) for key in stored):
        raise AssertionError(f"{name}: the TSK changed on restart")
    results.append((name, "TSK_0 at both levels bit for bit, no SST"))
    return results


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exe", type=Path)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()

    bases = {}
    for base, (checkpoint, wrfinput, origin, shape) in BASES.items():
        original = run_dir / checkpoint
        if (original / "surface_temperature_contract").read_text().split() != ["2"]:
            raise AssertionError(f"{checkpoint} is not a version 2 checkpoint")
        with nc.Dataset(run_dir / wrfinput) as ds:
            psfc_full = np.asarray(ds["PSFC"][0], dtype=np.float64)
        i0, j0 = origin
        nx, ny = shape
        psfc = psfc_full[j0:j0 + ny, i0:i0 + nx]
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
    try:
        print(f"nest_regrid_v1_restart: {run_nest_regrid_pressure_case(args.exe, run_dir, env)}")
    except AssertionError as error:
        print(f"FAILED {error}")
        failures += 1
    try:
        print(f"nest_aloft_metgrid_v1_restart: {run_nest_aloft_metgrid_case(args.exe, run_dir, env)}")
    except AssertionError as error:
        print(f"FAILED {error}")
        failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
