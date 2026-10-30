#!/usr/bin/env python3
"""Create the hill case's land setup and stage dependencies in an empty directory.
Requires numpy and netCDF4. Existing run files are never overwritten.
"""
import argparse
from pathlib import Path
import re
import shutil

import numpy as np
from netCDF4 import Dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    case = Path(__file__).resolve().parent
    repo = case.parents[3]
    run = args.directory.resolve()
    run.mkdir(parents=True, exist_ok=True)
    if any(run.iterdir()):
        parser.error(f"working directory must be empty: {run}")
    nx, ny, dx = 64, 32, 100.0
    # The established NoahMP_Ideal fixture supplies all uniform land properties.
    with Dataset(case.parent / "wrfinput_d01") as src, Dataset(run / "wrfinput_d01", "w", format="NETCDF3_CLASSIC") as dst:
        dst.setncatts(src.__dict__)
        dst.setncattr("WEST-EAST_GRID_DIMENSION", nx + 1)
        dst.setncattr("SOUTH-NORTH_GRID_DIMENSION", ny + 1)
        dst.setncattr("DX", np.float32(dx))
        dst.setncattr("DY", np.float32(dx))
        for name, dim in src.dimensions.items():
            dst.createDimension(name, None if dim.isunlimited() else
                                {"west_east": nx, "south_north": ny}.get(name, len(dim)))
        for name, var in src.variables.items():
            new = dst.createVariable(name, var.datatype, var.dimensions)
            new.setncatts(var.__dict__)
            data = var[:]
            if "west_east" in var.dimensions:
                # These fixture fields are spatially uniform, including each soil layer.
                new[:] = np.broadcast_to(data[..., :1, :1], data.shape[:-2] + (ny, nx))
            else:
                new[:] = data
        dst.setncattr("SIMULATION_START_DATE", "2024-08-05_13:00:00")
        dst["Times"][0, :] = np.array(list("2024-08-05_13:00:00"), dtype="S1")
        x, y = np.meshgrid(np.arange(nx + 1)*dx - nx*dx/2,
                           np.arange(ny + 1)*dx - ny*dx/2)
        h = 250.0/(1.0 + (x*x + y*y)/250.0**2)
        dst["HGT"][:] = 0.25*(h[:-1, :-1] + h[1:, :-1] + h[:-1, 1:] + h[1:, 1:])
    namelist = (case.parent / "namelist.erf").read_text()
    for key in ("FORCING_TIMESTEP", "NOAH_TIMESTEP", "OUTPUT_TIMESTEP"):
        namelist = re.sub(rf"({key}\s*=\s*)3600", r"\g<1>1", namelist)
    (run / "namelist.erf").write_text(namelist)
    shutil.copy2(case.parent / "NoahmpTable.TBL", run)
    shutil.copy2(case / "inputs", run)
    # A quiescent, moist, stably stratified atmosphere extending above the model top.
    # The sounding uses potential temperature [K] and mixing ratio [g/kg].
    z = np.arange(0.0, 30001.0, 250.0)
    temp = np.maximum(300.0 - 0.0065*z, 222.0)
    pressure = 1.e5*(temp/300.0)**(9.81/(287.0*0.0065))
    high = z > 12000
    pressure[high] = 1.e5*(222/300)**(9.81/(287*0.0065))*np.exp(-9.81*(z[high]-12000)/(287*222))
    theta = temp*(1.e5/pressure)**(287.0/1004.5)
    with (run / "input_sounding").open("w") as f:
        f.write("1000.0 300.0 6.0\n")
        for zz, th in zip(z, theta):
            f.write(f"{zz:.1f} {th:.8f} {6*np.exp(-zz/2000):.8f} 0.0 0.0\n")
    data_dir = run / "radiation_data"
    data_dir.mkdir()
    rrtmgp = repo / "Submodules/RRTMGP"
    files = list((rrtmgp / "rrtmgp/data").glob("rrtmgp-data-*-2018-12-04.nc"))
    files += list((rrtmgp / "extensions/cloud_optics").glob("rrtmgp-cloud-optics-coeffs-*.nc"))
    required = {"rrtmgp-data-sw-g224-2018-12-04.nc", "rrtmgp-data-lw-g256-2018-12-04.nc",
                "rrtmgp-cloud-optics-coeffs-sw.nc", "rrtmgp-cloud-optics-coeffs-lw.nc"}
    for path in files:
        if path.name in required:
            (data_dir / path.name).symlink_to(path)
            required.remove(path.name)
    if required:
        raise FileNotFoundError(f"RRTMGP submodule data missing: {sorted(required)}")
    print(f"Prepared {run}; run erf_exec inputs in that directory.")


if __name__ == "__main__":
    main()
