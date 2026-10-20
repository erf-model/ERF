#!/usr/bin/env python3
"""Check that ERF's copies of Noah-MP's soil and vegetation tables hold Noah-MP's values.

The two-stream surface energy balance takes soil parameters by Noah-MP soil category
(erf.radiation.seb_soil_type) from Source/Radiation/TwoStream/ERF_NoahMPSoilTable.H, a copy
of the STAS block (&noahmp_soil_stas_parameters) of NoahmpTable.TBL, and canopy parameters
by land-use category (erf.radiation.seb_vegetation_type) from ERF_NoahMPVegetationTable.H,
a copy of the MODIS block (&noahmp_modis_parameters). This compares every copied value,
CSOIL_DATA, RSURF_EXP and Z0SOIL with the table, so a change to the submodule's table cannot leave
the copies behind silently.

Usage: check_noahmp_soil_table.py --table NoahmpTable.TBL --header ERF_NoahMPSoilTable.H
                                  --vegetation-header ERF_NoahMPVegetationTable.H
"""

import argparse
import re
import sys

# Header struct order -> NoahmpTable.TBL names.
COLUMNS = ['BB', 'DRYSMC', 'MAXSMC', 'REFSMC', 'SATPSI', 'SATDK', 'WLTSMC', 'QTZ']


VEGETATION_COLUMNS = ['RS', 'RGL', 'HS', 'TOPT', 'RSMAX', 'Z0MVT'] + [
    'LAI_' + m for m in ('JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
                         'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC')]


def block_values(text, name):
    """{key: [values]} of the namelist group &name."""
    start = text.index('&' + name)
    block = text[start:text.index('\n/', start)]
    values, current = {}, None
    for line in block.splitlines():
        line = line.split('!')[0].rstrip()
        match = re.match(r'\s*([A-Z0-9_]+)\s*=\s*(.*)', line)
        if match:
            current = match.group(1)
            values[current] = match.group(2)
        elif current and line.strip():
            values[current] += ' ' + line.strip()
    parsed = {}
    for k, v in values.items():
        try:
            parsed[k] = [float(x) for x in re.split(r'[,\s]+', v.strip().rstrip(',')) if x]
        except ValueError:
            continue  # a string-valued entry
    return parsed


def scalar(text, name):
    match = re.search(r'^\s*' + name + r'\s*=\s*([0-9.Ee+-]+)', text, re.M)
    return float(match.group(1)) if match else None


def table_values(path):
    text = open(path).read()
    return (block_values(text, 'noahmp_soil_stas_parameters'),
            block_values(text, 'noahmp_modis_parameters'),
            scalar(text, 'CSOIL_DATA'), scalar(text, 'RSURF_EXP'), scalar(text, 'Z0SOIL'))


def table_body(text):
    """The numbers of the static table in a header, in order."""
    start = text.index('table[')
    body = text[start:text.index('};', start)]
    return [float(x) for x in re.findall(r'amrex::Real\(([^)]*)\)', body)]


def header_constant(text, name):
    match = re.search(name + r'\s*=\s*amrex::Real\(([^)]*)\)', text)
    return float(match.group(1)) if match else None


def header_values(path):
    text = open(path).read()
    numbers = table_body(text)
    n = len(COLUMNS)
    rows = [numbers[i:i + n] for i in range(0, len(numbers), n)]
    return rows, header_constant(text, 'noahmp_soil_heat_capacity'), \
        header_constant(text, 'noahmp_soil_resistance_exponent'), \
        header_constant(text, 'noahmp_soil_roughness')


def compare(failures, what, table, rows, columns, count):
    if len(rows) != count:
        failures.append(f"{what}: the header has {len(rows)} categories, not {count}")
    for c, name in enumerate(columns):
        if name not in table:
            failures.append(f"{what}: {name} is missing from the table")
            continue
        if len(table[name]) != len(rows):
            failures.append(f"{what}: {name} has {len(table[name])} values in the table, "
                            f"{len(rows)} in the header")
            continue
        for cat, (t, row) in enumerate(zip(table[name], rows), start=1):
            if abs(row[c] - t) > 1.0e-12 * max(abs(t), 1.0):
                failures.append(f"{what} category {cat} {name}: header {row[c]}, table {t}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--table', required=True)
    parser.add_argument('--header', required=True)
    parser.add_argument('--vegetation-header', required=True)
    args = parser.parse_args()

    soil, modis, table_csoil, table_rsurf_exp, table_z0soil = table_values(args.table)
    rows, header_csoil, header_rsurf_exp, header_z0soil = header_values(args.header)
    numbers = table_body(open(args.vegetation_header).read())
    n = len(VEGETATION_COLUMNS)
    vegetation_rows = [numbers[i:i + n] for i in range(0, len(numbers), n)]

    failures = []
    compare(failures, 'soil', soil, rows, COLUMNS, 19)
    compare(failures, 'vegetation', modis, vegetation_rows, VEGETATION_COLUMNS, 20)
    if table_csoil is None or header_csoil != table_csoil:
        failures.append(f"CSOIL_DATA: header {header_csoil}, table {table_csoil}")
    if table_rsurf_exp is None or header_rsurf_exp != table_rsurf_exp:
        failures.append(f"RSURF_EXP: header {header_rsurf_exp}, table {table_rsurf_exp}")
    if table_z0soil is None or header_z0soil != table_z0soil:
        failures.append(f"Z0SOIL: header {header_z0soil}, table {table_z0soil}")

    if failures:
        for message in failures:
            print(f"FAIL: {message}")
        return 1
    print(f"PASS: {len(rows)} soil categories x {len(COLUMNS)} parameters, "
          f"{len(vegetation_rows)} land-use categories x {n} parameters, CSOIL_DATA, "
          f"RSURF_EXP and Z0SOIL match")
    return 0


if __name__ == '__main__':
    sys.exit(main())
