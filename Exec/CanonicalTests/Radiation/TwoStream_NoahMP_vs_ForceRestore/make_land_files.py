#!/usr/bin/env python3
"""Write the land setup files of the comparison as CDL and, with ncgen, as NetCDF.

Noah-MP reads its land state from wrfinput-format files: wrfinput_d01 for level 0 and,
since every level here runs the land model itself, wrfinput_d02 for level 1 (named by
ERF_SETUP_FILE_01/02 in namelist.erf). Both are the grassland of
Exec/RegTests/NoahMP_Ideal (IVGTYP 10, silty clay loam, 290 K soil at 0.25 m3/m3) on this
case's grids:

  level 0  8 x 8 cells of 250 m (amr.n_cell 8 8, geometry.prob_extent 2000 2000);
  level 1  8 x 8 cells of 125 m over the middle 4 x 4 of level 0 (erf.centre.in_box
           500 to 1500 m): parent cell (3, 3), 1-based, ratio 2.

With --barren the land is instead bare (IVGTYP 16, barren or sparsely vegetated; no
vegetation, no leaves) on the same soil at its wilting point (0.12 m3/m3), for the
dry-surface variant of the comparison (README.md). --smois sets the soil water instead
(0.25 for the moist bare variant).

Usage: python3 make_land_files.py [--ncgen PATH] [--barren] [--smois Q]
"""
import argparse
import subprocess

FIELDS_2D = [  # name, type, units, value
    ('XLAT', 'float', 'degree_north', '40.0'), ('XLONG', 'float', 'degree_east', '-100.0'),
    ('XLAND', 'float', '1', '1.0'), ('IVGTYP', 'int', None, '10'), ('ISLTYP', 'int', None, '8'),
    ('TSK', 'float', 'K', '300.0'), ('TMN', 'float', 'K', '287.0'), ('HGT', 'float', 'm', '0.0'),
    ('SEAICE', 'float', '1', '0.0'), ('SNOW', 'float', 'kg m-2', '0.0'), ('SNOWH', 'float', 'm', '0.0'),
    ('SNOWC', 'float', '1', '0.0'), ('CANWAT', 'float', 'kg m-2', '0.0'), ('VEGFRA', 'float', '1', '50.0'),
    ('SHDMAX', 'float', '1', '80.0'), ('SHDMIN', 'float', '1', '10.0'), ('LAI', 'float', 'm2 m-2', '2.0'),
    ('MAPFAC_MX', 'float', '1', '1.0'), ('MAPFAC_MY', 'float', '1', '1.0')]
FIELDS_SOIL = [('TSLB', 'K', '290.0'), ('SMOIS', 'm3 m-3', '0.25')]


def cdl(name, n, dx, grid_id, ratio, parent_start, overrides):
    dims = '(Time, south_north, west_east)'
    lines = [f'netcdf {name} {{', 'dimensions:', '    Time = UNLIMITED ;', '    DateStrLen = 19 ;',
             f'    west_east = {n} ;', f'    south_north = {n} ;', '    soil_layers_stag = 4 ;',
             'variables:', '    char Times(Time, DateStrLen) ;']
    for var, typ, units, _ in FIELDS_2D:
        lines.append(f'    {typ} {var}{dims} ;' + (f' {var}:units = "{units}" ;' if units else ''))
    for var, units, _ in FIELDS_SOIL:
        lines.append(f'    float {var}(Time, soil_layers_stag, south_north, west_east) ; {var}:units = "{units}" ;')
    lines.append('    float DZS(Time, soil_layers_stag) ; DZS:units = "m" ;')
    attrs = [('TITLE', '"SYNTHETIC IDEAL WRFINPUT FOR TwoStream_NoahMP_vs_ForceRestore"'),
             ('SIMULATION_START_DATE', '"2024-08-05_15:00:00"'),
             ('WEST-EAST_GRID_DIMENSION', n + 1), ('SOUTH-NORTH_GRID_DIMENSION', n + 1),
             ('BOTTOM-TOP_GRID_DIMENSION', 2), ('DX', f'{dx}f'), ('DY', f'{dx}f'),
             ('GRID_ID', grid_id), ('grid_id', grid_id), ('PARENT_GRID_RATIO', ratio),
             ('I_PARENT_START', parent_start), ('J_PARENT_START', parent_start), ('MAP_PROJ', 0),
             ('MMINLU', '"MODIFIED_IGBP_MODIS_NOAH"'), ('ISWATER', 17), ('ISLAKE', 21), ('ISICE', 15),
             ('ISURBAN', 13), ('TRUELAT1', '30.f'), ('TRUELAT2', '60.f'), ('STAND_LON', '-100.f')]
    lines.append('')
    lines += [f'    :{k} = {v} ;' for k, v in attrs]
    lines += ['data:', '    Times = "2024-08-05_15:00:00" ;']
    for var, _, _, val in FIELDS_2D:
        val = overrides.get(var, val)
        lines.append(f'    {var} = ' + ', '.join([val] * (n * n)) + ' ;')
    for var, _, val in FIELDS_SOIL:
        val = overrides.get(var, val)
        lines.append(f'    {var} = ' + ', '.join([val] * (4 * n * n)) + ' ;')
    lines += ['    DZS = 0.1, 0.3, 0.6, 1.0 ;', '}']
    return '\n'.join(lines) + '\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--ncgen', default='ncgen')
    parser.add_argument('--barren', action='store_true',
                        help='bare land on the same soil at its wilting point')
    parser.add_argument('--smois', type=float, default=None,
                        help='soil water content of every layer [m3/m3] (default: 0.25, or '
                             '0.12 with --barren)')
    args = parser.parse_args()
    overrides = {}
    if args.barren:
        overrides = {'IVGTYP': '16', 'VEGFRA': '0.0', 'SHDMAX': '0.0', 'SHDMIN': '0.0',
                     'LAI': '0.0', 'SMOIS': '0.12'}
    if args.smois is not None:
        if not 0.0 < args.smois < 1.0:
            parser.error(f'--smois must be in (0, 1), got {args.smois}')
        overrides['SMOIS'] = repr(args.smois)
    for name, n, dx, grid_id, ratio, start in (('wrfinput_d01', 8, 250.0, 1, 1, 1),
                                               ('wrfinput_d02', 8, 125.0, 2, 2, 3)):
        with open(f'{name}.cdl', 'w') as handle:
            handle.write(cdl(name, n, dx, grid_id, ratio, start, overrides))
        subprocess.run([args.ncgen, '-o', name, f'{name}.cdl'], check=True)
        print(f'wrote {name}.cdl and {name}')


if __name__ == '__main__':
    main()
