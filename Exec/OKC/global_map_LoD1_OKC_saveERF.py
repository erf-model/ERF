"""
Convert global map LoD1 shapefile data for Oklahoma City to ERF format.
Processes building heights from shapefile and outputs grid data for ERF simulations.
"""

import shapefile
from sys import exit
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from matplotlib import cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pyproj import Transformer


# ---DETERMINE IF A POINT IS WITHIN A POLYGON-------------------------------
# NOTE: this code is based off of [Franklin, 2000]
def inpolygon(xy, p):
    n = 0  # number of times we cross
    if (p[-1] != p[0]):
        p = tuple(p[:]) + (p[0],)  # copies the first point at the end
    for i in range(len(p) - 1):  # loop over edges
        if (((p[i][1] <= xy[1]) and (p[i + 1][1] > xy[1])) or
            ((p[i][1] > xy[1]) and (p[i + 1][1] <= xy[1]))):
            t = (xy[1] - p[i][1]) / float(p[i + 1][1] - p[i][1])
            if (xy[0] < p[i][0] + t * (p[i + 1][0] - p[i][0])):
                n += 1
    # if the number of crossings is even then the point is outside the polygon
    # if the number of crossings is odd then it is inside the polygon
    return n % 2


def main():
    # Set up the transformer for WGS84 to UTM Zone 14N (Oklahoma City)
    transformer_check = Transformer.from_crs("epsg:4326", "epsg:32614", always_xy=True)

    # Oklahoma City approximate coordinates
    lon, lat = -97.5, 35.5

    # Transform to UTM
    easting, northing = transformer_check.transform(lon, lat)

    print("Easting:", easting)
    print("Northing:", northing)

    # ---OTHER INPUTS-----------------------------------------------------------
    dpi = 300                # Resolution for saved figures (dots per inch)
    cmap = cm.CMRmap_r       # Colormap for height visualization
    showTF = False           # Show plots interactively (False = save only)
    debug = True             # Print debug information during processing

    # ---GRID PARAMETERS--------------------------------------------------------
    dx = 2.5                 # Grid spacing in x-direction [m]
    dy = 2.5                 # Grid spacing in y-direction [m]
    e_we = 100              # Number of grid cells in x-direction (west-east)
    e_sn = 100              # Number of grid cells in y-direction (south-north)
    bbxs = 633000 + 1600          # Bounding box x start (UTM easting) [m]
    bbys = 3923800 + 2450    # Bounding box y start (UTM northing) [m]
    
    # ---GLOBAL MAP LOD1 SHAPE FILE---------------------------------------------
    output_file = "GBA_SK_2.5m_100.txt"

    pn_shp = "/Users/kang18/2024_01_ERF-EB/discretization_urban/OKC_LoD1_v2"
    fn_shp = "GUF04_DLR_v02_w100_n40_w095_n35_OGR04_lod1.shp"

    # Set up the transformer for Web Mercator (3857) to UTM Zone 14N (32614)
    transformer = Transformer.from_crs(3857, 32614, always_xy=True)

    # ---READ SHAPES FROM THE SHAPEFILE-----------------------------------------
    sf = shapefile.Reader(pn_shp + "/" + fn_shp, encoding='latin1')

    # ---GET THE POLYGONS FROM THE SHAPEFILE DATA-------------------------------
    shapes = sf.shapes()
    if (debug):
        print("there are #%4d shapes in " % (len(shapes)) + pn_shp + fn_shp)

    # ---FIND THE SPATIAL EXTENTS OF THE DATA-----------------------------------
    xmin = 9999999999.9
    xmax = 0.0
    ymin = 9999999999.9
    ymax = 0.0
    for shape in shapes:
        for point in shape.points:
            if (point[0] < xmin):
                xmin = point[0]
            if (point[0] > xmax):
                xmax = point[0]
            if (point[1] < ymin):
                ymin = point[1]
            if (point[1] > ymax):
                ymax = point[1]
    if (debug):
        print("[xmin, xmax] = [%9.1f, %9.1f]" % (xmin, xmax))
        print("[ymin, ymax] = [%9.1f, %9.1f]" % (ymin, ymax))

    # ---MODIFY THE BOUNDING BOX------------------------------------------------
    if (bbxs == 0.0):
        bbxs = xmin
    if (bbxs < xmin):
        print("WARNING[process_shp]: bbxs < xmin.")
    xmin = bbxs
    if (bbys == 0.0):
        bbys = ymin
    if (bbys < ymin):
        print("WARNING[process_shp]: bbys < ymin.")
    ymin = bbys
    if (e_we > 0):
        if (np.floor(xmin) + (e_we - 1) * dx > xmax):
            print("WARNING[process_shp]: bbxe > xmax, decrease e_we.")
        xmax = np.floor(xmin) + (e_we - 1) * dx
    elif (e_we != -1):
        print("ERROR[process_shp]: e_we (%3d) is out of bounds." % (e_we))
        exit(1)
    if (e_sn > 0):
        if (np.floor(ymin) + (e_sn - 1) * dy > ymax):
            print("WARNING[process_shp]: bbxe > ymax, decrease e_sn.")
        ymax = np.floor(ymin) + (e_sn - 1) * dy
    elif (e_sn != -1):
        print("ERROR[process_shp]: e_sn (%3d) is out of bounds." % (e_sn))
        exit(1)


    # ---CREATE THE GRID--------------------------------------------------------
    nx = int(np.ceil((xmax - np.floor(xmin)) / dx)) + 1
    ny = int(np.ceil((ymax - np.floor(ymin)) / dy)) + 1
    if (nx % 2 != 1):
        nx += 1
        xmax += dx
    if (ny % 2 != 1):
        ny += 1
        ymax += dy
    if (debug):
        print("[dx, dy] = [%6.2f, %6.2f]" % (dx, dy))
    if (debug):
        print("[nx, ny] = [%6d, %6d]" % (nx, ny))
    xs = np.linspace(np.floor(xmin), np.floor(xmin) + (nx - 1) * dx, nx)
    ys = np.linspace(np.floor(ymin), np.floor(ymin) + (ny - 1) * dy, ny)
    x, y = np.meshgrid(xs, ys)

    # ---READ HEIGHTS AND OTHER ATTRIBUTES FROM THE SHAPEFILE-------------------
    records = sf.records()

    # ---DETERMINE WHICH RECORD IS THE AVERAGE HEIGHT---------------------------
    fields = sf.fields
    if (debug):
        for i in range(0, len(fields)):
            print("field-%2d: " % (i) + fields[i][0])
    i = 0
    while i < len(fields):
        if (fields[i][0] == "height"):
            ir = i
            if (debug):
                print("record 'AVGHT_M' found at ir=%d" % (i))
            break
        elif (fields[i][0] == "ELEVATION"):
            ir = i - 1
            if (debug):
                print("record 'ELEVATION' found at ir=%d" % (i))
            break
        else:
            i += 1
            if (i == len(fields)):
                print("ERROR[process_shp]: AVGHT_M was not found in the record.")
                exit(1)
            else:
                continue

    ir = ir - 1

    # ---CALCULATE GRID HEIGHTS-------------------------------------------------
    v = []  # holds the index for shapes that are in the grid
    h = np.zeros((ny, nx), dtype=float)

    # Pre-filter shapes: only process those that could intersect the grid
    print("Total shapes in shapefile: %d" % len(shapes))
    shapes_to_process = []
    for k in range(len(shapes)):
        # Get shape bounding box in original coordinates and transform to UTM
        shape_bbox = shapes[k].bbox  # [minx, miny, maxx, maxy] in Web Mercator
        bbox_min_utm = transformer.transform(shape_bbox[0], shape_bbox[1])
        bbox_max_utm = transformer.transform(shape_bbox[2], shape_bbox[3])

        shape_xmin = min(bbox_min_utm[0], bbox_max_utm[0])
        shape_xmax = max(bbox_min_utm[0], bbox_max_utm[0])
        shape_ymin = min(bbox_min_utm[1], bbox_max_utm[1])
        shape_ymax = max(bbox_min_utm[1], bbox_max_utm[1])

        # Check if shape bbox intersects with grid bbox
        if (shape_xmax >= xmin and shape_xmin <= xmax and
            shape_ymax >= ymin and shape_ymin <= ymax):
            shapes_to_process.append(k)

    print("Shapes intersecting grid region: %d" % len(shapes_to_process))
    print("Filtering reduced processing by %.1f%%" % (100.0 * (1.0 - len(shapes_to_process) / len(shapes))))

    for k in shapes_to_process:
        if k % 1000 == 0:
            print("Processing shape %d..." % k)
        addTF = True

        ista = nx
        iend = 0
        jsta = ny
        jend = 0
        for point in shapes[k].points:
            point_utm = transformer.transform(point[0], point[1])
            point = point_utm

            if (np.floor((point[0] - x[0, 0]) / dx) < ista):
                ista = max(int(np.floor((point[0] - x[0, 0]) / dx)), 1) - 1
            if (np.ceil((point[0] - x[0, 0]) / dx) > iend):
                iend = min(int(np.ceil((point[0] - x[0, 0]) / dx)), nx - 1) + 1
            if (np.floor((point[1] - y[0, 0]) / dy) < jsta):
                jsta = max(int(np.floor((point[1] - y[0, 0]) / dy)), 1) - 1
            if (np.ceil((point[1] - y[0, 0]) / dy) > jend):
                jend = min(int(np.ceil((point[1] - y[0, 0]) / dy)), ny - 1) + 1

        for i in range(ista, iend):
            for j in range(jsta, jend):
                points_utm = [transformer.transform(point[0], point[1]) for point in shapes[k].points]
                if (inpolygon([x[j, i], y[j, i]], points_utm)):
                    h[j, i] += records[k][ir - 1]
                    if (addTF):
                        v.append(k)
                        addTF = False

    # ---PLOT THE UNMODIFIED GRID HEIGHTS---------------------------------------
    if (debug):
        print("plotting unmodified grid heights")
    fig, ax = plt.subplots(figsize=(12, 12))
    norm = colors.Normalize(vmin=0.0, vmax=150.0)
    im = plt.imshow(h, cmap=cmap, interpolation='none', aspect='equal',
                    origin='lower', norm=norm,
                    extent=[xmin, xmax, ymin, ymax])
    print(xmin, xmax, ymin, ymax)
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%6d'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%7d'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=250.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=250.0))
    ax.set_xlabel('UTM easting')
    ax.set_ylabel('UTM northing')
    cbar = plt.colorbar()
    cbar.set_label("height AGL [m]")
    ax.set_xlim(634000, 635000)
    ax.set_ylim(3925000, 3927000)
    plt.savefig('heights_unmodified.png', dpi=dpi)
    if showTF:
        plt.show()
    plt.close()

    # ---REMOVE WALKWAYS--------------------------------------------------------
    if (debug):
        print("plotting grid heights")
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    div = 5.0
    norm = colors.Normalize(vmin=0, vmax=12)

    x_plot = np.arange(0 + dx / 2, e_we * dx + 1.5 * dx / 2, dx)
    y_plot = np.arange(0 + dy / 2, e_sn * dy + 1.5 * dy / 2, dy)

    dlevel = 10
    levels = np.arange(5, np.max(h) + dlevel, dlevel)

    # Remove walkways
    h_nowalkways = h
    xlo_i = int(1847.5 / dx)
    xhi_i = int(1850 / dy)
    ylo_i = int(1975 / dx)
    yhi_i = int(2007.5 / dy)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    xlo_i = int(1732.5 / dx)
    xhi_i = int(1737.5 / dy)
    ylo_i = int(2105 / dx)
    yhi_i = int(2140 / dy)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    xlo_i = int(1395 / dx)
    xhi_i = int(1400 / dy)
    ylo_i = int(2307.5 / dx)
    yhi_i = int(2322.5 / dy)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    ylo_i = int(2107.5 / dx)
    yhi_i = int(2135 / dx)
    xlo_i = int(1875 / dx)
    xhi_i = int(1880 / dx)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    ylo_i = int(2177.5 / dx)
    yhi_i = int(2202.5 / dx)
    xlo_i = int(1490 / dx)
    xhi_i = int(1497.5 / dx)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    ylo_i = int(2180 / dx)
    yhi_i = int(2187.5 / dx)
    xlo_i = int(1547.5 / dx)
    xhi_i = int(1562.5 / dx)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    ylo_i = int(2187.5 / dx)
    yhi_i = int(2197.5 / dx)
    xlo_i = int(1547.5 / dx)
    xhi_i = int(1560 / dx)
    h_nowalkways[ylo_i:yhi_i, xlo_i:xhi_i] = 0

    im = ax.pcolormesh(x_plot, y_plot, h_nowalkways, cmap=plt.get_cmap('jet'), vmax=100, vmin=0)

    ax.set_aspect('equal', adjustable='box')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.25)
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label('Height [m]')

    ax.set_xlim(1000, 2000)
    ax.set_ylim(2000, 2750)
    plt.savefig('heights_no_walkways.png', dpi=dpi)
    if showTF:
        plt.show()
    plt.close()

    # ---PLOT MODIFIED HEIGHTS--------------------------------------------------
    if (debug):
        print("plotting unmodified grid heights")
    fig, ax = plt.subplots(figsize=(12, 12))
    norm = colors.Normalize(vmin=0.0, vmax=150.0)
    im = plt.imshow(h_nowalkways, cmap=cmap, interpolation='none', aspect='equal',
                    origin='lower', norm=norm,
                    extent=[xmin, xmax, ymin, ymax])
    print(xmin, xmax, ymin, ymax)
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%6d'))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%7d'))
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=250.0))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=250.0))
    ax.set_xlabel('UTM easting')
    ax.set_ylabel('UTM northing')
    cbar = plt.colorbar()
    cbar.set_label("height AGL [m]")
    ax.set_xlim(634000, 635000)
    ax.set_ylim(3925000, 3927000)
    plt.savefig('heights_modified.png', dpi=dpi)
    if showTF:
        plt.show()
    plt.close()

    # ---SAVE FOR ERF-----------------------------------------------------------
    x_nocc = np.arange(0, e_we * dx + dx, dx)
    y_nocc = np.arange(0, e_sn * dy + dy, dy)

    write = np.ones((len(x_plot) * len(y_plot), 3))

    count = 0
    for ii in range(len(x_plot) - 1):
        for jj in range(len(y_plot) - 1):
            write[count, 0] = x_nocc[ii]
            write[count, 1] = y_nocc[jj]
            write[count, 2] = h_nowalkways[jj, ii]
            count += 1

    new_format_write = []
    new_format_write.extend(x_nocc)
    new_format_write.extend(y_nocc)
    new_format_write.extend(write[:, 2])
    print(len(x_plot), len(y_plot))

    header = f"{e_we}\n{e_sn}"
    np.savetxt(output_file, new_format_write, fmt='%.3f', header=header, comments='')
    print("Output saved to ", output_file)


if __name__ == "__main__":
    main()
