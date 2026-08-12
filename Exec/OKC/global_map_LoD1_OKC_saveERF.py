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
from scipy.ndimage import gaussian_filter, laplace, median_filter


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


def spike_targeted_smoothing(h,
                            valley_threshold=2.0, valley_sigma=1.5,
                            peak_threshold=None, peak_sigma=None):
    """
    Apply smoothing to high-curvature regions with separate control for valleys and peaks.
    Preserves building edges and corners.

    Parameters:
    -----------
    h : ndarray
        Height field to smooth
    valley_threshold : float
        Laplacian threshold for detecting valleys (concave) [m]
        Lower = more aggressive, Higher = only extreme valleys
    valley_sigma : float
        Gaussian sigma for smoothing valleys (1.0-3.0)
    peak_threshold : float or None
        Laplacian threshold for detecting peaks (convex) [m]
        If None, peaks are not smoothed
    peak_sigma : float or None
        Gaussian sigma for smoothing peaks (1.0-3.0)
        If None, same as valley_sigma

    Returns:
    --------
    h_result : ndarray
        Smoothed height field
    valley_mask : ndarray (bool)
        Mask showing where valleys were smoothed
    peak_mask : ndarray (bool)
        Mask showing where peaks were smoothed
    """
    # Calculate curvature (Laplacian = second derivative)
    laplacian = laplace(h.astype(float))

    # Detect valleys (positive Laplacian = concave)
    valley_mask = laplacian > valley_threshold

    # Detect peaks (negative Laplacian = convex)
    if peak_threshold is not None:
        peak_mask = laplacian < -peak_threshold
    else:
        peak_mask = np.zeros_like(valley_mask, dtype=bool)

    # Start with original height field
    h_result = h.copy()

    # Smooth valleys
    if valley_sigma > 0 and np.any(valley_mask):
        h_valley_smoothed = gaussian_filter(h, sigma=valley_sigma)
        h_result[valley_mask] = h_valley_smoothed[valley_mask]

    # Smooth peaks (if enabled)
    if peak_threshold is not None and np.any(peak_mask):
        sigma_peak = peak_sigma if peak_sigma is not None else valley_sigma
        h_peak_smoothed = gaussian_filter(h, sigma=sigma_peak)
        h_result[peak_mask] = h_peak_smoothed[peak_mask]

    return h_result, valley_mask, peak_mask


def slope_sign_change_smoothing(h, dx, dy, slope_threshold=3.0, smoothing_sigma=1.5):
    """
    Detect and smooth points where slope sign changes sharply between neighbors.
    Targets V-shaped valleys (slope: -5 to +5) and sharp peaks (slope: +5 to -5).

    Parameters:
    -----------
    h : ndarray
        Height field to smooth
    dx, dy : float
        Grid spacing in x and y directions [m]
    slope_threshold : float
        Minimum absolute slope magnitude to consider (e.g., 3.0 means both slopes must be > 3.0)
    smoothing_sigma : float
        Gaussian sigma for smoothing detected points (1.0-2.0)

    Returns:
    --------
    h_smoothed : ndarray
        Smoothed height field
    mask : ndarray (bool)
        Mask showing where slope-sign-changes were detected and smoothed
    """
    h_smoothed = h.copy().astype(float)
    mask = np.zeros_like(h, dtype=bool)

    ny, nx = h.shape

    # Compute slopes in x-direction (forward difference)
    slope_x = np.zeros_like(h)
    slope_x[:, :-1] = (h[:, 1:] - h[:, :-1]) / dx

    # Compute slopes in y-direction (forward difference)
    slope_y = np.zeros_like(h)
    slope_y[:-1, :] = (h[1:, :] - h[:-1, :]) / dy

    # Detect slope sign changes in x-direction
    for j in range(ny):
        for i in range(1, nx-1):
            slope_left = slope_x[j, i-1]
            slope_right = slope_x[j, i]

            # Check for sign change with large magnitude
            if (slope_left * slope_right < 0 and
                abs(slope_left) > slope_threshold and
                abs(slope_right) > slope_threshold):
                mask[j, i] = True

    # Detect slope sign changes in y-direction
    for j in range(1, ny-1):
        for i in range(nx):
            slope_bottom = slope_y[j-1, i]
            slope_top = slope_y[j, i]

            # Check for sign change with large magnitude
            if (slope_bottom * slope_top < 0 and
                abs(slope_bottom) > slope_threshold and
                abs(slope_top) > slope_threshold):
                mask[j, i] = True

    # Apply smoothing only to detected points
    if np.any(mask):
        h_temp = gaussian_filter(h, sigma=smoothing_sigma)
        h_smoothed[mask] = h_temp[mask]

    return h_smoothed, mask


def limit_gradients(h, dx, dy, max_slope=1.0, max_iterations=50):
    """
    Enforce maximum slope constraint between adjacent cells.
    Iteratively adjusts heights to ensure no slope exceeds max_slope.

    Parameters:
    -----------
    h : ndarray
        Height field to process
    dx, dy : float
        Grid spacing in x and y directions [m]
    max_slope : float
        Maximum allowed dh/dx (e.g., 0.5 = 26.6°, 1.0 = 45°, 1.5 = 56.3°)
    max_iterations : int
        Maximum number of iterations (usually converges in 5-20)

    Returns:
    --------
    h_limited : ndarray
        Height field with slope constraints enforced
    num_iterations : int
        Number of iterations needed for convergence
    """
    h_limited = h.copy().astype(float)
    max_dh = max_slope * min(dx, dy)

    for iteration in range(max_iterations):
        changed = False
        h_new = h_limited.copy()

        for j in range(1, h.shape[0]-1):
            for i in range(1, h.shape[1]-1):
                # Check all 8 neighbors (including diagonals)
                neighbors = [
                    (h_limited[j-1, i], 1.0),       # N
                    (h_limited[j+1, i], 1.0),       # S
                    (h_limited[j, i-1], 1.0),       # W
                    (h_limited[j, i+1], 1.0),       # E
                    (h_limited[j-1, i-1], 1.414),   # NW (sqrt(2) distance)
                    (h_limited[j-1, i+1], 1.414),   # NE
                    (h_limited[j+1, i-1], 1.414),   # SW
                    (h_limited[j+1, i+1], 1.414),   # SE
                ]

                for neighbor_h, dist in neighbors:
                    max_dh_neighbor = max_dh * dist

                    # If current cell is too high relative to neighbor, lower it
                    if h_limited[j, i] - neighbor_h > max_dh_neighbor:
                        new_height = neighbor_h + max_dh_neighbor
                        if new_height < h_new[j, i]:
                            h_new[j, i] = new_height
                            changed = True

        h_limited = h_new

        if not changed:
            return h_limited, iteration + 1

    print(f"WARNING: Gradient limiting did not converge after {max_iterations} iterations")
    return h_limited, max_iterations


def main():
    # ---SMOOTHING CONFIGURATION------------------------------------------------
    # SMOOTHING_METHOD: 'none'|'gaussian'|'median'|'laplacian'|'gradient_limit'|'slope_sign'|'combined' (combined = laplacian + gradient limit)
    # VALLEY_THRESHOLD: Laplacian threshold [m] for valley detection (1.0-2.0=aggressive, 5.0-10.0=conservative)
    # VALLEY_SIGMA: Gaussian sigma for valleys (1.0-1.5=light, 2.0-2.5=moderate, 3.0+=heavy)
    # PEAK_THRESHOLD: Laplacian threshold [m] for peak detection (None=don't smooth peaks, recommended)
    # PEAK_SIGMA: Gaussian sigma for peaks (None=use VALLEY_SIGMA)
    # MAX_SLOPE: Maximum slope (0.5=26.6°, 1.0=45°, 1.5=56.3°)
    # APPLY_GRADIENT_LIMIT_AFTER: Apply gradient limiting after primary smoothing (True=safety net)
    # GAUSSIAN_SIGMA: Gaussian sigma (only for 'gaussian' method)
    # MEDIAN_SIZE: Median filter kernel size (3=3x3 light, 5=5x5 moderate, 7=7x7 heavy, only for 'median' method)
    # SLOPE_THRESHOLD: Minimum slope magnitude for sign-change detection (3.0-5.0=moderate, 8.0-10.0=only extreme)
    # SLOPE_SIGN_SIGMA: Gaussian sigma for slope-sign-change smoothing (1.0-2.0)
    #
    # --- PARAMETER VALUES (select and copy this block) -----------------------
    SMOOTHING_METHOD = 'slope_sign'
    VALLEY_THRESHOLD = 10.0
    VALLEY_SIGMA = 1.5
    PEAK_THRESHOLD = 30.0
    PEAK_SIGMA = 1.0
    MAX_SLOPE = 2.0
    APPLY_GRADIENT_LIMIT_AFTER = False
    GAUSSIAN_SIGMA = 1.5
    MEDIAN_SIZE = 3
    SLOPE_THRESHOLD = 3.0
    SLOPE_SIGN_SIGMA = 2.0
    # --------------------------------------------------------------------------

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
    e_we = 800              # Number of grid cells in x-direction (west-east)
    e_sn = 1680               # Number of grid cells in y-direction (south-north)
    bbxs = 633000 + 500          # Bounding box x start (UTM easting) [m]
    bbys = 3923800 + 400    # Bounding box y start (UTM northing) [m]
    
    # ---GLOBAL MAP LOD1 SHAPE FILE---------------------------------------------
    output_file = f"GBA_SK_2.5m_800_{SMOOTHING_METHOD}.txt"

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

    # ---APPLY SMOOTHING--------------------------------------------------------
    if SMOOTHING_METHOD == 'gaussian':
        if debug:
            print(f"Applying Gaussian smoothing with sigma={GAUSSIAN_SIGMA}")
        h_nowalkways = gaussian_filter(h_nowalkways, sigma=GAUSSIAN_SIGMA)
    elif SMOOTHING_METHOD == 'median':
        if debug:
            print(f"Applying median filter with kernel size={MEDIAN_SIZE}x{MEDIAN_SIZE}")
        h_nowalkways = median_filter(h_nowalkways, size=MEDIAN_SIZE)
        if debug:
            print(f"  Median filter applied to entire grid")
    elif SMOOTHING_METHOD == 'laplacian':
        if debug:
            print(f"Applying Laplacian-weighted smoothing (spike-targeted)")
            print(f"  Valleys: threshold={VALLEY_THRESHOLD}, sigma={VALLEY_SIGMA}")
            if PEAK_THRESHOLD is not None:
                peak_sig = PEAK_SIGMA if PEAK_SIGMA is not None else VALLEY_SIGMA
                print(f"  Peaks:   threshold={PEAK_THRESHOLD}, sigma={peak_sig}")
            else:
                print(f"  Peaks:   not smoothed")

        h_nowalkways, valley_mask, peak_mask = spike_targeted_smoothing(
            h_nowalkways,
            valley_threshold=VALLEY_THRESHOLD,
            valley_sigma=VALLEY_SIGMA,
            peak_threshold=PEAK_THRESHOLD,
            peak_sigma=PEAK_SIGMA
        )

        if debug:
            num_valley_cells = np.sum(valley_mask)
            num_peak_cells = np.sum(peak_mask)
            total_cells = valley_mask.size
            valley_percent = 100.0 * num_valley_cells / total_cells
            peak_percent = 100.0 * num_peak_cells / total_cells
            print(f"  Smoothed {num_valley_cells} valley cells ({valley_percent:.1f}% of grid)")
            print(f"  Smoothed {num_peak_cells} peak cells ({peak_percent:.1f}% of grid)")

            # Save spike mask visualization
            fig, ax = plt.subplots(figsize=(12, 12))

            # Create combined mask with different colors
            # 0 = no smoothing, 1 = valley, 2 = peak
            combined_mask = np.zeros_like(valley_mask, dtype=float)
            combined_mask[valley_mask] = 1.0
            combined_mask[peak_mask] = 2.0

            im = ax.imshow(combined_mask, cmap='RdYlBu_r', alpha=0.7, vmin=0, vmax=2,
                          extent=[xmin, xmax, ymin, ymax], origin='lower')
            ax.set_xlabel('UTM easting')
            ax.set_ylabel('UTM northing')
            ax.set_title('Detected regions (red=valleys, blue=peaks)')
            ax.set_xlim(634000, 635000)
            ax.set_ylim(3925000, 3927000)
            cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2])
            cbar.ax.set_yticklabels(['None', 'Valley', 'Peak'])
            plt.savefig(f'spike_mask_{SMOOTHING_METHOD}.png', dpi=dpi)
            if showTF:
                plt.show()
            plt.close()
    elif SMOOTHING_METHOD == 'slope_sign':
        if debug:
            print(f"Applying slope-sign-change smoothing")
            print(f"  Slope threshold={SLOPE_THRESHOLD}, sigma={SLOPE_SIGN_SIGMA}")
        h_nowalkways, slope_mask = slope_sign_change_smoothing(
            h_nowalkways, dx, dy,
            slope_threshold=SLOPE_THRESHOLD,
            smoothing_sigma=SLOPE_SIGN_SIGMA
        )
        if debug:
            num_smoothed = np.sum(slope_mask)
            total_cells = slope_mask.size
            smoothed_percent = 100.0 * num_smoothed / total_cells
            print(f"  Smoothed {num_smoothed} cells with slope-sign-changes ({smoothed_percent:.1f}% of grid)")

            # Save mask visualization
            fig, ax = plt.subplots(figsize=(12, 12))
            im = ax.imshow(slope_mask.astype(float), cmap='Reds', alpha=0.7,
                          extent=[xmin, xmax, ymin, ymax], origin='lower')
            ax.set_xlabel('UTM easting')
            ax.set_ylabel('UTM northing')
            ax.set_title('Detected slope-sign-change points (red=smoothed)')
            ax.set_xlim(634000, 635000)
            ax.set_ylim(3925000, 3927000)
            cbar = plt.colorbar(im, ax=ax)
            plt.savefig(f'slope_sign_mask_{SMOOTHING_METHOD}.png', dpi=dpi)
            if showTF:
                plt.show()
            plt.close()
    elif SMOOTHING_METHOD == 'gradient_limit':
        if debug:
            print(f"Applying gradient limiting with max_slope={MAX_SLOPE} ({np.degrees(np.arctan(MAX_SLOPE)):.1f}°)")
        h_nowalkways, num_iters = limit_gradients(h_nowalkways, dx, dy, max_slope=MAX_SLOPE)
        if debug:
            print(f"  Converged in {num_iters} iterations")
    elif SMOOTHING_METHOD == 'combined':
        if debug:
            print(f"Applying combined smoothing (Laplacian + Gradient Limiting)")
            print(f"  Step 1 - Laplacian-weighted smoothing:")
            print(f"    Valleys: threshold={VALLEY_THRESHOLD}, sigma={VALLEY_SIGMA}")
            if PEAK_THRESHOLD is not None:
                peak_sig = PEAK_SIGMA if PEAK_SIGMA is not None else VALLEY_SIGMA
                print(f"    Peaks:   threshold={PEAK_THRESHOLD}, sigma={peak_sig}")
            else:
                print(f"    Peaks:   not smoothed")

        # Step 1: Laplacian smoothing
        h_nowalkways, valley_mask, peak_mask = spike_targeted_smoothing(
            h_nowalkways,
            valley_threshold=VALLEY_THRESHOLD,
            valley_sigma=VALLEY_SIGMA,
            peak_threshold=PEAK_THRESHOLD,
            peak_sigma=PEAK_SIGMA
        )

        if debug:
            num_valley_cells = np.sum(valley_mask)
            num_peak_cells = np.sum(peak_mask)
            total_cells = valley_mask.size
            valley_percent = 100.0 * num_valley_cells / total_cells
            peak_percent = 100.0 * num_peak_cells / total_cells
            print(f"    Smoothed {num_valley_cells} valley cells ({valley_percent:.1f}% of grid)")
            print(f"    Smoothed {num_peak_cells} peak cells ({peak_percent:.1f}% of grid)")

        # Step 2: Gradient limiting
        if debug:
            print(f"  Step 2 - Gradient limiting: max_slope={MAX_SLOPE} ({np.degrees(np.arctan(MAX_SLOPE)):.1f}°)")
        h_nowalkways, num_iters = limit_gradients(h_nowalkways, dx, dy, max_slope=MAX_SLOPE)
        if debug:
            print(f"    Converged in {num_iters} iterations")

        # Save spike mask visualization
        if debug:
            fig, ax = plt.subplots(figsize=(12, 12))
            combined_mask = np.zeros_like(valley_mask, dtype=float)
            combined_mask[valley_mask] = 1.0
            combined_mask[peak_mask] = 2.0
            im = ax.imshow(combined_mask, cmap='RdYlBu_r', alpha=0.7, vmin=0, vmax=2,
                          extent=[xmin, xmax, ymin, ymax], origin='lower')
            ax.set_xlabel('UTM easting')
            ax.set_ylabel('UTM northing')
            ax.set_title('Detected regions before gradient limiting (red=valleys, blue=peaks)')
            ax.set_xlim(634000, 635000)
            ax.set_ylim(3925000, 3927000)
            cbar = plt.colorbar(im, ax=ax, ticks=[0, 1, 2])
            cbar.ax.set_yticklabels(['None', 'Valley', 'Peak'])
            plt.savefig(f'spike_mask_{SMOOTHING_METHOD}.png', dpi=dpi)
            if showTF:
                plt.show()
            plt.close()
    elif SMOOTHING_METHOD == 'none':
        if debug:
            print("No smoothing applied")
    else:
        print(f"WARNING: Unknown smoothing method '{SMOOTHING_METHOD}', skipping smoothing")

    # ---OPTIONAL POST-PROCESSING: GRADIENT LIMITING----------------------------
    if SMOOTHING_METHOD not in ['none', 'gradient_limit', 'combined'] and APPLY_GRADIENT_LIMIT_AFTER:
        if debug:
            print(f"Applying gradient limiting post-process: max_slope={MAX_SLOPE} ({np.degrees(np.arctan(MAX_SLOPE)):.1f}°)")
        h_nowalkways, num_iters = limit_gradients(h_nowalkways, dx, dy, max_slope=MAX_SLOPE)
        if debug:
            print(f"  Converged in {num_iters} iterations")

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
    plt.savefig(f'heights_modified_{SMOOTHING_METHOD}.png', dpi=dpi)
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
