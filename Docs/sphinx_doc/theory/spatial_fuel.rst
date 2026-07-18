
 .. role:: cpp(code)
    :language: c++

.. _sec:SpatialFuel:

Spatial Fuel Map and Firebreak Barriers
========================================

Physical Background
--------------------

Real wildfire simulations require spatially varying fuel properties because fuel models control fire spread rate, intensity, and flame characteristics. Anderson's 13 fuel model classification captures key differences in fuel bed structure, moisture behavior, and fire dynamics across vegetation types (grass, shrub, timber litter, etc.). A spatially continuous representation of fuel variation improves accuracy of simulated fire perimeter growth and heat flux distribution, enabling better prediction of fire behavior in heterogeneous landscapes. Firebreak barriers represent suppression infrastructure such as fuel breaks, cleared roads, and engineered barriers that prevent fire propagation.

Fuel Map File Formats
---------------------

Two file formats are supported for loading spatially varying fuel model codes:

**ESRI ASCII Grid Format**

This is a simple, human-readable raster format. File structure:

.. code-block:: text

    ncols <N>
    nrows <N>
    xllcorner <x>
    yllcorner <y>
    cellsize <dx>
    nodata_value <val>
    <row 1 of integer codes>
    <row 2 of integer codes>
    ...
    <row N of integer codes>

Header lines specify grid dimensions, geospatial coordinates, cell size, and the sentinel value for missing/invalid cells. Data rows are ordered from north (highest y-coordinate) to south (lowest y-coordinate), following the ESRI standard convention. The ERF fire reader reverses this row order internally to match fire domain coordinates (south to north), ensuring that cell (i, j) in the fire grid is correctly mapped to its corresponding fuel model code.

**FARSITE LCP Binary Format**

FARSITE landscape files (.lcp) are binary rasters used in operational fire simulation. The reader extracts the fuel model layer (layer 2 in the LCP specification) and maps it to the fire grid. Cell dimensions must match the fire grid exactly.

Reference: Finney, M.A. (1998/2004). FARSITE: Fire Area Simulator. RMRS-RP-4.

Fuel Boundary Blending
----------------------

At boundaries between different fuel model zones, ROS discontinuities can produce unphysical artifacts. To mitigate this, ROS values at zone boundaries are blended toward neighboring cells with different fuel codes. The blending formula is:

.. math::

    \text{ROS}_{\text{blended}}(i,j) = (1-f) \, \text{ROS}_{\text{cell}}(i,j) + 
    f \, \frac{1}{n_{\text{diff}}} \sum_{\text{neighbors with different fuel code}} \text{ROS}_{\text{neighbor}}

where :math:`f` is the blending fraction (default 0.0, range [0, 1]), :math:`n_{\text{diff}}` is the count of neighbors with fuel codes differing from the central cell, and neighbors are the 4-cell von Neumann stencil (±i, ±j directions). A value :math:`f=0` disables blending (sharp boundaries); :math:`f=1` fully replaces the cell ROS with the mean neighbor ROS (maximum smoothing).

When a cell has no neighbors with different fuel codes, it is not modified. Blending is applied after all ROS calculations but before fire propagation, allowing smooth transitions in fire front shape and intensity across zone boundaries.

Firebreak Barriers
------------------

Firebreak barriers are permanent non-burnable zones established at fire initialization. Two barrier types are supported:

**Rectangular Barriers**

A rectangular barrier is defined by x-extent [:math:`x_{\text{lo}}, x_{\text{hi}}`] and y-extent [:math:`y_{\text{lo}}, y_{\text{hi}}`]. A cell is contained within the barrier if its centre coordinates (x, y) satisfy:

.. math::

    x_{\text{lo}} \le x \le x_{\text{hi}} \quad \text{and} \quad y_{\text{lo}} \le y \le y_{\text{hi}}

**Circular Barriers**

A circular barrier is defined by centre :math:`(c_x, c_y)` and radius :math:`r`. A cell is contained within the barrier if its centre lies within the radius:

.. math::

    (x - c_x)^2 + (y - c_y)^2 \le r^2

**Sentinel Value and Permanence**

All cells within any barrier have their level-set field ``fire_phi`` stamped to the constant value ``FIREBREAK_PHI_SENTINEL = 1.0e6`` at initialization, after ignition. This sentinel is strictly greater than ``farsite_phi_threshold`` (default 0.1) and any unburned cell phi value. Once set, firebreak cells remain at the sentinel value throughout the simulation, permanently preventing fire arrival and propagation. Firebreak barriers cannot be burned through or overcome during the fire simulation.

Reference: Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4.

Input Parameters
----------------

Phase 10 adds fuel map and firebreak configuration parameters to the fire module. All parameters use the ``erf.fire.*`` prefix. Consult ``Source/Fire/ERF_FireParams.H`` for default values and code documentation.

.. list-table:: Spatial Fuel Map Parameters
   :widths: 35 15 15 35
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``fuel_map.file``
     - string
     - ""
     - Path to fuel map file (ASCII or LCP); empty = uniform fuel (fallback to ``fuel_model_id``)
   * - ``fuel_map.format``
     - string
     - "ascii"
     - File format: "ascii" (ESRI ASCII grid) or "lcp" (FARSITE binary)
   * - ``fuel_map.blending_fraction``
     - Real
     - 0.0
     - ROS blending weight at fuel boundaries [0–1]

.. list-table:: Firebreak Barrier Parameters (Indexed)
   :widths: 35 15 15 35
   :header-rows: 1

   * - Parameter
     - Type
     - Default
     - Description
   * - ``firebreak.N.type``
     - string
     - —
     - Barrier type: "rect" (rectangular) or "circle" (circular)
   * - ``firebreak.N.x_lo``
     - Real
     - —
     - Rectangle: west edge x-coordinate [m]
   * - ``firebreak.N.x_hi``
     - Real
     - —
     - Rectangle: east edge x-coordinate [m]
   * - ``firebreak.N.y_lo``
     - Real
     - —
     - Rectangle: south edge y-coordinate [m]
   * - ``firebreak.N.y_hi``
     - Real
     - —
     - Rectangle: north edge y-coordinate [m]
   * - ``firebreak.N.cx``
     - Real
     - —
     - Circle: centre x-coordinate [m]
   * - ``firebreak.N.cy``
     - Real
     - —
     - Circle: centre y-coordinate [m]
   * - ``firebreak.N.radius``
     - Real
     - —
     - Circle: radius [m]

Barriers are indexed N = 0, 1, 2, ... and defined by setting at least the ``type`` key. The configuration loop stops at the first missing barrier index, allowing up to 10 barriers. Unused parameters for a given barrier type (e.g., ``cx``, ``cy``, ``radius`` for a rectangular barrier) are ignored.

Output Variables
-----------------

Phase 10 adds a single new plotfile variable:

- ``fire_fuel_model_code`` — per-cell fuel model code (Anderson FBFM13, 1–13) as read from the fuel map or defaulting to ``fuel_model_id``. Integer values stored as Real; includes code 0 for cells with nodata values in the input file.

This variable enables verification that fuel maps were loaded correctly and visualization of fuel model spatial distribution.

Limitations
-----------

- **Exact dimension match**: Fuel map grid dimensions must equal the fire grid dimensions exactly. No interpolation or regridding is performed.
- **Permanent barriers**: Firebreak cells are stamped at initialization and remain non-burnable for the duration of the simulation. Barriers cannot evolve or be modified during the run.
- **LCP layer extraction**: The FARSITE LCP reader extracts only the fuel model layer (layer 2). Elevation, aspect, slope, fuel moisture, and other LCP layers are not currently loaded.
- **Valid fuel codes**: FBFM13 recognizes burnable fuel codes 1–13. Code 0 (nodata), code 14 and higher (out-of-range), and code -9999 (sentinel) are treated as non-burnable or invalid and default to the fallback ``fuel_model_id``.

References
----------

Anderson, H.E. (1982). Aids to determining fuel models for estimating fire behavior. General Technical Report INT-122. USDA Forest Service, Intermountain Forest and Range Experiment Station.

Finney, M.A. (1998). FARSITE: Fire Area Simulator—Model development and evaluation. Research Paper RMRS-RP-4, Revised 2004. USDA Forest Service, Rocky Mountain Research Station.

Scott, J.H., Burgan, R.E. (2005). Standard fire behavior fuel models: A comprehensive set for use with ArcGIS. General Technical Report RMRS-GTR-153. USDA Forest Service, Rocky Mountain Research Station.
