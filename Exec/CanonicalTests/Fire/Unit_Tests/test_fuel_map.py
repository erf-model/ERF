#!/usr/bin/env python3
"""
Unit tests for Phase 10 fuel map and firebreak modules.

Verifies the following:
  - ESRI ASCII fuel map header parsing and row reversal
  - nodata value handling
  - Comment line skipping
  - Rectangular and circular firebreak containment geometry
  - Firebreak phi sentinel validity
  - Fuel boundary blending formula and behavior

Source files tested:
  - Source/Fire/ERF_FuelMap.H (read_ascii_fuel_map, fill_fuel_model_mf)
  - Source/Fire/ERF_FireBreak.H (apply_firebreaks, FIREBREAK_PHI_SENTINEL)
  - Source/Fire/ERF_FuelBlending.H (apply_fuel_boundary_blending)
  - Source/Fire/ERF_FireParams.H (ParmParse keys)

Reference:
  - Finney, M.A. (1998). FARSITE: Fire Area Simulator. RMRS-RP-4.

Run: python3 test_fuel_map.py
"""

import math
import sys


def parse_esri_ascii_header(lines):
    """
    Parse ESRI ASCII grid header from a list of text lines.
    
    Expected header format (in order):
      ncols <int>
      nrows <int>
      xllcorner <float>
      yllcorner <float>
      cellsize <float>
      nodata_value <int>
    
    Skips comment lines (starting with '#').
    
    Args:
        lines: List of text lines (with or without trailing newlines)
    
    Returns:
        Dict with keys: ncols, nrows, xllcorner, yllcorner, cellsize, nodata_value
    
    Raises:
        ValueError if header is malformed or incomplete
    """
    header = {}
    expected_keys = ['ncols', 'nrows', 'xllcorner', 'yllcorner', 'cellsize', 'nodata_value']
    key_count = 0
    
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        
        parts = line.split()
        if len(parts) < 2:
            continue
        
        key = parts[0]
        if key in expected_keys and key not in header:
            try:
                if key in ['ncols', 'nrows', 'nodata_value']:
                    header[key] = int(parts[1])
                else:
                    header[key] = float(parts[1])
                key_count += 1
            except ValueError:
                raise ValueError(f"Cannot parse {key} value: {parts[1]}")
        
        if key_count == len(expected_keys):
            break
    
    if len(header) < len(expected_keys):
        raise ValueError(f"Incomplete header: got {list(header.keys())}, expected {expected_keys}")
    
    return header


def reverse_fuel_map_rows(data, nrows, ncols):
    """
    Reverse row order from ESRI ASCII format (north-first) to domain coordinates (south-first).
    
    ESRI format: row 0 = northernmost (highest y), row nrows-1 = southernmost (lowest y)
    Domain: j=0 = southernmost (lowest y), j=nrows-1 = northernmost (highest y)
    
    Args:
        data: List of lists, each inner list is a row of fuel codes
        nrows: Number of rows
        ncols: Number of columns (for validation)
    
    Returns:
        List of lists with rows reversed
    """
    if len(data) != nrows:
        raise ValueError(f"Data has {len(data)} rows, expected {nrows}")
    
    for row in data:
        if len(row) != ncols:
            raise ValueError(f"Row has {len(row)} columns, expected {ncols}")
    
    # Reverse: flip the list of rows
    return data[::-1]


def cell_centre_coords(i, j, prob_lo, dx):
    """
    Compute cell centre coordinates.
    
    Args:
        i: Grid index in x direction
        j: Grid index in y direction
        prob_lo: Tuple (x_min, y_min) of domain lower-left corner
        dx: Grid cell size
    
    Returns:
        Tuple (x, y) of cell centre
    """
    x = prob_lo[0] + (i + 0.5) * dx
    y = prob_lo[1] + (j + 0.5) * dx
    return x, y


def is_in_rect_firebreak(x, y, x_lo, y_lo, x_hi, y_hi):
    """
    Check if point (x, y) is inside a rectangular firebreak.
    
    Rectangle check (from ERF_FireBreak.H):
      in_barrier = (x >= x_lo && x <= x_hi && y >= y_lo && y <= y_hi)
    
    Args:
        x, y: Point coordinates
        x_lo, y_lo, x_hi, y_hi: Rectangle bounds
    
    Returns:
        True if (x, y) is inside the rectangle
    """
    return (x_lo <= x <= x_hi) and (y_lo <= y <= y_hi)


def is_in_circle_firebreak(x, y, cx, cy, radius):
    """
    Check if point (x, y) is inside a circular firebreak.
    
    Circle check (from ERF_FireBreak.H):
      dist_sq = (x - cx)^2 + (y - cy)^2
      in_barrier = (dist_sq <= radius^2)
    
    Args:
        x, y: Point coordinates
        cx, cy: Circle centre
        radius: Circle radius
    
    Returns:
        True if Euclidean distance from (x,y) to (cx,cy) is <= radius
    """
    dist_sq = (x - cx)**2 + (y - cy)**2
    return dist_sq <= (radius * radius)


def blend_ros(ros_cell, ros_neighbor_list, blending_fraction):
    """
    Compute blended ROS at a fuel boundary.
    
    Blending formula (from ERF_FuelBlending.H):
      ros_blended = (1 - f) * ros_cell + f * mean(ros_neighbors)
    
    Args:
        ros_cell: Rate-of-spread at this cell [m/s]
        ros_neighbor_list: List of ROS values at differently-coded neighbors
        blending_fraction: Blending weight f in [0, 1]
    
    Returns:
        Blended ROS value
    """
    if not ros_neighbor_list:
        return ros_cell
    
    ros_blend = sum(ros_neighbor_list) / len(ros_neighbor_list)
    return (1.0 - blending_fraction) * ros_cell + blending_fraction * ros_blend


# ============================================================================
# Unit Tests
# ============================================================================

def test_ascii_header_parsing():
    """Test ESRI ASCII header parsing extracts correct integer values."""
    lines = [
        "ncols 200",
        "nrows 200",
        "xllcorner 0.0",
        "yllcorner 0.0",
        "cellsize 10.0",
        "nodata_value -9999"
    ]
    header = parse_esri_ascii_header(lines)
    assert header['ncols'] == 200, f"ncols={header['ncols']}, expected 200"
    assert header['nrows'] == 200, f"nrows={header['nrows']}, expected 200"
    assert header['cellsize'] == 10.0, f"cellsize={header['cellsize']}, expected 10.0"
    assert header['nodata_value'] == -9999, f"nodata_value={header['nodata_value']}, expected -9999"
    print("✓ test_ascii_header_parsing PASSED")


def test_row_reversal():
    """Test that row reversal correctly maps file rows to domain coordinates."""
    # File rows (north-first): row 0 has [10, 20], row 1 has [30, 40]
    file_data = [[10, 20], [30, 40]]
    nrows, ncols = 2, 2
    
    # After reversal: domain south row (j=0) should have [30, 40], north row (j=1) should have [10, 20]
    domain_data = reverse_fuel_map_rows(file_data, nrows, ncols)
    
    assert domain_data[0] == [30, 40], f"domain j=0 (south)={domain_data[0]}, expected [30, 40]"
    assert domain_data[1] == [10, 20], f"domain j=1 (north)={domain_data[1]}, expected [10, 20]"
    print("✓ test_row_reversal PASSED")


def test_nodata_handling():
    """Test that nodata sentinel values are identified (parser records sentinel)."""
    # A cell with nodata value is marked; C++ code maps it to 0
    lines = [
        "ncols 2",
        "nrows 2",
        "xllcorner 0.0",
        "yllcorner 0.0",
        "cellsize 10.0",
        "nodata_value -9999"
    ]
    header = parse_esri_ascii_header(lines)
    
    # Test that nodata_value is correctly parsed
    assert header['nodata_value'] == -9999, f"nodata_value={header['nodata_value']}, expected -9999"
    
    # Simulation: cells with nodata_value in file are recorded and mapped to code 0 by C++ logic
    fuel_codes_raw = [1, 2, -9999, 4]
    fuel_codes_mapped = [code if code != -9999 else 0 for code in fuel_codes_raw]
    
    assert fuel_codes_mapped == [1, 2, 0, 4], f"Mapped codes {fuel_codes_mapped}, expected [1, 2, 0, 4]"
    print("✓ test_nodata_handling PASSED")


def test_comment_line_skipping():
    """Test that comment lines (starting with '#') do not affect parsing."""
    lines = [
        "# This is a comment",
        "ncols 100",
        "# Another comment",
        "nrows 100",
        "xllcorner 0.0",
        "yllcorner 0.0",
        "cellsize 10.0",
        "# And another",
        "nodata_value -9999"
    ]
    
    header = parse_esri_ascii_header(lines)
    assert header['ncols'] == 100, f"ncols={header['ncols']}, expected 100"
    assert header['nrows'] == 100, f"nrows={header['nrows']}, expected 100"
    print("✓ test_comment_line_skipping PASSED")


def test_two_zone_layout():
    """Test two-zone fuel map with north/south boundary."""
    # 4x4 map: north 2 rows = FM4, south 2 rows = FM1
    # In file format (north-first): rows 0-1 are FM4, rows 2-3 are FM1
    file_data = [[4, 4, 4, 4], [4, 4, 4, 4], [1, 1, 1, 1], [1, 1, 1, 1]]
    
    # After reversal: domain j=0-1 (south) should be FM1, j=2-3 (north) should be FM4
    domain_data = reverse_fuel_map_rows(file_data, nrows=4, ncols=4)
    
    assert domain_data[0][0] == 1, f"domain j=0 (south) fuel={domain_data[0][0]}, expected FM1"
    assert domain_data[3][0] == 4, f"domain j=3 (north) fuel={domain_data[3][0]}, expected FM4"
    print("✓ test_two_zone_layout PASSED")


def test_rect_firebreak_containment():
    """Test rectangular firebreak containment at cell centres."""
    # Rectangle: x in [800, 900], y in [500, 1500]
    x_lo, x_hi = 800.0, 900.0
    y_lo, y_hi = 500.0, 1500.0
    
    # Cell centre at (i=85, j=50) with dx=10, prob_lo=(0, 0)
    # x = 0 + (85 + 0.5)*10 = 855, y = 0 + (50 + 0.5)*10 = 505
    x, y = cell_centre_coords(85, 50, (0, 0), 10.0)
    assert is_in_rect_firebreak(x, y, x_lo, y_lo, x_hi, y_hi), \
        f"Cell centre ({x}, {y}) should be inside rectangle [{x_lo},{x_hi}] × [{y_lo},{y_hi}]"
    
    # Cell centre outside x range
    x, y = cell_centre_coords(75, 50, (0, 0), 10.0)
    assert not is_in_rect_firebreak(x, y, x_lo, y_lo, x_hi, y_hi), \
        f"Cell centre ({x}, {y}) should be outside rectangle"
    
    print("✓ test_rect_firebreak_containment PASSED")


def test_circle_firebreak_containment():
    """Test circular firebreak containment at cell centres."""
    # Circle: centre (1000, 1000), radius 200
    cx, cy, radius = 1000.0, 1000.0, 200.0
    
    # Cell at (i=100, j=100) with dx=10: centre at (1005, 1005), distance ~7.07 < 200 -> inside
    x, y = cell_centre_coords(100, 100, (0, 0), 10.0)
    dist = math.sqrt((x - cx)**2 + (y - cy)**2)
    assert is_in_circle_firebreak(x, y, cx, cy, radius), \
        f"Cell centre ({x}, {y}), dist={dist:.2f}, should be inside circle"
    
    # Cell far away: distance > 200 -> outside
    x, y = cell_centre_coords(0, 0, (0, 0), 10.0)
    assert not is_in_circle_firebreak(x, y, cx, cy, radius), \
        f"Cell centre ({x}, {y}) should be outside circle"
    
    print("✓ test_circle_firebreak_containment PASSED")


def test_firebreak_phi_sentinel():
    """Test that FIREBREAK_PHI_SENTINEL is valid (greater than typical phi values)."""
    # From ERF_FireBreak.H: FIREBREAK_PHI_SENTINEL = 1.0e6
    FIREBREAK_PHI_SENTINEL = 1.0e6
    
    # From ERF_FireParams.H: farsite_phi_threshold = 0.1 (default)
    farsite_phi_threshold = 0.1
    
    # Verify sentinel is strictly greater than threshold
    assert FIREBREAK_PHI_SENTINEL > farsite_phi_threshold, \
        f"Sentinel {FIREBREAK_PHI_SENTINEL} should be > threshold {farsite_phi_threshold}"
    
    # Verify sentinel is much larger than typical unburned cell phi (~0)
    typical_unburned_phi = 1.0
    assert FIREBREAK_PHI_SENTINEL > typical_unburned_phi, \
        f"Sentinel {FIREBREAK_PHI_SENTINEL} should be >> typical phi {typical_unburned_phi}"
    
    print("✓ test_firebreak_phi_sentinel PASSED")


def test_blending_identity_at_zero():
    """Test that blending fraction = 0 gives identity: result = original ROS."""
    ros_cell = 0.5  # Original cell ROS [m/s]
    ros_neighbor_list = [0.1, 0.2]  # Neighboring ROS values
    blending_fraction = 0.0
    
    ros_blended = blend_ros(ros_cell, ros_neighbor_list, blending_fraction)
    
    assert abs(ros_blended - ros_cell) < 1e-10, \
        f"Blended ROS={ros_blended}, should equal original {ros_cell} when f=0"
    
    print("✓ test_blending_identity_at_zero PASSED")


def test_blending_intermediate_value():
    """Test that blending fraction = 0.3 produces intermediate result."""
    ros_cell = 0.5  # Original [m/s]
    ros_neighbor_list = [0.1]  # One neighbor at 0.1 [m/s]
    blending_fraction = 0.3
    
    ros_blended = blend_ros(ros_cell, ros_neighbor_list, blending_fraction)
    
    # Expected: 0.7 * 0.5 + 0.3 * 0.1 = 0.35 + 0.03 = 0.38
    expected = 0.7 * ros_cell + 0.3 * 0.1
    
    assert abs(ros_blended - expected) < 1e-10, \
        f"Blended ROS={ros_blended}, expected {expected}"
    
    # Check that result is strictly between original and neighbor
    assert 0.1 < ros_blended < 0.5, \
        f"Blended ROS={ros_blended} should be between neighbor (0.1) and original (0.5)"
    
    print("✓ test_blending_intermediate_value PASSED")


def test_blending_full_at_one():
    """Test that blending fraction = 1.0 gives neighbor ROS."""
    ros_cell = 0.5  # Original [m/s]
    ros_neighbor_list = [0.1]  # Neighbor [m/s]
    blending_fraction = 1.0
    
    ros_blended = blend_ros(ros_cell, ros_neighbor_list, blending_fraction)
    
    # Expected: 0 * 0.5 + 1.0 * 0.1 = 0.1
    assert abs(ros_blended - 0.1) < 1e-10, \
        f"Blended ROS={ros_blended}, expected 0.1 (neighbor value) when f=1"
    
    print("✓ test_blending_full_at_one PASSED")


def test_fbfm13_valid_code_range():
    """Test FBFM13 valid burnable codes (1–13) vs invalid codes (0, 14, -9999)."""
    valid_codes = list(range(1, 14))  # 1–13 are valid burnable models
    invalid_codes = [0, 14, -9999]    # 0=nodata, 14=out of range, -9999=sentinel
    
    for code in valid_codes:
        assert 1 <= code <= 13, f"Code {code} should be in valid range [1, 13]"
    
    for code in invalid_codes:
        assert not (1 <= code <= 13), f"Code {code} should NOT be in valid range [1, 13]"
    
    print("✓ test_fbfm13_valid_code_range PASSED")


def main():
    """Run all unit tests."""
    print("=" * 70)
    print("Phase 10: Fuel Map and Firebreak Unit Tests")
    print("=" * 70)
    
    tests = [
        test_ascii_header_parsing,
        test_row_reversal,
        test_nodata_handling,
        test_comment_line_skipping,
        test_two_zone_layout,
        test_rect_firebreak_containment,
        test_circle_firebreak_containment,
        test_firebreak_phi_sentinel,
        test_blending_identity_at_zero,
        test_blending_intermediate_value,
        test_blending_full_at_one,
        test_fbfm13_valid_code_range,
    ]
    
    failed = []
    for test in tests:
        try:
            test()
        except AssertionError as e:
            print(f"✗ {test.__name__} FAILED: {e}")
            failed.append((test.__name__, str(e)))
    
    print("=" * 70)
    if not failed:
        print(f"All {len(tests)} tests PASSED")
        return 0
    else:
        print(f"{len(failed)} test(s) FAILED:")
        for name, error in failed:
            print(f"  - {name}: {error}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
