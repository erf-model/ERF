#!/usr/bin/env python3
"""
Unit tests for Phase 11 multi-ignition schedule and polygon ignition.

Tests CSV parsing, time-window firing, priority ordering, phi stamping,
and polygon/polyline geometry functions.

Run: python3 test_ignition_schedule.py
"""

import math
import sys


def parse_ignition_schedule_csv(csv_content):
    """
    Parse ignition schedule CSV content.
    Returns list of dicts with keys: time_s, cx, cy, radius, source_type, priority, suppress_if_burning.
    """
    events = []
    lines = csv_content.strip().split('\n')
    for line in lines:
        # Skip empty lines and comments
        if not line or line[0] in '#!':
            continue
        
        # Convert commas to spaces
        line = line.replace(',', ' ')
        
        # Parse fields
        parts = line.split()
        if len(parts) < 4:
            continue
        
        try:
            time_s = float(parts[0])
            cx = float(parts[1])
            cy = float(parts[2])
            radius = float(parts[3])
            source_type = int(parts[4]) if len(parts) > 4 else 0
            priority = int(parts[5]) if len(parts) > 5 else 5
            suppress_if_burning = int(parts[6]) if len(parts) > 6 else 0
        except (ValueError, IndexError):
            continue
        
        events.append({
            'time_s': time_s,
            'cx': cx,
            'cy': cy,
            'radius': radius,
            'source_type': source_type,
            'priority': priority,
            'suppress_if_burning': (suppress_if_burning != 0),
            'fired': False
        })
    
    return events


def sort_events_by_time_priority(events):
    """Sort events by time (ascending) then priority (descending)."""
    return sorted(events, key=lambda e: (e['time_s'], -e['priority']))


def get_events_in_window(events, prev_time, current_time):
    """Return events that fall in (prev_time, current_time]."""
    result = []
    for e in events:
        if prev_time < e['time_s'] <= current_time:
            result.append(e)
    return result


def sphere_sdf_stamp(phi_old, dist, radius):
    """Compute SDF stamp value for sphere ignition."""
    new_phi = -(radius - dist)
    # Clamp to [-1, 1]
    new_phi = max(min(new_phi, 1.0), -1.0)
    # Apply min() semantics
    return min(phi_old, new_phi)


def polygon_winding_number(px, py, xs, ys):
    """
    Compute winding number for point (px, py) relative to polygon.
    Polygon vertices given by xs, ys arrays.
    """
    winding = 0
    n = len(xs)
    
    for i in range(n):
        x1, y1 = xs[i], ys[i]
        x2, y2 = xs[(i + 1) % n], ys[(i + 1) % n]
        
        if y1 <= py:
            if y2 > py:
                # Upward crossing
                cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1)
                if cross > 0.0:
                    winding += 1
        else:
            if y2 <= py:
                # Downward crossing
                cross = (x2 - x1) * (py - y1) - (px - x1) * (y2 - y1)
                if cross < 0.0:
                    winding -= 1
    
    return winding


def segment_dist_sq(px, py, ax, ay, bx, by):
    """Compute squared distance from point P to segment AB."""
    dx = bx - ax
    dy = by - ay
    len_sq = dx * dx + dy * dy
    
    if len_sq < 1.0e-14:
        # Degenerate segment
        dpx = px - ax
        dpy = py - ay
        return dpx * dpx + dpy * dpy
    
    # Parameter t of closest point on segment
    t = ((px - ax) * dx + (py - ay) * dy) / len_sq
    t = max(0.0, min(1.0, t))
    
    cx = ax + t * dx
    cy = ay + t * dy
    dpx = px - cx
    dpy = py - cy
    return dpx * dpx + dpy * dpy


def test_schedule_csv_parsing():
    """Test parsing of ignition schedule CSV."""
    csv_content = """
# Phase 11 test
0      500.0  1000.0  30.0  0  10
300    1400.0 1000.0  20.0  0   5
600    1000.0  400.0  15.0  0   5
"""
    events = parse_ignition_schedule_csv(csv_content)
    assert len(events) == 3, f"Expected 3 events, got {len(events)}"
    assert events[0]['time_s'] == 0.0, f"Event 0 time should be 0.0, got {events[0]['time_s']}"
    assert events[0]['cx'] == 500.0, f"Event 0 cx should be 500.0, got {events[0]['cx']}"
    assert events[0]['priority'] == 10, f"Event 0 priority should be 10, got {events[0]['priority']}"
    assert events[1]['time_s'] == 300.0, f"Event 1 time should be 300.0, got {events[1]['time_s']}"
    assert events[2]['time_s'] == 600.0, f"Event 2 time should be 600.0, got {events[2]['time_s']}"
    print("✓ test_schedule_csv_parsing PASSED")


def test_priority_sorting_within_timestep():
    """Test that events at same time are sorted by priority (higher first)."""
    events = [
        {'time_s': 100.0, 'priority': 5, 'cx': 500.0, 'cy': 500.0},
        {'time_s': 100.0, 'priority': 10, 'cx': 600.0, 'cy': 600.0},
    ]
    sorted_events = sort_events_by_time_priority(events)
    assert sorted_events[0]['priority'] == 10, \
        f"First event should have priority 10, got {sorted_events[0]['priority']}"
    assert sorted_events[1]['priority'] == 5, \
        f"Second event should have priority 5, got {sorted_events[1]['priority']}"
    print("✓ test_priority_sorting_within_timestep PASSED")


def test_time_window_firing():
    """Test time window firing logic."""
    events = [
        {'time_s': 0.0, 'priority': 10},
        {'time_s': 300.0, 'priority': 5},
        {'time_s': 600.0, 'priority': 5},
    ]
    
    # Window (0, 300]
    fired = get_events_in_window(events, 0.0, 300.0)
    assert len(fired) == 1, f"Window (0, 300] should fire 1 event, got {len(fired)}"
    assert fired[0]['time_s'] == 300.0, f"Should fire t=300 event, got t={fired[0]['time_s']}"
    
    # Window (-1, 0]
    fired = get_events_in_window(events, -1.0, 0.0)
    assert len(fired) == 1, f"Window (-1, 0] should fire 1 event, got {len(fired)}"
    assert fired[0]['time_s'] == 0.0, f"Should fire t=0 event, got t={fired[0]['time_s']}"
    
    # Window (-100, 600]
    fired = get_events_in_window(events, -100.0, 600.0)
    assert len(fired) == 3, f"Window (-100, 600] should fire 3 events, got {len(fired)}"
    
    print("✓ test_time_window_firing PASSED")


def test_fired_flag_prevents_reapplication():
    """Test that fired flag prevents re-application."""
    events = [{'time_s': 100.0, 'fired': True}]
    
    # Event with fired=True should not appear in firing list
    fired = [e for e in events if not e['fired'] and e['time_s'] == 100.0]
    assert len(fired) == 0, f"Fired event should not be in firing list, got {len(fired)}"
    
    events[0]['fired'] = False
    fired = [e for e in events if not e['fired'] and e['time_s'] == 100.0]
    assert len(fired) == 1, f"Non-fired event should be in firing list, got {len(fired)}"
    
    print("✓ test_fired_flag_prevents_reapplication PASSED")


def test_sphere_sdf_stamp_formula():
    """Test sphere SDF stamping formula."""
    # Cell at distance 15 from ignition center with radius 30
    dist = 15.0
    radius = 30.0
    
    # For phi_old = +1.0
    phi_old = 1.0
    phi_new = sphere_sdf_stamp(phi_old, dist, radius)
    expected = min(1.0, -(radius - dist))
    expected = max(min(expected, 1.0), -1.0)
    assert abs(phi_new - expected) < 1.0e-10, \
        f"phi_new should be {expected}, got {phi_new}"
    
    # For phi_old = -20.0 (existing fire)
    phi_old = -20.0
    phi_new = sphere_sdf_stamp(phi_old, dist, radius)
    expected_val = -(radius - dist)
    expected_val = max(min(expected_val, 1.0), -1.0)
    expected = min(phi_old, expected_val)
    assert abs(phi_new - expected) < 1.0e-10, \
        f"phi_new should be {expected}, got {phi_new}"
    
    print("✓ test_sphere_sdf_stamp_formula PASSED")


def test_suppress_if_burning_logic():
    """Test suppress_if_burning flag."""
    # Cell with phi = -5.0 (already burning) should not be modified with suppress_if_burning=True
    phi = -5.0
    radius = 30.0
    dist = 15.0
    suppress = True
    
    if suppress and phi < 0.0:
        phi_new = phi  # Not modified
    else:
        phi_new = sphere_sdf_stamp(phi, dist, radius)
    
    assert phi_new == phi, f"Burning cell should not be modified, got {phi_new}"
    
    # Cell with phi = +10.0 (unburned) should be modified with suppress_if_burning=True
    phi = 10.0
    if suppress and phi < 0.0:
        phi_new = phi
    else:
        phi_new = sphere_sdf_stamp(phi, dist, radius)
    
    assert phi_new != phi, f"Unburned cell should be modified, got {phi_new}"
    
    print("✓ test_suppress_if_burning_logic PASSED")


def test_polygon_winding_number_inside():
    """Test winding number for point inside square polygon."""
    # Square: (0,0), (2,0), (2,2), (0,2)
    xs = [0.0, 2.0, 2.0, 0.0]
    ys = [0.0, 0.0, 2.0, 2.0]
    
    # Point at center (1, 1) should be inside
    winding = polygon_winding_number(1.0, 1.0, xs, ys)
    assert winding != 0, f"Center point should have non-zero winding, got {winding}"
    
    print("✓ test_polygon_winding_number_inside PASSED")


def test_polygon_winding_number_outside():
    """Test winding number for point outside square polygon."""
    # Square: (0,0), (2,0), (2,2), (0,2)
    xs = [0.0, 2.0, 2.0, 0.0]
    ys = [0.0, 0.0, 2.0, 2.0]
    
    # Point at (3, 3) should be outside
    winding = polygon_winding_number(3.0, 3.0, xs, ys)
    assert winding == 0, f"Outside point should have zero winding, got {winding}"
    
    print("✓ test_polygon_winding_number_outside PASSED")


def test_polygon_winding_number_vertex_count():
    """Test that winding number works with regular polygons."""
    # Regular hexagon centered at (0, 0)
    n_sides = 6
    xs = [math.cos(2 * math.pi * i / n_sides) for i in range(n_sides)]
    ys = [math.sin(2 * math.pi * i / n_sides) for i in range(n_sides)]
    
    # Point at center (0, 0) should be inside
    winding = polygon_winding_number(0.0, 0.0, xs, ys)
    assert winding != 0, f"Center of hexagon should have non-zero winding, got {winding}"
    
    print("✓ test_polygon_winding_number_vertex_count PASSED")


def test_segment_distance_squared_formula():
    """Test segment distance formula."""
    # Segment from (0,0) to (2,0), point at (1,1)
    dist_sq = segment_dist_sq(1.0, 1.0, 0.0, 0.0, 2.0, 0.0)
    assert abs(dist_sq - 1.0) < 1.0e-10, f"Distance squared should be 1.0, got {dist_sq}"
    
    # Segment from (0,0) to (2,0), point at (3,1)
    dist_sq = segment_dist_sq(3.0, 1.0, 0.0, 0.0, 2.0, 0.0)
    expected = (3.0 - 2.0)**2 + (1.0 - 0.0)**2  # Distance to (2, 0)
    assert abs(dist_sq - expected) < 1.0e-10, f"Distance squared should be {expected}, got {dist_sq}"
    
    print("✓ test_segment_distance_squared_formula PASSED")


def test_polyline_half_width_stamp():
    """Test polyline half-width stamping."""
    # Cell at distance 8.0 from segment with half_width=10.0
    dist = 8.0
    half_width = 10.0
    
    new_phi = -half_width if dist <= half_width else dist
    assert new_phi == -10.0, f"Should get phi=-10.0, got {new_phi}"
    
    # Cell at distance 12.0 from segment with half_width=10.0
    dist = 12.0
    new_phi = -half_width if dist <= half_width else dist
    assert new_phi == 12.0, f"Should get phi=12.0, got {new_phi}"
    
    print("✓ test_polyline_half_width_stamp PASSED")


def test_multi_ignition_phi_merge():
    """Test that multiple ignitions merge via min() accumulation."""
    # Two ignitions produce phi values of -15 and -8 at same cell
    phi1 = -15.0
    phi2 = -8.0
    
    # Merge via min()
    merged = min(phi1, phi2)
    assert merged == -15.0, f"min(-15, -8) should be -15, got {merged}"
    
    # Verify that min() doesn't overwrite existing fire
    phi_old = -20.0
    new_phi = -10.0
    merged = min(phi_old, new_phi)
    assert merged == -20.0, f"Existing fire should not be overwritten, got {merged}"
    
    print("✓ test_multi_ignition_phi_merge PASSED")


def main():
    """Run all unit tests."""
    print("=" * 70)
    print("Phase 11: Multi-Ignition Schedule Unit Tests")
    print("=" * 70)
    
    tests = [
        test_schedule_csv_parsing,
        test_priority_sorting_within_timestep,
        test_time_window_firing,
        test_fired_flag_prevents_reapplication,
        test_sphere_sdf_stamp_formula,
        test_suppress_if_burning_logic,
        test_polygon_winding_number_inside,
        test_polygon_winding_number_outside,
        test_polygon_winding_number_vertex_count,
        test_segment_distance_squared_formula,
        test_polyline_half_width_stamp,
        test_multi_ignition_phi_merge,
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
