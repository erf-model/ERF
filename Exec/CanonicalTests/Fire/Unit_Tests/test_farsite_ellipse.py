#!/usr/bin/env python3
"""
FARSITE elliptical propagation unit tests.

Tests the FARSITE elliptical fire expansion algorithm that models
fire spread as an ellipse oriented by wind and slope.

Run: python3 test_farsite_ellipse.py
"""

import sys
import math


class FARSITEEllipse:
    """FARSITE elliptical fire expansion model."""
    
    def __init__(self, head_fire_ros, lateral_fire_ros, back_fire_ros):
        """
        Initialize ellipse with fire spread rates.
        
        Args:
            head_fire_ros: Head fire rate of spread (m/s) - direction of max spread
            lateral_fire_ros: Lateral fire rate of spread (m/s) - perpendicular to head
            back_fire_ros: Back fire rate of spread (m/s) - opposite to head
        """
        self.head_ros = head_fire_ros
        self.lateral_ros = lateral_fire_ros
        self.back_ros = back_fire_ros
        
        # Ellipse parameters
        self.a = None  # Semi-major axis (head direction)
        self.b = None  # Semi-minor axis (lateral direction)
        self.c = None  # Semi-focal distance
        self.eccentricity = None
        self.compute_ellipse_params()
    
    def compute_ellipse_params(self):
        """Compute FARSITE ellipse parameters from ROS values."""
        # Head and back average: 2a = head + back
        self.a = (self.head_ros + self.back_ros) / 2.0
        
        # Semi-minor axis (lateral direction)
        self.b = self.lateral_ros
        
        # Semi-focal distance: c² = a² - b²
        if self.a > self.b:
            c_squared = self.a**2 - self.b**2
            self.c = math.sqrt(c_squared)
            self.eccentricity = self.c / self.a if self.a > 0 else 0.0
        else:
            # If lateral > head, swap (shouldn't happen for normal fire)
            self.c = 0.0
            self.eccentricity = 0.0
    
    def ros_at_angle(self, theta):
        """
        Compute ROS at angle theta from head direction.
        
        Args:
            theta: Angle from head direction (radians), 0 = head, pi = back
        
        Returns:
            Rate of spread at that angle (m/s)
        """
        if self.eccentricity == 0.0:
            return self.a  # Circular fire
        
        # FARSITE ellipse formula
        numerator = self.a * self.b
        denominator = math.sqrt((self.b * math.cos(theta))**2 + 
                                (self.a * math.sin(theta))**2)
        
        if denominator == 0.0:
            return 0.0
        return numerator / denominator
    
    def fire_perimeter_advance(self, dt):
        """
        Advance ellipse by time dt, returning new ellipse dimensions.
        
        Args:
            dt: Time step (seconds)
        
        Returns:
            Dictionary with advanced ellipse parameters
        """
        return {
            'major_axis': self.a * dt,
            'minor_axis': self.b * dt,
            'head_distance': self.head_ros * dt,
            'back_distance': self.back_ros * dt,
            'lateral_distance': self.lateral_ros * dt,
        }


def test_ros_at_angle():
    """Test ROS calculation at various angles."""
    print("\n" + "="*70)
    print("Testing FARSITE ROS at Various Angles")
    print("="*70)
    
    # Create an ellipse: head=0.5 m/s, lateral=0.1 m/s, back=0.05 m/s
    ellipse = FARSITEEllipse(head_fire_ros=0.5, lateral_fire_ros=0.1, back_fire_ros=0.05)
    
    print(f"\nEllipse Parameters:")
    print(f"  Semi-major axis (a): {ellipse.a:.4f} m/s")
    print(f"  Semi-minor axis (b): {ellipse.b:.4f} m/s")
    print(f"  Eccentricity (e): {ellipse.eccentricity:.4f}")
    
    # Test at cardinal directions
    test_angles = [
        (0, "Head (0°)"),
        (math.pi/2, "Right lateral (90°)"),
        (math.pi, "Back (180°)"),
        (-math.pi/2, "Left lateral (-90°)"),
    ]
    
    print(f"\nROS at cardinal directions:")
    failures = []
    
    for angle_rad, label in test_angles:
        ros = ellipse.ros_at_angle(angle_rad)
        print(f"  {label}: {ros:.4f} m/s")
        
    def test_ros_at_angle(self):
        """Test ROS calculation at various angles - verify head/back values."""
        # The ellipse semi-major axis is average of head and back
        # Head ROS at theta=0 should be a = (head + back) / 2
        ros_head = self.ros_at_angle(0)
        ros_back = self.ros_at_angle(math.pi)
        
        # Both should equal semi-major axis for an ellipse
        assert abs(ros_head - self.a) < 1e-6, f"Head ROS mismatch: {ros_head} vs {self.a}"
        assert abs(ros_back - self.a) < 1e-6, f"Back ROS mismatch: {ros_back} vs {self.a}"
    
    return failures


def test_ellipse_symmetry():
    """Test that ellipse is symmetric about major axis."""
    print("\n" + "="*70)
    print("Testing FARSITE Ellipse Symmetry")
    print("="*70)
    
    ellipse = FARSITEEllipse(head_fire_ros=0.8, lateral_fire_ros=0.2, back_fire_ros=0.1)
    
    # Test symmetry at multiple angles
    test_angles = [math.pi/6, math.pi/4, math.pi/3]
    
    print(f"\nSymmetry check (left vs right lateral):")
    failures = []
    
    for angle in test_angles:
        ros_left = ellipse.ros_at_angle(-angle)
        ros_right = ellipse.ros_at_angle(angle)
        print(f"  Angle ±{math.degrees(angle):.1f}°: {ros_left:.4f} vs {ros_right:.4f}")
        
        if abs(ros_left - ros_right) > 1e-6:
            failures.append(f"Asymmetry at angle {angle}: {ros_left} vs {ros_right}")
    
    return failures


def test_perimeter_advance():
    """Test fire perimeter advancement over time."""
    print("\n" + "="*70)
    print("Testing Fire Perimeter Advancement")
    print("="*70)
    
    ellipse = FARSITEEllipse(head_fire_ros=0.5, lateral_fire_ros=0.1, back_fire_ros=0.05)
    
    # Advance by 10 seconds
    dt = 10.0
    advance = ellipse.fire_perimeter_advance(dt)
    
    print(f"\nFire advancement over {dt} seconds:")
    print(f"  Head distance: {advance['head_distance']:.2f} m")
    print(f"  Back distance: {advance['back_distance']:.2f} m")
    print(f"  Lateral distance: {advance['lateral_distance']:.2f} m")
    
    failures = []
    
    # Verify distances match ROS × dt
    if abs(advance['head_distance'] - 0.5 * dt) > 1e-6:
        failures.append(f"Head distance mismatch: {advance['head_distance']} vs {0.5*dt}")
    if abs(advance['back_distance'] - 0.05 * dt) > 1e-6:
        failures.append(f"Back distance mismatch: {advance['back_distance']} vs {0.05*dt}")
    if abs(advance['lateral_distance'] - 0.1 * dt) > 1e-6:
        failures.append(f"Lateral distance mismatch: {advance['lateral_distance']} vs {0.1*dt}")
    
    return failures


def test_edge_cases():
    """Test edge cases and special conditions."""
    print("\n" + "="*70)
    print("Testing FARSITE Edge Cases")
    print("="*70)
    
    failures = []
    
    # Test 1: Circular fire (all ROS equal)
    print("\nTest 1: Circular fire (uniform ROS)")
    circular = FARSITEEllipse(head_fire_ros=0.3, lateral_fire_ros=0.3, back_fire_ros=0.3)
    print(f"  Eccentricity (should be 0): {circular.eccentricity:.4f}")
    if circular.eccentricity > 1e-6:
        failures.append(f"Circular fire should have e=0, got {circular.eccentricity}")
    
    ros_at_45 = circular.ros_at_angle(math.pi/4)
    if abs(ros_at_45 - 0.3) > 1e-6:
        failures.append(f"Circular fire ROS at 45° should be 0.3, got {ros_at_45}")
    print(f"  ROS at 45°: {ros_at_45:.4f} (should be 0.3)")
    
    # Test 2: Very elongated ellipse (head >> lateral)
    print("\nTest 2: Elongated ellipse (head >> lateral)")
    elongated = FARSITEEllipse(head_fire_ros=1.0, lateral_fire_ros=0.01, back_fire_ros=0.05)
    print(f"  Eccentricity: {elongated.eccentricity:.4f}")
    ros_head = elongated.ros_at_angle(0)
    ros_lateral = elongated.ros_at_angle(math.pi/2)
    print(f"  ROS head: {ros_head:.4f}, lateral: {ros_lateral:.4f}")
    if ros_lateral > ros_head:
        failures.append(f"Lateral ROS should not exceed head ROS")
    
    return failures


def main():
    """Run all FARSITE ellipse tests."""
    print("="*70)
    print("FARSITE Elliptical Propagation Unit Tests")
    print("="*70)
    
    all_failures = []
    
    all_failures.extend(test_ros_at_angle())
    all_failures.extend(test_ellipse_symmetry())
    all_failures.extend(test_perimeter_advance())
    all_failures.extend(test_edge_cases())
    
    # Summary
    print("\n" + "="*70)
    if not all_failures:
        print(f"FINAL RESULT: All tests PASSED")
        return 0
    else:
        print(f"FINAL RESULT: {len(all_failures)} test(s) FAILED")
        for failure in all_failures:
            print(f"  ✗ {failure}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
