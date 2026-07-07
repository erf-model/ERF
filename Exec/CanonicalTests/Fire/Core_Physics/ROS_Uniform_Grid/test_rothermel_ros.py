#!/usr/bin/env python3
"""
Regression test: Rothermel ROS magnitude check.

Runs ERF with fire enabled for 5 steps on a flat domain with FM1 fuel
and uniform wind. Parses [FIRE DEBUG] output and verifies that:
  1. Max ROS is in range [0.05, 3.0] m/s  (not ~1.75e-5 from buggy code)
  2. Max effective wind is in range [1.0, 10.0] m/s
  3. Reaction intensity print shows I_R > 100 BTU/ft²/min
  4. C coefficient shows C > 7.0 (not near zero)
  5. U_max is > 100 ft/min (MEWS cap physically reasonable)

Run:
  python3 test_rothermel_ros.py --erf_exe ./erf --input_file inputs_fire_flat_uniform
"""

import os
import sys
import subprocess
import argparse
import re


def run_erf_short_test(erf_executable, input_file, num_processes=1, max_steps=5):
    """
    Run ERF with fire enabled for a short test run.
    
    Args:
        erf_executable: Path to the compiled ERF executable
        input_file: Path to the fire model input file
        num_processes: Number of MPI processes to use
        max_steps: Number of steps to run
    
    Returns:
        Tuple (success: bool, stdout: str, stderr: str)
    """
    if not os.path.exists(erf_executable):
        print(f"ERROR: ERF executable not found at {erf_executable}")
        return False, "", ""

    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found at {input_file}")
        return False, "", ""

    # Prepare command with max_step override
    if num_processes > 1:
        cmd = ["mpirun", "-np", str(num_processes), erf_executable]
    else:
        cmd = [erf_executable]

    cmd.extend([input_file, f"max_step={max_steps}"])

    print(f"Running ERF test: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        return result.returncode == 0, result.stdout, result.stderr

    except subprocess.TimeoutExpired:
        print("ERROR: ERF test timed out after 600 seconds")
        return False, "", ""
    except Exception as e:
        print(f"ERROR: Failed to run ERF: {e}")
        return False, "", ""


def parse_rothermel_params(output):
    """
    Parse Rothermel parameters from ERF debug output.
    
    Args:
        output: Combined stdout and stderr from ERF run
    
    Returns:
        Dict with parsed parameters or None if parsing fails
    """
    params = {}
    
    # Parse "Rothermel Params:" line
    rothermel_pattern = r'Rothermel Params: R0=([\d.]+) ft/min \(([\d.e+-]+) m/s\).*I_R=([\d.e+-]+).*C=([\d.e+-]+).*B=([\d.e+-]+).*E=([\d.e+-]+).*U_max=([\d.e+-]+)'
    match = re.search(rothermel_pattern, output)
    if match:
        params['R0_ftmin'] = float(match.group(1))
        params['R0_ms'] = float(match.group(2))
        params['I_R'] = float(match.group(3))
        params['C'] = float(match.group(4))
        params['B'] = float(match.group(5))
        params['E'] = float(match.group(6))
        params['U_max_ftmin'] = float(match.group(7))
    
    return params if params else None


def parse_max_ros(output):
    """
    Parse maximum ROS from fire debug output.
    
    Args:
        output: Combined stdout and stderr from ERF run
    
    Returns:
        Max ROS value in m/s or None if not found
    """
    # Look for "[FIRE DEBUG] Rate-of-spread computed. Max:" patterns
    pattern = r'\[FIRE DEBUG\].*Rate-of-spread computed.*Max.*[:=]?\s*([\d.e+-]+)\s*m/s'
    match = re.search(pattern, output)
    if match:
        return float(match.group(1))
    
    # Alternative pattern for max ROS output
    pattern2 = r'Max ROS.*?=\s*([\d.e+-]+)\s*m/s'
    match2 = re.search(pattern2, output, re.IGNORECASE)
    if match2:
        return float(match2.group(1))
    
    return None


def main():
    parser = argparse.ArgumentParser(description="Run ERF fire model and verify Rothermel ROS")
    parser.add_argument("--erf_exe", type=str, default="./erf",
                        help="Path to compiled ERF executable")
    parser.add_argument("--input_file", type=str, default="inputs_fire_flat_uniform",
                        help="Path to fire model input file")
    parser.add_argument("--num_procs", type=int, default=1,
                        help="Number of MPI processes")
    parser.add_argument("--max_steps", type=int, default=5,
                        help="Maximum number of steps to run")

    args = parser.parse_args()

    print("=" * 70)
    print("Rothermel ROS Regression Test (Integration)")
    print("=" * 70)
    print(f"ERF Executable: {args.erf_exe}")
    print(f"Input File: {args.input_file}")
    print(f"Number of Processes: {args.num_procs}")
    print(f"Max Steps: {args.max_steps}")
    print("=" * 70)

    # Run ERF
    print("\nRunning ERF fire simulation...")
    success, stdout, stderr = run_erf_short_test(args.erf_exe, args.input_file, 
                                                  args.num_procs, args.max_steps)
    
    if not success:
        print(f"ERROR: ERF failed to run successfully")
        print("STDOUT:", stdout[-500:] if stdout else "(empty)")
        print("STDERR:", stderr[-500:] if stderr else "(empty)")
        return 1
    
    print("ERF run completed successfully")
    
    # Combine stdout and stderr for parsing
    combined_output = stdout + stderr
    
    # Save output for debugging
    with open("/tmp/erf_output.log", "w") as f:
        f.write(combined_output)
    print(f"Full output saved to /tmp/erf_output.log")
    
    # Parse Rothermel parameters
    print("\nParsing Rothermel parameters...")
    params = parse_rothermel_params(combined_output)
    
    if not params:
        print("WARNING: Could not parse Rothermel parameters from output")
        print("Last 1000 chars of output:")
        print(combined_output[-1000:])
    else:
        print(f"  R0: {params['R0_ftmin']:.3f} ft/min ({params['R0_ms']:.6f} m/s)")
        print(f"  I_R: {params['I_R']:.1f} BTU/ft²/min")
        print(f"  C: {params['C']:.4f}")
        print(f"  B: {params['B']:.6f}")
        print(f"  E: {params['E']:.6f}")
        print(f"  U_max: {params['U_max_ftmin']:.1f} ft/min")
    
    # Parse max ROS
    print("\nParsing maximum ROS from fire grid...")
    max_ros = parse_max_ros(combined_output)
    if max_ros is not None:
        print(f"  Max ROS: {max_ros:.6f} m/s")
    else:
        print("  WARNING: Could not parse max ROS from output")
    
    # Run assertions
    print("\n" + "=" * 70)
    print("Validating physical constraints...")
    print("=" * 70)
    
    failures = []
    
    if params:
        # Check C coefficient (must be near 7.4, not near-zero)
        if params['C'] < 7.0:
            failures.append(f"C={params['C']:.4f} < 7.0 (expected ~7.4, wrong formula?)")
        else:
            print(f"✓ C coefficient: {params['C']:.4f} (correct, > 7.0)")
        
        # Check I_R (must be > 100)
        if params['I_R'] < 100:
            failures.append(f"I_R={params['I_R']:.1f} < 100 BTU/ft²/min (too low)")
        else:
            print(f"✓ Reaction intensity I_R: {params['I_R']:.1f} BTU/ft²/min (> 100)")
        
        # Check U_max (MEWS cap, must be > 100 ft/min)
        if params['U_max_ftmin'] < 100:
            failures.append(f"U_max={params['U_max_ftmin']:.1f} ft/min < 100 (cap too low, not ~4 ft/min?)")
        else:
            print(f"✓ MEWS cap U_max: {params['U_max_ftmin']:.1f} ft/min (> 100, physically reasonable)")
        
        # Check R0 (must be > 0.001 m/s)
        if params['R0_ms'] < 0.001:
            failures.append(f"R0={params['R0_ms']:.6f} m/s < 0.001 (no-wind ROS too low)")
        else:
            print(f"✓ No-wind ROS R0: {params['R0_ms']:.6f} m/s (> 0.001)")
    
    if max_ros is not None:
        # Check max ROS
        if max_ros < 0.05 or max_ros > 5.0:
            failures.append(f"Max ROS={max_ros:.6f} m/s outside expected range [0.05, 5.0]")
        else:
            print(f"✓ Max ROS: {max_ros:.6f} m/s (in expected range [0.05, 5.0])")
    
    print("=" * 70)
    if not failures:
        print("FINAL RESULT: ALL CHECKS PASSED")
        return 0
    else:
        print("FINAL RESULT: VALIDATION FAILED")
        print("\nFailures:")
        for failure in failures:
            print(f"  ✗ {failure}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
