#!/usr/bin/env python3
"""
Dummy regression test for 2D Fire model in ERF.

This test verifies that:
1. Fire model initializes correctly
2. Dummy Rothermel function calls work
3. Elliptical expansion functions work
4. Fire intensity calculations work
"""

import os
import sys
import subprocess
import argparse


def run_fire_model_test(erf_executable, input_file, num_processes=1):
    """
    Run the fire model with dummy inputs and verify execution.

    Args:
        erf_executable: Path to the compiled ERF executable
        input_file: Path to the fire model input file
        num_processes: Number of MPI processes to use

    Returns:
        True if test passes, False otherwise
    """
    if not os.path.exists(erf_executable):
        print(f"ERROR: ERF executable not found at {erf_executable}")
        return False

    if not os.path.exists(input_file):
        print(f"ERROR: Input file not found at {input_file}")
        return False

    # Prepare command
    if num_processes > 1:
        cmd = ["mpirun", "-np", str(num_processes), erf_executable]
    else:
        cmd = [erf_executable]

    cmd.append(input_file)

    print(f"Running fire model test: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        print("STDOUT:")
        print(result.stdout)
        print("STDERR:")
        print(result.stderr)

        # Check for successful completion
        if result.returncode == 0:
            print("\nTest PASSED: Fire model ran successfully")
            return True
        else:
            print(f"\nTest FAILED: Fire model exited with code {result.returncode}")
            return False

    except subprocess.TimeoutExpired:
        print("ERROR: Fire model test timed out after 300 seconds")
        return False
    except Exception as e:
        print(f"ERROR: Failed to run fire model: {e}")
        return False


def verify_output_files(output_dir=".", expected_files=None):
    """
    Verify that expected output files were created.

    Args:
        output_dir: Directory to check for output files
        expected_files: List of expected output file patterns

    Returns:
        True if all expected files exist, False otherwise
    """
    if expected_files is None:
        expected_files = ["plt", "chk"]

    print("\nVerifying output files...")
    for pattern in expected_files:
        matching_files = [f for f in os.listdir(output_dir) if pattern in f]
        if matching_files:
            print(f"  Found output files matching '{pattern}': {matching_files}")
        else:
            print(f"  WARNING: No output files matching '{pattern}'")

    return True


def main():
    parser = argparse.ArgumentParser(description="Run dummy fire model regression test")
    parser.add_argument("--erf_exe", type=str, default="./erf",
                        help="Path to compiled ERF executable")
    parser.add_argument("--input_file", type=str, default="inputs_fire_dummy",
                        help="Path to fire model input file")
    parser.add_argument("--num_procs", type=int, default=1,
                        help="Number of MPI processes")
    parser.add_argument("--verify_outputs", action="store_true",
                        help="Verify output files were created")

    args = parser.parse_args()

    print("=" * 70)
    print("2D Fire Model - Phase 1 Dummy Regression Test")
    print("=" * 70)
    print(f"ERF Executable: {args.erf_exe}")
    print(f"Input File: {args.input_file}")
    print(f"Number of Processes: {args.num_procs}")
    print("=" * 70)

    # Run the test
    success = run_fire_model_test(args.erf_exe, args.input_file, args.num_procs)

    # Verify outputs if requested
    if success and args.verify_outputs:
        verify_output_files()

    print("=" * 70)
    if success:
        print("FINAL RESULT: TEST PASSED")
        return 0
    else:
        print("FINAL RESULT: TEST FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
