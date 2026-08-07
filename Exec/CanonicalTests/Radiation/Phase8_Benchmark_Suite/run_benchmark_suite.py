#!/usr/bin/env python3
"""
Phase 8 Benchmark Suite: Main Driver Script

Orchestrates the execution of all benchmark cases, extracts metrics,
validates against tolerances, and generates reports.

Usage:
    python3 run_benchmark_suite.py [--verbose] [--no-run]

Options:
    --verbose:  Print detailed output during validation
    --no-run:   Skip test execution; only validate existing outputs
"""

import sys
import json
import argparse
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any

import benchmark_tolerances as tol
from benchmark_config import CASES, list_all_cases, BenchmarkCase
from check_benchmark_metrics import validate_case


class BenchmarkSuite:
    """Orchestrates benchmark execution and reporting."""
    
    def __init__(self, base_dir: Path, verbose: bool = False):
        """
        Initialize benchmark suite.
        
        Args:
            base_dir: Base directory for the suite (Phase8_Benchmark_Suite/)
            verbose: Enable verbose output
        """
        self.base_dir = base_dir
        self.verbose = verbose
        self.results: List[Dict[str, Any]] = []
        self.timestamp = datetime.now().isoformat()
    
    def print_header(self):
        """Print suite header."""
        print("\n" + "=" * 80)
        print("PHASE 8: VALIDATION & BENCHMARKING SUITE FOR TWOSTREAM RADIATION")
        print("=" * 80)
        print(f"Timestamp: {self.timestamp}")
        print(f"Base directory: {self.base_dir}")
        print()
    
    def print_case_header(self, case: BenchmarkCase):
        """Print case header."""
        print("-" * 80)
        print(f"CASE: {case.short_name.upper()}")
        print(f"  Name: {case.name}")
        print(f"  Physics: {case.physics_description}")
        print(f"  Diagnostics mode: {case.diag_callsite_mode}")
        print(f"  Expected rows: {case.expected_diag_rows} (mode={case.diag_callsite_mode})")
        print()
    
    def validate_case_output(self, case: BenchmarkCase) -> bool:
        """
        Validate a single case's output.
        
        Args:
            case: BenchmarkCase definition
        
        Returns:
            True if case passes, False otherwise
        """
        self.print_case_header(case)
        
        # Resolve case directory
        case_dir = self.base_dir.parent / case.base_case_dir
        if not case_dir.exists():
            print(f"  [ERROR] Case directory not found: {case_dir}")
            self.results.append({
                "case": case.short_name,
                "case_name": case.name,
                "is_pass": False,
                "errors": [f"Case directory not found: {case_dir}"],
                "warnings": [],
                "metrics": {},
            })
            return False
        
        if self.verbose:
            print(f"  Case directory: {case_dir}")
        
        # Validate case output
        is_pass, summary = validate_case(case, case_dir)
        self.results.append(summary)
        
        # Print results
        if is_pass:
            print(f"  [PASS] Case validation succeeded")
        else:
            print(f"  [FAIL] Case validation failed")
            if summary.get("errors"):
                for error in summary["errors"]:
                    print(f"    - {error}")
        
        if summary.get("warnings") and self.verbose:
            print(f"  [WARNINGS]")
            for warning in summary["warnings"]:
                print(f"    - {warning}")
        
        if self.verbose and summary.get("metrics"):
            print(f"  [METRICS]")
            for key, val in summary["metrics"].items():
                if isinstance(val, float):
                    print(f"    {key}: {val:.6e}")
                else:
                    print(f"    {key}: {val}")
        
        print()
        return is_pass
    
    def run_suite(self) -> bool:
        """
        Run the entire benchmark suite validation.
        
        Returns:
            True if all cases pass, False if any case fails
        """
        self.print_header()
        tol.print_tolerance_summary()
        
        all_pass = True
        for short_name in list_all_cases():
            case = CASES[short_name]
            case_pass = self.validate_case_output(case)
            if not case_pass:
                all_pass = False
        
        return all_pass
    
    def generate_json_report(self, output_path: Path) -> None:
        """
        Generate JSON report of benchmark results.
        
        Args:
            output_path: Path to output JSON file
        """
        report = {
            "timestamp": self.timestamp,
            "suite": "Phase8_TwoStream_Radiation_Benchmark",
            "overall_pass": all(r.get("is_pass", False) for r in self.results),
            "total_cases": len(self.results),
            "passed_cases": sum(1 for r in self.results if r.get("is_pass", False)),
            "failed_cases": sum(1 for r in self.results if not r.get("is_pass", False)),
            "cases": self.results,
        }
        
        with output_path.open("w") as f:
            json.dump(report, f, indent=2)
        
        if self.verbose:
            print(f"JSON report saved to: {output_path}")
    
    def generate_markdown_report(self, output_path: Path) -> None:
        """
        Generate Markdown report of benchmark results.
        
        Args:
            output_path: Path to output Markdown file
        """
        lines = []
        lines.append("# Phase 8 Benchmark Suite Results\n")
        lines.append(f"**Timestamp:** {self.timestamp}\n")
        
        # Summary
        total = len(self.results)
        passed = sum(1 for r in self.results if r.get("is_pass", False))
        failed = total - passed
        
        lines.append("## Summary\n")
        lines.append(f"- **Total Cases:** {total}")
        lines.append(f"- **Passed:** {passed}")
        lines.append(f"- **Failed:** {failed}")
        lines.append(f"- **Overall Status:** {'✅ PASS' if failed == 0 else '❌ FAIL'}\n")
        
        # Case table
        lines.append("## Case Results\n")
        lines.append("| Case | Name | Status | Errors |\n")
        lines.append("|------|------|--------|--------|\n")
        
        for result in self.results:
            case_name = result.get("case", "unknown")
            display_name = result.get("case_name", case_name)
            status = "✅ PASS" if result.get("is_pass", False) else "❌ FAIL"
            error_count = len(result.get("errors", []))
            error_summary = f"{error_count} error(s)" if error_count > 0 else "None"
            
            lines.append(f"| `{case_name}` | {display_name} | {status} | {error_summary} |\n")
        
        lines.append()
        
        # Detailed results
        lines.append("## Detailed Results\n")
        
        for result in self.results:
            case_name = result.get("case", "unknown")
            display_name = result.get("case_name", case_name)
            status = "✅ PASS" if result.get("is_pass", False) else "❌ FAIL"
            
            lines.append(f"### {case_name}: {display_name} [{status}]\n")
            
            if result.get("errors"):
                lines.append("**Errors:**\n")
                for error in result["errors"]:
                    lines.append(f"- {error}\n")
            else:
                lines.append("**Status:** No errors\n")
            
            if result.get("metrics"):
                lines.append("\n**Key Metrics:**\n")
                for key, val in result["metrics"].items():
                    if isinstance(val, bool):
                        lines.append(f"- `{key}`: {val}\n")
                    elif isinstance(val, float):
                        lines.append(f"- `{key}`: {val:.6e}\n")
                    else:
                        lines.append(f"- `{key}`: {val}\n")
            
            lines.append()
        
        # Tolerance reference
        lines.append("## Tolerance Configuration\n")
        lines.append(f"- SW_TOA relative error: {tol.SW_TOA_RELATIVE_TOL_PCT}%\n")
        lines.append(f"- SW_surface relative error: {tol.SW_SURFACE_RELATIVE_TOL_PCT}%\n")
        lines.append(f"- LW_net_surface relative error: {tol.LW_NET_SURFACE_RELATIVE_TOL_PCT}%\n")
        lines.append(f"- Heating CV upper bound: {tol.HEATING_CV_UPPER_BOUND * 100}%\n")
        lines.append(f"- Row count tolerance: ±{tol.ROW_COUNT_ABS_TOL}\n")
        
        with output_path.open("w") as f:
            f.writelines(lines)
        
        if self.verbose:
            print(f"Markdown report saved to: {output_path}")
    
    def finalize(self) -> int:
        """
        Finalize suite: generate reports and return exit code.
        
        Returns:
            0 if all cases pass, 1 if any case fails
        """
        all_pass = all(r.get("is_pass", False) for r in self.results)
        
        # Generate reports
        json_report = self.base_dir / "benchmark_summary.json"
        md_report = self.base_dir / "benchmark_summary.md"
        
        self.generate_json_report(json_report)
        self.generate_markdown_report(md_report)
        
        print("=" * 80)
        print("BENCHMARK SUITE SUMMARY")
        print("=" * 80)
        
        total = len(self.results)
        passed = sum(1 for r in self.results if r.get("is_pass", False))
        failed = total - passed
        
        print(f"Total cases: {total}")
        print(f"Passed: {passed}")
        print(f"Failed: {failed}")
        print()
        print(f"Reports generated:")
        print(f"  - JSON: {json_report}")
        print(f"  - Markdown: {md_report}")
        print()
        
        if all_pass:
            print("✅ ALL BENCHMARKS PASSED")
            print("=" * 80)
            return 0
        else:
            print("❌ SOME BENCHMARKS FAILED")
            print("=" * 80)
            return 1


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Phase 8 Benchmark Suite for TwoStream Radiation"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose output"
    )
    parser.add_argument(
        "--no-run",
        action="store_true",
        help="Skip test execution; only validate existing outputs"
    )
    
    args = parser.parse_args()
    
    # Base directory is the directory containing this script
    base_dir = Path(__file__).parent.resolve()
    
    # Create and run suite
    suite = BenchmarkSuite(base_dir, verbose=args.verbose)
    suite.run_suite()
    
    # Finalize and exit
    exit_code = suite.finalize()
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
