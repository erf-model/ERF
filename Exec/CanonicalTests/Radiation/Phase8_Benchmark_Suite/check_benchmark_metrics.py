#!/usr/bin/env python3
"""
Phase 8 Benchmark Suite: Metrics Extraction and Validation

This module reads radiation diagnostic CSV files, computes aggregate metrics,
and validates them against configured tolerances. Returns pass/fail results
for each case.
"""

import csv
import math
import statistics
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from collections import Counter

import benchmark_tolerances as tol
from benchmark_config import BenchmarkCase


class MetricsExtractor:
    """Extracts and validates metrics from a diagnostic CSV file."""
    
    def __init__(self, case: BenchmarkCase, csv_path: Path):
        """
        Initialize metrics extractor for a case.
        
        Args:
            case: BenchmarkCase definition
            csv_path: Path to the diagnostic CSV file
        """
        self.case = case
        self.csv_path = csv_path
        self.rows: List[Dict[str, str]] = []
        self.metrics: Dict[str, Any] = {}
        self.errors: List[str] = []
        self.warnings: List[str] = []
    
    def read_csv(self) -> bool:
        """
        Read and parse the diagnostic CSV file.
        
        Returns:
            True if successful, False if file not found or parse error
        """
        if not self.csv_path.exists():
            self.errors.append(f"Diagnostic CSV not found: {self.csv_path}")
            return False
        
        try:
            with self.csv_path.open("r", newline="") as f:
                reader = csv.DictReader(f)
                self.rows = list(reader)
        except Exception as e:
            self.errors.append(f"Failed to read CSV: {e}")
            return False
        
        if not self.rows:
            self.errors.append("Diagnostic CSV is empty")
            return False
        
        return True
    
    def validate_schema(self) -> bool:
        """
        Validate that CSV has all required columns.
        
        Returns:
            True if schema is valid, False otherwise
        """
        if not self.rows:
            return False
        
        missing = [c for c in self.case.required_metrics if c not in self.rows[0]]
        if missing:
            self.errors.append(f"Missing required columns: {missing}")
            return False
        
        return True
    
    def to_float(self, row: Dict, key: str) -> Optional[float]:
        """Safely convert row value to float."""
        try:
            return float(row[key])
        except Exception as e:
            self.errors.append(f"Could not parse float '{key}' from row: {row}")
            return None
    
    def to_int(self, row: Dict, key: str) -> Optional[int]:
        """Safely convert row value to int."""
        try:
            return int(float(row[key]))
        except Exception as e:
            self.errors.append(f"Could not parse int '{key}' from row: {row}")
            return None
    
    def extract_metrics(self) -> bool:
        """
        Extract and compute aggregate metrics from CSV.
        
        Returns:
            True if extraction successful, False if any error occurred
        """
        if not self.rows:
            return False
        
        nrows = len(self.rows)
        self.metrics["row_count"] = nrows
        
        # Row count check
        expected = self.case.expected_diag_rows
        if abs(nrows - expected) <= tol.ROW_COUNT_ABS_TOL:
            self.metrics["row_count_pass"] = True
        else:
            self.metrics["row_count_pass"] = False
            self.errors.append(
                f"Row count mismatch: got {nrows}, expected ~{expected} "
                f"(±{tol.ROW_COUNT_ABS_TOL})"
            )
        
        # Per-step multiplicity check
        steps = []
        for row in self.rows:
            s = self.to_int(row, "step")
            if s is None:
                return False
            steps.append(s)
        
        step_counts = Counter(steps)
        bad_counts = {s: c for s, c in sorted(step_counts.items()) 
                     if c != self.case.rows_per_step}
        
        if bad_counts:
            self.metrics["step_multiplicity_pass"] = False
            self.errors.append(
                f"Step multiplicity mismatch: expected {self.case.rows_per_step} "
                f"rows/step, but got: " + 
                ", ".join(f"step {s}: {c}" for s, c in bad_counts.items())
            )
        else:
            self.metrics["step_multiplicity_pass"] = True
        
        # Call-site validation
        call_sites = [row.get("call_site", "") for row in self.rows]
        has_pre = any("pre" in cs.lower() for cs in call_sites)
        has_post = any("post" in cs.lower() for cs in call_sites)
        
        if self.case.diag_callsite_mode == "pre_only":
            if has_pre and not has_post:
                self.metrics["callsite_mode_pass"] = True
            else:
                self.metrics["callsite_mode_pass"] = False
                self.errors.append(
                    f"Expected pre_only but found: pre={has_pre}, post={has_post}"
                )
        elif self.case.diag_callsite_mode == "post_only":
            if has_post and not has_pre:
                self.metrics["callsite_mode_pass"] = True
            else:
                self.metrics["callsite_mode_pass"] = False
                self.errors.append(
                    f"Expected post_only but found: pre={has_pre}, post={has_post}"
                )
        else:  # "both"
            if has_pre and has_post:
                self.metrics["callsite_mode_pass"] = True
            else:
                self.metrics["callsite_mode_pass"] = False
                self.errors.append(
                    f"Expected both pre and post but found: pre={has_pre}, post={has_post}"
                )
        
        # Extract flux metrics
        for metric_name in self.case.flux_metrics:
            values = []
            for row in self.rows:
                val = self.to_float(row, metric_name)
                if val is None:
                    return False
                values.append(val)
            
            if values:
                self.metrics[f"{metric_name}_mean"] = statistics.mean(values)
                self.metrics[f"{metric_name}_final"] = values[-1]
                self.metrics[f"{metric_name}_max"] = max(values)
                self.metrics[f"{metric_name}_min"] = min(values)
                
                # Check for NaN/Inf
                if any(not math.isfinite(v) for v in values):
                    self.errors.append(f"{metric_name} contains NaN or Inf")
                    return False
        
        # Extract heating metrics
        for metric_name in self.case.heating_metrics:
            values = []
            for row in self.rows:
                val = self.to_float(row, metric_name)
                if val is None:
                    return False
                values.append(val)
            
            if values:
                self.metrics[f"{metric_name}_mean"] = statistics.mean(values)
                self.metrics[f"{metric_name}_final"] = values[-1]
                self.metrics[f"{metric_name}_max"] = max(values)
                self.metrics[f"{metric_name}_min"] = min(values)
                
                # Check for NaN/Inf
                if any(not math.isfinite(v) for v in values):
                    self.errors.append(f"{metric_name} contains NaN or Inf")
                    return False
                
                # Check nonzero (where applicable)
                if not any(abs(v) > tol.HEATING_NONZERO_TOL for v in values):
                    self.warnings.append(
                        f"{metric_name} is effectively zero (all < {tol.HEATING_NONZERO_TOL})"
                    )
                
                # Coefficient of variation check
                if len(values) > 1:
                    mean = statistics.mean(values)
                    stdev = statistics.pstdev(values)
                    cv = abs(stdev / mean) if abs(mean) > 0 else 0.0
                    self.metrics[f"{metric_name}_cv"] = cv
                    
                    if cv > tol.HEATING_CV_UPPER_BOUND:
                        self.errors.append(
                            f"{metric_name} CV too high: {cv:.6f} "
                            f"(threshold: {tol.HEATING_CV_UPPER_BOUND:.6f})"
                        )
        
        return True
    
    def validate_all(self) -> Tuple[bool, Dict[str, Any]]:
        """
        Run all validation checks.
        
        Returns:
            Tuple (is_pass, metrics_dict)
        """
        if not self.read_csv():
            return False, self.metrics
        
        if not self.validate_schema():
            return False, self.metrics
        
        if not self.extract_metrics():
            return False, self.metrics
        
        # Overall pass/fail: no errors means pass
        is_pass = len(self.errors) == 0
        return is_pass, self.metrics
    
    def get_summary(self) -> Dict[str, Any]:
        """
        Get a summary of validation results.
        
        Returns:
            Dictionary with pass/fail status, errors, and warnings
        """
        return {
            "case": self.case.short_name,
            "case_name": self.case.name,
            "is_pass": len(self.errors) == 0,
            "errors": self.errors,
            "warnings": self.warnings,
            "metrics": self.metrics,
        }


def validate_case(case: BenchmarkCase, case_dir: Path) -> Tuple[bool, Dict[str, Any]]:
    """
    Validate a single benchmark case.
    
    Args:
        case: BenchmarkCase definition
        case_dir: Path to the case directory
    
    Returns:
        Tuple (is_pass, summary_dict)
    """
    csv_path = case_dir / case.diag_file
    extractor = MetricsExtractor(case, csv_path)
    is_pass, metrics = extractor.validate_all()
    summary = extractor.get_summary()
    return is_pass, summary


if __name__ == "__main__":
    print("Benchmark Metrics Extractor Module")
    print("Use validate_case() to check a benchmark case")
