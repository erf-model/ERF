#!/usr/bin/env python3
"""
Phase 8 Benchmark Suite: Case Definitions

Defines the 5 benchmark test cases:
1. Clear-sky SW baseline
2. LW isothermal baseline
3. Cloud-layer absorption
4. Cloud scattering
5. Coupled SW+LW non-isothermal time-integration (Phase 6/7 style)

Each case references its input file and specifies expected diagnostics behavior,
metrics to validate, and pass/fail thresholds.
"""

from dataclasses import dataclass
from typing import Dict, Any, Optional, List
from pathlib import Path

# =========================================================================
# CASE DEFINITIONS
# =========================================================================

@dataclass
class BenchmarkCase:
    """Definition of a single benchmark test case."""
    
    # Case identification
    name: str                           # "Clear-sky SW baseline", etc.
    short_name: str                     # "sw_clearsky", "lw_isothermal", etc.
    description: str                    # Physics description
    
    # Input configuration
    base_case_dir: str                  # Relative path to existing test (or symlink target)
    inputs_file: str                    # "inputs" or specific config file name
    diag_file: str                      # Expected diagnostic CSV output filename
    
    # Simulation parameters
    dt: float                           # Fixed timestep [s]
    stop_time: float                    # Simulation stop time [s]
    expected_steps: int                 # Expected number of coarse steps
    
    # Diagnostics configuration
    diag_callsite_mode: str            # "both", "pre_only", or "post_only"
    expected_diag_rows: int            # Expected total rows in diagnostic CSV
    rows_per_step: int                 # Expected rows per timestep (depends on mode)
    
    # Metrics to validate
    required_metrics: List[str]        # Columns to check exist in CSV
    flux_metrics: Optional[Dict[str, str]] = None  # {"SW_TOA": "W/m^2", ...}
    heating_metrics: Optional[Dict[str, str]] = None  # {"heating_rate_max": "K/s", ...}
    
    # Physics-specific checks
    physics_description: str = ""      # Brief physics note
    
    def __post_init__(self):
        if self.flux_metrics is None:
            self.flux_metrics = {}
        if self.heating_metrics is None:
            self.heating_metrics = {}


# =========================================================================
# BENCHMARK CASE SUITE
# =========================================================================

CASES: Dict[str, BenchmarkCase] = {
    
    # =====================================================================
    # CASE 1: Clear-sky SW baseline
    # =====================================================================
    "sw_clearsky": BenchmarkCase(
        name="Clear-sky SW baseline",
        short_name="sw_clearsky",
        description="Beer-Lambert shortwave direct-beam solar radiation in clear sky",
        base_case_dir="SW_ClearSky_Analytical",
        inputs_file="inputs",
        diag_file="radiation_sw_diag.dat",
        dt=0.5,
        stop_time=0.25,
        expected_steps=1,  # stop_time / dt = 0.25 / 0.5 = 0.5, rounds to 1 step
        diag_callsite_mode="both",
        expected_diag_rows=2,  # 2 rows/step × 1 step
        rows_per_step=2,
        required_metrics=["step", "time", "call_site", "SW_TOA", "SW_surface", "heating_rate_max"],
        flux_metrics={"SW_TOA": "W/m^2", "SW_surface": "W/m^2"},
        heating_metrics={"heating_rate_max": "K/s"},
        physics_description="S0=1361 W/m^2, zenith=60°, tau=0.05/layer, Beer-Lambert",
    ),
    
    # =====================================================================
    # CASE 2: LW isothermal baseline
    # =====================================================================
    "lw_isothermal": BenchmarkCase(
        name="LW isothermal baseline",
        short_name="lw_isothermal",
        description="Gray-gas longwave radiation in isothermal mode (energy balance check)",
        base_case_dir="LW_Isothermal",
        inputs_file="inputs",
        diag_file="radiation_lw_diag.dat",
        dt=0.5,
        stop_time=0.25,
        expected_steps=1,
        diag_callsite_mode="both",
        expected_diag_rows=2,
        rows_per_step=2,
        required_metrics=["step", "time", "call_site", "LW_surface", "LW_TOA", "heating_rate_max"],
        flux_metrics={"LW_surface": "W/m^2", "LW_TOA": "W/m^2"},
        heating_metrics={"heating_rate_max": "K/s"},
        physics_description="T_iso=288.15K, tau_lw=1.0/layer, expect heating=0 (energy balance)",
    ),
    
    # =====================================================================
    # CASE 3: Cloud-layer absorption
    # =====================================================================
    "sw_cloud_layer": BenchmarkCase(
        name="Cloud-layer absorption",
        short_name="sw_cloud_layer",
        description="Shortwave radiation with cloud-layer absorption",
        base_case_dir="SW_Cloud_Layer",
        inputs_file="inputs",
        diag_file="radiation_sw_diag.dat",
        dt=0.5,
        stop_time=0.25,
        expected_steps=1,
        diag_callsite_mode="both",
        expected_diag_rows=2,
        rows_per_step=2,
        required_metrics=["step", "time", "call_site", "SW_TOA", "SW_surface", "heating_rate_max"],
        flux_metrics={"SW_TOA": "W/m^2", "SW_surface": "W/m^2"},
        heating_metrics={"heating_rate_max": "K/s"},
        physics_description="Cloud layer with height-varying optical depth, S0=1361",
    ),
    
    # =====================================================================
    # CASE 4: Cloud scattering
    # =====================================================================
    "sw_scattering": BenchmarkCase(
        name="Cloud scattering",
        short_name="sw_scattering",
        description="Shortwave radiation with cloud scattering via two-stream approximation",
        base_case_dir="SW_Scattering_Cloud",
        inputs_file="inputs",
        diag_file="radiation_sw_diag.dat",
        dt=0.5,
        stop_time=0.25,
        expected_steps=1,
        diag_callsite_mode="both",
        expected_diag_rows=2,
        rows_per_step=2,
        required_metrics=["step", "time", "call_site", "SW_TOA", "SW_surface", "heating_rate_max"],
        flux_metrics={"SW_TOA": "W/m^2", "SW_surface": "W/m^2"},
        heating_metrics={"heating_rate_max": "K/s"},
        physics_description="Cloud layer with Meador-Weaver two-stream scattering",
    ),
    
    # =====================================================================
    # CASE 5: Coupled SW+LW non-isothermal time-integration
    # =====================================================================
    "phase6_timing": BenchmarkCase(
        name="Coupled SW+LW non-isothermal time-integration",
        short_name="phase6_timing",
        description="Phase 6/7 style: coupled SW+LW with time-stepping and call-site diagnostics",
        base_case_dir="./cases/phase6_timing",  # Local case directory with pre_only mode inputs
        inputs_file="inputs",
        diag_file="radiation_phase6_timing_diag.dat",
        dt=0.5,
        stop_time=5.0,
        expected_steps=10,  # stop_time / dt = 5.0 / 0.5 = 10 steps
        diag_callsite_mode="pre_only",  # Test single call-site mode (pre_only)
        expected_diag_rows=10,  # 1 row/step × 10 steps
        rows_per_step=1,  # Only pre_dycore, not post_dycore
        required_metrics=["step", "time", "call_site", "SW_TOA", "LW_surface", "heating_rate_max"],
        flux_metrics={"SW_TOA": "W/m^2", "LW_surface": "W/m^2"},
        heating_metrics={"heating_rate_max": "K/s"},
        physics_description="SW+LW coupling over 10 time steps, pre_only mode, RhoTheta coupling",
    ),
}

# =========================================================================
# HELPER FUNCTIONS
# =========================================================================

def get_case(short_name: str) -> BenchmarkCase:
    """Get a benchmark case by short name."""
    if short_name not in CASES:
        raise KeyError(f"Unknown benchmark case: {short_name}")
    return CASES[short_name]


def list_all_cases() -> List[str]:
    """Return list of all available case short names in deterministic order."""
    return sorted(CASES.keys())


def get_all_cases() -> Dict[str, BenchmarkCase]:
    """Return all cases as a dictionary."""
    return CASES


def print_case_matrix():
    """Print a summary table of all benchmark cases."""
    print("\n" + "=" * 100)
    print("PHASE 8 BENCHMARK SUITE: CASE MATRIX")
    print("=" * 100)
    print(f"{'Case':<20} {'Physics':<30} {'Steps':<8} {'Mode':<12} {'Expected Rows':<15}")
    print("-" * 100)
    for short_name in list_all_cases():
        case = CASES[short_name]
        print(
            f"{short_name:<20} {case.physics_description:<30} "
            f"{case.expected_steps:<8} {case.diag_callsite_mode:<12} {case.expected_diag_rows:<15}"
        )
    print("=" * 100 + "\n")


if __name__ == "__main__":
    print_case_matrix()
