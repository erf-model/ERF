#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Optional


BASE_DIR = Path(__file__).resolve().parent
CASES_DIR = BASE_DIR / "cases"
RUNS_DIR = BASE_DIR / "_runs"
ERF_EXE = Path(os.environ.get("ERF_EXE", "../../../build/Exec/erf_exec")).resolve()

TOL = {
    "sw_toa_rel_pct": 0.1,
    "row_count_abs": 2,
    "heating_nonzero_min": 1e-12,
    "heating_cv_max": 0.05,
}

CASE_PATHS = {
    "LW_ISOTHERMAL": CASES_DIR / "lw_isothermal",
    "PHASE6_TIMING": CASES_DIR / "phase6_timing",
    "SW_CLEARSKY": CASES_DIR / "sw_clearsky",
    "SW_CLOUD_LAYER": CASES_DIR / "sw_cloud_layer",
    "SW_SCATTERING": CASES_DIR / "sw_scattering",
}


@dataclass
class CaseSpec:
    key: str
    name: str
    case_dir: Path
    input_file: str
    diag_file: str
    diag_mode: str
    expected_steps: int
    expected_rows: int
    expected_sw_toa: float
    check_sw_toa: bool = True
    check_heating_nonzero: bool = True


@dataclass
class CaseResult:
    key: str
    name: str
    status: str
    reason: str
    rows_found: int = 0
    expected_rows: int = 0
    sw_toa_mean: float = float("nan")
    sw_toa_final: float = float("nan")
    heating_mean: float = float("nan")
    heating_final: float = float("nan")
    heating_cv: float = float("nan")
    call_sites: Dict[str, int] = None
    work_dir: str = ""


def strip_comment(line: str) -> str:
    return line.split("#", 1)[0].strip()


def extract_value(text: str, key: str) -> Optional[str]:
    pat = re.compile(rf"^\s*{re.escape(key)}\s*=\s*(.*?)\s*$", re.MULTILINE)
    m = pat.search(text)
    if not m:
        return None
    return strip_comment(m.group(1))


def first_numeric_token(raw: Optional[str]) -> Optional[float]:
    if raw is None:
        return None
    s = raw.strip().strip('"').strip("'")
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def parse_str(text: str, keys: List[str], default: Optional[str] = None) -> Optional[str]:
    for k in keys:
        v = extract_value(text, k)
        if v is not None:
            return v.strip().strip('"').strip("'")
    return default


def parse_float(text: str, keys: List[str], default: Optional[float] = None) -> Optional[float]:
    for k in keys:
        raw = extract_value(text, k)
        val = first_numeric_token(raw)
        if val is not None:
            return val
    return default


def parse_int(text: str, keys: List[str], default: Optional[int] = None) -> Optional[int]:
    v = parse_float(text, keys, None)
    if v is None:
        return default
    return int(round(v))


def load_case_spec(key: str, case_dir: Path) -> CaseSpec:
    candidates = sorted(case_dir.glob("inputs*"))
    if not candidates:
        raise FileNotFoundError(f"No inputs* file found in {case_dir}")
    input_path = candidates[0]
    txt = input_path.read_text()

    stop_time = parse_float(txt, ["stop_time", "amr.stop_time"], None)
    fixed_dt = parse_float(txt, ["fixed_dt", "erf.fixed_dt", "erf.fixed_dt[0]"], None)
    max_step = parse_int(txt, ["max_step", "amr.max_step"], None)

    # Fallback derivations
    if fixed_dt is None and stop_time is not None and max_step is not None and max_step > 0:
        fixed_dt = stop_time / float(max_step)

    if stop_time is None and fixed_dt is not None and max_step is not None and max_step >= 0:
        stop_time = fixed_dt * float(max_step)

    if fixed_dt is None:
        raise ValueError(
            f"{key}: missing fixed_dt (tried fixed_dt / erf.fixed_dt / erf.fixed_dt[0]) in {input_path}"
        )

    # Determine steps
    if stop_time is not None:
        steps = int(round(stop_time / fixed_dt))
    elif max_step is not None:
        steps = max_step
    else:
        steps = 0  # snapshot-style case

    diag_file = parse_str(txt, ["erf.radiation.diag_file"], "radiation_diagnostics.dat")
    mode = parse_str(txt, ["erf.radiation.diag_callsite_mode"], "both")
    diag_enable_s = parse_str(txt, ["erf.radiation.diag_enable"], None)
    if diag_enable_s is not None and diag_enable_s.lower() == "false":
        mode = "off"

    s0 = parse_float(txt, ["erf.radiation.S0"], 1361.0)
    zen = parse_float(txt, ["erf.radiation.solar_zenith"], 60.0)
    expected_sw_toa = s0 * math.cos(math.radians(zen))

    # Expected rows:
    # If steps>0: cadence scales with steps
    # If steps==0: treat as one-shot diagnostics event
    if mode == "both":
        expected_rows = 2 * steps if steps > 0 else 2
    elif mode in ("pre_only", "post_only"):
        expected_rows = steps if steps > 0 else 1
    elif mode == "off":
        expected_rows = 0
    else:
        raise ValueError(f"{key}: invalid diag_callsite_mode='{mode}'")

    check_sw_toa = True
    check_heating_nonzero = True
    if key == "LW_ISOTHERMAL":
        check_sw_toa = False
        check_heating_nonzero = False

    return CaseSpec(
        key=key,
        name=key.replace("_", " ").title(),
        case_dir=case_dir,
        input_file=input_path.name,
        diag_file=diag_file,
        diag_mode=mode,
        expected_steps=steps,
        expected_rows=expected_rows,
        expected_sw_toa=expected_sw_toa,
        check_sw_toa=check_sw_toa,
        check_heating_nonzero=check_heating_nonzero,
    )


def run_case(spec: CaseSpec) -> CaseResult:
    work_dir = RUNS_DIR / spec.key.lower()
    if work_dir.exists():
        shutil.rmtree(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    for item in spec.case_dir.iterdir():
        dest = work_dir / item.name
        if item.is_dir():
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)

    diag_path = work_dir / spec.diag_file
    if diag_path.exists():
        diag_path.unlink()

    if not ERF_EXE.exists():
        return CaseResult(spec.key, spec.name, "FAIL", f"ERF executable not found: {ERF_EXE}",
                          expected_rows=spec.expected_rows, call_sites={}, work_dir=str(work_dir))

    log_file = work_dir / "run.log"
    cmd = [str(ERF_EXE), spec.input_file]
    proc = subprocess.run(cmd, cwd=work_dir, stdout=log_file.open("w"), stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        return CaseResult(spec.key, spec.name, "FAIL",
                          f"Simulation failed (exit {proc.returncode}); see {log_file}",
                          expected_rows=spec.expected_rows, call_sites={}, work_dir=str(work_dir))

    if spec.diag_mode == "off" and not diag_path.exists():
        return CaseResult(spec.key, spec.name, "PASS",
                          "Diagnostics disabled and no file produced (expected)",
                          0, 0, call_sites={}, work_dir=str(work_dir))

    if not diag_path.exists():
        return CaseResult(spec.key, spec.name, "FAIL", f"Diagnostics file missing: {diag_path.name}",
                          expected_rows=spec.expected_rows, call_sites={}, work_dir=str(work_dir))

    rows = list(csv.DictReader(diag_path.open()))
    rows_found = len(rows)

    call_sites: Dict[str, int] = {}
    sw_toa_vals: List[float] = []
    heat_vals: List[float] = []

    for r in rows:
        cs = (r.get("call_site") or "").strip()
        call_sites[cs] = call_sites.get(cs, 0) + 1
        try:
            sw_toa_vals.append(float(r["SW_TOA"]))
            heat_vals.append(float(r["heating_rate_max"]))
        except Exception:
            return CaseResult(spec.key, spec.name, "FAIL",
                              "Could not parse numeric diagnostics fields",
                              rows_found, spec.expected_rows, call_sites=call_sites, work_dir=str(work_dir))

    def mean(xs: List[float]) -> float:
        return sum(xs) / len(xs) if xs else float("nan")

    def cv(xs: List[float]) -> float:
        if not xs:
            return float("nan")
        m = mean(xs)
        if abs(m) < 1e-30:
            return 0.0
        var = sum((x - m) ** 2 for x in xs) / len(xs)
        return math.sqrt(var) / abs(m)

    sw_toa_mean = mean(sw_toa_vals)
    sw_toa_final = sw_toa_vals[-1] if sw_toa_vals else float("nan")
    heating_mean = mean(heat_vals)
    heating_final = heat_vals[-1] if heat_vals else float("nan")
    heating_cv = cv(heat_vals)

    reasons: List[str] = []

    if abs(rows_found - spec.expected_rows) > TOL["row_count_abs"]:
        reasons.append(f"row_count {rows_found} != expected {spec.expected_rows}±{TOL['row_count_abs']}")

    if spec.diag_mode == "both":
        if "pre_dycore" not in call_sites or "post_dycore" not in call_sites:
            reasons.append("both mode requires pre_dycore and post_dycore")
    elif spec.diag_mode == "pre_only":
        if any(k and k != "pre_dycore" for k in call_sites.keys()):
            reasons.append(f"pre_only found unexpected call_site(s): {list(call_sites.keys())}")
    elif spec.diag_mode == "post_only":
        if any(k and k != "post_dycore" for k in call_sites.keys()):
            reasons.append(f"post_only found unexpected call_site(s): {list(call_sites.keys())}")

    if spec.check_sw_toa:
        rel_sw = abs(sw_toa_mean - spec.expected_sw_toa) / max(abs(spec.expected_sw_toa), 1e-12) * 100.0
        if rel_sw > TOL["sw_toa_rel_pct"]:
            reasons.append(f"SW_TOA rel error {rel_sw:.4f}% > {TOL['sw_toa_rel_pct']}%")

    if any(not math.isfinite(x) for x in heat_vals):
        reasons.append("heating_rate_max contains NaN/Inf")

    if spec.diag_mode != "off" and spec.check_heating_nonzero:
        if abs(heating_mean) < TOL["heating_nonzero_min"]:
            reasons.append(f"heating mean too small: {heating_mean:.3e}")
        if heating_cv > TOL["heating_cv_max"]:
            reasons.append(f"heating CV too large: {heating_cv:.4f} > {TOL['heating_cv_max']:.4f}")

    status = "PASS" if not reasons else "FAIL"
    return CaseResult(
        key=spec.key,
        name=spec.name,
        status=status,
        reason="; ".join(reasons) if reasons else "Case validation succeeded",
        rows_found=rows_found,
        expected_rows=spec.expected_rows,
        sw_toa_mean=sw_toa_mean,
        sw_toa_final=sw_toa_final,
        heating_mean=heating_mean,
        heating_final=heating_final,
        heating_cv=heating_cv,
        call_sites=call_sites,
        work_dir=str(work_dir),
    )


def write_reports(results: List[CaseResult]) -> Tuple[Path, Path]:
    js_path = BASE_DIR / "benchmark_summary.json"
    md_path = BASE_DIR / "benchmark_summary.md"

    payload = {
        "timestamp": datetime.now().isoformat(),
        "base_dir": str(BASE_DIR),
        "erf_exe": str(ERF_EXE),
        "tolerances": TOL,
        "results": [asdict(r) for r in results],
        "totals": {
            "cases": len(results),
            "passed": sum(1 for r in results if r.status == "PASS"),
            "failed": sum(1 for r in results if r.status != "PASS"),
        },
    }
    js_path.write_text(json.dumps(payload, indent=2))

    lines = [
        "# Phase 8 Benchmark Summary",
        "",
        f"- Timestamp: `{payload['timestamp']}`",
        f"- ERF_EXE: `{payload['erf_exe']}`",
        "",
        "| Case | Status | Rows | Expected | SW_TOA(mean) | Heating CV | Reason |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for r in results:
        lines.append(
            f"| {r.key} | {r.status} | {r.rows_found} | {r.expected_rows} | "
            f"{r.sw_toa_mean:.6f} | {r.heating_cv:.6f} | {r.reason} |"
        )
    md_path.write_text("\n".join(lines) + "\n")
    return js_path, md_path


def main() -> int:
    print("=" * 80)
    print("PHASE 8: VALIDATION & BENCHMARKING SUITE FOR TWOSTREAM RADIATION")
    print("=" * 80)
    print(f"Timestamp: {datetime.now().isoformat()}")
    print(f"Base directory: {BASE_DIR}")
    print(f"Cases directory: {CASES_DIR}")
    print(f"ERF executable: {ERF_EXE}")
    print()

    if not CASES_DIR.exists():
        print(f"[FATAL] Cases directory not found: {CASES_DIR}")
        return 2

    RUNS_DIR.mkdir(parents=True, exist_ok=True)

    specs: List[CaseSpec] = []
    for key, path in CASE_PATHS.items():
        if not path.exists():
            print(f"[ERROR] Case directory not found: {path}")
            continue
        try:
            spec = load_case_spec(key, path)
            specs.append(spec)
        except Exception as e:
            print(f"[ERROR] Failed to load case spec for {key}: {e}")

    results: List[CaseResult] = []
    for spec in specs:
        print("-" * 80)
        print(f"CASE: {spec.key}")
        print(f"  Dir: {spec.case_dir}")
        print(f"  Input: {spec.input_file}")
        print(f"  Diagnostics mode: {spec.diag_mode}")
        print(f"  Expected steps: {spec.expected_steps}")
        print(f"  Expected rows: {spec.expected_rows}")
        print(f"  Expected SW_TOA: {spec.expected_sw_toa:.6f}")
        print(f"  SW_TOA check enabled: {spec.check_sw_toa}")
        print(f"  Heating nonzero check enabled: {spec.check_heating_nonzero}")
        res = run_case(spec)
        results.append(res)
        print(f"  [{res.status}] {res.reason}")

    js, md = write_reports(results)

    passed = sum(1 for r in results if r.status == "PASS")
    failed = sum(1 for r in results if r.status != "PASS")

    print("=" * 80)
    print("BENCHMARK SUITE SUMMARY")
    print("=" * 80)
    print(f"Total cases run: {len(results)}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    print(f"Reports:\n  - {js}\n  - {md}")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
