#!/usr/bin/env python3
"""Write the inputs schema for one ERF build: every literal ParmParse query in ERF and its AMReX, typed, with the
build requirements each key or enum value has, and the canonical form of each type.

Offline and read-only: source is read by `git archive` from a local clone, so no checkout changes and no build runs.

    python3 Exec/Tools/build_erf_schema.py --repo . --erf HEAD --cmake="-DERF_ENABLE_NETCDF=ON" --out schema.json
    python3 Exec/Tools/build_erf_schema.py --erf-dir <erf tree> --amrex-dir <amrex tree> --out schema.json

AMReX is the commit ERF's Submodules/AMReX points at, unless --amrex is given. Build options come from --cmake=... (a -D
list, repeatable; the = form, since the value starts with a dash) and/or --cmake-cache (a CMakeCache.txt). Facts that are the site's, not the build's (for example
whether the MPI is GPU-aware), come from --site-facts. Check a deck against the result with
`check_erf_inputs.py --schema schema.json inputs`.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path

TOOLS = Path(__file__).resolve().parent
sys.path.insert(0, str(TOOLS))

from erf_inputs_schema import cmake, schema, source  # noqa: E402

GENERATOR = "Exec/Tools/build_erf_schema.py"


def parse_args(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--repo", type=Path, help="an ERF git clone (git mode)")
    p.add_argument("--erf", help="the ERF commit to read (git mode)")
    p.add_argument("--amrex-repo", type=Path, help="AMReX clone (default: <repo>/Submodules/AMReX)")
    p.add_argument("--amrex", help="AMReX commit (default: ERF's gitlink at --erf)")
    p.add_argument("--erf-dir", type=Path, help="an ERF tree with Source/, CMake/, CMakeLists.txt (directory mode)")
    p.add_argument("--amrex-dir", type=Path, help="an AMReX tree with Src/ (directory mode)")
    p.add_argument("--cmake", action="append", default=[], help="cmake -D options, given as --cmake='-DERF_ENABLE_NETCDF=ON ...'")
    p.add_argument("--cmake-cache", type=Path, help="a CMakeCache.txt to take ERF_* options from")
    p.add_argument("--site-facts", type=Path, help="a JSON list of site facts (see erf_inputs_schema/README)")
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    git_mode, dir_mode = (a.repo, a.erf), (a.erf_dir, a.amrex_dir)
    if any(git_mode) == any(dir_mode) or not all(git_mode if any(git_mode) else dir_mode):
        p.error("give either --repo and --erf, or --erf-dir and --amrex-dir")
    return a


def build_options(a):
    text = " ".join(a.cmake)
    opts = cmake.parse_d_options(a.cmake_cache.read_text()) if a.cmake_cache else {}
    opts.update(cmake.parse_d_options(text))
    return dict(sorted(opts.items()))


def facts_from(a, site):
    """Which kind of build facts the schema came from: a CMakeCache.txt holds every ERF option, a -D list only those
    passed (the rest take ERF's CMakeLists defaults)."""
    kind = "+".join(k for k, given in (("cmake_cache", a.cmake_cache), ("d_list", a.cmake)) if given) or "none"
    return {"cmake": kind, "site_facts": site}


def write(result, out):
    out.write_text(json.dumps(result, indent=1, sort_keys=True, ensure_ascii=False) + "\n")


def main(argv=None):
    a = parse_args(argv)
    site = json.loads(a.site_facts.read_text()) if a.site_facts else None
    meta = {"generator": GENERATOR, "generator_commit": source.generator_commit(TOOLS)}
    given = (build_options(a), site, meta, facts_from(a, site))
    if a.erf_dir:
        meta.update(erf_commit=None, amrex_commit=None)
        result = schema.generate(a.erf_dir, a.amrex_dir, *given)
    else:
        with tempfile.TemporaryDirectory() as td:
            erf_root, amrex_root, erf_c, amrex_c = source.checkout(a.repo, a.erf, td, a.amrex_repo, a.amrex)
            meta.update(erf_commit=erf_c, amrex_commit=amrex_c)
            result = schema.generate(erf_root, amrex_root, *given)
    write(result, a.out)
    print(json.dumps({"keys": len(result["keys"]), "patterns": len(result["patterns"]),
                      "unresolved": len(result["unresolved"]),
                      "features_on": [k for k, v in result["features"].items() if v]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
