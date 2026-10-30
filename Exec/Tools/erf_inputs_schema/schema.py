"""Assemble the inputs schema from an ERF tree and its AMReX tree."""
from pathlib import Path

from . import canonical, cmake, scan

SCHEMA_VERSION = 1


def _files(root, rel_base):
    return [(f"{rel_base}/{p.relative_to(root)}", scan.strip_comments(p.read_text(errors="replace")))
            for p in sorted(root.rglob("*")) if p.suffix in (".H", ".cpp")]


def _decls(files, enum_names):
    out = {}
    for _, text in files:
        for name, kinds in scan.declarations(text, enum_names).items():
            out.setdefault(name, set()).update(kinds)
    return out


def layer_of(source):
    return "AMReX" if source.startswith("AMReX/") else "ERF"


def merge(old, new):
    """The first read's fact, its type refined by later reads; a real disagreement is recorded, never hidden.
    AMReX's reads come first, so a shared key's fact (and its source) is AMReX's: layer names it, and layers keeps
    each layer's first read."""
    lay = layer_of(new["source"])
    if old is None:
        return dict(new, layer=lay, layers={lay: new["source"]})
    out = dict(old, layers=dict(old["layers"]))
    out["layers"].setdefault(lay, new["source"])
    a, b = old["type"], new["type"]
    if a == "unknown" or b == a + "_list":
        out["type"] = b
    elif b != a and b != "unknown" and a != b + "_list":
        out.setdefault("type_conflict", sorted({a, b}))
    if old.get("per_level") or new.get("per_level"):
        out["per_level"] = True
    return out


def _requirements(hits):
    """key -> {macro: evidence line} for the macros every read of the key sits under."""
    needs = {}
    for key, fact, conds, _, _ in hits:
        here = {m: f'{fact["source"].rsplit(":", 1)[0]}:{line}' for m, positive, line in conds if positive}
        if key not in needs:
            needs[key] = here
        else:
            needs[key] = {m: e for m, e in needs[key].items() if m in here}
    return needs


def _enum_facts(table, enum_names, values):
    for fact in table.values():
        if not fact["type"].startswith("enum"):
            continue
        name = fact["type"].split(":", 1)[1].replace("_list", "")
        fact["type"] = "enum_list" if fact["type"].startswith("enum_list") else "enum"
        fact["values"] = enum_names[name]
        fact.setdefault("case_insensitive", False)
        vr = {v: dict(values[(name, v)]) for v in fact["values"] if (name, v) in values}
        if vr:
            fact["value_requires"] = vr


def generate(erf_root, amrex_root, cmake_options=None, site_facts=None, meta=None):
    """The schema for the ERF tree at erf_root (Source/, CMakeLists.txt, CMake/) and AMReX at amrex_root (Src/)."""
    erf_root, amrex_root = Path(erf_root), Path(amrex_root)
    cmake_options = dict(cmake_options or {})
    erf_files = _files(erf_root / "Source", "Source")
    enum_names = scan.enums(erf_files)
    hits, unresolved = scan.scan(erf_root / "Source", "Source", enum_names, _decls(erf_files, enum_names), True,
                                 erf_files)
    amrex_files = _files(amrex_root / "Src", "AMReX/Src")
    amrex_hits, _ = scan.scan(amrex_root / "Src", "AMReX/Src", {}, _decls(amrex_files, {}), False, amrex_files)
    features, feature_sources = cmake.cmake_features(erf_root, cmake_options)
    for macro in sorted({c[0] for _, _, conds, _, _ in hits for c in conds}):
        features.setdefault(macro, False)
    values, candidates = scan.value_requirements(erf_files, enum_names)

    keys, patterns = {}, {}
    for key, fact, _, _, _ in amrex_hits + hits:
        table = patterns if "<" in key else keys
        table[key] = merge(table.get(key), fact)
    needs = _requirements(amrex_hits + hits)
    for table in (keys, patterns):
        _enum_facts(table, enum_names, values)
        for key, fact in table.items():
            if needs.get(key):
                fact["requires"] = sorted(needs[key])
                fact["requires_evidence"] = dict(sorted(needs[key].items()))
    return {"schema_version": SCHEMA_VERSION, **(meta or {}), "build": {"cmake": cmake_options},
            "canonical_rules": canonical.RULES, "features": dict(sorted(features.items())),
            "feature_sources": dict(sorted(feature_sources.items())), "enums": dict(sorted(enum_names.items())),
            "keys": dict(sorted(keys.items())), "patterns": dict(sorted(patterns.items())),
            "unresolved": sorted(unresolved),
            "value_requires_candidates": [dict(enum=e, value=v, **f) for (e, v), f in sorted(candidates.items())]}
