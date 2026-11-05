"""Assemble the inputs schema from an ERF tree and its AMReX tree."""
from pathlib import Path

from . import canonical, cmake, scan

SCHEMA_VERSION = 1


def _files(root, rel_base):
    return [(f"{rel_base}/{scan.posix_rel(p, root)}", scan.strip_comments(p.read_text(errors="replace")))
            for p in sorted(root.rglob("*")) if p.suffix in (".H", ".cpp")]


def _decls(files, enum_names):
    out = {}
    for _, text in files:
        for name, kinds in scan.declarations(text, enum_names).items():
            out.setdefault(name, set()).update(kinds)
    return out


def layer_of(source):
    return "AMReX" if source.startswith("AMReX/") else "ERF"


def _schema_type(kind):
    """A scan type as the schema writes it: enum:<Name> -> enum, enum:<Name>_list -> enum_list."""
    if kind.startswith("enum:"):
        return "enum_list" if kind.endswith("_list") else "enum"
    return kind


PROBLEM_CODE = "Source/Prob/"      # compiled into one problem's executable only (ERF_CloudChamber.H)


def _merge_removed(old, new, lay):
    """removed wins: ERF stops the run if the key is present. A typed read elsewhere is kept in shadowed_reads, so
    the conflict stays visible. A refusal in problem-specific code (Source/Prob/) over a typed read elsewhere is the
    other way round: the read's fact is kept and the refusal is listed in refused_by, since only that problem stops
    (ERF_CloudChamber.H refuses <face>.velocity, which ERF_InitBCs.cpp reads for every other problem)."""
    keep, other = (old, new) if old["type"] == "removed" else (new, old)
    if keep["source"].startswith(PROBLEM_CODE) and other["type"] not in ("removed", "unknown"):
        out = dict(other, layer=layer_of(other["source"]), layers=dict(old["layers"]))
        out["layers"].setdefault(lay, new["source"])
        out["refused_by"] = sorted(set(other.get("refused_by", [])) | {keep["removed_evidence"]})
        return out
    keep_lay = layer_of(keep["source"])
    out = dict(keep, layer=keep_lay, layers=dict(old["layers"]))
    out["layers"].setdefault(lay, new["source"])
    out["layers"][keep_lay] = keep["source"]         # the record names the removal, whichever file sorted first
    shadowed = list(keep.get("shadowed_reads", [])) + list(other.get("shadowed_reads", []))
    if other["type"] not in ("removed", "unknown"):
        shadowed.append(f'{other["source"]} {_schema_type(other["type"])}')
    if shadowed:
        out["shadowed_reads"] = sorted(set(shadowed))
    return out


def _scalar(kind):
    return kind[:-5] if kind.endswith("_list") else kind


def _resolve_types(a, b):
    """(type, conflict) of two reads of one key: unknown yields to anything; a list read of the same scalar makes the
    key a list (ERF reads it both ways, as erf.vert_implicit_fac); an untyped list read (unknown_list) makes the other
    read's scalar a list; an enum read wins over a std::string probe of the same key (erf.init_type, erf.moisture_model,
    erf.terrain_type read the string first), with the conflict recorded; anything else keeps a and records the pair."""
    if a == b or b == "unknown":
        return a, None
    if a == "unknown":
        return b, None
    sa, sb = _scalar(a), _scalar(b)
    if sa == sb:
        return sa + "_list", None
    if sb == "unknown":
        return sa + "_list", None
    if sa == "unknown":
        return sb + "_list", None
    if sa.startswith("enum:") and sb == "string":
        return a, sorted({a, b})
    if sb.startswith("enum:") and sa == "string":
        return b, sorted({a, b})
    return a, sorted({a, b})


DEFAULT_FIELDS = ("default", "default_origin", "default_source")


def merge(old, new):
    """The first read's fact, its type refined by later reads (_resolve_types); a real disagreement is recorded,
    never hidden. AMReX's reads come first, so a shared key's fact (and its source) is AMReX's: layer names it, and
    layers keeps each layer's first read. A removed key stays removed (_merge_removed). A default is taken only from
    a read of the type the key ends with."""
    lay = layer_of(new["source"])
    if old is None:
        return dict(new, layer=lay, layers={lay: new["source"]})
    if "removed" in (old["type"], new["type"]):
        return _merge_removed(old, new, lay)
    out = dict(old, layers=dict(old["layers"]))
    out["layers"].setdefault(lay, new["source"])
    a, b = old["type"], new["type"]
    kind, conflict = _resolve_types(a, b)
    out["type"] = kind
    if conflict:
        out.setdefault("type_conflict", conflict)
    if kind != a and "case_insensitive" in new:      # a presence read carries none; the typed read's counts
        out["case_insensitive"] = new["case_insensitive"]
    if old.get("per_level") or new.get("per_level"):
        out["per_level"] = True
    out["required"] = bool(old.get("required") or new.get("required"))
    for flag in ("line", "table", "parser", "as_double"):
        if new.get(flag):
            out[flag] = True
    if "ignores" in new:
        out.setdefault("ignores", new["ignores"])
    if _scalar(kind) != _scalar(a) and _scalar(kind) == _scalar(b):      # the type changed to the new read's
        for k in DEFAULT_FIELDS:
            out.pop(k, None)
            if k in new:
                out[k] = new[k]
    elif (old.get("default_origin") == "none" and new.get("default_origin") != "none"
          and _scalar(b) == _scalar(kind)):
        for k in DEFAULT_FIELDS:
            if k in new:
                out[k] = new[k]
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
        fact["type"] = "enum_list" if fact["type"].endswith("_list") else "enum"   # scan writes enum:<Name>_list
        fact["values"] = enum_names[name]
        fact.setdefault("case_insensitive", False)
        if fact.get("default") is not None and fact["type"] == "enum":
            try:
                fact["default"] = canonical.canonical_value("enum", fact["default"], fact)
            except ValueError:
                fact.update(default=None, default_origin="expression")
        vr = {v: dict(values[(name, v)]) for v in fact["values"] if (name, v) in values}
        if vr:
            fact["value_requires"] = vr


def _apply_adds(adds, keys):
    """add_par's pp.add("k", literal) is k's default: it runs before the inputs are read, which override it."""
    for key, raw, source in adds:
        fact = keys.get(key)
        if fact is None or fact["type"] in ("removed", "unknown") or fact["type"].endswith("_list"):
            continue
        try:
            value = canonical.canonical_value(fact["type"], scan._literal(raw), fact)
        except ValueError:
            continue
        fact.update(default=value, default_origin="add", default_source=source)


def _apply_site_facts(site_facts, keys, features, feature_sources):
    """A site fact is a feature the build options cannot show (e.g. whether the MPI is GPU-aware): it sets the
    feature and, when it names a key and value, that value's requirement, with the key's read line as evidence."""
    for f in site_facts:
        features[f["feature"]] = bool(f["on"])
        feature_sources[f["feature"]] = "site"
        if f.get("key") in keys and f.get("value") is not None:
            fact = keys[f["key"]]
            where = fact["source"]
            value = str(f["value"])
            try:                                     # checks/schema.py looks the value up in canonical form
                value = canonical.canonical_value(_schema_type(fact["type"]).replace("_list", ""), value, fact)
            except ValueError:
                pass
            keys[f["key"]].setdefault("value_requires", {})[value] = {
                "feature": f["feature"], "evidence": f'{f.get("evidence", "site")}; {where}'}


def generate(erf_root, amrex_root, cmake_options=None, site_facts=None, meta=None, facts_from=None):
    """The schema for the ERF tree at erf_root (Source/, CMakeLists.txt, CMake/) and AMReX at amrex_root (Src/)."""
    erf_root, amrex_root = Path(erf_root), Path(amrex_root)
    cmake_options = dict(cmake_options or {})
    erf_files = _files(erf_root / "Source", "Source")
    enum_names = scan.enums(erf_files)
    hits, unresolved, families_hit, adds = scan.scan(erf_root / "Source", "Source", enum_names,
                                                     _decls(erf_files, enum_names), True, erf_files)
    present = {rel for rel, _ in erf_files}
    for (rel, expr) in sorted(scan.FAMILIES):
        if rel in present and (rel, expr) not in families_hit:
            raise ValueError(f"scan.FAMILIES entry {rel} {expr!r} matched no ParmParse declaration: the source "
                             "changed, so update the entry (its keys would otherwise drop out of the schema)")
    amrex_files = _files(amrex_root / "Src", "AMReX/Src")
    amrex_hits, _, _, _ = scan.scan(amrex_root / "Src", "AMReX/Src", {}, _decls(amrex_files, {}), False, amrex_files)
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
    _apply_adds(adds, keys)
    _apply_site_facts(site_facts or [], keys, features, feature_sources)
    build = {"cmake": cmake_options,
             "facts_from": facts_from or {"cmake": "given" if cmake_options else "none", "site_facts": site_facts}}
    return {"schema_version": SCHEMA_VERSION, **(meta or {}), "build": build,
            "canonical_rules": canonical.RULES, "features": dict(sorted(features.items())),
            "feature_sources": dict(sorted(feature_sources.items())), "enums": dict(sorted(enum_names.items())),
            "keys": dict(sorted(keys.items())), "patterns": dict(sorted(patterns.items())),
            "unresolved": sorted(unresolved),
            "value_requires_candidates": [dict(enum=e, value=v, **f) for (e, v), f in sorted(candidates.items())]}
