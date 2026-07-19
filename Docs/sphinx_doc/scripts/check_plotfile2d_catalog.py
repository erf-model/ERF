#!/usr/bin/env python3
"""Check the documented fixed 2-D diagnostic catalog against the C++ catalog."""

from pathlib import Path
import re
import sys


BEGIN = ".. BEGIN ERF BUILT-IN 2D DIAGNOSTIC CATALOG"
END = ".. END ERF BUILT-IN 2D DIAGNOSTIC CATALOG"

DESCRIPTOR_RE = re.compile(
    r'\{(?:DiagnosticID::\w+),\s*"([^"]+)",\s*"([^"]+)",\s*'
    r'"([^"]+)",\s*DiagnosticCategory::(\w+),\s*'
    r'MissingPolicy::(\w+)\}',
    re.DOTALL,
)
ROW_RE = re.compile(r'^\s*\* - ``([^`]+)``\s*$')


def fail(message: str) -> None:
    print(f"plotfile2d catalog check: {message}", file=sys.stderr)
    raise SystemExit(1)


def parse_cpp_catalog(path: Path):
    source = path.read_text(encoding="utf-8")
    start = source.find("catalog_storage ()")
    end = source.find("return catalog;", start)
    if start < 0 or end < 0:
        fail(f"could not locate catalog initializer in {path}")
    rows = DESCRIPTOR_RE.findall(source[start:end])
    if not rows:
        fail(f"could not parse descriptors in {path}")
    return [
        {
            "name": name,
            "category": category,
            "units": units,
            "policy": policy,
            "long_name": long_name,
        }
        for name, long_name, units, category, policy in rows
    ]


def strip_markup(value: str) -> str:
    value = value.strip()
    if value.startswith("``") and value.endswith("``"):
        return value[2:-2]
    return value


def parse_rst_catalog(path: Path):
    lines = path.read_text(encoding="utf-8").splitlines()
    try:
        begin = lines.index(BEGIN)
        end = lines.index(END, begin + 1)
    except ValueError:
        fail(f"missing or duplicated catalog markers in {path}")

    rows = []
    i = begin + 1
    while i < end:
        match = ROW_RE.match(lines[i])
        if not match:
            i += 1
            continue
        values = [match.group(1)]
        for offset in range(1, 5):
            if i + offset >= end:
                fail(f"truncated RST row for {values[0]!r}")
            line = lines[i + offset]
            if not re.match(r'^\s*-\s+', line):
                fail(f"unparseable RST row for {values[0]!r}: {line.strip()}")
            values.append(strip_markup(re.sub(r'^\s*-\s+', '', line)))
        rows.append(dict(zip(("name", "category", "units", "policy", "long_name"), values)))
        i += 5
    return rows, "\n".join(lines[begin : end + 1])


def main() -> int:
    repo = Path(__file__).resolve().parents[3]
    cpp = repo / "Source/IO/ERF_Plotfile2DCatalog.cpp"
    rst = repo / "Docs/sphinx_doc/Plotfiles.rst"
    expected = parse_cpp_catalog(cpp)
    documented, marked_region = parse_rst_catalog(rst)

    if len(expected) != 64:
        fail(f"expected 64 C++ descriptors, found {len(expected)}")
    if len(documented) != len(expected):
        fail(f"expected {len(expected)} documented rows, found {len(documented)}")

    seen = set()
    for index, (source, docs) in enumerate(zip(expected, documented), start=1):
        name = docs["name"]
        if name in seen:
            fail(f"duplicate documented row {name!r} at position {index}")
        seen.add(name)
        if source != docs:
            for field in source:
                if source[field] != docs[field]:
                    fail(
                        f"row {index} {name!r} field {field}: "
                        f"expected {source[field]!r}, documented {docs[field]!r}"
                    )
            fail(f"row {index} {name!r} does not match the C++ descriptor")

    for literal in ("smois_<layer>", "sh2o_<layer>", "tslb_<layer>"):
        if literal not in rst.read_text(encoding="utf-8"):
            fail(f"missing dynamic soil pattern {literal!r}")

    print(f"plotfile2d catalog check: {len(documented)} fixed descriptors match")
    print("plotfile2d catalog check: dynamic soil families match required patterns")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
