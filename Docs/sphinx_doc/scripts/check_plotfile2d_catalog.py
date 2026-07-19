#!/usr/bin/env python3
"""Check the public fixed 2-D diagnostic catalog against its C++ source.

ERF presents fixed 2-D diagnostic metadata in two places: the authoritative
``Source/IO/ERF_Plotfile2DCatalog.cpp`` file and the marked table in
``Docs/sphinx_doc/Plotfiles.rst``. The repetition helps users, but it can
drift. This checker keeps the public documentation and C++ catalog aligned.

``Docs/BuildDocs.sh`` runs this script before Doxygen and Sphinx. The GitHub
documentation workflow runs ``Docs/BuildDocs.sh``, so a checker failure stops
documentation CI before generation. Developers can run it from the repository
root with::

    python3 Docs/sphinx_doc/scripts/check_plotfile2d_catalog.py

The checker compares public name, canonical order, category, exact unit
string, missing-value policy, and long name. It also checks the three dynamic
soil family patterns in their separate marked RST region. The standard-library
tests run with::

    python3 -m unittest discover -s Docs/sphinx_doc/scripts \
        -p 'test_check_plotfile2d_catalog.py'

This checker does not prove configuration-level selectability, runtime field
values, numerical formulas, source-selection behavior, correctness of
``2DMetadata.json``, or full dynamic-soil metadata parity. Source review,
C++ tests, and Sphinx review cover those contracts.

The parser is deliberately narrow and fail-closed. A changed C++ initializer
or RST table format requires an intentional parser and test update; unfamiliar
entries are never silently omitted. Expected validation failures print one
prefixed error on stderr and exit with status 1. ``Docs/BuildDocs.sh`` uses
``set -e``, so that status stops the documentation build.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys
from typing import List, Optional, Sequence, Tuple


FIXED_CATALOG_BEGIN = ".. BEGIN ERF BUILT-IN 2D DIAGNOSTIC CATALOG"
FIXED_CATALOG_END = ".. END ERF BUILT-IN 2D DIAGNOSTIC CATALOG"
DYNAMIC_SOIL_BEGIN = ".. BEGIN ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES"
DYNAMIC_SOIL_END = ".. END ERF DYNAMIC SOIL DIAGNOSTIC FAMILIES"
DYNAMIC_SOIL_PATTERNS = ("smois_<layer>", "sh2o_<layer>", "tslb_<layer>")

CATALOG_FUNCTION_RE = re.compile(
    r"(?ms)^[ \t]*const[ \t]+amrex::Vector<DiagnosticDescriptor>&"
    r"\s*catalog_storage\s*\(\s*\)\s*\{"
)
CATALOG_INITIALIZER_RE = re.compile(
    r"(?m)^[ \t]*static[ \t]+const[ \t]+"
    r"amrex::Vector<DiagnosticDescriptor>[ \t]+catalog[ \t]*\{"
)
CATALOG_RETURN_RE = re.compile(r"(?m)^[ \t]*return[ \t]+catalog[ \t]*;")
CANDIDATE_RE = re.compile(
    r"(?m)^[ \t]*\{[ \t]*DiagnosticID::(?P<diagnostic_id>[A-Za-z_]\w*)\b"
)

# This expression handles only the fixed catalog's six-field initializer
# layout. It is not intended to parse arbitrary C++.
DESCRIPTOR_RE = re.compile(
    r"(?m)^[ \t]*\{[ \t]*DiagnosticID::(?P<diagnostic_id>[A-Za-z_]\w*)"
    r"\s*,\s*\"(?P<name>[^\"]*)\"\s*,\s*"
    r"\"(?P<long_name>[^\"]*)\"\s*,\s*"
    r"\"(?P<units>[^\"]*)\"\s*,\s*"
    r"DiagnosticCategory::(?P<category>[A-Za-z_]\w*)\s*,\s*"
    r"MissingPolicy::(?P<policy>[A-Za-z_]\w*)\s*\}[ \t]*,?"
)
RST_ROW_RE = re.compile(r"^[ \t]*\* - ``(?P<name>[^`]+)``[ \t]*$")
RST_CELL_RE = re.compile(r"^[ \t]*- (?P<value>.*)$")


class CatalogCheckError(Exception):
    """A source or documentation contract failed validation."""


@dataclass(frozen=True)
class Descriptor:
    """One fixed descriptor and its source line."""

    diagnostic_id: str
    name: str
    long_name: str
    units: str
    category: str
    policy: str
    source_line: int


@dataclass(frozen=True)
class MarkedRegion:
    """The lines and bounds of one uniquely marked source region."""

    lines: Tuple[str, ...]
    begin_line: int
    end_line: int


def read_utf8(path: Path) -> str:
    """Read ``path`` as UTF-8 and convert read failures to check errors."""

    try:
        return path.read_text(encoding="utf-8")
    except OSError as exc:
        raise CatalogCheckError(f"could not read {path}: {exc}") from exc


def extract_marked_region(
    lines: Sequence[str], begin_marker: str, end_marker: str, path: Path
) -> MarkedRegion:
    """Return the unique ordered region between exact marker lines."""

    begin_positions = [i for i, line in enumerate(lines) if line == begin_marker]
    end_positions = [i for i, line in enumerate(lines) if line == end_marker]
    if len(begin_positions) != 1:
        raise CatalogCheckError(
            f"expected exactly one begin marker in {path}; found {len(begin_positions)}"
        )
    if len(end_positions) != 1:
        raise CatalogCheckError(
            f"expected exactly one end marker in {path}; found {len(end_positions)}"
        )
    begin = begin_positions[0]
    end = end_positions[0]
    if begin > end:
        raise CatalogCheckError(
            f"reversed markers in {path}: begin line {begin + 1}, end line {end + 1}"
        )
    return MarkedRegion(tuple(lines[begin + 1 : end]), begin + 1, end + 1)


def source_line(text: str, offset: int) -> int:
    """Return the one-based source line containing ``offset``."""

    return text.count("\n", 0, offset) + 1


def parse_cpp_catalog(path: Path) -> List[Descriptor]:
    """Parse only the fixed initializer inside ``catalog_storage()``."""

    source = read_utf8(path)
    function = CATALOG_FUNCTION_RE.search(source)
    if function is None:
        raise CatalogCheckError(f"could not locate catalog_storage() in {path}")

    initializer = CATALOG_INITIALIZER_RE.search(source, function.end())
    returned = CATALOG_RETURN_RE.search(source, function.end())
    if initializer is None or returned is None or initializer.start() > returned.start():
        raise CatalogCheckError(f"could not locate fixed catalog initializer in {path}")

    body = source[initializer.end() : returned.start()]
    candidates = list(CANDIDATE_RE.finditer(body))
    parsed = list(DESCRIPTOR_RE.finditer(body))
    candidate_starts = [match.start() for match in candidates]
    parsed_starts = [match.start() for match in parsed]
    if len(candidates) != len(parsed) or candidate_starts != parsed_starts:
        raise CatalogCheckError(
            f"found {len(candidates)} candidate C++ descriptors but parsed "
            f"{len(parsed)} in {path}"
        )
    if not candidates:
        raise CatalogCheckError(f"fixed catalog contains no descriptors in {path}")

    descriptors = []
    for match in parsed:
        values = match.groupdict()
        descriptors.append(
            Descriptor(
                diagnostic_id=values["diagnostic_id"],
                name=values["name"],
                long_name=values["long_name"],
                units=values["units"],
                category=values["category"],
                policy=values["policy"],
                source_line=source_line(source, initializer.end() + match.start()),
            )
        )

    ids = [descriptor.diagnostic_id for descriptor in descriptors]
    duplicate_ids = sorted({item for item in ids if ids.count(item) > 1})
    if duplicate_ids:
        raise CatalogCheckError(f"duplicate C++ diagnostic ID(s): {', '.join(duplicate_ids)}")
    names = [descriptor.name for descriptor in descriptors]
    duplicate_names = sorted({item for item in names if names.count(item) > 1})
    if duplicate_names:
        raise CatalogCheckError(f"duplicate C++ public name(s): {', '.join(duplicate_names)}")
    return descriptors


def strip_markup(value: str) -> str:
    """Remove only the RST literal markup used by catalog cells."""

    value = value.strip()
    if value.startswith("``") and value.endswith("``"):
        return value[2:-2]
    return value


def parse_rst_catalog(text: str, path: Path) -> List[Descriptor]:
    """Parse five-cell descriptor rows from the unique fixed-catalog region."""

    region = extract_marked_region(
        text.splitlines(), FIXED_CATALOG_BEGIN, FIXED_CATALOG_END, path
    )
    rows = []
    i = 0
    while i < len(region.lines):
        match = RST_ROW_RE.match(region.lines[i])
        if match is None:
            i += 1
            continue
        values = [match.group("name")]
        for offset in range(1, 5):
            row_line = i + offset
            if row_line >= len(region.lines):
                raise CatalogCheckError(
                    f"{path} line {region.begin_line + i}: "
                    f"truncated RST row for {values[0]!r}"
                )
            cell = RST_CELL_RE.match(region.lines[row_line])
            if cell is None:
                raise CatalogCheckError(
                    f"{path} line {region.begin_line + row_line}: "
                    f"unparsable RST row for {values[0]!r}"
                )
            values.append(strip_markup(cell.group("value")))
        next_line = i + 5
        if next_line < len(region.lines) and RST_CELL_RE.match(region.lines[next_line]):
            raise CatalogCheckError(
                f"{path} line {region.begin_line + next_line}: "
                f"too many cells in RST row for {values[0]!r}"
            )
        rows.append(
            Descriptor(
                diagnostic_id="",
                name=values[0],
                category=values[1],
                units=values[2],
                policy=values[3],
                long_name=values[4],
                source_line=region.begin_line + i,
            )
        )
        i += 5
    return rows


def parse_dynamic_soil_patterns(text: str, path: Path) -> List[str]:
    """Extract first-column code literals from the marked soil-family table."""

    region = extract_marked_region(
        text.splitlines(), DYNAMIC_SOIL_BEGIN, DYNAMIC_SOIL_END, path
    )
    return [
        match.group("name")
        for line in region.lines
        if (match := RST_ROW_RE.match(line)) is not None
    ]


def validate_catalog(expected: Sequence[Descriptor], documented: Sequence[Descriptor]) -> None:
    """Require exact names, order, and metadata for the fixed catalog."""

    expected_names = [descriptor.name for descriptor in expected]
    documented_names = [descriptor.name for descriptor in documented]
    duplicate_names = sorted(
        {item for item in documented_names if documented_names.count(item) > 1}
    )
    if duplicate_names:
        raise CatalogCheckError(
            f"duplicate documented name(s): {', '.join(repr(item) for item in duplicate_names)}"
        )

    missing = [name for name in expected_names if name not in documented_names]
    extra = [name for name in documented_names if name not in expected_names]
    if missing or extra:
        details = []
        if missing:
            details.append(f"missing {', '.join(repr(item) for item in missing)}")
        if extra:
            details.append(f"extra {', '.join(repr(item) for item in extra)}")
        raise CatalogCheckError("documented catalog names differ: " + "; ".join(details))

    if expected_names != documented_names:
        for position, (expected_name, documented_name) in enumerate(
            zip(expected_names, documented_names), start=1
        ):
            if expected_name != documented_name:
                raise CatalogCheckError(
                    f"documented catalog order differs at position {position}: "
                    f"expected {expected_name!r}, found {documented_name!r}"
                )

    documented_by_name = {descriptor.name: descriptor for descriptor in documented}
    errors = []
    for expected_descriptor in expected:
        documented_descriptor = documented_by_name[expected_descriptor.name]
        for field in ("category", "units", "policy", "long_name"):
            expected_value = getattr(expected_descriptor, field)
            documented_value = getattr(documented_descriptor, field)
            if expected_value != documented_value:
                errors.append(
                    f"{documented_descriptor.source_line}: {expected_descriptor.name!r} "
                    f"field {field!r}: expected {expected_value!r}, "
                    f"documented {documented_value!r}"
                )
    if errors:
        raise CatalogCheckError("; ".join(errors))


def validate_dynamic_soil_patterns(patterns: Sequence[str]) -> None:
    """Require the exact ordered dynamic-soil family sequence."""

    if list(patterns) != list(DYNAMIC_SOIL_PATTERNS):
        raise CatalogCheckError(
            "dynamic soil families differ: "
            f"expected {list(DYNAMIC_SOIL_PATTERNS)!r}, found {list(patterns)!r}"
        )


def check_catalog(cpp_path: Path, rst_path: Path) -> Tuple[int, int]:
    """Read both sources once and validate fixed and dynamic contracts."""

    expected = parse_cpp_catalog(cpp_path)
    rst_text = read_utf8(rst_path)
    documented = parse_rst_catalog(rst_text, rst_path)
    dynamic_patterns = parse_dynamic_soil_patterns(rst_text, rst_path)
    validate_catalog(expected, documented)
    validate_dynamic_soil_patterns(dynamic_patterns)
    return len(expected), len(dynamic_patterns)


def main() -> int:
    """Run the repository check and convert expected failures to status 1."""

    repo = Path(__file__).resolve().parents[3]
    cpp = repo / "Source/IO/ERF_Plotfile2DCatalog.cpp"
    rst = repo / "Docs/sphinx_doc/Plotfiles.rst"
    try:
        fixed_count, dynamic_count = check_catalog(cpp, rst)
    except CatalogCheckError as exc:
        print(f"plotfile2d catalog check: error: {exc}", file=sys.stderr)
        return 1
    print(f"plotfile2d catalog check: {fixed_count} fixed descriptors match")
    print(f"plotfile2d catalog check: {dynamic_count} dynamic soil families match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
