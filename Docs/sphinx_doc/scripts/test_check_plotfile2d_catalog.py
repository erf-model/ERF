#!/usr/bin/env python3
"""Regression tests for the 2-D catalog drift checker."""

from dataclasses import replace
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import check_plotfile2d_catalog as checker  # noqa: E402


CPP_HEADER = """\
const amrex::Vector<DiagnosticDescriptor>&
catalog_storage ()
{
    static const amrex::Vector<DiagnosticDescriptor> catalog{
"""
CPP_FOOTER = """\
    };

    return catalog;
}
"""


def cpp_entry(
    diagnostic_id: str,
    name: str,
    long_name: str,
    units: str,
    category: str = "Geometry",
    policy: str = "AlwaysAvailable",
) -> str:
    """Build one fixture descriptor in the production initializer layout."""

    return (
        f'        {{DiagnosticID::{diagnostic_id}, "{name}", "{long_name}", '
        f'"{units}", DiagnosticCategory::{category}, MissingPolicy::{policy}}},\n'
    )


def cpp_fixture() -> str:
    """Return a minimal exact C++ catalog fixture."""

    return CPP_HEADER + "".join(
        (
            cpp_entry("Alpha", "alpha", "Alpha field", "K"),
            cpp_entry("Beta", "beta", "Beta field", "m", "SurfaceLayer", "FillMinus999WhenUnavailable"),
        )
    ) + CPP_FOOTER


def rst_fixture(
    names=("alpha", "beta"),
    category="Geometry",
    units=("K", "m"),
    policy=("AlwaysAvailable", "FillMinus999WhenUnavailable"),
    long_names=("Alpha field", "Beta field"),
    dynamic=checker.DYNAMIC_SOIL_PATTERNS,
) -> str:
    """Return a short fixed-catalog and dynamic-soil RST fixture."""

    lines = [
        checker.FIXED_CATALOG_BEGIN,
        "",
        ".. list-table::",
        "   :header-rows: 1",
        "",
        "   * - Variable",
        "     - Category",
        "     - Units",
        "     - Metadata policy",
        "     - Description",
    ]
    categories = (category, "SurfaceLayer")
    for index, name in enumerate(names):
        lines.extend(
            (
                f"   * - ``{name}``",
                f"     - ``{categories[index]}``",
                f"     - ``{units[index]}``",
                f"     - ``{policy[index]}``",
                f"     - {long_names[index]}",
            )
        )
    lines.extend(
        (
            checker.FIXED_CATALOG_END,
            "",
            checker.DYNAMIC_SOIL_BEGIN,
            "",
            ".. list-table::",
            "   :header-rows: 1",
            "",
            "   * - Name family",
            "     - Meaning",
        )
    )
    for name in dynamic:
        lines.extend((f"   * - ``{name}``", "     - Test family"))
    lines.extend((checker.DYNAMIC_SOIL_END, ""))
    return "\n".join(lines)


class CheckerTests(unittest.TestCase):
    def write_pair(self, cpp_text=None, rst_text=None):
        """Write fixture sources and return their temporary paths."""

        directory = tempfile.TemporaryDirectory()
        root = Path(directory.name)
        cpp = root / "catalog.cpp"
        rst = root / "Plotfiles.rst"
        cpp.write_text(cpp_text or cpp_fixture(), encoding="utf-8")
        rst.write_text(rst_text or rst_fixture(), encoding="utf-8")
        self.addCleanup(directory.cleanup)
        return cpp, rst

    def test_minimal_exact_pair_passes(self):
        cpp, rst = self.write_pair()
        self.assertEqual(checker.check_catalog(cpp, rst), (2, 3))

    def test_duplicate_fixed_begin_marker_fails(self):
        lines = [checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END]
        with self.assertRaisesRegex(checker.CatalogCheckError, "exactly one begin marker"):
            checker.extract_marked_region(lines, checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END, Path("Plotfiles.rst"))

    def test_duplicate_fixed_end_marker_fails(self):
        lines = [checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END, checker.FIXED_CATALOG_END]
        with self.assertRaisesRegex(checker.CatalogCheckError, "exactly one end marker"):
            checker.extract_marked_region(lines, checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END, Path("Plotfiles.rst"))

    def test_missing_or_reversed_markers_fail(self):
        path = Path("Plotfiles.rst")
        with self.assertRaisesRegex(checker.CatalogCheckError, "begin marker"):
            checker.extract_marked_region([checker.FIXED_CATALOG_END], checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END, path)
        with self.assertRaisesRegex(checker.CatalogCheckError, "reversed markers"):
            checker.extract_marked_region([checker.FIXED_CATALOG_END, checker.FIXED_CATALOG_BEGIN], checker.FIXED_CATALOG_BEGIN, checker.FIXED_CATALOG_END, path)

    def test_unparsed_cpp_candidate_fails_closed(self):
        malformed = cpp_entry("Beta", "beta", "Beta field", "m").replace(
            "MissingPolicy::AlwaysAvailable", "MissingPolicy::"
        )
        with self.assertRaisesRegex(checker.CatalogCheckError, "candidate C\+\+ descriptors"):
            checker.parse_cpp_catalog(
                self.write_pair(cpp_text=CPP_HEADER + cpp_entry("Alpha", "alpha", "Alpha field", "K") + malformed + CPP_FOOTER)[0]
            )

    def test_duplicate_cpp_ids_fail(self):
        text = CPP_HEADER + cpp_entry("Alpha", "alpha", "Alpha field", "K") + cpp_entry("Alpha", "beta", "Beta field", "m") + CPP_FOOTER
        with self.assertRaisesRegex(checker.CatalogCheckError, "duplicate C\+\+ diagnostic ID"):
            checker.parse_cpp_catalog(self.write_pair(cpp_text=text)[0])

    def test_duplicate_cpp_names_fail(self):
        text = CPP_HEADER + cpp_entry("Alpha", "alpha", "Alpha field", "K") + cpp_entry("Beta", "alpha", "Beta field", "m") + CPP_FOOTER
        with self.assertRaisesRegex(checker.CatalogCheckError, "duplicate C\+\+ public name"):
            checker.parse_cpp_catalog(self.write_pair(cpp_text=text)[0])

    def parsed_fixture(self):
        cpp, rst = self.write_pair()
        return checker.parse_cpp_catalog(cpp), checker.parse_rst_catalog(rst.read_text(encoding="utf-8"), rst)

    def test_duplicate_rst_names_fail(self):
        expected, documented = self.parsed_fixture()
        duplicate = replace(documented[1], name=documented[0].name)
        with self.assertRaisesRegex(checker.CatalogCheckError, "duplicate documented name"):
            checker.validate_catalog(expected, [documented[0], duplicate])

    def test_missing_and_extra_names_fail(self):
        expected, documented = self.parsed_fixture()
        changed = replace(documented[1], name="gamma")
        with self.assertRaisesRegex(checker.CatalogCheckError, "missing 'beta'.*extra 'gamma'"):
            checker.validate_catalog(expected, [documented[0], changed])

    def test_reordered_names_fail(self):
        expected, documented = self.parsed_fixture()
        with self.assertRaisesRegex(checker.CatalogCheckError, "order differs at position 1"):
            checker.validate_catalog(expected, list(reversed(documented)))

    def test_metadata_mismatches_are_reported(self):
        expected, documented = self.parsed_fixture()
        for field, value in (
            ("category", "Radiation"),
            ("units", "degC"),
            ("policy", "FillZeroWhenUnavailable"),
            ("long_name", "Wrong field"),
        ):
            with self.subTest(field=field):
                changed = replace(documented[0], **{field: value})
                with self.assertRaisesRegex(checker.CatalogCheckError, f"field '{field}'"):
                    checker.validate_catalog(expected, [changed, documented[1]])

    def test_dynamic_text_outside_region_does_not_satisfy_check(self):
        text = rst_fixture(dynamic=checker.DYNAMIC_SOIL_PATTERNS) + "\nsmois_<layer>\n"
        path = self.write_pair(rst_text=text)[1]
        patterns = checker.parse_dynamic_soil_patterns(path.read_text(encoding="utf-8"), path)
        checker.validate_dynamic_soil_patterns(patterns)

    def test_dynamic_patterns_must_match_exact_order(self):
        for patterns in (
            checker.DYNAMIC_SOIL_PATTERNS[:2],
            ("smois_<layer>", "smois_<layer>", "tslb_<layer>"),
            (*checker.DYNAMIC_SOIL_PATTERNS, "extra_<layer>"),
            ("sh2o_<layer>", "smois_<layer>", "tslb_<layer>"),
        ):
            with self.subTest(patterns=patterns):
                with self.assertRaisesRegex(checker.CatalogCheckError, "dynamic soil families differ"):
                    checker.validate_dynamic_soil_patterns(patterns)

    def test_repository_catalog_passes(self):
        repo = SCRIPT_DIR.parents[2]
        fixed_count, dynamic_count = checker.check_catalog(
            repo / "Source/IO/ERF_Plotfile2DCatalog.cpp",
            repo / "Docs/sphinx_doc/Plotfiles.rst",
        )
        self.assertGreater(fixed_count, 0)
        self.assertEqual(dynamic_count, len(checker.DYNAMIC_SOIL_PATTERNS))


if __name__ == "__main__":
    unittest.main()
