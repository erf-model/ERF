"""Tests for the per-build ERF inputs-schema generator (Exec/Tools/build_erf_schema.py).

Written before the implementation; each test fails until its capability exists. Run from the ERF root:
    python3 -m unittest discover -s Exec/Tools/erf_inputs_schema/tests
"""
import sys
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2]
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from erf_inputs_schema import canonical  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"


class CanonicalValueTests(unittest.TestCase):
    def test_real_spellings_normalize_to_one_form(self):
        forms = {canonical.canonical_value("real", w, {}) for w in ("1", "1.0", "1e0", "1.", "+1.0E+00")}
        self.assertEqual(forms, {"1.0"})

    def test_int_bool_and_enum_forms(self):
        self.assertEqual(canonical.canonical_value("int", "+03", {}), "3")
        self.assertEqual(canonical.canonical_value("bool", "T", {}), "true")
        self.assertEqual(canonical.canonical_value("bool", "0", {}), "false")
        fact = {"values": ["None", "Kessler"], "case_insensitive": True}
        self.assertEqual(canonical.canonical_value("enum", "kessler", fact), "Kessler")
        self.assertEqual(canonical.canonical_value("string", "a b", {}), "a b")

    def test_a_word_that_is_not_of_the_type_is_refused(self):
        with self.assertRaises(ValueError):
            canonical.canonical_value("real", "1.0.0", {})
        with self.assertRaises(ValueError):
            canonical.canonical_value("enum", "Morrison", {"values": ["None"], "case_insensitive": False})


if __name__ == "__main__":
    unittest.main()
