"""Tests for the per-build ERF inputs-schema generator (Exec/Tools/build_erf_schema.py).

Written before the implementation; each test fails until its capability exists. Run from the ERF root:
    python3 -m unittest discover -s Exec/Tools/erf_inputs_schema/tests
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

TOOLS = Path(__file__).resolve().parents[2]
if str(TOOLS) not in sys.path:
    sys.path.insert(0, str(TOOLS))

from erf_inputs_schema import canonical  # noqa: E402

FIXTURES = Path(__file__).resolve().parent / "fixtures"
SYNTHETIC = FIXTURES / "synthetic"


def generate(**kwargs):
    from erf_inputs_schema import schema
    return schema.generate(SYNTHETIC / "erf", SYNTHETIC / "AMReX", **kwargs)


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


class BuildRequirementTests(unittest.TestCase):
    def test_a_key_read_only_under_ifdef_needs_that_feature(self):
        fact = generate()["keys"]["erf.widget_on"]
        self.assertEqual(fact["source"], "Source/ERF_Synthetic.cpp:12")
        self.assertEqual(fact["requires"], ["ERF_USE_WIDGET"])
        self.assertEqual(fact["requires_evidence"], {"ERF_USE_WIDGET": "Source/ERF_Synthetic.cpp:10"})

    def test_an_unconditional_key_needs_nothing(self):
        self.assertNotIn("requires", generate()["keys"]["erf.cfl"])

    def test_the_feature_follows_the_build_options(self):
        off = generate()
        on = generate(cmake_options={"ERF_ENABLE_WIDGET": True})
        self.assertIs(off["features"]["ERF_USE_WIDGET"], False)
        self.assertIs(on["features"]["ERF_USE_WIDGET"], True)
        self.assertEqual(on["feature_sources"]["ERF_USE_WIDGET"], "CMakeLists.txt:3")


class LayerTests(unittest.TestCase):
    def test_a_shared_key_records_both_layers_and_which_won(self):
        fact = generate()["keys"]["amr.max_level"]
        self.assertEqual(fact["layers"], {"AMReX": "AMReX/Src/Base/AMReX_ParallelDescriptor.cpp:9",
                                          "ERF": "Source/ERF_Synthetic.cpp:19"})
        self.assertEqual(fact["layer"], "AMReX")
        self.assertEqual(fact["source"], fact["layers"]["AMReX"])

    def test_a_key_of_one_layer_names_that_layer(self):
        keys = generate()["keys"]
        self.assertEqual(keys["erf.cfl"]["layer"], "ERF")
        self.assertEqual(keys["erf.cfl"]["layers"], {"ERF": "Source/ERF_Synthetic.cpp:6"})
        self.assertEqual(keys["amrex.use_gpu_aware_mpi"]["layer"], "AMReX")


CLI = TOOLS / "build_erf_schema.py"
GIT_ENV = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid", "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@example.invalid", "GIT_AUTHOR_DATE": "2026-10-29T00:00:00+00:00",
           "GIT_COMMITTER_DATE": "2026-10-29T00:00:00+00:00"}


def git(repo, *args):
    return subprocess.run(["git", "-c", "core.hooksPath=/dev/null", "-C", str(repo), *args], check=True,
                          capture_output=True, text=True, env={**os.environ, **GIT_ENV}).stdout.strip()


def git_repo(root, tree):
    """A one-commit repo holding a copy of tree; returns its commit."""
    shutil.copytree(tree, root)
    git(root, "init", "-q")
    git(root, "add", ".")
    git(root, "commit", "-q", "-m", "fixture")
    return git(root, "rev-parse", "HEAD")


def run_cli(*args):
    return subprocess.run([sys.executable, str(CLI), *map(str, args)], check=True, capture_output=True, text=True)


class DeterminismTests(unittest.TestCase):
    def test_two_runs_give_identical_bytes(self):
        with tempfile.TemporaryDirectory() as td:
            outs = [Path(td) / f"s{i}.json" for i in (1, 2)]
            for out in outs:
                run_cli("--erf-dir", SYNTHETIC / "erf", "--amrex-dir", SYNTHETIC / "AMReX", "--out", out)
            self.assertEqual(outs[0].read_bytes(), outs[1].read_bytes())
            self.assertTrue(outs[0].read_bytes().endswith(b"}\n"))

    def test_git_mode_reads_by_archive_and_takes_amrex_from_the_gitlink(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            amrex = git_repo(td / "amrex", SYNTHETIC / "AMReX")
            erf = td / "erf"
            git_repo(erf, SYNTHETIC / "erf")
            git(erf, "update-index", "--add", "--cacheinfo", f"160000,{amrex},Submodules/AMReX")
            git(erf, "commit", "-q", "-m", "gitlink")
            commit = git(erf, "rev-parse", "HEAD")
            out = td / "s.json"
            run_cli("--repo", erf, "--erf", commit, "--amrex-repo", td / "amrex", "--out", out)
            import json
            got = json.loads(out.read_text())
            self.assertEqual((got["erf_commit"], got["amrex_commit"]), (commit, amrex))
            self.assertEqual(got["keys"], generate()["keys"])


class ParityTests(unittest.TestCase):
    """The synthetic tree's schema, key for key: Build B's fields (type, source, per_level, enum values,
    case_insensitive, requires) plus this tool's (layer, layers, default, default_origin, required)."""

    def test_every_key_matches_the_expected_record(self):
        import json
        expected = json.loads((FIXTURES / "synthetic_expected_keys.json").read_text())
        got = generate()["keys"]
        self.assertEqual(sorted(got), sorted(expected))
        for key, record in expected.items():
            with self.subTest(key=key):
                self.assertEqual(got[key], record)


class DeclarationTests(unittest.TestCase):
    def test_every_name_on_a_multi_declarator_line_is_typed(self):
        keys = generate()["keys"]
        self.assertEqual(keys["erf.eye_latitude"]["type"], "real")
        self.assertEqual(keys["erf.eye_longitude"]["type"], "real")
        self.assertEqual(keys["erf.eye_longitude"]["default"], "-10000000000.0")


class RemovedKeyTests(unittest.TestCase):
    def test_a_key_whose_presence_aborts_as_removed_is_marked_removed(self):
        fact = generate()["keys"]["erf.old_knob"]
        self.assertEqual(fact["type"], "removed")
        self.assertEqual(fact["removed_evidence"], "Source/ERF_Removed.cpp:5")

    def test_contains_without_an_abort_is_not_removed(self):
        fact = generate()["keys"]["erf.probe_only"]
        self.assertEqual(fact["type"], "unknown")
        self.assertNotIn("removed_evidence", fact)

    def test_query_greater_than_zero_then_error_is_removed(self):
        fact = generate()["keys"]["erf.use_terrain"]
        self.assertEqual((fact["type"], fact["removed_shape"]), ("removed", "query_gt0_error"))
        self.assertEqual(fact["removed_evidence"], "Source/ERF_Removed.cpp:17")

    def test_countval_greater_than_zero_then_error_is_removed(self):
        fact = generate()["keys"]["erf.no_substepping"]
        self.assertEqual((fact["type"], fact["removed_shape"]), ("removed", "countval_gt0_error"))

    def test_contains_guarding_a_typed_read_is_not_removed(self):
        fact = generate()["keys"]["erf.also_typed"]
        self.assertEqual(fact["type"], "int")
        self.assertNotIn("removed_shape", fact)

    def test_removed_wins_over_a_typed_read_elsewhere_and_keeps_both_lines(self):
        fact = generate()["keys"]["erf.old_knob"]
        self.assertEqual(fact["type"], "removed")
        self.assertEqual(fact["removed_evidence"], "Source/ERF_Removed.cpp:5")
        self.assertEqual(fact["shadowed_reads"], ["Source/ERF_Shadow.cpp:5 int"])


class SeveralKindsTests(unittest.TestCase):
    """One name with several kinds in a file (the owner's rule, 10-30): keep the kinds the read call takes, then the
    nearest declaration before the read, otherwise unknown."""

    def test_a_list_read_takes_the_list_kind(self):
        self.assertEqual(generate()["keys"]["xlo.velocity"]["type"], "real_list")

    def test_a_scalar_read_takes_the_scalar_kind(self):
        self.assertEqual(generate()["keys"]["xlo.speed"]["type"], "real")

    def test_the_nearest_declaration_before_the_read_decides(self):
        keys = generate()["keys"]
        self.assertEqual((keys["xlo.k_int"]["type"], keys["xlo.k_real"]["type"]), ("int", "real"))


CHECKER = TOOLS / "check_erf_inputs.py"
DECK = """erf.old_knob = 1
erf.cfl = fast
erf.widget_on = true
erf.notakey = 1
erf.nsteps = 20
"""


def check_deck(schema_args):
    import json
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        (td / "s.json").write_text(json.dumps(generate()))
        (td / "inputs").write_text(DECK)
        args = [sys.executable, str(CHECKER), "--json", *[a.replace("S", str(td / "s.json")) for a in schema_args],
                str(td / "inputs")]
        out = subprocess.run(args, capture_output=True, text=True).stdout
    return {(f["checker"], f["key"]): f["severity"] for f in json.loads(out)["findings"]}


class CheckerSchemaTests(unittest.TestCase):
    def test_the_checker_uses_a_given_schema(self):
        got = {k: v for (c, k), v in check_deck(["--schema", "S"]).items() if c == "schema"}
        self.assertEqual(got, {"erf.old_knob": "ERROR", "erf.cfl": "ERROR", "erf.widget_on": "ERROR",
                               "erf.notakey": "WARNING"})

    def test_without_a_schema_the_schema_check_does_not_run(self):
        self.assertFalse([k for (c, k) in check_deck([]) if c == "schema"])

    def test_the_source_derived_removed_set_is_cross_checked_against_retired_keys(self):
        from erf_input_checker.checks import schema as schema_check
        curated_only, source_only = schema_check.cross_check(generate(), {"erf.old_knob": "x", "erf.gone": "y"})
        self.assertEqual(curated_only, ["erf.gone"])
        self.assertEqual(source_only, [("erf.no_substepping", "Source/ERF_Removed.cpp:20"),
                                       ("erf.use_terrain", "Source/ERF_Removed.cpp:17")])


SITE_FACTS = [{"feature": "mpi_gpu_aware", "on": False, "key": "amrex.use_gpu_aware_mpi", "value": "true",
               "evidence": "site: the MPI is configured without CUDA"}]


class SiteFactTests(unittest.TestCase):
    def test_a_site_fact_sets_a_feature_and_a_value_requirement(self):
        s = generate(site_facts=SITE_FACTS)
        self.assertIs(s["features"]["mpi_gpu_aware"], False)
        need = s["keys"]["amrex.use_gpu_aware_mpi"]["value_requires"]["true"]
        self.assertEqual(need["feature"], "mpi_gpu_aware")
        self.assertEqual(need["evidence"], "site: the MPI is configured without CUDA; "
                                           "AMReX/Src/Base/AMReX_ParallelDescriptor.cpp:6")

    def test_the_schema_records_where_its_build_facts_came_from(self):
        with tempfile.TemporaryDirectory() as td:
            out, facts = Path(td) / "s.json", Path(td) / "site.json"
            import json
            facts.write_text(json.dumps(SITE_FACTS))
            run_cli("--erf-dir", SYNTHETIC / "erf", "--amrex-dir", SYNTHETIC / "AMReX", "--out", out,
                    "--cmake=-DERF_ENABLE_WIDGET=ON", "--site-facts", facts)
            got = json.loads(out.read_text())["build"]
        self.assertEqual(got["cmake"], {"ERF_ENABLE_WIDGET": True})
        self.assertEqual(got["facts_from"], {"cmake": "d_list", "site_facts": SITE_FACTS})


class MoreShapesTests(unittest.TestCase):
    """Shapes the 87f45a6e parity run found (ERF_DataStruct.H, AMReX_Amr.cpp, ERF_WriteBndryPlanes.cpp)."""

    def test_a_count_kept_in_a_variable_then_error_is_removed(self):
        fact = generate()["keys"]["erf.old_flag"]
        self.assertEqual((fact["type"], fact["removed_shape"]), ("removed", "countval_var_gt0_error"))
        self.assertEqual(fact["removed_evidence"], "Source/ERF_MoreShapes.cpp:6")

    def test_contains_whose_body_reads_the_key_is_not_removed_despite_an_error(self):
        fact = generate()["keys"]["erf.legacy"]
        self.assertEqual(fact["type"], "int")
        self.assertNotIn("removed_shape", fact)

    def test_constructor_style_declarators_are_typed(self):
        keys = generate()["keys"]
        self.assertEqual((keys["erf.chv"]["type"], keys["erf.chv"]["default"]), ("int", "1"))
        self.assertEqual(keys["erf.box_hi"]["type"], "real_list")


if __name__ == "__main__":
    unittest.main()
