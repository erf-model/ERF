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


def cmake_features(cmakelists, options=None):
    """The features a one-file CMake tree gives, with these -D options."""
    from erf_inputs_schema import cmake
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        (root / "CMake").mkdir()
        (root / "CMakeLists.txt").write_text(cmakelists)
        return cmake.cmake_features(root, options or {})[0]


CHAIN = """option(ERF_A "a" {a})
option(ERF_B "b" {b})
option(ERF_C "c" {c})
if(ERF_A)
  target_compile_definitions(t PUBLIC ERF_USE_A)
elseif(ERF_B)
  target_compile_definitions(t PUBLIC ERF_USE_B)
elseif(ERF_C)
  target_compile_definitions(t PUBLIC ERF_USE_C)
else()
  target_compile_definitions(t PUBLIC ERF_USE_D)
endif()
"""


class CMakeBranchTests(unittest.TestCase):
    """Only one branch of an if/elseif/else chain is taken, as in CMake."""

    def test_after_a_taken_branch_no_later_branch_is_on(self):
        got = cmake_features(CHAIN.format(a="ON", b="ON", c="ON"))
        self.assertEqual(got, {"ERF_USE_A": True, "ERF_USE_B": False, "ERF_USE_C": False, "ERF_USE_D": False})

    def test_else_is_on_only_when_no_branch_was_taken(self):
        self.assertEqual(cmake_features(CHAIN.format(a="OFF", b="OFF", c="OFF"))["ERF_USE_D"], True)
        self.assertEqual(cmake_features(CHAIN.format(a="OFF", b="ON", c="OFF"))["ERF_USE_D"], False)

    def test_a_trailing_comment_does_not_hide_an_if(self):
        got = cmake_features('option(ERF_X "x" OFF)\nif(ERF_X) # optional\n'
                             '  target_compile_definitions(t PUBLIC ERF_USE_X)\nendif()\n')
        self.assertIs(got["ERF_USE_X"], False)

    def test_an_option_default_naming_another_option_follows_its_d_override(self):
        # ERF's own CMakeLists: option(ERF_ENABLE_HDF5 "Enable HDF5 IO" ${ERF_ENABLE_NETCDF})
        text = ('option(ERF_ENABLE_NETCDF "n" OFF)\noption(ERF_ENABLE_HDF5 "h" ${ERF_ENABLE_NETCDF})\n'
                'if(ERF_ENABLE_HDF5)\n  target_compile_definitions(t PUBLIC ERF_USE_HDF5)\nendif()\n')
        self.assertIs(cmake_features(text, {"ERF_ENABLE_NETCDF": True})["ERF_USE_HDF5"], True)
        self.assertIs(cmake_features(text)["ERF_USE_HDF5"], False)

    def test_any_cmake_true_spelling_of_a_default_is_on(self):
        for word in ("ON", "TRUE", "YES", "1"):
            got = cmake_features(f'option(ERF_T "t" {word})\nif(ERF_T)\n  target_compile_definitions(t PUBLIC ERF_USE_T)\nendif()\n')
            self.assertIs(got["ERF_USE_T"], True, word)


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
    """The synthetic tree's schema, key for key: the reference schema's fields (type, source, per_level, enum values,
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


def tree_keys(source):
    """The keys a one-file ERF tree (Source/ERF_T.cpp) gives, with an empty AMReX."""
    from erf_inputs_schema import schema
    with tempfile.TemporaryDirectory() as td:
        erf, amrex = Path(td) / "erf", Path(td) / "AMReX"
        (erf / "Source").mkdir(parents=True)
        (erf / "CMake").mkdir()
        (amrex / "Src").mkdir(parents=True)
        (erf / "CMakeLists.txt").write_text("")
        (erf / "Source" / "ERF_T.cpp").write_text(source)
        return schema.generate(erf, amrex)["keys"]


class DefaultTests(unittest.TestCase):
    """A default is a literal of the key's type on the declaration the read fills; anything else is not a default."""

    def test_a_string_initialized_from_an_expression_has_no_default(self):
        keys = tree_keys('void f ()\n{\n    ParmParse pp("erf");\n    std::string pfx = default_prefix;\n'
                         '    pp.query("pfx", pfx);\n    std::string lit = "plain";\n    pp.query("lit", lit);\n}\n')
        self.assertEqual((keys["erf.pfx"]["default"], keys["erf.pfx"]["default_origin"]), (None, "expression"))
        self.assertEqual((keys["erf.lit"]["default"], keys["erf.lit"]["default_origin"]), ("plain", "initializer"))

    def test_an_initializer_of_an_earlier_declaration_is_not_the_default(self):
        keys = tree_keys('int k = 1;\nvoid f ()\n{\n    ParmParse pp("erf");\n    amrex::Real k;\n'
                         '    pp.query("k", k);\n}\n')
        self.assertEqual(keys["erf.k"]["type"], "real")
        self.assertEqual((keys["erf.k"]["default"], keys["erf.k"]["default_origin"]), (None, "none"))

    def test_a_comparison_is_not_an_initializer(self):
        # ERF_FireCouplingParams.H: std::string coupling_type = "lagged"; ... return coupling_type == "synchronous";
        keys = tree_keys('struct P {\n    std::string mode = "lagged";\n'
                         '    bool sync () const { return mode == "synchronous"; }\n};\n'
                         'void f (P& p)\n{\n    ParmParse pp("erf");\n    pp.query("mode", p.mode);\n}\n')
        self.assertEqual((keys["erf.mode"]["default"], keys["erf.mode"].get("default_source")), ("lagged", "Source/ERF_T.cpp:2"))


class EnumListTests(unittest.TestCase):
    """A Vector of an enum read with queryarr is a list of enum values, not one value."""

    def test_an_enum_list_read_keeps_the_list_type(self):
        fact = generate()["keys"]["erf.moisture_species"]
        self.assertEqual((fact["type"], fact["values"]), ("enum_list", ["None", "Kessler", "SAM"]))


CHECKER = TOOLS / "check_erf_inputs.py"
DECK = """erf.old_knob = 1
erf.cfl = fast
erf.widget_on = true
erf.notakey = 1
erf.nsteps = 20
erf.moisture_species = Kessler SAM
erf.probe_only = 1 2
erf.eye_latitude = 2*0.5
erf.no_substepping = true
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
                               "erf.notakey": "WARNING", "erf.eye_latitude": "WARNING"})

    def test_an_unknown_type_takes_any_number_of_words(self):
        self.assertNotIn(("schema", "erf.probe_only"), check_deck(["--schema", "S"]))

    def test_an_expression_is_a_warning_since_parmparse_evaluates_it(self):
        from erf_input_checker.checks import schema as schema_check
        real = {"type": "real"}
        self.assertIsNone(schema_check.value_problem(real, ["2*0.5"], 1))
        self.assertEqual(schema_check.expression_words(real, ["2*0.5", "1.0", "(a+b)", "${x}"]), ["2*0.5", "(a+b)", "${x}"])
        self.assertEqual(schema_check.value_problem(real, ["fast"], 1), "not a real literal: 'fast'")

    def test_a_retired_key_gets_one_error_not_two(self):
        found = check_deck(["--schema", "S"])
        self.assertEqual(found.get(("keys.retired", "erf.no_substepping")), "ERROR")
        self.assertNotIn(("schema", "erf.no_substepping"), found)

    def test_without_a_schema_the_schema_check_does_not_run(self):
        self.assertFalse([k for (c, k) in check_deck([]) if c == "schema"])

    def test_the_source_derived_removed_set_is_cross_checked_against_retired_keys(self):
        from erf_input_checker.checks import schema as schema_check
        curated_only, source_only = schema_check.cross_check(generate(), {"erf.old_knob": "x", "erf.gone": "y"})
        self.assertEqual(curated_only, ["erf.gone"])
        self.assertEqual(source_only, [("erf.no_substepping", "Source/ERF_Removed.cpp:20"),
                                       ("erf.old_flag", "Source/ERF_MoreShapes.cpp:6"),
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

    def test_a_site_fact_value_is_stored_in_canonical_form(self):
        for spelling in ("1", "TRUE", True):
            facts = [dict(SITE_FACTS[0], value=spelling)]
            got = generate(site_facts=facts)["keys"]["amrex.use_gpu_aware_mpi"]["value_requires"]
            self.assertEqual(list(got), ["true"], spelling)


class MergeOrderTests(unittest.TestCase):
    """A removed key's record names the removal, whichever file sorts first; types in shadowed_reads are schema types."""

    def test_a_removal_after_a_typed_read_names_the_removal_and_its_layer(self):
        from erf_inputs_schema import schema
        typed = {"type": "real", "source": "AMReX/Src/A.cpp:5"}
        removed = {"type": "removed", "source": "Source/ERF_Removed.cpp:4", "removed_evidence": "Source/ERF_Removed.cpp:5"}
        got = schema.merge(schema.merge(None, typed), removed)
        self.assertEqual((got["source"], got["layer"]), ("Source/ERF_Removed.cpp:4", "ERF"))
        self.assertEqual(got["layers"], {"AMReX": "AMReX/Src/A.cpp:5", "ERF": "Source/ERF_Removed.cpp:4"})
        self.assertEqual(got["shadowed_reads"], ["AMReX/Src/A.cpp:5 real"])

    def test_a_shadowed_enum_read_shows_the_schema_type(self):
        from erf_inputs_schema import schema
        removed = {"type": "removed", "source": "Source/ERF_Removed.cpp:4", "removed_evidence": "Source/ERF_Removed.cpp:5"}
        got = schema.merge(schema.merge(None, removed), {"type": "enum:MoistureType_list", "source": "Source/A.cpp:3"})
        self.assertEqual(got["shadowed_reads"], ["Source/A.cpp:3 enum_list"])

    def test_a_presence_read_then_a_case_insensitive_enum_read_keeps_case_insensitive(self):
        from erf_inputs_schema import schema
        unknown = {"type": "unknown", "source": "Source/A.cpp:2"}
        enum = {"type": "enum:MoistureType", "source": "Source/A.cpp:9", "case_insensitive": True}
        self.assertIs(schema.merge(schema.merge(None, unknown), enum).get("case_insensitive"), True)

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

    def test_a_size_check_inside_the_body_is_not_a_removal(self):
        fact = generate()["keys"]["erf.z_levels"]
        self.assertEqual(fact["type"], "real_list")
        self.assertNotIn("removed_shape", fact)

    def test_constructor_style_declarators_are_typed(self):
        keys = generate()["keys"]
        self.assertEqual((keys["erf.chv"]["type"], keys["erf.chv"]["default"]), ("int", "1"))
        self.assertEqual(keys["erf.box_hi"]["type"], "real_list")


def run_tool(script, *args):
    return subprocess.run([sys.executable, str(TOOLS / script), *map(str, args)], capture_output=True, text=True)


class CommandLineTests(unittest.TestCase):
    """Bad arguments give a usage error, not a traceback; the help names the whole tool."""

    def test_the_help_gives_the_whole_first_paragraph(self):
        out = run_tool("build_erf_schema.py", "--help").stdout
        self.assertIn("canonical form of each type", " ".join(out.split()))

    def test_mixed_git_and_directory_flags_are_a_usage_error(self):
        with tempfile.TemporaryDirectory() as td:
            r = run_tool("build_erf_schema.py", "--repo", ".", "--erf", "HEAD", "--erf-dir", td, "--out", Path(td) / "o")
        self.assertEqual(r.returncode, 2)
        self.assertNotIn("Traceback", r.stderr)

    def test_an_unreadable_schema_is_a_usage_error(self):
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            (td / "inputs").write_text("erf.nsteps = 1\n")
            (td / "bad.json").write_text("{not json")
            (td / "nokeys.json").write_text("{}")
            for schema_file in ("missing.json", "bad.json", "nokeys.json"):
                r = run_tool("check_erf_inputs.py", "--schema", td / schema_file, td / "inputs")
                self.assertEqual(r.returncode, 2, schema_file)
                self.assertNotIn("Traceback", r.stderr, schema_file)

    def test_min_severity_also_hides_the_info_cross_check(self):
        import json
        with tempfile.TemporaryDirectory() as td:
            td = Path(td)
            (td / "s.json").write_text(json.dumps(generate()))
            (td / "inputs").write_text("erf.nsteps = 1\n")
            out = run_tool("check_erf_inputs.py", "--schema", td / "s.json", "--min-severity", "ERROR", td / "inputs").stdout
        self.assertNotIn("[INFO]", out)


class SourceTests(unittest.TestCase):
    def test_an_uninitialized_amrex_submodule_is_named(self):
        from erf_inputs_schema import source
        with tempfile.TemporaryDirectory() as td:
            repo = Path(td) / "erf"
            (repo / "Submodules" / "AMReX").mkdir(parents=True)
            (repo / "README").write_text("x\n")
            git(repo, "init", "-q")
            git(repo, "add", "README")
            git(repo, "commit", "-q", "-m", "x")
            with self.assertRaisesRegex(ValueError, "Submodules/AMReX is not initialized"):
                source.checkout(repo, "HEAD", Path(td) / "out")

    def test_extraction_works_without_tarfile_data_filter(self):
        """A Python older than 3.12 (and the 3.8-3.11 patch releases) has no data_filter and no filter= keyword."""
        import tarfile as real
        from types import SimpleNamespace
        from unittest import mock
        from erf_inputs_schema import source

        class OldTar:
            def __init__(self, tar):
                self.tar = tar

            def __enter__(self):
                return self

            def __exit__(self, *exc):
                self.tar.close()

            def extractall(self, path, members=None, **kw):
                if kw:
                    raise TypeError(f"extractall() got an unexpected keyword argument {next(iter(kw))!r}")
                self.tar.extractall(path, members, filter="data")

        old = SimpleNamespace(open=lambda *a, **k: OldTar(real.open(*a, **k)))
        with tempfile.TemporaryDirectory() as td:
            commit = git_repo(Path(td) / "erf", SYNTHETIC / "erf")
            with mock.patch.object(source, "tarfile", old):
                source.extract(Path(td) / "erf", commit, ["CMakeLists.txt"], Path(td) / "out")
            self.assertTrue((Path(td) / "out" / "CMakeLists.txt").exists())

    def test_relative_paths_are_posix_on_any_platform(self):
        from pathlib import PureWindowsPath
        from erf_inputs_schema import scan
        self.assertEqual(scan.posix_rel(PureWindowsPath("C:/erf/Source/Init/A.cpp"), PureWindowsPath("C:/erf")),
                         "Source/Init/A.cpp")


if __name__ == "__main__":
    unittest.main()
