"""Mutate real workflows to prove cache contracts catch wiring drift."""

import copy
import unittest

import cache_workflows as contracts


class MatrixTests(unittest.TestCase):
    def test_cartesian_exclude_and_include(self):
        # Added fields can be overwritten; original axis values cannot.
        # Extra include rows do not receive later include merges.
        matrix = {
            "fruit": ["apple", "pear"], "animal": ["cat", "dog"],
            "exclude": [{"fruit": "pear", "animal": "dog"}],
            "include": [
                {"color": "green"}, {"color": "pink", "animal": "cat"},
                {"fruit": "apple", "shape": "circle"},
                {"fruit": "banana"}, {"fruit": "banana", "animal": "cat"},
            ],
        }
        self.assertEqual(contracts.expand_matrix(matrix, "push"), [
            {"fruit": "apple", "animal": "cat", "color": "pink", "shape": "circle"},
            {"fruit": "apple", "animal": "dog", "color": "green", "shape": "circle"},
            {"fruit": "pear", "animal": "cat", "color": "pink"},
            {"fruit": "banana"}, {"fruit": "banana", "animal": "cat"},
        ])

    def test_include_only(self):
        rows = [{"os": "ubuntu-22.04"}, {"os": "macos-latest"}]
        self.assertEqual(contracts.expand_matrix({"include": rows}, "push"), rows)

    def test_include_can_restore_excluded_row(self):
        matrix = {"x": ["a", "b"], "exclude": [{"x": "b"}], "include": [{"x": "b"}]}
        self.assertEqual(contracts.expand_matrix(matrix, "push"), [{"x": "a"}, {"x": "b"}])

    def test_pr_dependent_axis(self):
        axis = "${{ fromJSON(github.event_name == 'pull_request' && '[\"ON\"]' || '[\"OFF\", \"ON\"]') }}"
        self.assertEqual(contracts.expand_matrix({"x": axis}, "pull_request"), [{"x": "ON"}])
        for event in ("push", "workflow_dispatch"):
            self.assertEqual(contracts.expand_matrix({"x": axis}, event), [{"x": "OFF"}, {"x": "ON"}])

    def test_unknown_and_empty_matrices_fail(self):
        for matrix in ("${{ needs.setup.outputs.matrix }}", {"x": "${{ inputs.axes }}"},
                       {"x": []}, {"x": ["a"], "exclude": [{"x": "a"}]}):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                contracts.expand_matrix(matrix, "push")


class WorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original, cls.defaults = contracts.repository_inputs()

    def setUp(self):
        self.workflows = copy.deepcopy(self.original)

    def job(self, name="ci.yml"):
        return self.workflows[name]["jobs"][contracts.REQUIRED_JOBS[name]]

    def calls(self, name="ci.yml"):
        return [s for s in self.job(name)["steps"] if contracts.is_cache(s)]

    def fallback(self):
        return next(s for s in self.job()["steps"] if s.get("uses", "").startswith("actions/cache/restore@"))

    def check(self):
        return contracts.check_workflows(self.workflows, self.defaults)

    def test_repository_contracts(self):
        counts, tags = self.check()
        # Inventory assertions, not limits imposed by the checker.
        self.assertEqual(counts, {"push": 21, "pull_request": 9, "workflow_dispatch": 21})
        self.assertEqual(len(tags), 23)
        self.assertIn("ubuntu-22.04-x64-gcc-mpich-double-particles-off", tags)
        for arch in ("x64", "arm64"):
            self.assertIn(f"macos-latest-{arch}-appleclang-openmpi-double-particles-on", tags)
        for mpi in ("msmpi", "mpi-off"):
            self.assertIn(f"windows-latest-x64-msvc-{mpi}-double-particles-on", tags)

    def test_lost_key_dimension(self):
        self.job()["env"]["CCACHE_TRUNK_KEY"] = "ubuntu-22.04-x64-gcc-mpich-double"
        with self.assertRaisesRegex(ValueError, "ci.yml:.*cache tag collision"):
            self.check()

    def test_hip_version_change_updates_installer_and_tags(self):
        _, old_tags = self.check()
        job = self.job("hip.yml")
        job["strategy"]["matrix"]["rocm_version"] = ["6.4.0"]
        _, new_tags = self.check()
        expected = {tag.replace("-rocm-6.3.2-", "-rocm-6.4.0-") for tag in old_tags}
        self.assertEqual(set(new_tags), expected)
        installer = next(s["run"] for s in job["steps"] if "dependencies_hip.sh" in s.get("run", ""))
        for event in ("push", "pull_request", "workflow_dispatch"):
            for row in contracts.expand_matrix(job["strategy"]["matrix"], event):
                command = contracts.render(installer, row, {}, {})
                self.assertEqual(command, 'Submodules/AMReX/.github/workflows/dependencies/dependencies_hip.sh "6.4.0"')

    def test_hip_installer_cannot_fall_back_to_default_or_another_version(self):
        path = "Submodules/AMReX/.github/workflows/dependencies/dependencies_hip.sh"
        for argument in ("", ' "6.3.2"', ' "${{ matrix.mesh_precision }}"'):
            with self.subTest(argument=argument):
                self.setUp()
                step = next(s for s in self.job("hip.yml")["steps"] if "dependencies_hip.sh" in s.get("run", ""))
                step["run"] = path + argument
                with self.assertRaisesRegex(ValueError, "HIP installer must take"):
                    self.check()

    def test_hip_key_cannot_hardcode_version(self):
        for version in ("6.3.2", "6.4.0"):
            with self.subTest(version=version):
                self.setUp()
                env = self.job("hip.yml")["env"]
                env["CCACHE_TRUNK_KEY"] = env["CCACHE_TRUNK_KEY"].replace("${{ matrix.rocm_version }}", version)
                with self.assertRaisesRegex(ValueError, "HIP cache key must use matrix.rocm_version"):
                    self.check()

    def test_hip_step_key_override_must_match_selected_version(self):
        job = self.job("hip.yml")
        wrong_key = job["env"]["CCACHE_TRUNK_KEY"].replace("${{ matrix.rocm_version }}", "6.4.0")
        for step in self.calls("hip.yml"):
            step["with"]["key"] = wrong_key
        with self.assertRaisesRegex(ValueError, "HIP cache key differs from selected ROCm version"):
            self.check()

    def test_hip_installer_must_appear_once(self):
        for count in (0, 2):
            with self.subTest(count=count):
                self.setUp()
                steps = self.job("hip.yml")["steps"]
                installer = next(s for s in steps if "dependencies_hip.sh" in s.get("run", ""))
                steps.remove(installer)
                steps.extend([copy.deepcopy(installer) for _ in range(count)])
                with self.assertRaisesRegex(ValueError, "HIP installer must take"):
                    self.check()

    def test_expanding_gpu_particles_requires_key_dimension(self):
        self.job("cuda-ci.yml")["strategy"]["matrix"]["particles"] = ["OFF", "ON"]
        with self.assertRaisesRegex(ValueError, "cuda-ci.yml:.*cache tag collision"):
            self.check()

    def test_case_and_sanitization_collisions(self):
        for variants in (["ON", "on"], ["a/b", "a b"], ["a?b", "a-b"]):
            with self.subTest(variants=variants):
                self.job()["strategy"]["matrix"]["particles"] = variants
                with self.assertRaisesRegex(ValueError, "cache tag collision"):
                    self.check()

    def test_cross_workflow_collision_and_new_job_discovery(self):
        duplicate = copy.deepcopy(self.workflows["ci.yml"])
        self.workflows["new-build.yml"] = duplicate
        with self.assertRaisesRegex(ValueError, "new-build.yml:.*cache tag collision"):
            self.check()
        duplicate["jobs"]["Build-And-Test"]["env"]["CCACHE_TRUNK_KEY"] += "-new"
        counts, tags = self.check()
        self.assertEqual(counts["push"], 23)
        self.assertEqual(len(tags), 25)

    def test_restore_save_settings_match(self):
        for field, value in (("key", "different"), ("path", "~/different"), ("transport", "oras")):
            with self.subTest(field=field):
                self.setUp()
                self.calls()[1]["with"][field] = value
                with self.assertRaisesRegex(ValueError, "restore/save .* mismatch"):
                    self.check()

    def test_step_env_override_cannot_hide_mismatch(self):
        self.calls()[1]["env"] = {"CCACHE_TRUNK_KEY": "different"}
        with self.assertRaisesRegex(ValueError, "restore/save .* mismatch"):
            self.check()

    def test_effective_transport_and_action_defaults(self):
        restore, save = self.calls()
        del restore["with"]["path"]
        save["with"]["transport"] = "docker"
        self.check()

    def test_invalid_runtime_tag_and_transport(self):
        for field, value, error in (("key", "-invalid", "valid GHCR tag"),
                                    ("key", "x" * 129, "valid GHCR tag"),
                                    ("transport", "unknown", "transport must be")):
            with self.subTest(field=field, value=value):
                self.setUp()
                for step in self.calls():
                    step["with"][field] = value
                with self.assertRaisesRegex(ValueError, error):
                    self.check()
        self.setUp()
        for step in self.calls("windows.yml"):
            step["with"]["transport"] = "docker"
        with self.assertRaisesRegex(ValueError, "Docker cache transport requires a Linux runner"):
            self.check()

    def test_save_guard_cannot_be_weakened(self):
        guards = ("", "always()", "github.event_name != 'pull_request'",
                  "github.ref_name == github.event.repository.default_branch",
                  contracts.SAVE_GUARD + " || failure()",
                  contracts.SAVE_GUARD.replace("&&", "||"),
                  contracts.SAVE_GUARD.replace("pull_request", "pull_ request"))
        for guard in guards:
            with self.subTest(guard=guard):
                self.calls()[1]["if"] = guard
                with self.assertRaisesRegex(ValueError, "trusted publication guard"):
                    self.check()

    def test_guard_formatting_is_not_significant(self):
        self.calls()[1]["if"] = "${{\n" + contracts.SAVE_GUARD.replace(" && ", "&&\n") + " }}"
        self.check()

    def test_untrusted_trigger_fails(self):
        for event in ("pull_request_target", "workflow_run", "workflow_call"):
            with self.subTest(event=event):
                self.setUp()
                self.workflows["ci.yml"]["on"][event] = {}
                with self.assertRaisesRegex(ValueError, "unsupported cache workflow trigger"):
                    self.check()

    def test_publication_requires_permissions_and_github_token(self):
        self.job()["permissions"]["packages"] = "read"
        with self.assertRaisesRegex(ValueError, "packages: write"):
            self.check()
        self.setUp()
        for token in ("", "secrets.GITHUB_TOKEN", "${{ secrets.OTHER_TOKEN }}"):
            with self.subTest(token=token):
                self.calls()[1]["with"]["token"] = token
                with self.assertRaisesRegex(ValueError, "save must use secrets.GITHUB_TOKEN"):
                    self.check()

    def test_pair_required_even_when_both_calls_removed(self):
        for mode in ("restore", "save", "both"):
            with self.subTest(mode=mode):
                self.setUp()
                steps = self.job()["steps"]
                steps[:] = [s for s in steps if not (contracts.is_cache(s)
                            and (mode == "both" or s["with"]["mode"] == mode))]
                with self.assertRaisesRegex(ValueError, "exactly one restore followed by one save"):
                    self.check()

    def test_wrong_pair_order_and_duplicate_call(self):
        self.job()["steps"].reverse()
        with self.assertRaisesRegex(ValueError, "one restore followed by one save"):
            self.check()
        self.setUp()
        self.job()["steps"].append(copy.deepcopy(self.calls()[1]))
        with self.assertRaisesRegex(ValueError, "exactly one restore"):
            self.check()

    def test_known_workflow_cannot_disappear(self):
        del self.workflows["windows.yml"]
        with self.assertRaisesRegex(ValueError, "missing required cache job: windows.yml"):
            self.check()

    def test_legacy_fallback_contract(self):
        for field, value, error in (("if", "steps.wrong.outputs.hit != 'true'", "hit output"),
                                    ("uses", "actions/cache@v6", "read-only"),
                                    ("uses", "actions/cache/save@v6", "read-only")):
            with self.subTest(field=field, value=value):
                self.setUp()
                self.fallback()[field] = value
                with self.assertRaisesRegex(ValueError, error):
                    self.check()
        self.setUp()
        self.fallback()["with"]["path"] = "~/other"
        with self.assertRaisesRegex(ValueError, "legacy fallback path"):
            self.check()

    def test_fallback_tracks_actual_restore_id(self):
        self.calls()[0]["id"] = "renamed-cache"
        self.fallback()["if"] = "steps.renamed-cache.outputs.hit != 'true'"
        self.check()

    def test_native_architecture_must_be_in_key(self):
        setup = next(s for s in self.job("macos.yml")["steps"] if "TRUNK_KEY" in s.get("env", {}))
        setup["env"]["TRUNK_KEY"] = setup["env"]["TRUNK_KEY"].replace("${{ runner.arch }}", "x64")
        with self.assertRaisesRegex(ValueError, "macos.yml:.*cache tag collision"):
            self.check()

    def test_native_setup_is_explicit_not_shell_interpretation(self):
        setup = next(s for s in self.job("macos.yml")["steps"] if "TRUNK_KEY" in s.get("env", {}))
        setup["run"] += "\necho CCACHE_TRUNK_KEY=other >> $GITHUB_ENV"
        with self.assertRaisesRegex(ValueError, "unsupported CCACHE_TRUNK_KEY setup"):
            self.check()

    def test_unknown_expression_and_runner_fail_loudly(self):
        self.job()["env"]["CCACHE_TRUNK_KEY"] = "${{ format('cache-{0}', matrix.particles) }}"
        with self.assertRaisesRegex(ValueError, "unsupported cache expression"):
            self.check()
        self.setUp()
        self.job()["runs-on"] = "future-runner"
        with self.assertRaisesRegex(ValueError, "unsupported runner: future-runner"):
            self.check()

    def test_cyclic_env_fails(self):
        self.job()["env"]["CCACHE_TRUNK_KEY"] = "${{ env.CCACHE_TRUNK_KEY }}"
        with self.assertRaisesRegex(ValueError, "cyclic env reference"):
            self.check()


if __name__ == "__main__":
    unittest.main()
