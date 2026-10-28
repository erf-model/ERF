#!/usr/bin/env python3
"""Exercise cache failures without Docker, credentials, or network access."""

import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


ACTION = Path(__file__).resolve().parent
MOCK = r'''#!/usr/bin/env python3
import json, os, pathlib, subprocess, sys
tool = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
calls = pathlib.Path(os.environ['MOCK_CALLS'])
previous = [json.loads(line) for line in calls.read_text().splitlines()] if calls.exists() else []
with calls.open('a') as out:
    out.write(json.dumps([tool, args, os.environ.get('DOCKER_CONFIG')]) + '\n')
scenario = os.environ.get('MOCK_SCENARIO', '')
if tool == 'sleep':
    sys.exit(0)
if tool == 'timeout':
    assert args[0] == '--kill-after=5s' and args[1].endswith('s'), args
    command = args[2:]
    if scenario == 'pull_timeout' and command[:2] == ['docker', 'pull']:
        sys.exit(124)
    if scenario == 'prune_timeout' and command[0] == 'bash':
        sys.exit(124)
    sys.exit(subprocess.call(command))
if tool == 'mv':
    source, target = args[-2:]
    if scenario == 'promotion_failure' and '.ghcr-' in source and not source.endswith('.previous'):
        print('simulated rename failure', file=sys.stderr)
        sys.exit(1)
    sys.exit(subprocess.call(['/bin/mv', *args]))
if tool == 'docker':
    operation = args[0]
    count = sum(t == 'docker' and a[0] == operation for t, a, _ in previous)
    if operation == 'login':
        assert sys.stdin.read() == 'synthetic-test-token'
        pathlib.Path(os.environ['DOCKER_CONFIG'], 'config.json').write_text('synthetic credentials')
    fail = scenario == operation + '_failure'
    fail |= scenario == operation + '_transient' and count == 0
    if operation == 'pull' and scenario == 'missing':
        print('manifest unknown', file=sys.stderr)
        sys.exit(1)
    if operation == 'cp':
        destination = pathlib.Path(args[-1])
        (destination / 'restored-object').write_text('complete object')
        if scenario == 'partial_copy':
            (destination / 'partial-object').write_text('incomplete')
            fail = True
    if fail:
        print('503 registry unavailable / simulated transfer failure', file=sys.stderr)
        sys.exit(1)
    if operation == 'create':
        print('synthetic-container')
    sys.exit(0)
if tool == 'gh':
    if '--method' not in args:
        if scenario == 'prune_failure':
            print('gh: Forbidden (HTTP 403)', file=sys.stderr)
            sys.exit(1)
        print('1\t2026-10-27T01:00:00Z\tubuntu-22.04-x64-gcc-mpich-double-particles-on')
        print('2\t2026-10-27T00:00:00Z\t')
        print('3\t2026-10-27T02:00:00Z\tubuntu-22.04-x64-gcc-mpich-double-particles-off')
        print('4\t2026-10-26T00:00:00Z\tsha-old')
        if scenario == 'retention':
            print('5\t2026-10-27T03:00:00Z\tsha-recent')
            print('6\t2026-10-27T04:00:00Z\tsha-newest')
            print('7\t2026-10-25T00:00:00Z\tsha-older,release')
    elif scenario == 'delete_race' and args[-1].endswith('/2'):
        print('gh: Not Found (HTTP 404)', file=sys.stderr)
        sys.exit(1)
    sys.exit(0)
raise AssertionError((tool, args))
'''


class CacheTransferTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="erf-cache-test-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bin = self.root / "bin"
        self.bin.mkdir()
        for name in ("docker", "gh", "timeout", "sleep", "mv"):
            script = self.bin / name
            script.write_text(MOCK)
            script.chmod(0o755)
        self.cache = self.root / "cache with spaces"
        self.cache.mkdir()
        (self.cache / "local-object").write_text("original")
        self.runner = self.root / "runner"
        self.runner.mkdir()
        self.shared_config = self.root / "shared-docker-config"
        self.shared_config.mkdir()
        (self.shared_config / "config.json").write_text("existing credentials")
        self.output = self.root / "output"
        self.call_log = self.root / "calls.jsonl"
        self.env = dict(os.environ, PATH=f"{self.bin}{os.pathsep}{os.environ['PATH']}",
                        CACHE_MODE="restore", CACHE_KEY="Ubuntu-22.04 X64 GCC MPICH DOUBLE Particles ON",
                        CACHE_PATH=str(self.cache), GHCR_TOKEN="synthetic-test-token",
                        RUNNER_OS="Linux", RUNNER_TEMP=str(self.runner),
                        GITHUB_OUTPUT=str(self.output), GITHUB_REPOSITORY_OWNER="ExampleOwner",
                        GITHUB_ACTOR="test-actor", GITHUB_SHA="0" * 40,
                        CACHE_SOURCE="https://example.invalid/repo", CACHE_OWNER_TYPE="User",
                        DOCKER_CONFIG=str(self.shared_config), MOCK_CALLS=str(self.call_log))

    def run_action(self, mode="restore", scenario="", **overrides):
        env = dict(self.env, CACHE_MODE=mode, MOCK_SCENARIO=scenario, **overrides)
        result = subprocess.run(["bash", str(ACTION / "cache.sh")], env=env,
                                capture_output=True, text=True, timeout=20)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(list(self.runner.iterdir()), [], "temporary credentials leaked")
        self.assertEqual(list(self.root.glob("*.ghcr-*")), [], "staging directory leaked")
        self.assertEqual((self.shared_config / "config.json").read_text(), "existing credentials")
        return result

    def calls(self, tool):
        if not self.call_log.exists():
            return []
        return [args for name, args, _ in map(json.loads, self.call_log.read_text().splitlines())
                if name == tool]

    def assert_local_preserved(self):
        self.assertEqual([p.name for p in self.cache.iterdir()], ["local-object"])
        self.assertEqual((self.cache / "local-object").read_text(), "original")
        self.assertEqual(self.output.read_text(), "hit=false\n")

    def test_complete_restore_replaces_cache(self):
        self.run_action()
        self.assertEqual([p.name for p in self.cache.iterdir()], ["restored-object"])
        self.assertTrue(self.output.read_text().endswith("hit=true\n"))
        self.assertTrue(any(a[0] == "rm" for a in self.calls("docker")))

    def test_failed_restores_preserve_cache_and_allow_build(self):
        for scenario in ("login_failure", "pull_failure", "pull_timeout", "create_failure",
                         "partial_copy", "promotion_failure", "missing"):
            with self.subTest(scenario=scenario):
                self.output.unlink(missing_ok=True)
                self.call_log.unlink(missing_ok=True)
                result = self.run_action(scenario=scenario)
                self.assert_local_preserved()
                if scenario == "missing":
                    self.assertIn("cold build", result.stdout)
                else:
                    self.assertIn("::warning::", result.stderr)
                if scenario in ("login_failure", "pull_failure"):
                    operation = scenario.split("_")[0]
                    self.assertEqual(sum(a[0] == operation for a in self.calls("docker")), 2)

    def test_restore_without_credentials_or_existing_directory(self):
        (self.cache / "local-object").unlink()
        self.cache.rmdir()
        self.run_action(GHCR_TOKEN="")
        self.assertEqual((self.cache / "restored-object").read_text(), "complete object")
        self.assertFalse(any(a[0] == "login" for a in self.calls("docker")))

    def test_transient_pull_recovers_with_two_attempts(self):
        self.run_action(scenario="pull_transient")
        self.assertEqual(sum(a[0] == "pull" for a in self.calls("docker")), 2)
        self.assertTrue(self.output.read_text().endswith("hit=true\n"))

    def test_save_failure_never_prunes_and_allows_tests(self):
        for scenario in ("login_failure", "build_failure", "push_failure"):
            with self.subTest(scenario=scenario):
                self.call_log.unlink(missing_ok=True)
                result = self.run_action("save", scenario)
                self.assertIn("::warning::", result.stderr)
                self.assertEqual(self.calls("gh"), [])

    def test_save_retries_and_prunes_only_superseded_versions(self):
        self.run_action("save", "push_transient", CACHE_OWNER_TYPE="Organization")
        self.assertEqual(sum(a[0] == "push" for a in self.calls("docker")), 2)
        deleted = [a[-1] for a in self.calls("gh") if "DELETE" in a]
        self.assertEqual(deleted, [f"orgs/ExampleOwner/packages/container/erf-ccache/versions/{i}"
                                   for i in (2, 4)])

    def test_prune_failure_timeout_and_delete_race_are_nonfatal(self):
        for scenario in ("prune_failure", "prune_timeout", "delete_race"):
            with self.subTest(scenario=scenario):
                self.call_log.unlink(missing_ok=True)
                result = self.run_action("save", scenario)
                self.assertIn(
                    "saved ghcr.io/exampleowner/erf-ccache:ubuntu-22.04-x64-gcc-mpich-double-particles-on",
                    result.stdout)
                if scenario == "delete_race":
                    deleted = [a[-1] for a in self.calls("gh") if "DELETE" in a]
                    self.assertEqual(deleted, [f"users/ExampleOwner/packages/container/erf-ccache/versions/{i}"
                                               for i in (2, 4)])
                else:
                    self.assertIn("::warning::", result.stderr)

    def test_prune_retains_newest_sha_versions_and_all_named_tags(self):
        cases = (
            ("erf-ccache:0", {2, 4, 5, 6}),
            ("erf-ccache:1", {2, 4, 5}),
            ("erf-ccache", {2, 4}),
        )
        for package_spec, expected_deleted in cases:
            with self.subTest(package_spec=package_spec):
                self.call_log.unlink(missing_ok=True)
                env = dict(self.env, MOCK_SCENARIO="retention", KEEP="2",
                           OWNER="ExampleOwner", OWNER_TYPE="User",
                           GH_TOKEN="synthetic-test-token")
                result = subprocess.run(["bash", str(ACTION / "prune.sh"), package_spec],
                                        env=env, capture_output=True, text=True, timeout=20)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                deleted = {int(args[-1].rsplit("/", 1)[1]) for args in self.calls("gh")
                           if "DELETE" in args}
                self.assertEqual(deleted, expected_deleted)

    def test_empty_save_does_not_login_or_prune(self):
        (self.cache / "local-object").unlink()
        self.run_action("save")
        self.assertFalse(self.call_log.exists())

    def test_invalid_configuration_fails_before_transport(self):
        for override in ({"CACHE_MODE": "invalid"}, {"CACHE_KEY": ""},
                         {"CACHE_MODE": "save", "GHCR_TOKEN": ""},
                         {"CACHE_PATH": "/"}, {"CACHE_PATH": os.environ["HOME"]}):
            with self.subTest(override=override):
                result = subprocess.run(["bash", str(ACTION / "cache.sh")],
                                        env=dict(self.env, **override), capture_output=True, text=True)
                self.assertEqual(result.returncode, 2)
                self.assertIn("::error::", result.stderr)
                self.assertFalse(self.call_log.exists())


if __name__ == "__main__":
    unittest.main()
