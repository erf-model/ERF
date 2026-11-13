"""Test cleanup safeguards with mocked APIs; never contact GitHub."""

import contextlib
import copy
import importlib.util
import io
import os
import subprocess
import unittest
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "cleanup", ROOT / ".github/scripts/cleanup_closed_pr_ci.py"
)
cleanup = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cleanup)


def run_record(**changes):
    return {
        "id": 1,
        "event": "pull_request",
        "status": "queued",
        "run_attempt": 1,
        "created_at": "2026-11-11T12:00:00Z",
        "head_sha": "a" * 40,
        "head_branch": "feature",
        "path": ".github/workflows/gcc.yml",
        "head_repository": {"id": 42, "owner": {"login": "contributor"}},
        "pull_requests": [{"number": 10}],
        **changes,
    }


def pr_record(**changes):
    return {
        "number": 10,
        "state": "closed",
        "merged": True,
        "base": {"repo": {"full_name": "erf-model/ERF"}},
        "head": {"repo": {"id": 42}, "ref": "feature", "sha": "a" * 40},
        **changes,
    }


class CleanupTests(unittest.TestCase):
    def setUp(self):
        self.listed = [run_record()]
        self.current = {1: run_record()}
        self.prs = {10: pr_record()}
        self.associated, self.open_prs, self.writes = [], [], []
        self.errors = {}
        self.stuck = False
        self.sleep = mock.Mock()
        patcher = mock.patch.object(cleanup, "api", side_effect=self.api)
        self.client = patcher.start()
        self.addCleanup(patcher.stop)

    def api(self, method, path, paginate=False):
        endpoint = urlsplit(path).path
        if (method, endpoint) in self.errors:
            raise cleanup.APIError(self.errors[method, endpoint], "request failed")
        if method == "POST":
            self.writes.append(path)
            if not self.stuck:
                self.current[int(endpoint.split("/")[-2])]["status"] = "completed"
            return None
        if endpoint.endswith("/actions/runs/9000"):
            return {"created_at": "2026-11-12T12:00:00Z"}
        if endpoint.endswith("/actions/runs"):
            status = parse_qs(urlsplit(path).query)["status"][0]
            runs = [run for run in self.listed if run["status"] == status]
            # Model multiple API pages, not just the first page.
            return copy.deepcopy([{"workflow_runs": runs[:1]}, {"workflow_runs": runs[1:]}])
        if "/actions/runs/" in endpoint:
            return copy.deepcopy(self.current[int(endpoint.rsplit("/", 1)[1])])
        if "/commits/" in endpoint:
            return [copy.deepcopy(self.associated)]
        if endpoint.endswith("/pulls"):
            return [copy.deepcopy(self.open_prs)]
        return copy.deepcopy(self.prs[int(endpoint.rsplit("/", 1)[1])])

    def perform(self, **kwargs):
        self.output = io.StringIO()
        with contextlib.redirect_stdout(self.output):
            return cleanup.cleanup("erf-model/ERF", 9000, sleep=self.sleep, **kwargs)

    def test_default_cancels_closed_and_merged_prs_in_all_unfinished_states(self):
        for merged in (False, True):
            for status in cleanup.UNFINISHED_STATUSES:
                with self.subTest(merged=merged, status=status):
                    self.listed = [run_record(status=status)]
                    self.current[1] = run_record(status=status)
                    self.prs[10] = pr_record(merged=merged)
                    self.writes.clear()
                    self.assertEqual(self.perform(), 0)
                    self.assertEqual(self.writes, ["repos/erf-model/ERF/actions/runs/1/cancel"])

    def test_dry_run_lists_proposals_without_writes_or_sleep(self):
        self.assertEqual(self.perform(dry_run=True), 0)
        self.assertEqual(self.writes, [])
        self.sleep.assert_not_called()
        self.assertIn("Would cancel: run 1", self.output.getvalue())

    def test_open_pr_and_reused_branch_are_protected(self):
        for number in (10, 11):
            with self.subTest(number=number):
                self.prs[10] = pr_record(state="open" if number == 10 else "closed")
                self.open_prs = [pr_record(number=number, state="open")]
                self.perform()
                self.assertEqual(self.writes, [])

    def test_open_pr_skip_identifies_pr_and_current_or_different_commit(self):
        for sha, message in (
            ("a" * 40, "matches its current commit"),
            ("b" * 40, "is from a different commit"),
        ):
            with self.subTest(sha=sha):
                self.prs[10] = pr_record(state="open")
                self.prs[10]["head"]["sha"] = sha
                self.perform(dry_run=True)
                self.assertIn(
                    "Skipping run 1: PR #10 is open; run " + message, self.output.getvalue()
                )
                self.assertEqual(self.writes, [])

    def test_missing_commit_association_can_identify_open_pr_but_not_authorize_cancel(self):
        self.listed = [run_record(pull_requests=[])]
        self.current[1] = run_record(pull_requests=[])
        self.open_prs = [pr_record(number=11, state="open")]
        self.perform()
        self.assertIn("PR #11 is open; run matches its current commit", self.output.getvalue())
        self.assertEqual(self.writes, [])
        self.open_prs = []
        self.perform()
        self.assertIn("no closed PR association", self.output.getvalue())
        self.assertEqual(self.writes, [])

    def test_other_fork_same_branch_is_not_mistaken_for_source_pr(self):
        self.open_prs = [pr_record(state="open", head={"repo": {"id": 99}, "ref": "feature"})]
        self.perform()
        self.assertEqual(len(self.writes), 1)
        self.writes.clear()
        self.current[1] = run_record()
        self.prs[10] = pr_record(head={"repo": {"id": 99}, "ref": "feature"})
        self.perform()
        self.assertEqual(self.writes, [])

    def test_missing_fork_association_resolves_by_commit_or_skips(self):
        self.listed = [run_record(pull_requests=[])]
        self.current[1] = run_record(pull_requests=[])
        self.perform()
        self.assertEqual(self.writes, [])
        self.associated = [{"number": 10}]
        self.perform()
        self.assertEqual(len(self.writes), 1)
        self.assertTrue(
            any("/commits/" + "a" * 40 in call.args[1] for call in self.client.call_args_list)
        )

    def test_protected_events_workflows_and_uncertain_metadata_are_skipped(self):
        for record in (
            run_record(event="push"),
            run_record(event="workflow_dispatch"),
            run_record(status="completed"),
            run_record(head_repository=None),
            run_record(created_at="invalid"),
            run_record(run_attempt=None),
            run_record(created_at="2026-11-12T12:00:01Z"),
            run_record(path=".github/workflows/spack-containers.yml"),
        ):
            with self.subTest(run=record):
                self.listed = [record]
                self.perform()
                self.assertEqual(self.writes, [])

    def test_force_cancel_rechecks_closed_state_and_run_attempt(self):
        self.stuck = True
        self.perform()
        self.sleep.assert_called_once_with(60)
        self.assertTrue(self.writes[-1].endswith("/force-cancel"))
        for change in (
            lambda: self.prs[10].update(state="open"),
            lambda: self.current[1].update(run_attempt=2),
        ):
            self.current[1] = run_record()
            self.prs[10] = pr_record()
            self.writes.clear()
            self.sleep.side_effect = lambda seconds, change=change: change()
            self.perform()
            self.assertEqual(len(self.writes), 1)

    def test_multiple_pages_and_duplicate_runs(self):
        self.listed = [run_record(), run_record(id=2), run_record()]
        self.current[2] = run_record(id=2)
        self.perform()
        self.assertEqual(len(self.writes), 2)
        listings = [call for call in self.client.call_args_list if "?event=" in call.args[1]]
        self.assertEqual(len(listings), 5)
        self.assertTrue(all(call.kwargs["paginate"] for call in listings))

    def test_finished_runs_and_new_attempts_are_not_cancelled(self):
        for changes in ({"status": "completed"}, {"run_attempt": 2}, {"id": 2}):
            self.current[1] = run_record(**changes)
            self.perform()
            self.assertEqual(self.writes, [])

    def test_api_errors_fail_visibly_except_expected_disappearance_or_cancel_races(self):
        for method, endpoint in (
            ("GET", "repos/erf-model/ERF/pulls/10"),
            ("POST", "repos/erf-model/ERF/actions/runs/1/cancel"),
        ):
            for status in (404, 409, 403, 500):
                with self.subTest(method=method, status=status):
                    self.errors = {(method, endpoint): status}
                    benign = status == 404 or (method == "POST" and status == 409)
                    self.assertEqual(self.perform(), 0 if benign else 1)
                    self.assertEqual(self.writes, [])

    def test_manual_invocation_uses_start_time_without_an_actions_run(self):
        with (
            mock.patch.object(
                cleanup.time, "time", return_value=cleanup.timestamp("2026-11-12T12:00:00Z")
            ),
            contextlib.redirect_stdout(io.StringIO()),
        ):
            cleanup.cleanup("erf-model/ERF", dry_run=True, sleep=self.sleep)
        self.assertFalse(
            any(call.args[1].endswith("/runs/None") for call in self.client.call_args_list)
        )
        self.assertEqual(self.writes, [])


class EntrypointTests(unittest.TestCase):
    def test_cli_defaults_to_cancellation_and_dry_run_is_explicit(self):
        for flags, expected in (([], False), (["--dry-run"], True)):
            with (
                mock.patch.dict(os.environ, {}, clear=True),
                mock.patch("sys.argv", ["cleanup", *flags]),
                mock.patch.object(cleanup, "cleanup", return_value=0) as sweep,
            ):
                self.assertEqual(cleanup.main(), 0)
                sweep.assert_called_once_with("erf-model/ERF", None, dry_run=expected)

    def test_gh_transport_uses_argument_arrays_pagination_and_http_status(self):
        with mock.patch.object(cleanup.subprocess, "run") as command:
            command.return_value = subprocess.CompletedProcess([], 0, "[[]]", "")
            self.assertEqual(cleanup.api("GET", "repos/erf-model/ERF/pulls", paginate=True), [[]])
            self.assertIn("--paginate", command.call_args.args[0])
            self.assertNotIn("shell", command.call_args.kwargs)
            command.return_value = subprocess.CompletedProcess([], 1, "", "gh: details (HTTP 409)")
            with self.assertRaises(cleanup.APIError) as caught:
                cleanup.api("POST", "repos/erf-model/ERF/actions/runs/1/cancel")
            self.assertEqual(caught.exception.status, 409)

    def test_workflow_is_independent_scoped_and_has_no_path_filters(self):
        workflow = yaml.load(
            (ROOT / ".github/workflows/cleanup-merged-pr-ci.yml").read_text(),
            Loader=yaml.BaseLoader,
        )
        self.assertEqual(workflow["on"], {"push": {"branches": ["development"]}})
        self.assertEqual(workflow["concurrency"]["cancel-in-progress"], "false")
        job = workflow["jobs"]["cleanup"]
        self.assertNotIn("needs", job)
        self.assertEqual(job["if"], "github.repository == 'erf-model/ERF'")
        self.assertEqual(
            job["permissions"], {"contents": "read", "actions": "write", "pull-requests": "read"}
        )
        self.assertEqual(job["steps"][-1]["run"], "python3 .github/scripts/cleanup_closed_pr_ci.py")
        self.assertEqual(job["steps"][-1]["env"]["GH_TOKEN"], "${{ secrets.GITHUB_TOKEN }}")


if __name__ == "__main__":
    unittest.main()
