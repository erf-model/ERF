"""Test cancellation policy and gh transport without contacting GitHub."""

import contextlib
import copy
import importlib.util
import io
import subprocess
import unittest
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlsplit

import yaml

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "closed_pr_cleanup", ROOT / ".github/scripts/cleanup_closed_pr_ci.py"
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
        "head": {"repo": {"id": 42}, "ref": "feature"},
        **changes,
    }


class CleanupPolicyTests(unittest.TestCase):
    def setUp(self):
        self.api = mock.create_autospec(cleanup.GitHubAPI, instance=True)
        self.listed = [run_record()]
        self.current = {1: run_record()}
        self.prs = {10: pr_record()}
        self.cancelled = []
        self.forced = []
        self.logs = []
        self.sleeps = []
        self.stuck = False
        self.on_sleep = lambda: None
        self.api.list_pr_runs.side_effect = lambda status: copy.deepcopy(
            [run for run in self.listed if run["status"] == status]
        )
        self.api.get_run.side_effect = self.get_run
        self.api.get_pr.side_effect = lambda number: copy.deepcopy(self.prs[number])
        self.api.prs_for_commit.return_value = []
        self.api.open_prs_for_branch.return_value = []
        self.api.cancel_run.side_effect = self.cancel_run

    def get_run(self, run_id):
        if run_id == 9000:
            return {"created_at": "2026-11-12T12:00:00Z"}
        return copy.deepcopy(self.current[run_id])

    def cancel_run(self, run_id, force=False):
        if force:
            self.forced.append(run_id)
        else:
            self.cancelled.append(run_id)
            if not self.stuck:
                self.current[run_id]["status"] = "completed"

    def sleep(self, seconds):
        self.sleeps.append(seconds)
        self.on_sleep()

    def perform(self, **kwargs):
        return cleanup.cleanup(
            self.api, "erf-model/ERF", 9000, sleep=self.sleep, log=self.logs.append, **kwargs
        )

    def test_closed_and_merged_prs_across_unfinished_states(self):
        for merged in (False, True):
            for status in cleanup.UNFINISHED_STATUSES:
                with self.subTest(merged=merged, status=status):
                    self.setUp()
                    self.listed = [run_record(status=status)]
                    self.current[1] = run_record(status=status)
                    self.prs[10] = pr_record(merged=merged)
                    self.assertEqual(self.perform(apply=True), 0)
                    self.assertEqual(self.cancelled, [1])
                    self.assertEqual(self.forced, [])

    def test_default_is_read_only_without_sleeping(self):
        self.assertEqual(self.perform(), 0)
        self.api.cancel_run.assert_not_called()
        self.assertEqual(self.sleeps, [])
        self.assertTrue(any("Would cancel" in line for line in self.logs))

    def test_open_pr_is_not_cancelled(self):
        self.prs[10] = pr_record(state="open", merged=False)
        self.assertEqual(self.perform(apply=True), 0)
        self.api.cancel_run.assert_not_called()
        self.assertTrue(any("PR #10 is not closed" in line for line in self.logs))

    def test_reopened_pr_is_not_force_cancelled(self):
        self.stuck = True
        self.on_sleep = lambda: self.prs[10].update(state="open")
        self.perform(apply=True)
        self.assertEqual(self.cancelled, [1])
        self.assertEqual(self.forced, [])

    def test_branch_reuse_blocks_cleanup_but_other_fork_does_not(self):
        for source_id, expected in ((42, []), (99, [1])):
            with self.subTest(source_id=source_id):
                self.setUp()
                self.api.open_prs_for_branch.return_value = [
                    pr_record(
                        number=11, state="open", head={"repo": {"id": source_id}, "ref": "feature"}
                    )
                ]
                self.perform(apply=True)
                self.assertEqual(self.cancelled, expected)

    def test_wrong_source_or_base_pr_cannot_authorize_cancellation(self):
        for pr in (
            pr_record(head={"repo": {"id": 99}, "ref": "feature"}),
            pr_record(head={"repo": {"id": 42}, "ref": "another"}),
            pr_record(head={"repo": None, "ref": "feature"}),
            pr_record(base={"repo": {"full_name": "other/ERF"}}),
        ):
            with self.subTest(pr=pr):
                self.setUp()
                self.prs[10] = pr
                self.perform(apply=True)
                self.api.cancel_run.assert_not_called()

    def test_missing_fork_association_is_resolved_by_commit(self):
        self.listed = [run_record(pull_requests=[])]
        self.current[1] = run_record(pull_requests=[])
        self.api.prs_for_commit.return_value = [{"number": 10}]
        self.perform(apply=True)
        self.api.prs_for_commit.assert_called_once_with("a" * 40)
        self.assertEqual(self.cancelled, [1])

    def test_missing_or_malformed_association_is_skipped(self):
        for associations in ([], [{}]):
            with self.subTest(associations=associations):
                self.setUp()
                self.listed = [run_record(pull_requests=associations)]
                self.current[1] = run_record(pull_requests=associations)
                self.perform(apply=True)
                self.api.cancel_run.assert_not_called()

    def test_open_association_blocks_a_closed_one(self):
        record = run_record(pull_requests=[{"number": 10}, {"number": 11}])
        self.listed = [record]
        self.current[1] = record
        self.prs[11] = pr_record(number=11, state="open")
        self.perform(apply=True)
        self.api.cancel_run.assert_not_called()

    def test_protected_or_uncertain_runs_are_untouched(self):
        cases = [
            run_record(event="push"),
            run_record(event="workflow_dispatch"),
            run_record(event="pull_request_target"),
            run_record(status="completed"),
            run_record(created_at="2026-11-12T12:00:01Z"),
            run_record(created_at="invalid"),
            run_record(created_at="2026-11-11T12:00:00"),
            run_record(head_repository=None),
            run_record(head_sha=None),
            run_record(head_branch=None),
            run_record(run_attempt=None),
            run_record(id=None),
            run_record(path=".github/workflows/spack-containers.yml"),
            run_record(path=".github/workflows/cleanup-merged-pr-ci.yml"),
        ]
        for record in cases:
            with self.subTest(run=record):
                self.setUp()
                self.listed = [record]
                self.perform(apply=True)
                self.api.cancel_run.assert_not_called()

    def test_changed_run_state_identity_or_attempt_is_skipped(self):
        for changes in ({"status": "completed"}, {"id": 2}, {"run_attempt": 2}):
            with self.subTest(changes=changes):
                self.setUp()
                self.current[1].update(changes)
                self.perform(apply=True)
                self.api.cancel_run.assert_not_called()

    def test_force_cancel_waits_and_rechecks_attempt(self):
        self.stuck = True
        self.perform(apply=True)
        self.assertEqual(self.sleeps, [60])
        self.assertEqual(self.forced, [1])
        self.setUp()
        self.stuck = True
        self.on_sleep = lambda: self.current[1].update(run_attempt=2)
        self.perform(apply=True)
        self.assertEqual(self.forced, [])

    def test_duplicate_listings_are_deduplicated(self):
        self.api.list_pr_runs.side_effect = lambda status: [run_record()]
        self.perform(apply=True)
        self.assertEqual(self.cancelled, [1])
        self.assertEqual(
            self.api.list_pr_runs.call_args_list,
            [mock.call(status) for status in cleanup.UNFINISHED_STATUSES],
        )

    def test_cancel_completion_races_are_harmless_but_real_errors_fail(self):
        for status in (404, 409, 403, 500):
            for force in (False, True):
                with self.subTest(status=status, force=force):
                    self.setUp()
                    self.stuck = force

                    def post(run_id, force=False, error_status=status):
                        if force == self.stuck:
                            raise cleanup.GitHubAPIError(error_status, "request failed")
                        self.cancel_run(run_id, force)

                    self.api.cancel_run.side_effect = post
                    self.assertEqual(self.perform(apply=True), 0 if status in (404, 409) else 1)

    def test_read_errors_never_authorize_cancellation(self):
        for status in (404, 409, 403, 500):
            with self.subTest(status=status):
                self.setUp()
                self.api.get_pr.side_effect = cleanup.GitHubAPIError(status, "lookup failed")
                self.assertEqual(self.perform(apply=True), 0 if status == 404 else 1)
                self.api.cancel_run.assert_not_called()

    def test_one_run_error_does_not_prevent_other_cleanup(self):
        self.listed.append(run_record(id=2))
        self.current[2] = run_record(id=2)

        def post(run_id, force=False):
            if run_id == 1:
                raise cleanup.GitHubAPIError(403, "not permitted")
            self.cancel_run(run_id, force)

        self.api.cancel_run.side_effect = post
        self.assertEqual(self.perform(apply=True), 1)
        self.assertEqual(self.cancelled, [2])

    def test_branch_metacharacters_remain_api_data(self):
        branch = "feature;$(do-not-execute)"
        self.listed = [run_record(head_branch=branch)]
        self.current[1] = run_record(head_branch=branch)
        self.prs[10] = pr_record(head={"repo": {"id": 42}, "ref": branch})
        self.perform(apply=True)
        self.api.open_prs_for_branch.assert_called_once_with("contributor", branch)
        self.assertEqual(self.cancelled, [1])

    def test_invalid_sweep_time_stops_before_listing_or_cancellation(self):
        self.api.get_run.return_value = {"created_at": None}
        self.api.get_run.side_effect = None
        with self.assertRaisesRegex(ValueError, "creation time"):
            self.perform(apply=True)
        self.api.list_pr_runs.assert_not_called()
        self.api.cancel_run.assert_not_called()


class GitHubTransportTests(unittest.TestCase):
    def setUp(self):
        self.api = cleanup.GitHubAPI("erf-model/ERF")
        patcher = mock.patch.object(cleanup.subprocess, "run")
        self.command = patcher.start()
        self.addCleanup(patcher.stop)

    def response(self, output, error="", code=0):
        self.command.return_value = subprocess.CompletedProcess([], code, output, error)

    def test_paginated_run_objects_and_pr_arrays_are_flattened(self):
        self.response('[{"workflow_runs": [{"id": 1}]}, {"workflow_runs": [{"id": 2}]}]')
        self.assertEqual(self.api.list_pr_runs("queued"), [{"id": 1}, {"id": 2}])
        command = self.command.call_args.args[0]
        self.assertEqual(
            command[:7], ["gh", "api", "--method", "GET", "--paginate", "--slurp", command[-1]]
        )
        query = parse_qs(urlsplit(command[-1]).query)
        self.assertEqual(
            query, {"event": ["pull_request"], "status": ["queued"], "per_page": ["100"]}
        )
        self.response('[[{"number": 10}], [{"number": 11}]]')
        self.assertEqual(self.api.prs_for_commit("a" * 40), [{"number": 10}, {"number": 11}])

    def test_branch_names_are_encoded_and_no_shell_is_used(self):
        self.response("[[]]")
        branch = "feature&state=closed;$(do-not-execute)"
        self.api.open_prs_for_branch("contributor", branch)
        command = self.command.call_args.args[0]
        self.assertEqual(parse_qs(urlsplit(command[-1]).query)["head"], ["contributor:" + branch])
        self.assertNotIn("shell", self.command.call_args.kwargs)
        self.assertEqual(self.command.call_args.kwargs["timeout"], 60)

    def test_only_explicit_cancel_method_uses_post(self):
        self.response("")
        self.api.cancel_run(1)
        self.assertEqual(
            self.command.call_args.args[0],
            ["gh", "api", "--method", "POST", "repos/erf-model/ERF/actions/runs/1/cancel"],
        )
        self.api.cancel_run(1, force=True)
        self.assertTrue(self.command.call_args.args[0][-1].endswith("/force-cancel"))

    def test_http_status_is_preserved_without_echoing_raw_errors(self):
        self.response("", "gh: untrusted details (HTTP 409)", code=1)
        with self.assertRaises(cleanup.GitHubAPIError) as caught:
            self.api.cancel_run(1)
        self.assertEqual(caught.exception.status, 409)
        self.assertNotIn("untrusted details", str(caught.exception))

    def test_transport_timeout_is_an_error(self):
        self.command.side_effect = subprocess.TimeoutExpired("gh", 60)
        with self.assertRaises(cleanup.GitHubAPIError):
            self.api.get_run(1)

    def test_malformed_responses_fail_closed(self):
        for output in ("", "not JSON", "null", "[]"):
            with self.subTest(output=output):
                self.response(output)
                with self.assertRaises(cleanup.GitHubAPIError):
                    self.api.get_run(1)
        for output in ("{}", '[{"workflow_runs": null}]', '[{"workflow_runs": [1]}]'):
            with self.subTest(output=output):
                self.response(output)
                with self.assertRaises(cleanup.GitHubAPIError):
                    self.api.list_pr_runs("queued")


class CleanupEntrypointTests(unittest.TestCase):
    def test_cli_is_read_only_unless_apply_is_explicit(self):
        for extra, expected in (([], False), (["--apply"], True)):
            with (
                self.subTest(extra=extra),
                mock.patch(
                    "sys.argv", ["cleanup", "--repo", "erf-model/ERF", "--run-id", "9000", *extra]
                ),
                mock.patch.object(cleanup, "cleanup", return_value=0) as sweep,
            ):
                self.assertEqual(cleanup.main(), 0)
                self.assertEqual(sweep.call_args.kwargs["apply"], expected)

    def test_workflow_uses_scoped_token_and_python_on_every_development_push(self):
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
        checkout, execute = job["steps"]
        self.assertNotIn("ref", checkout.get("with", {}))
        self.assertEqual(checkout["with"]["persist-credentials"], "false")
        self.assertEqual(execute["env"]["GH_TOKEN"], "${{ secrets.GITHUB_TOKEN }}")
        self.assertEqual(execute["run"], "python3 .github/scripts/cleanup_closed_pr_ci.py --apply")

    def test_cli_reports_failure_and_rejects_invalid_scope(self):
        with (
            mock.patch("sys.argv", ["cleanup", "--repo", "erf-model/ERF", "--run-id", "9000"]),
            mock.patch.object(cleanup, "cleanup", return_value=256),
        ):
            self.assertEqual(cleanup.main(), 1)
        with (
            mock.patch("sys.argv", ["cleanup", "--repo", "bad/repo/path", "--run-id", "9000"]),
            contextlib.redirect_stderr(io.StringIO()),
            self.assertRaises(SystemExit) as caught,
        ):
            cleanup.main()
        self.assertEqual(caught.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
