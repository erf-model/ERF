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
        "workflow_id": 7,
        "run_number": 1,
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
        self.history = []
        self.stuck = False
        self.force_stuck = False
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
            stuck = self.force_stuck if endpoint.endswith("/force-cancel") else self.stuck
            if not stuck:
                self.current[int(endpoint.split("/")[-2])]["status"] = "completed"
            return None
        if endpoint.endswith("/actions/runs/9000"):
            return {"created_at": "2026-11-12T12:00:00Z"}
        if endpoint.endswith("/actions/runs"):
            status = parse_qs(urlsplit(path).query)["status"][0]
            runs = [run for run in self.listed if run["status"] == status]
            # Model multiple API pages, not just the first page.
            return copy.deepcopy([{"workflow_runs": runs[:1]}, {"workflow_runs": runs[1:]}])
        if "/actions/workflows/" in endpoint:
            return copy.deepcopy(
                [{"workflow_runs": self.history[:1]}, {"workflow_runs": self.history[1:]}]
            )
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

    def add_replacement(self, **changes):
        self.prs[10] = pr_record(state="open")
        newer = run_record(id=2, run_number=2, head_sha="b" * 40, **changes)
        self.history = [newer]
        self.current[2] = copy.deepcopy(newer)

    def test_superseded_open_pr_run_is_cancelled_with_live_or_completed_replacement(self):
        for status, conclusion in (
            ("queued", None),
            ("in_progress", None),
            ("completed", "success"),
            ("completed", "failure"),
            ("completed", "timed_out"),
        ):
            with self.subTest(status=status, conclusion=conclusion):
                self.current[1] = run_record()
                self.writes.clear()
                self.add_replacement(status=status, conclusion=conclusion)
                self.assertEqual(self.perform(), 0)
                self.assertEqual(self.writes, ["repos/erf-model/ERF/actions/runs/1/cancel"])
                self.assertIn("superseded by run 2", self.output.getvalue())
                self.assertEqual(self.current[2]["status"], status)

    def test_supersession_does_not_require_fork_commit_pr_associations(self):
        self.add_replacement(pull_requests=[])
        self.listed[0]["pull_requests"] = []
        self.current[1]["pull_requests"] = []
        self.assertEqual(self.perform(), 0)
        self.assertEqual(len(self.writes), 1)

    def test_dry_run_of_superseded_open_pr_run_has_no_writes_or_sleep(self):
        self.add_replacement()
        self.assertEqual(self.perform(dry_run=True), 0)
        self.assertIn("superseded by run 2", self.output.getvalue())
        self.assertEqual(self.writes, [])
        self.sleep.assert_not_called()

    def test_cancelled_skipped_or_unusable_replacements_do_not_supersede(self):
        for conclusion in ("cancelled", "skipped", "action_required", "startup_failure", None):
            with self.subTest(conclusion=conclusion):
                self.add_replacement(status="completed", conclusion=conclusion)
                self.perform()
                self.assertEqual(self.writes, [])

    def test_replacements_must_match_source_workflow_and_known_pr(self):
        for changes in (
            {"workflow_id": 8},
            {"head_branch": "other"},
            {"head_repository": {"id": 99}},
            {"pull_requests": [{"number": 11}]},
            {"event": "push"},
            {"path": ".github/workflows/windows.yml"},
        ):
            with self.subTest(changes=changes):
                self.add_replacement(**changes)
                self.perform()
                self.assertEqual(self.writes, [])

    def test_latest_run_or_missing_replacement_metadata_is_kept(self):
        self.add_replacement()
        for changes in ({"run_number": 1}, {"run_number": None}, {"workflow_id": None}):
            with self.subTest(changes=changes):
                self.current[1] = run_record(**changes)
                if changes == {"run_number": 1}:
                    self.history[0]["run_number"] = 1
                self.perform()
                self.assertEqual(self.writes, [])

    def test_missing_or_null_replacement_number_keeps_current_checks(self):
        for missing in (False, True):
            with self.subTest(missing=missing):
                self.add_replacement()
                for record in (self.history[0], self.current[2]):
                    if missing:
                        record.pop("run_number")
                    else:
                        record["run_number"] = None
                self.assertEqual(self.perform(), 0)
                self.assertEqual(self.writes, [])

    def test_replacement_is_refetched_and_cancelled_replacement_cannot_authorize_cancel(self):
        self.add_replacement()
        self.current[2].update(status="completed", conclusion="cancelled")
        self.perform()
        self.assertEqual(self.writes, [])

    def test_replacement_pages_and_query_include_completed_runs(self):
        self.add_replacement(status="completed", conclusion="success")
        self.history.insert(0, run_record(id=3, run_number=3, head_repository={"id": 99}))
        self.perform()
        self.assertIn("superseded by run 2", self.output.getvalue())
        listing = next(
            call
            for call in self.client.call_args_list
            if "/actions/workflows/7/runs?" in call.args[1]
        )
        query = parse_qs(urlsplit(listing.args[1]).query)
        self.assertEqual(query["event"], ["pull_request"])
        self.assertNotIn("status", query)
        self.assertEqual(query["created"], [">=2026-11-11T12:00:00Z"])
        self.assertTrue(listing.kwargs["paginate"])

    def test_force_rechecks_replacement_after_normal_request(self):
        self.add_replacement()
        self.stuck = True
        self.sleep.side_effect = lambda seconds: self.current[2].update(
            status="completed", conclusion="cancelled"
        )
        self.assertEqual(self.perform(), 0)
        self.assertEqual(self.writes, ["repos/erf-model/ERF/actions/runs/1/cancel"])
        self.sleep.assert_called_once_with(60)

    def test_normal_409_on_persistent_run_reaches_force_fallback(self):
        self.errors = {("POST", "repos/erf-model/ERF/actions/runs/1/cancel"): 409}
        self.assertEqual(self.perform(), 0)
        self.assertEqual(self.writes, ["repos/erf-model/ERF/actions/runs/1/force-cancel"])
        self.assertIn("cancel returned HTTP 409", self.output.getvalue())
        self.assertIn("cancellation requests: 1", self.output.getvalue())
        self.assertEqual(self.sleep.call_args_list, [mock.call(60), mock.call(5)])

    def test_normal_409_completion_race_does_not_force_or_sleep(self):
        def finish_during_cancel(method, path, paginate=False):
            if method == "POST":
                self.current[1]["status"] = "completed"
                raise cleanup.APIError(409, "completion race")
            return self.api(method, path, paginate=paginate)

        self.client.side_effect = finish_during_cancel
        self.assertEqual(self.perform(), 0)
        self.assertEqual(self.writes, [])
        self.sleep.assert_not_called()

    def test_normal_409_does_not_bypass_force_eligibility(self):
        self.errors = {("POST", "repos/erf-model/ERF/actions/runs/1/cancel"): 409}
        for change in (
            lambda: self.prs[10].update(state="open"),
            lambda: self.current[1].update(run_attempt=2),
            lambda: self.current[1].update(event="push"),
        ):
            self.prs[10] = pr_record()
            self.current[1] = run_record()
            self.sleep.side_effect = lambda seconds, change=change: change()
            self.assertEqual(self.perform(), 0)
            self.assertEqual(self.writes, [])

    def test_force_409_completion_race_is_not_an_unresolved_error(self):
        self.stuck = True

        def finish_during_force(method, path, paginate=False):
            if method == "POST" and path.endswith("/force-cancel"):
                self.current[1]["status"] = "completed"
                raise cleanup.APIError(409, "completion race")
            return self.api(method, path, paginate=paginate)

        self.client.side_effect = finish_during_force
        self.assertEqual(self.perform(), 0)
        self.assertNotIn("Unresolved", self.output.getvalue())
        self.sleep.assert_called_once_with(60)

    def test_force_409_remaining_unfinished_is_an_unresolved_error(self):
        self.stuck = True
        self.errors = {("POST", "repos/erf-model/ERF/actions/runs/1/force-cancel"): 409}
        self.assertEqual(self.perform(), 1)
        self.assertIn(
            "Unresolved run 1: force-cancel returned HTTP 409; still queued", self.output.getvalue()
        )
        self.sleep.assert_called_once_with(60)

    def test_force_acceptance_is_verified_and_persistent_run_is_unresolved(self):
        self.stuck = self.force_stuck = True
        self.assertEqual(self.perform(), 1)
        self.assertEqual(len(self.writes), 2)
        self.assertIn("cancellation not confirmed; still queued", self.output.getvalue())
        self.assertIn("cancellation requests: 2; errors: 1", self.output.getvalue())
        self.assertEqual(self.sleep.call_args_list, [mock.call(60), mock.call(5)])

    def test_protected_events_workflows_and_uncertain_metadata_are_skipped(self):
        for record in (
            run_record(event="push"),
            run_record(event="workflow_dispatch"),
            run_record(event="schedule"),
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
        self.assertEqual(self.sleep.call_args_list, [mock.call(60), mock.call(5)])
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
        listings = [call for call in self.client.call_args_list if "/actions/runs?" in call.args[1]]
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
                    expected = (
                        ["repos/erf-model/ERF/actions/runs/1/force-cancel"]
                        if method == "POST" and status == 409
                        else []
                    )
                    self.assertEqual(self.writes, expected)
                    self.writes.clear()
                    self.current[1] = run_record()

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

    def test_development_builds_keep_cache_producing_runs(self):
        for name in ("gcc", "cuda-ci", "ci"):
            with self.subTest(workflow=name):
                workflow = yaml.load(
                    (ROOT / f".github/workflows/{name}.yml").read_text(), Loader=yaml.BaseLoader
                )
                self.assertEqual(
                    workflow["concurrency"]["cancel-in-progress"],
                    "${{ github.ref != 'refs/heads/development' }}",
                )
                self.assertIn("latest", workflow["concurrency"]["group"])

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
