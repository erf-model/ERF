"""Exercise the trusted push cleanup with mocked GitHub APIs; never contact GitHub."""

import json
from pathlib import Path
import subprocess
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / ".github/scripts/cleanup-closed-pr-ci.cjs"

HARNESS = r"""
const fixture = JSON.parse(require('fs').readFileSync(0, 'utf8'));
const result = {cancelled: [], forced: [], failed: [], calls: [], sleeps: []};
const endpoints = {runs: 'runs', commits: 'commits', open: 'open'};
const context = {repo: {owner: 'erf-model', repo: 'ERF'}, runId: 9000};
const prReads = {};
const github = {
  rest: {
    actions: {
      listWorkflowRunsForRepo: endpoints.runs,
      getWorkflowRun: async p => {
        if (p.run_id === 9000) return {data: {created_at: '2026-11-12T12:00:00Z'}};
        const run = fixture.runs.find(r => r.id === p.run_id);
        if (result.cancelled.includes(run.id)) {
          return {data: {...run, status: fixture.stuck ? run.status : 'completed',
            run_attempt: fixture.restarted ? 2 : run.run_attempt}};
        }
        return {data: {...run, status: fixture.finished ? 'completed' : run.status}};
      },
      cancelWorkflowRun: async p => {
        if (fixture.cancelError) throw {status: fixture.cancelError, message: 'cancel error'};
        result.cancelled.push(p.run_id);
      },
    },
    repos: {listPullRequestsAssociatedWithCommit: endpoints.commits},
    pulls: {
      list: endpoints.open,
      get: async p => {
        const entries = fixture.prs[p.pull_number];
        if (!entries) throw {status: 404};
        const index = prReads[p.pull_number] || 0;
        prReads[p.pull_number] = index + 1;
        return {data: entries[Math.min(index, entries.length - 1)]};
      },
    },
  },
  paginate: async (method, p) => {
    result.calls.push({method, ...p});
    if (method === endpoints.runs) {
      return fixture.runs.filter(r => r.status === p.status || fixture.duplicate);
    }
    if (method === endpoints.commits) {
      if (fixture.associationError) throw {status: fixture.associationError, message: 'lookup error'};
      return fixture.associated || [];
    }
    if (method === endpoints.open) return fixture.open || [];
    throw new Error('Unexpected paginated API');
  },
  request: async (route, p) => {
    if (!route.endsWith('/force-cancel')) throw new Error('Unexpected write');
    if (fixture.forceError) throw {status: fixture.forceError, message: 'force error'};
    result.forced.push(p.run_id);
  },
};
const core = {info: () => {}, warning: () => {}, setFailed: s => result.failed.push(s)};
require(process.argv[1])({github, context, core, sleep: async ms => result.sleeps.push(ms)})
  .then(() => process.stdout.write(JSON.stringify(result)))
  .catch(error => {console.error(error); process.exitCode = 1;});
"""


def run_record(**changes):
    return dict({"id": 1, "event": "pull_request", "status": "queued", "run_attempt": 1,
                 "created_at": "2026-11-11T12:00:00Z", "head_sha": "abc",
                 "head_branch": "feature", "path": ".github/workflows/gcc.yml",
                 "head_repository": {"id": 42, "owner": {"login": "contributor"}},
                 "pull_requests": [{"number": 10}]}, **changes)


def pr_record(**changes):
    return dict({"number": 10, "state": "closed", "merged": True,
                 "base": {"repo": {"full_name": "erf-model/ERF"}},
                 "head": {"repo": {"id": 42}, "ref": "feature"}}, **changes)


class ClosedPrCleanupTests(unittest.TestCase):
    def run_cleanup(self, runs=None, prs=None, **options):
        fixture = {"runs": runs if runs is not None else [run_record()],
                   "prs": prs if prs is not None else {10: [pr_record()]}, **options}
        process = subprocess.run(["node", "-e", HARNESS, str(SCRIPT)],
                                 input=json.dumps(fixture), capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        return json.loads(process.stdout)

    def test_closed_and_merged_prs_across_unfinished_states(self):
        for merged in (False, True):
            for status in ("queued", "pending", "in_progress", "waiting", "requested"):
                with self.subTest(merged=merged, status=status):
                    result = self.run_cleanup([run_record(status=status)], {10: [pr_record(merged=merged)]})
                    self.assertEqual(result["cancelled"], [1])
                    self.assertEqual(result["forced"], [])
                    self.assertEqual(result["failed"], [])

    def test_open_or_reopened_pr_is_not_cancelled(self):
        result = self.run_cleanup(prs={10: [pr_record(state="open", merged=False)]})
        self.assertEqual(result["cancelled"], [])
        result = self.run_cleanup(stuck=True, prs={10: [pr_record(), pr_record(state="open")]})
        self.assertEqual(result["cancelled"], [1])
        self.assertEqual(result["forced"], [])

    def test_reused_source_branch_is_not_cancelled(self):
        result = self.run_cleanup(open=[pr_record(number=11, state="open")])
        self.assertEqual(result["cancelled"], [])

    def test_other_fork_with_same_branch_does_not_block_cleanup(self):
        other = pr_record(state="open", head={"repo": {"id": 99}, "ref": "feature"})
        self.assertEqual(self.run_cleanup(open=[other])["cancelled"], [1])

    def test_wrong_source_or_base_pr_cannot_authorize_cancellation(self):
        for pr in (pr_record(head={"repo": {"id": 99}, "ref": "feature"}),
                   pr_record(head={"repo": {"id": 42}, "ref": "another"}),
                   pr_record(head={"repo": None, "ref": "feature"}),
                   pr_record(base={"repo": {"full_name": "other/ERF"}})):
            with self.subTest(pr=pr):
                self.assertEqual(self.run_cleanup(prs={10: [pr]})["cancelled"], [])

    def test_missing_fork_association_is_resolved_by_commit(self):
        result = self.run_cleanup([run_record(pull_requests=[])], associated=[{"number": 10}])
        self.assertEqual(result["cancelled"], [1])
        lookup = next(c for c in result["calls"] if c["method"] == "commits")
        self.assertEqual(lookup["commit_sha"], "abc")
        self.assertEqual(lookup["per_page"], 100)
        result = self.run_cleanup([run_record(pull_requests=[])])
        self.assertEqual(result["cancelled"], [])

    def test_open_pr_association_blocks_a_closed_one(self):
        record = run_record(pull_requests=[{"number": 10}, {"number": 11}])
        result = self.run_cleanup([record], {10: [pr_record()], 11: [pr_record(number=11, state="open")]})
        self.assertEqual(result["cancelled"], [])

    def test_unavailable_association_never_authorizes_cancellation(self):
        for error in (404, 403, 500):
            with self.subTest(error=error):
                result = self.run_cleanup([run_record(pull_requests=[])], associationError=error)
                self.assertEqual(result["cancelled"], [])
                self.assertEqual(bool(result["failed"]), error != 404)

    def test_missing_pr_or_source_metadata_is_skipped(self):
        self.assertEqual(self.run_cleanup(prs={})["cancelled"], [])
        self.assertEqual(self.run_cleanup([run_record(head_repository=None)])["cancelled"], [])

    def test_non_pr_manual_completed_newer_and_unlisted_runs_are_untouched(self):
        cases = [run_record(event="push"), run_record(event="workflow_dispatch"),
                 run_record(event="pull_request_target"), run_record(status="completed"),
                 run_record(created_at="2026-11-12T12:00:01Z"),
                 run_record(created_at="invalid"),
                 run_record(path=".github/workflows/spack-containers.yml"),
                 run_record(path=".github/workflows/cleanup-merged-pr-ci.yml")]
        for record in cases:
            with self.subTest(run=record):
                result = self.run_cleanup([record])
                self.assertEqual(result["cancelled"], [])
                self.assertEqual(result["forced"], [])

    def test_finished_run_race_and_duplicate_listings(self):
        self.assertEqual(self.run_cleanup(finished=True)["cancelled"], [])
        result = self.run_cleanup(duplicate=True)
        self.assertEqual(result["cancelled"], [1])
        for call in result["calls"]:
            self.assertEqual(call["per_page"], 100)
        listings = [c for c in result["calls"] if c["method"] == "runs"]
        self.assertEqual(len(listings), 5)
        self.assertTrue(all(c["event"] == "pull_request" for c in listings))

    def test_force_cancel_only_after_grace_period_and_state_recheck(self):
        result = self.run_cleanup(stuck=True)
        self.assertEqual(result["sleeps"], [60000])
        self.assertEqual(result["forced"], [1])
        self.assertEqual(self.run_cleanup(stuck=True, restarted=True)["forced"], [])

    def test_api_errors_fail_visibly_but_completion_races_are_harmless(self):
        for error in (404, 409):
            self.assertEqual(self.run_cleanup(cancelError=error)["failed"], [])
            self.assertEqual(self.run_cleanup(stuck=True, forceError=error)["failed"], [])
        self.assertTrue(self.run_cleanup(cancelError=403)["failed"])
        self.assertTrue(self.run_cleanup(stuck=True, forceError=500)["failed"])

    def test_branch_metacharacters_remain_api_data(self):
        branch = "feature;$(do-not-execute)"
        result = self.run_cleanup([run_record(head_branch=branch)],
                                  {10: [pr_record(head={"repo": {"id": 42}, "ref": branch})]})
        self.assertEqual(result["cancelled"], [1])
        lookup = next(c for c in result["calls"] if c["method"] == "open")
        self.assertEqual(lookup["head"], "contributor:" + branch)

    def test_workflow_is_independent_trusted_push_without_path_filters(self):
        workflow = yaml.load((ROOT / ".github/workflows/cleanup-merged-pr-ci.yml").read_text(),
                             Loader=yaml.BaseLoader)
        self.assertEqual(workflow["on"], {"push": {"branches": ["development"]}})
        self.assertEqual(workflow["concurrency"]["cancel-in-progress"], "false")
        job = workflow["jobs"]["cleanup"]
        self.assertNotIn("needs", job)
        self.assertEqual(job["if"], "github.repository == 'erf-model/ERF'")
        self.assertEqual(job["permissions"], {"contents": "read", "actions": "write", "pull-requests": "read"})
        checkout = job["steps"][0]
        self.assertNotIn("ref", checkout.get("with", {}))
        self.assertEqual(checkout["with"]["persist-credentials"], "false")


if __name__ == "__main__":
    unittest.main()
