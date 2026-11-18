"""Check CI metrics collection and incomplete-data handling without GitHub."""

import contextlib
import importlib.util
import io
from datetime import datetime, timedelta, timezone
from pathlib import Path
import tempfile
import unittest
from unittest import mock
from urllib.parse import parse_qs, urlsplit

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ci_metrics", ROOT / ".github/workflows/metrics/ci_metrics.py")
metrics = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(metrics)
BASE = datetime(2026, 11, 18, tzinfo=timezone.utc)


def time_at(seconds):
    return (BASE + timedelta(seconds=seconds)).strftime("%Y-%m-%dT%H:%M:%SZ")


def run_record(identifier=1, **changes):
    return {"id": identifier, "name": "Linux GCC", "head_sha": "a" * 40,
            "status": "completed", "conclusion": "success", "run_attempt": 1,
            "created_at": time_at(0), "updated_at": time_at(120), **changes}


def job_record(identifier=1, **changes):
    return {"id": identifier, "status": "completed", "conclusion": "success", "runner_id": 9,
            "created_at": time_at(0), "started_at": time_at(180), "completed_at": time_at(540),
            **changes}


class RunCollectionTests(unittest.TestCase):
    def test_partition_before_cap_and_include_boundary_records(self):
        records = [run_record(1, created_at=time_at(0)), run_record(2, created_at=time_at(1)),
                   run_record(3, created_at=time_at(2), name="ParentProject"),
                   run_record(4, created_at=time_at(3))]
        queries = []

        def api(path, paginate=True):
            query = parse_qs(urlsplit(path).query)
            queries.append(query)
            self.assertEqual(query["event"], ["push"])
            self.assertEqual(query["branch"], ["development"])
            left, right = query["created"][0].split("..")
            matched = [r for r in records if left <= r["created_at"] <= right]
            total = 1000 if left == time_at(0) and right == time_at(3) else len(matched)
            return [{"total_count": total, "workflow_runs": matched}]

        with mock.patch.object(metrics, "gh", side_effect=api):
            actual = metrics.runs("example/ERF", "push", BASE, BASE + timedelta(seconds=3), "development")
        self.assertEqual({r["id"] for r in actual}, {1, 2, 3, 4})
        self.assertTrue(any(q["created"] == [f"{time_at(0)}..{time_at(1)}"] for q in queries))
        self.assertTrue(any(q["created"] == [f"{time_at(2)}..{time_at(3)}"] for q in queries))
        self.assertFalse(any(q["created"] == [f"{time_at(0)}..{time_at(3)}"] and q["per_page"] == ["100"]
                             for q in queries))

    def test_oversubscribed_one_second_fails_closed(self):
        with mock.patch.object(metrics, "gh", return_value=[{"total_count": 1000, "workflow_runs": []}]):
            with self.assertRaisesRegex(RuntimeError, "one second; collection incomplete"):
                metrics.runs("example/ERF", "pull_request", BASE, BASE)

    def test_api_error_and_invalid_window_are_not_hidden(self):
        with mock.patch.object(metrics, "gh", side_effect=RuntimeError("API unavailable")):
            with self.assertRaisesRegex(RuntimeError, "API unavailable"):
                metrics.runs("example/ERF", "pull_request", BASE, BASE)
        with mock.patch.object(metrics, "gh") as api:
            with self.assertRaisesRegex(ValueError, "window starts after"):
                metrics.runs("example/ERF", "pull_request", BASE + timedelta(seconds=1), BASE)
            api.assert_not_called()

    def test_truncated_or_duplicate_pages_fail_closed(self):
        for records in ([run_record(1)], [run_record(1), run_record(1)]):
            with self.subTest(records=records), mock.patch.object(
                metrics, "gh", side_effect=[[{"total_count": 2}], [{"total_count": 2, "workflow_runs": records}]]
            ):
                with self.assertRaisesRegex(RuntimeError, "GitHub returned"):
                    metrics.runs("example/ERF", "pull_request", BASE, BASE + timedelta(seconds=1))


class SignalTests(unittest.TestCase):
    def test_observed_terminal_cohort_and_incomplete_cohorts(self):
        self.assertEqual(metrics.signal_minutes([run_record()]), 2.0)
        self.assertEqual(metrics.signal_minutes([run_record(), run_record(2, conclusion="failure")]), 2.0)
        for changes in ({"status": "queued", "conclusion": None},
                        {"status": "in_progress", "conclusion": None},
                        {"conclusion": "cancelled"}, {"conclusion": "skipped"}):
            with self.subTest(changes=changes):
                self.assertIsNone(metrics.signal_minutes([run_record(), run_record(2, **changes)]))
        self.assertIsNone(metrics.signal_minutes([]))
        self.assertIsNone(metrics.signal_minutes([run_record(updated_at=time_at(-1))]))


class JobSampleTests(unittest.TestCase):
    def test_completed_assigned_jobs_and_partial_groups(self):
        group = [run_record(), run_record(2, status="queued", conclusion=None)]
        jobs = [job_record(), job_record(2, conclusion="skipped", runner_id=None),
                job_record(3, conclusion="cancelled", runner_id=None)]
        with mock.patch.object(metrics, "gh", return_value=[{"total_count": 3, "jobs": jobs}]) as api:
            count, minutes, queue, partial = metrics.sample_jobs("example/ERF", group)
        self.assertIn("/runs/1/attempts/1/jobs?", api.call_args.args[0])
        self.assertEqual((count, minutes, queue, partial), (1, 6.0, [3.0], True))
        self.assertEqual(metrics.sample_value("6.0", partial, count > 0), "6.0 (partial)")

    def test_unknown_jobs_are_not_reported_as_zero(self):
        with mock.patch.object(metrics, "gh", return_value=[{"total_count": 0, "jobs": []}]):
            count, minutes, queue, partial = metrics.sample_jobs("example/ERF", [run_record()])
        self.assertEqual((count, minutes, queue, partial), (0, 0.0, [], True))
        self.assertEqual(metrics.sample_value(count, partial, count > 0), "n/a (partial)")

    def test_incomplete_job_pages_fail_closed(self):
        with mock.patch.object(metrics, "gh", return_value=[{"total_count": 2, "jobs": [job_record()]}]):
            with self.assertRaisesRegex(RuntimeError, "incomplete job list"):
                metrics.sample_jobs("example/ERF", [run_record()])

    def test_latest_attempt_or_cancellation_is_partial(self):
        for conclusion, attempt in (("success", 2), ("cancelled", 1)):
            with self.subTest(conclusion=conclusion, attempt=attempt):
                run = run_record(conclusion=conclusion, run_attempt=attempt)
                job = job_record(conclusion=conclusion)
                with mock.patch.object(metrics, "gh", return_value=[{"total_count": 1, "jobs": [job]}]) as api:
                    count, minutes, queue, partial = metrics.sample_jobs("example/ERF", [run])
                self.assertIn(f"/runs/1/attempts/{attempt}/jobs?", api.call_args.args[0])
                self.assertEqual((count, minutes, queue, partial), (1, 6.0, [3.0], True))

    def test_missing_attempt_and_inconsistent_times_are_partial(self):
        self.assertEqual(metrics.sample_jobs("example/ERF", [run_record(run_attempt=None)]), (0, 0.0, [], True))
        jobs = [job_record(completed_at=time_at(0))]
        with mock.patch.object(metrics, "gh", return_value=[{"total_count": 1, "jobs": jobs}]):
            self.assertEqual(metrics.sample_jobs("example/ERF", [run_record()]), (0, 0.0, [], True))


class ReportTests(unittest.TestCase):
    def test_sample_zero_and_honest_labels(self):
        records = [run_record(), run_record(2, head_sha="b" * 40, status="queued", conclusion=None),
                   run_record(3, name="ParentProject", head_sha="c" * 40)]
        with mock.patch.object(metrics, "runs", return_value=records), \
             mock.patch.object(metrics, "sample_jobs") as sample, \
             mock.patch("sys.argv", ["ci_metrics.py", "--days", "1", "--sample", "0"]), \
             mock.patch.dict("os.environ", {}, clear=True), contextlib.redirect_stdout(io.StringIO()) as output:
            self.assertEqual(metrics.main(), 0)
        sample.assert_not_called()
        text = output.getvalue()
        self.assertIn("| Development push workflow records returned by API | 3 |", text)
        self.assertIn("| Selected development push workflow records | 2 |", text)
        self.assertIn("Observed development SHA duration, median min (n=1)", text)
        self.assertIn("do not establish full required-check coverage", text)
        self.assertIn("latest 0 development SHAs", text)
        self.assertNotIn("Time to full signal", text)
        self.assertNotIn("per merge", text)

    def test_summary_appends_the_report(self):
        with tempfile.TemporaryDirectory() as directory:
            summary = Path(directory) / "summary.md"
            summary.write_text("existing\n")
            with mock.patch.object(metrics, "runs", return_value=[]), \
                 mock.patch("sys.argv", ["ci_metrics.py", "--days", "1", "--sample", "0"]), \
                 mock.patch.dict("os.environ", {"GITHUB_STEP_SUMMARY": str(summary)}), \
                 contextlib.redirect_stdout(io.StringIO()) as output:
                self.assertEqual(metrics.main(), 0)
            self.assertEqual(summary.read_text(), "existing\n" + output.getvalue().rstrip("\n") + "\n")


if __name__ == "__main__":
    unittest.main()
