#!/usr/bin/env python3
"""Summarize selected ERF CI workflow records from the GitHub API.

Durations cover observed SHA cohorts, not necessarily every required check.
Runs weekly from metrics.yml or locally with an authenticated gh CLI.
"""
import argparse, collections, json, re, statistics, subprocess, sys
from datetime import datetime, timedelta, timezone
from urllib.parse import urlencode

BUILD_WORKFLOWS = {"ERF CI", "Linux GCC", "Linux GCC NetCDF RRTMGP", "ERF CI (cuda)",
                   "ERF CI (hip)", "ERF CI (sycl)", "MacOS", "Windows", "Windows MPI",
                   "codespell", "Style"}


def gh(path, paginate=True):
    cmd = ["gh", "api"] + (["--paginate"] if paginate else []) + [path]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError as exc:
        match = re.search(r"\bHTTP ([1-5][0-9]{2})\b", exc.stderr or "")
        status = f" (HTTP {match.group(1)})" if match else ""
        raise RuntimeError(f"GitHub API request failed{status}; check permissions, rate limits and service status") from None
    except OSError:
        raise RuntimeError("Unable to run the GitHub CLI (gh)") from None
    # --paginate concatenates JSON objects; each has a workflow_runs / jobs list
    items = []
    dec = json.JSONDecoder(); i = 0
    while i < len(out):
        while i < len(out) and out[i].isspace(): i += 1
        if i >= len(out): break
        obj, j = dec.raw_decode(out, i); i = j
        items.append(obj)
    return items


def ts(s):
    return datetime.fromisoformat(s.replace("Z", "+00:00"))


def runs(repo, event, since, until, branch=None):
    """Split filtered queries before GitHub's 1,000-result limit."""
    first, last = int(since.timestamp()), int(until.timestamp())
    if first > last:
        raise ValueError("collection window starts after it ends")
    found = {}

    def collect(left, right):
        start = datetime.fromtimestamp(left, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        end = datetime.fromtimestamp(right, timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
        params = {"event": event, "created": f"{start}..{end}"}
        if branch:
            params["branch"] = branch
        endpoint = f"repos/{repo}/actions/runs?"
        count_page = gh(endpoint + urlencode({**params, "per_page": 1}), paginate=False)
        total = count_page[0]["total_count"]
        if total >= 1000:
            if left == right:
                raise RuntimeError(f"GitHub returned {total} runs in one second; collection incomplete")
            middle = (left + right) // 2
            collect(left, middle)
            collect(middle + 1, right)
            return
        pages = gh(endpoint + urlencode({**params, "per_page": 100}))
        records = [run for page in pages for run in page["workflow_runs"]]
        if len(records) != total:
            raise RuntimeError(f"GitHub returned {len(records)} of {total} runs for {start}..{end}")
        if len({run["id"] for run in records}) != len(records):
            raise RuntimeError(f"GitHub returned duplicate runs for {start}..{end}")
        for run in records:
            if run["id"] in found:
                raise RuntimeError(f"GitHub returned run {run['id']} in multiple time buckets")
            found[run["id"]] = run

    collect(first, last)
    return list(found.values())


def by_sha(rs):
    d = collections.defaultdict(list)
    for r in rs: d[r["head_sha"]].append(r)
    return d


def median(xs):
    xs = [x for x in xs if x is not None]
    return f"{statistics.median(xs):.0f}" if xs else "n/a"


def outcomes(records):
    counts = collections.Counter((run["conclusion"] or "unknown") if run["status"] == "completed"
                                 else "unfinished" for run in records)
    return ", ".join(f"{name}: {count}" for name, count in sorted(counts.items())) or "none"


def sample_jobs(repo, group):
    """Measure finished, assigned jobs; flag anything that makes totals partial."""
    count, minutes, queue = 0, 0.0, []
    partial = any(run["status"] != "completed" or run["conclusion"] not in ("success", "failure")
                  for run in group)
    for run in group:
        if run["status"] != "completed":
            continue
        attempt = run.get("run_attempt")
        if not isinstance(attempt, int) or attempt < 1:
            partial = True
            continue
        if attempt > 1:
            partial = True
        path = f"repos/{repo}/actions/runs/{run['id']}/attempts/{attempt}/jobs?per_page=100"
        pages = gh(path)
        jobs = [job for page in pages for job in page["jobs"]]
        if not pages or len(jobs) != pages[0]["total_count"] or len({job["id"] for job in jobs}) != len(jobs):
            raise RuntimeError(f"GitHub returned an incomplete job list for run {run['id']}")
        if not jobs and run["conclusion"] not in ("skipped", "cancelled"):
            partial = True
        for job in jobs:
            if job["status"] != "completed":
                partial = True
                continue
            if job["conclusion"] == "skipped":
                continue
            if not job.get("runner_id") or not job.get("started_at") or not job.get("completed_at"):
                partial = True
                continue
            try:
                started, completed = ts(job["started_at"]), ts(job["completed_at"])
                if started.tzinfo is None or completed.tzinfo is None or completed < started:
                    raise ValueError("invalid job timestamps")
            except (AttributeError, TypeError, ValueError):
                partial = True
                continue
            count += 1
            minutes += (completed - started).total_seconds() / 60
            if job.get("created_at"):
                try:
                    created = ts(job["created_at"])
                    if created.tzinfo is not None and created <= started:
                        queue.append((started - created).total_seconds() / 60)
                except (AttributeError, TypeError, ValueError):
                    pass
    return count, minutes, queue, partial


def sample_value(value, partial, has_data):
    if partial and not has_data:
        return "n/a (partial)"
    return f"{value} (partial)" if partial else str(value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="erf-model/ERF")
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--sample", type=int, default=3, help="latest development SHAs for job sampling; 0 disables it")
    a = ap.parse_args()
    if a.days <= 0 or a.sample < 0:
        ap.error("days must be positive and sample must be nonnegative")
    until = datetime.now(timezone.utc).replace(microsecond=0)
    since = until - timedelta(days=a.days)

    dev_records = runs(a.repo, "push", since, until, branch="development")
    pr_records = runs(a.repo, "pull_request", since, until)
    dev_runs = [run for run in dev_records if run["name"] in BUILD_WORKFLOWS]
    pr_runs = [run for run in pr_records if run["name"] in BUILD_WORKFLOWS]
    dev, prs = by_sha(dev_runs), by_sha(pr_runs)
    red = [group for group in dev.values() if any(run["conclusion"] == "failure" for run in group)]
    red_by = collections.Counter(run["name"] for group in red for run in group if run["conclusion"] == "failure")

    latest = sorted(dev.items(), key=lambda item: max(ts(run["created_at"]) for run in item[1]))
    latest = latest[-a.sample:] if a.sample else []
    jobs_per, minutes_per, queue = [], [], []
    queue_partial = False
    for sha, group in latest:
        count, minutes, waits, partial = sample_jobs(a.repo, group)
        jobs_per.append(f"{sha[:12]}: {sample_value(count, partial, count > 0)}")
        minutes_per.append(f"{sha[:12]}: {sample_value(f'{minutes:.1f}', partial, count > 0)}")
        queue.extend(waits)
        queue_partial = queue_partial or partial or len(waits) != count

    rows = [
        ("Development push workflow records returned by API", str(len(dev_records))),
        ("Selected development push workflow records", str(len(dev_runs))),
        ("  Workflow record outcomes", outcomes(dev_runs)),
        ("Selected development push SHAs", str(len(dev))),
        ("  SHAs with any observed workflow failure", str(len(red))),
        ("  Failed workflow records by name (top 4)",
         ", ".join(f"{name}: {count}" for name, count in red_by.most_common(4)) or "none"),
        ("Pull-request workflow records returned by API", str(len(pr_records))),
        ("Selected pull-request workflow records", str(len(pr_runs))),
        ("  Workflow record outcomes", outcomes(pr_runs)),
        ("Selected pull-request run SHAs", str(len(prs))),
        (f"Completed assigned-job created-to-start delay, median min (latest {len(latest)} development SHAs; n={len(queue)})",
         sample_value(median(queue), queue_partial, bool(queue))),
        (f"Observed completed runner jobs (latest {len(latest)} development SHAs)", " / ".join(jobs_per) or "n/a"),
        (f"Observed completed runner-minutes (latest {len(latest)} development SHAs)", " / ".join(minutes_per) or "n/a"),
    ]
    out = [f"### Selected CI records for {a.repo}",
           f"Window: {since:%Y-%m-%d %H:%M:%S} to {until:%Y-%m-%d %H:%M:%S} UTC.",
           f"Workflow scope: {', '.join(sorted(BUILD_WORKFLOWS))}.", "",
           "SHAs can span PRs or reruns. Job samples use completed latest attempts; partial values omit work. "
           "Created-to-start "
           "can include non-runner waiting; runner-minutes are not monetary cost.", "",
           "| Measure | Value |", "|---|---|"] + [f"| {k} | {v} |" for k, v in rows]
    print("\n".join(out))
    import os
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f: f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
