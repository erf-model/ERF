#!/usr/bin/env python3
"""Weekly CI numbers for ERF, from the GitHub API, as a Markdown table.

Four things decide whether the CI is doing its job, and none of them is visible
on the Actions tab:

  merges and red merges   how many pushes to development in the window, and how
                          many had at least one failed build workflow
  time to full signal     from the first workflow starting to the last one
                          finishing, per push to development and per PR commit
  queue wait              how long jobs sat queued before a runner picked them
                          up (merge tier, sampled); the org cap is 20 jobs
  runner-minutes per merge  what one merge costs, sampled from the latest pushes

Runs weekly from metrics.yml into the job summary; also runs by hand:
    GH_TOKEN=... python3 ci_metrics.py [--repo erf-model/ERF] [--days 7]
"""
import argparse, collections, json, statistics, subprocess, sys
from datetime import datetime, timedelta, timezone

BUILD_WORKFLOWS = {"ERF CI", "Linux GCC", "Linux GCC NetCDF RRTMGP", "ERF CI (cuda)",
                   "ERF CI (hip)", "ERF CI (sycl)", "MacOS", "Windows", "Windows MPI",
                   "codespell", "Style"}


def gh(path, paginate=True):
    cmd = ["gh", "api"] + (["--paginate"] if paginate else []) + [path]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
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


def runs(repo, event, since, branch=None):
    q = f"repos/{repo}/actions/runs?event={event}&created=>={since:%Y-%m-%dT%H:%M:%SZ}&per_page=100"
    if branch: q += f"&branch={branch}"
    out = []
    for page in gh(q):
        out += [r for r in page.get("workflow_runs", []) if r["name"] in BUILD_WORKFLOWS]
    return out


def by_sha(rs):
    d = collections.defaultdict(list)
    for r in rs: d[r["head_sha"]].append(r)
    return d


def signal_minutes(group):
    done = [r for r in group if r["status"] == "completed" and r["conclusion"] in ("success", "failure")]
    if len(done) < 4: return None
    t0 = min(ts(r["created_at"]) for r in done); t1 = max(ts(r["updated_at"]) for r in done)
    return (t1 - t0).total_seconds() / 60


def median(xs):
    xs = [x for x in xs if x is not None]
    return f"{statistics.median(xs):.0f}" if xs else "n/a"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", default="erf-model/ERF")
    ap.add_argument("--days", type=int, default=7)
    ap.add_argument("--sample", type=int, default=3, help="latest pushes to cost out via the jobs API")
    a = ap.parse_args()
    since = datetime.now(timezone.utc) - timedelta(days=a.days)

    dev = by_sha(runs(a.repo, "push", since, branch="development"))
    prs = by_sha(runs(a.repo, "pull_request", since))

    red = [s for s, g in dev.items() if any(r["conclusion"] == "failure" for r in g)]
    red_by = collections.Counter(r["name"] for s in red for r in dev[s] if r["conclusion"] == "failure")
    dev_signal = [signal_minutes(g) for g in dev.values()]
    pr_signal = [signal_minutes(g) for g in prs.values()]

    # cost and queue wait, sampled from the latest pushes
    latest = sorted(dev.items(), key=lambda kv: max(ts(r["created_at"]) for r in kv[1]))[-a.sample:]
    jobs_per, minutes_per, queue = [], [], []
    for sha, group in latest:
        nj = 0; mins = 0.0
        for r in group:
            js = [j for page in gh(f"repos/{a.repo}/actions/runs/{r['id']}/jobs?per_page=100")
                  for j in page.get("jobs", []) if j.get("started_at") and j.get("completed_at")]
            if not js: continue
            nj += len(js)
            mins += sum((ts(j["completed_at"]) - ts(j["started_at"])).total_seconds() / 60 for j in js)
            queue.append((min(ts(j["started_at"]) for j in js) - ts(r["created_at"])).total_seconds() / 60)
        jobs_per.append(nj); minutes_per.append(mins)

    rows = [
        ("Pushes to development", f"{len(dev)}"),
        ("  with a failed build workflow", f"{len(red)}" + (f" ({', '.join(f'{k} {v}' for k, v in red_by.most_common(4))})" if red_by else "")),
        ("PR commits built", f"{len(prs)}"),
        ("Time to full signal, merge, median min", median(dev_signal)),
        ("Time to full signal, PR commit, median min", median(pr_signal)),
        (f"Queue wait, merge tier, median min (last {a.sample} pushes)", median(queue)),
        (f"Jobs per merge (last {a.sample})", "/".join(str(x) for x in jobs_per) or "n/a"),
        (f"Runner-minutes per merge (last {a.sample})", "/".join(f"{x:.0f}" for x in minutes_per) or "n/a"),
    ]
    out = [f"### CI, last {a.days} days, {a.repo}, as of {datetime.now(timezone.utc):%Y-%m-%d %H:%M} UTC", "",
           "| Measure | Value |", "|---|---|"] + [f"| {k} | {v} |" for k, v in rows]
    print("\n".join(out))
    import os
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as f: f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
