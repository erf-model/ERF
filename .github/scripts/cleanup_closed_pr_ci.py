#!/usr/bin/env python3
"""Cancel obsolete ERF PR CI. Use --dry-run to preview without cancelling."""

import argparse
import json
import os
import re
import subprocess
import time
from datetime import datetime
from urllib.parse import urlencode

ALLOWED_WORKFLOWS = {
    f".github/workflows/{name}.yml"
    for name in (
        "ci",
        "cuda-ci",
        "docs",
        "gcc-rrtmgp",
        "gcc",
        "hip",
        "lint",
        "macos",
        "parent-project",
        "style",
        "sycl",
        "windows",
    )
}
UNFINISHED_STATUSES = ("queued", "pending", "in_progress", "waiting", "requested")


class APIError(RuntimeError):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


def api(method, path, paginate=False):
    """gh uses GH_TOKEN; no credentials or PR metadata are interpolated into shell."""
    command = ["gh", "api", "--method", method]
    if paginate:
        command += ["--paginate", "--slurp"]
    command.append(path)
    try:
        result = subprocess.run(command, capture_output=True, text=True, timeout=60, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise APIError(None, "gh api could not complete") from exc
    if result.returncode:
        match = re.search(r"\(HTTP (\d{3})\)", result.stderr)
        status = int(match[1]) if match else None
        # Do not echo raw response text or credentials.
        raise APIError(status, f"gh api {method} failed (HTTP {status or 'unknown'})")
    return json.loads(result.stdout) if method == "GET" else None


def timestamp(value):
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return parsed.timestamp() if parsed.tzinfo is not None else None
    except (AttributeError, ValueError):
        return None


def eligible(run, cutoff):
    source = run.get("head_repository") or {}
    created = timestamp(run.get("created_at"))
    return bool(
        run.get("event") == "pull_request"
        and run.get("path") in ALLOWED_WORKFLOWS
        and run.get("status") in UNFINISHED_STATUSES
        and created is not None
        and created <= cutoff
        and source.get("id")
        and (source.get("owner") or {}).get("login")
        and run.get("id")
        and run.get("run_attempt")
        and run.get("head_branch")
        and run.get("head_sha")
    )


def same_source(pr, run):
    head = pr.get("head") or {}
    return (head.get("repo") or {}).get("id") == run["head_repository"]["id"] and head.get(
        "ref"
    ) == run["head_branch"]


def open_pr_reason(pr, run):
    head_sha = (pr.get("head") or {}).get("sha")
    if head_sha == run["head_sha"]:
        commit = "run matches its current commit"
    elif head_sha:
        commit = "run is from a different commit"
    else:
        commit = "commit comparison unavailable"
    return f"PR #{pr['number']} is open; {commit}"


def pr_skip_reason(repo, run):
    base = f"repos/{repo}"
    associations = run.get("pull_requests") or []
    if not associations:
        # Fork runs often omit associations. Resolve by commit, never branch name alone.
        pages = api("GET", f"{base}/commits/{run['head_sha']}/pulls?per_page=100", paginate=True)
        associations = [pr for page in pages for pr in page]
    matched = False
    for number in sorted({pr["number"] for pr in associations}):
        pr = api("GET", f"{base}/pulls/{number}")  # Fresh state, not cached closed state.
        if ((pr.get("base") or {}).get("repo") or {}).get("full_name") != repo:
            continue
        if not same_source(pr, run):
            continue
        if pr.get("state") == "open":
            return open_pr_reason(pr, run)
        if pr.get("state") != "closed":
            return f"PR #{number} is not confirmed closed"
        matched = True
    # Also identify open PRs when the commit association is missing. A branch
    # match can explain a skip, but must never authorize cancellation.
    owner = run["head_repository"]["owner"]["login"]
    query = urlencode({"state": "open", "head": f"{owner}:{run['head_branch']}", "per_page": 100})
    pages = api("GET", f"{base}/pulls?{query}", paginate=True)
    for page in pages:
        for pr in page:
            if same_source(pr, run):
                return open_pr_reason(pr, run)
    if not matched:
        return "no closed PR association matches this source repository and branch"
    return None


def cancel_if_closed(repo, original, cutoff, force=False, dry_run=False):
    run = api("GET", f"repos/{repo}/actions/runs/{original['id']}")
    if (
        not eligible(run, cutoff)
        or run["id"] != original["id"]
        or run["run_attempt"] != original["run_attempt"]
    ):
        print(f"Skipping run {original['id']}: finished, changed, or restarted")
        return False
    reason = pr_skip_reason(repo, run)
    if reason is not None:
        print(f"Skipping run {run['id']}: {reason}")
        return False
    operation = "force-cancel" if force else "cancel"
    if dry_run:
        print(f"Would {operation}: run {run['id']} ({run['path']})")
        return False
    try:
        api("POST", f"repos/{repo}/actions/runs/{run['id']}/{operation}")
    except APIError as exc:
        if exc.status == 409:  # A race on cancellation, not on arbitrary API reads.
            print(f"Run {run['id']} is no longer cancellable")
            return False
        raise
    print(f"Requested {operation}: run {run['id']} ({run['path']})")
    return True


def cancel_batch(repo, runs, cutoff, force=False, dry_run=False):
    requested, failures = [], 0
    for run in runs:
        try:
            if cancel_if_closed(repo, run, cutoff, force=force, dry_run=dry_run):
                requested.append(run)
        except APIError as exc:
            if exc.status == 404:
                print(f"Skipping run {run['id']}: run or PR became unavailable")
            else:
                print(f"Could not inspect/cancel run {run['id']}: {exc}")
                failures += 1
    return requested, failures


def cleanup(repo, sweep_run_id=None, dry_run=False, sleep=time.sleep):
    # Actions uses the sweep's server creation time; manual runs use their start time.
    cutoff = time.time()
    if sweep_run_id is not None:
        cutoff = timestamp(api("GET", f"repos/{repo}/actions/runs/{sweep_run_id}")["created_at"])
        if cutoff is None:
            raise ValueError("Invalid cleanup run creation time")
    runs = {}
    for status in UNFINISHED_STATUSES:
        query = urlencode({"event": "pull_request", "status": status, "per_page": 100})
        pages = api("GET", f"repos/{repo}/actions/runs?{query}", paginate=True)
        for page in pages:
            for run in page["workflow_runs"]:
                if eligible(run, cutoff):
                    runs[run["id"]] = run  # Deduplicate state changes between API pages.
    requested, failures = cancel_batch(repo, runs.values(), cutoff, dry_run=dry_run)
    if requested:
        # Normal cancellation can leave deeply queued runs stuck. Recheck before forcing.
        sleep(60)
        _, force_failures = cancel_batch(repo, requested, cutoff, force=True)
        failures += force_failures
    mode = "Dry run" if dry_run else "Cleanup"
    print(
        f"{mode}: examined {len(runs)} candidates; cancellation requests: {len(requested)}; errors: {failures}"
    )
    return failures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true", help="list proposed cancellations only")
    args = parser.parse_args()
    repo = os.environ.get("GITHUB_REPOSITORY", "erf-model/ERF")
    run_id = os.environ.get("GITHUB_RUN_ID")
    try:
        failures = cleanup(repo, int(run_id) if run_id else None, dry_run=args.dry_run)
    except (APIError, ValueError) as exc:
        print(f"Cleanup failed: {exc}")
        return 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
