#!/usr/bin/env python3
"""Cancel closed or superseded ERF PR CI. Use --dry-run to preview."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from collections.abc import Iterable
from datetime import datetime
from typing import Callable, Literal, NamedTuple, Optional, TypedDict, cast
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


class Owner(TypedDict):
    login: str


class Repository(TypedDict, total=False):
    id: int
    full_name: str
    owner: Owner | None


class PRBranch(TypedDict, total=False):
    repo: Repository | None
    ref: str
    sha: str


class PullRequest(TypedDict, total=False):
    number: int
    state: str
    head: PRBranch | None
    base: PRBranch | None


class Run(TypedDict, total=False):
    id: int
    event: str
    path: str
    status: str
    conclusion: str | None
    created_at: str
    run_attempt: int
    workflow_id: int | None
    run_number: int | None
    head_sha: str
    head_branch: str
    head_repository: Repository | None
    pull_requests: list[PullRequest] | None


class RunPage(TypedDict):
    workflow_runs: list[Run]


CancelResult = Optional[Literal["requested", "conflict", "unresolved"]]


class BatchResult(NamedTuple):
    requested: list[Run]
    follow_up: list[Run]
    failures: int


class APIError(RuntimeError):
    def __init__(self, status: int | None, message: str) -> None:
        super().__init__(message)
        self.status = status


def api(method: Literal["GET", "POST"], path: str, paginate: bool = False) -> object:
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
    return cast(object, json.loads(result.stdout)) if method == "GET" else None


def read_run(repo: str, run_id: int) -> Run:
    # These casts describe GitHub's response shapes, not runtime schema validation.
    return cast(Run, api("GET", f"repos/{repo}/actions/runs/{run_id}"))


def timestamp(value: str | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
        return parsed.timestamp() if parsed.tzinfo is not None else None
    except (AttributeError, ValueError):
        return None


def eligible(run: Run, cutoff: float) -> bool:
    source: Repository = run.get("head_repository") or {}
    owner = source.get("owner")
    created = timestamp(run.get("created_at"))
    return bool(
        run.get("event") == "pull_request"
        and run.get("path") in ALLOWED_WORKFLOWS
        and run.get("status") in UNFINISHED_STATUSES
        and created is not None
        and created <= cutoff
        and source.get("id")
        and owner is not None
        and owner.get("login")
        and run.get("id")
        and run.get("run_attempt")
        and run.get("head_branch")
        and run.get("head_sha")
    )


def same_source(pr: PullRequest, run: Run) -> bool:
    source = run.get("head_repository")
    if source is None:
        return False
    head: PRBranch = pr.get("head") or {}
    return (head.get("repo") or {}).get("id") == source["id"] and head.get("ref") == run[
        "head_branch"
    ]


def open_pr_reason(pr: PullRequest, run: Run) -> str:
    head_sha = (pr.get("head") or {}).get("sha")
    if head_sha == run["head_sha"]:
        commit = "run matches its current commit"
    elif head_sha:
        commit = "run is from a different commit"
    else:
        commit = "commit comparison unavailable"
    return f"PR #{pr['number']} is open; {commit}"


def pr_skip_reason(repo: str, run: Run) -> str | None:
    base = f"repos/{repo}"
    associations = run.get("pull_requests") or []
    if not associations:
        # Fork runs often omit associations. Resolve by commit, never branch name alone.
        pages = cast(
            list[list[PullRequest]],
            api("GET", f"{base}/commits/{run['head_sha']}/pulls?per_page=100", paginate=True),
        )
        associations = [pr for page in pages for pr in page]
    matched = False
    for number in sorted({pr["number"] for pr in associations}):
        pr = cast(PullRequest, api("GET", f"{base}/pulls/{number}"))  # Fresh state.
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
    source = run.get("head_repository")
    owner = source.get("owner") if source else None
    if owner is None:
        return "source repository owner unavailable"
    query = urlencode(
        {"state": "open", "head": f"{owner['login']}:{run['head_branch']}", "per_page": 100}
    )
    pages = cast(list[list[PullRequest]], api("GET", f"{base}/pulls?{query}", paginate=True))
    for page in pages:
        for pr in page:
            if same_source(pr, run):
                return open_pr_reason(pr, run)
    if not matched:
        return "no closed PR association matches this source repository and branch"
    return None


def is_replacement(newer: Run, run: Run) -> bool:
    """A different fork, workflow, or known PR is not a replacement."""
    source = run.get("head_repository")
    if source is None:
        return False
    old_prs = {pr["number"] for pr in run.get("pull_requests") or []}
    new_prs = {pr["number"] for pr in newer.get("pull_requests") or []}
    status = newer.get("status")
    useful = status in UNFINISHED_STATUSES or (
        status == "completed" and newer.get("conclusion") in {"success", "failure", "timed_out"}
    )
    return bool(
        useful
        and newer.get("event") == "pull_request"
        and newer.get("workflow_id") == run["workflow_id"]
        and newer.get("path") == run["path"]
        and (newer.get("head_repository") or {}).get("id") == source["id"]
        and newer.get("head_branch") == run["head_branch"]
        and (newer.get("run_number") or 0) > (run.get("run_number") or 0)
        and newer.get("id") != run["id"]
        and (not old_prs or not new_prs or bool(old_prs & new_prs))
    )


def newer_replacement(repo: str, run: Run) -> Run | None:
    if not run.get("workflow_id") or not run.get("run_number"):
        return None
    # Include completed replacements, not just the unfinished sweep candidates.
    query = urlencode(
        {
            "event": "pull_request",
            "branch": run["head_branch"],
            "created": f">={run['created_at']}",
            "per_page": 100,
        }
    )
    pages = cast(
        list[RunPage],
        api(
            "GET",
            f"repos/{repo}/actions/workflows/{run['workflow_id']}/runs?{query}",
            paginate=True,
        ),
    )
    replacements = [
        newer for page in pages for newer in page["workflow_runs"] if is_replacement(newer, run)
    ]
    if not replacements:
        return None
    newest = max(replacements, key=lambda newer: newer.get("run_number") or 0)
    fresh = read_run(repo, newest["id"])
    return fresh if fresh["id"] == newest["id"] and is_replacement(fresh, run) else None


def unchanged_unfinished(run: Run, original: Run, cutoff: float) -> bool:
    return (
        eligible(run, cutoff)
        and run["id"] == original["id"]
        and run["run_attempt"] == original["run_attempt"]
    )


def cancel_if_obsolete(
    repo: str, original: Run, cutoff: float, force: bool = False, dry_run: bool = False
) -> CancelResult:
    run = read_run(repo, original["id"])
    if not unchanged_unfinished(run, original, cutoff):
        print(f"Skipping run {original['id']}: finished, changed, or restarted")
        return None
    replacement = newer_replacement(repo, run)
    if replacement:
        reason = f"superseded by run {replacement['id']}"
    else:
        skip = pr_skip_reason(repo, run)
        if skip is not None:
            print(f"Skipping run {run['id']}: {skip}")
            return None
        reason = "source PR is closed or merged"
    operation = "force-cancel" if force else "cancel"
    if dry_run:
        print(f"Would {operation}: run {run['id']} ({run['path']}): {reason}")
        return None
    try:
        api("POST", f"repos/{repo}/actions/runs/{run['id']}/{operation}")
    except APIError as exc:
        if exc.status != 409:
            raise
        fresh = read_run(repo, run["id"])
        if not unchanged_unfinished(fresh, original, cutoff):
            print(f"Skipping run {run['id']}: finished, changed, or restarted after HTTP 409")
            return None
        if force:
            print(
                f"Unresolved run {run['id']}: force-cancel returned HTTP 409; still {fresh['status']}"
            )
            return "unresolved"
        print(f"Run {run['id']}: cancel returned HTTP 409; will recheck before force-cancel")
        return "conflict"
    print(f"Requested {operation}: run {run['id']} ({run['path']}): {reason}")
    return "requested"


def cancel_batch(
    repo: str, runs: Iterable[Run], cutoff: float, force: bool = False, dry_run: bool = False
) -> BatchResult:
    requested: list[Run] = []
    follow_up: list[Run] = []
    failures = 0
    for run in runs:
        try:
            result = cancel_if_obsolete(repo, run, cutoff, force=force, dry_run=dry_run)
            if result == "requested":
                requested.append(run)
            if not force and result in ("requested", "conflict"):
                follow_up.append(run)
            if result == "unresolved":
                failures += 1
        except APIError as exc:
            if exc.status == 404:
                print(f"Skipping run {run['id']}: run or PR became unavailable")
            else:
                print(f"Could not inspect/cancel run {run['id']}: {exc}")
                failures += 1
    return BatchResult(requested, follow_up, failures)


def verify_cancellations(repo: str, runs: Iterable[Run], cutoff: float) -> int:
    failures = 0
    for original in runs:
        try:
            run = read_run(repo, original["id"])
            if unchanged_unfinished(run, original, cutoff):
                print(
                    f"Unresolved run {run['id']}: cancellation not confirmed; still {run['status']}"
                )
                failures += 1
        except APIError as exc:
            if exc.status != 404:
                print(f"Could not verify run {original['id']}: {exc}")
                failures += 1
    return failures


def cleanup(
    repo: str,
    sweep_run_id: int | None = None,
    dry_run: bool = False,
    sleep: Callable[[float], None] = time.sleep,
) -> int:
    # Actions uses the sweep's server creation time; manual runs use their start time.
    cutoff = time.time()
    if sweep_run_id is not None:
        created = timestamp(read_run(repo, sweep_run_id)["created_at"])
        if created is None:
            raise ValueError("Invalid cleanup run creation time")
        cutoff = created
    runs: dict[int, Run] = {}
    for status in UNFINISHED_STATUSES:
        query = urlencode({"event": "pull_request", "status": status, "per_page": 100})
        pages = cast(list[RunPage], api("GET", f"repos/{repo}/actions/runs?{query}", paginate=True))
        for page in pages:
            for run in page["workflow_runs"]:
                if eligible(run, cutoff):
                    runs[run["id"]] = run  # Deduplicate state changes between API pages.
    batch = cancel_batch(repo, runs.values(), cutoff, dry_run=dry_run)
    requested, failures = batch.requested, batch.failures
    if batch.follow_up:
        # An accepted request or HTTP 409 is not proof of completion. Recheck
        # eligibility, attempt, and obsolescence before the single force fallback.
        sleep(60)
        forced = cancel_batch(repo, batch.follow_up, cutoff, force=True)
        failures += forced.failures
        requested.extend(forced.requested)
        if forced.requested:
            sleep(5)
            failures += verify_cancellations(repo, forced.requested, cutoff)
    mode = "Dry run" if dry_run else "Cleanup"
    print(
        f"{mode}: examined {len(runs)} candidates; cancellation requests: {len(requested)}; errors: {failures}"
    )
    return failures


def main() -> int:
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
