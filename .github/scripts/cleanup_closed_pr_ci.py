#!/usr/bin/env python3
"""Clean up closed-PR CI runs. Read-only unless --apply is supplied."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import time
from collections.abc import Callable
from datetime import datetime
from typing import Any
from urllib.parse import quote, urlencode

JsonObject = dict[str, Any]
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
CANCELLATION_GRACE_SECONDS = 60


class GitHubAPIError(RuntimeError):
    def __init__(self, status: int | None, message: str):
        super().__init__(message)
        self.status = status


class GitHubAPI:
    """Use gh's GH_TOKEN authentication and pagination, without shell execution."""

    def __init__(self, repository: str):
        self.base = f"repos/{repository}"

    def _request(self, method: str, endpoint: str, paginate: bool = False) -> Any:
        command = ["gh", "api", "--method", method]
        if paginate:
            command.extend(["--paginate", "--slurp"])
        command.append(f"{self.base}/{endpoint}")
        try:
            result = subprocess.run(
                command, capture_output=True, text=True, timeout=60, check=False
            )
        except (OSError, subprocess.TimeoutExpired) as exc:
            raise GitHubAPIError(None, "Could not complete gh api request") from exc
        if result.returncode != 0:
            match = re.search(r"\(HTTP (\d{3})\)", result.stderr)
            status = int(match[1]) if match else None
            # Do not echo raw CLI output, credentials, or untrusted response text.
            raise GitHubAPIError(status, f"GitHub API {method} failed (HTTP {status or 'unknown'})")
        if method == "POST":
            return None  # Cancellation endpoints return an empty success response.
        try:
            return json.loads(result.stdout)
        except json.JSONDecodeError as exc:
            raise GitHubAPIError(None, "GitHub API returned invalid JSON") from exc

    def _list(self, endpoint: str, key: str | None = None, **params: Any) -> list[JsonObject]:
        pages = self._request("GET", endpoint + "?" + urlencode(params), paginate=True)
        if not isinstance(pages, list):
            raise GitHubAPIError(None, "Expected paginated GitHub API response")
        items = []
        for page in pages:
            if key is not None:
                if not isinstance(page, dict):
                    raise GitHubAPIError(None, "Expected a GitHub API page object")
                entries = page.get(key)
            else:
                entries = page
            if not isinstance(entries, list) or not all(isinstance(item, dict) for item in entries):
                raise GitHubAPIError(None, "Expected a list of GitHub API objects")
            items.extend(entries)
        return items

    def _get_object(self, endpoint: str) -> JsonObject:
        response = self._request("GET", endpoint)
        if not isinstance(response, dict):
            raise GitHubAPIError(None, "Expected a GitHub API object")
        return response

    def get_run(self, run_id: int) -> JsonObject:
        return self._get_object(f"actions/runs/{run_id}")

    def list_pr_runs(self, status: str) -> list[JsonObject]:
        return self._list(
            "actions/runs", key="workflow_runs", event="pull_request", status=status, per_page=100
        )

    def get_pr(self, number: int) -> JsonObject:
        return self._get_object(f"pulls/{number}")

    def prs_for_commit(self, sha: str) -> list[JsonObject]:
        return self._list(f"commits/{quote(sha, safe='')}/pulls", per_page=100)

    def open_prs_for_branch(self, owner: str, branch: str) -> list[JsonObject]:
        return self._list("pulls", state="open", head=f"{owner}:{branch}", per_page=100)

    def cancel_run(self, run_id: int, force: bool = False) -> None:
        operation = "force-cancel" if force else "cancel"
        self._request("POST", f"actions/runs/{run_id}/{operation}")


def parse_timestamp(value: Any) -> datetime | None:
    if not isinstance(value, str):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    return parsed if parsed.tzinfo is not None else None


def candidate_skip_reason(run: JsonObject, cutoff: datetime) -> str | None:
    """None means eligible for further checks, not permission to cancel."""
    if run.get("event") != "pull_request":
        return "not a pull_request run"
    if run.get("path") not in ALLOWED_WORKFLOWS:
        return "workflow is not allowlisted"
    if run.get("status") not in UNFINISHED_STATUSES:
        return "run is no longer unfinished"
    created = parse_timestamp(run.get("created_at"))
    if created is None or created > cutoff:
        return "creation time is missing, invalid, or newer than this sweep"
    source = run.get("head_repository") or {}
    if not source.get("id") or not (source.get("owner") or {}).get("login"):
        return "source repository identity is unavailable"
    if not run.get("head_branch") or not run.get("head_sha"):
        return "source branch or commit is unavailable"
    for field in ("id", "run_attempt"):
        value = run.get(field)
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            return f"{field} is missing or invalid"
    return None


def matches_source_branch(pr: JsonObject, run: JsonObject) -> bool:
    head = pr.get("head") or {}
    source = head.get("repo") or {}
    return bool(
        source.get("id") == run["head_repository"]["id"] and head.get("ref") == run["head_branch"]
    )


def pr_skip_reason(api: GitHubAPI, repository: str, run: JsonObject) -> str | None:
    associations = run.get("pull_requests") or []
    if not associations:
        # Fork runs often omit associations. A matching branch name alone is not proof.
        associations = api.prs_for_commit(run["head_sha"])
    if not associations:
        return "no PR association could be established"

    numbers = set()
    for association in associations:
        number = association.get("number")
        if not isinstance(number, int) or isinstance(number, bool) or number < 1:
            return "PR association has no valid number"
        numbers.add(number)
    matched = False
    for number in sorted(numbers):
        pr = api.get_pr(number)  # Re-read state; a previously closed PR can be reopened.
        base_repository = (pr.get("base") or {}).get("repo") or {}
        if base_repository.get("full_name") != repository or not matches_source_branch(pr, run):
            continue
        if pr.get("state") != "closed":
            return f"associated PR #{number} is not closed"
        matched = True
    if not matched:
        return "no associated PR matches this repository and source branch"

    # A closed PR's branch may now belong to a different open PR.
    source_owner = run["head_repository"]["owner"]["login"]
    open_prs = api.open_prs_for_branch(source_owner, run["head_branch"])
    if any(matches_source_branch(pr, run) for pr in open_prs):
        return "source branch is used by an open PR"
    return None


def cancel_if_obsolete(
    api: GitHubAPI,
    repository: str,
    original: JsonObject,
    cutoff: datetime,
    apply: bool,
    force: bool,
    log: Callable[[str], None],
) -> bool:
    """Recheck an immutable run ID and its PRs before requesting cancellation."""
    run = api.get_run(original["id"])
    reason = candidate_skip_reason(run, cutoff)
    if reason is None and run["id"] != original["id"]:
        reason = "run identity changed"
    if reason is None and run["run_attempt"] != original["run_attempt"]:
        reason = "run was restarted as a different attempt"
    if reason is None:
        reason = pr_skip_reason(api, repository, run)
    if reason is not None:
        log(f"Skipping run {original['id']}: {reason}")
        return False
    operation = "force-cancel" if force else "cancel"
    if not apply:
        log(f"Would {operation} run {run['id']} ({run['path']})")
        return False
    try:
        api.cancel_run(run["id"], force=force)
    except GitHubAPIError as exc:
        # These races are harmless only on the cancellation request, not arbitrary reads.
        if exc.status in (404, 409):
            log(f"Run {run['id']} is unavailable or no longer cancellable")
            return False
        raise
    log(f"Requested {operation}: run {run['id']} ({run['path']})")
    return True


def cleanup(
    api: GitHubAPI,
    repository: str,
    sweep_run_id: int,
    apply: bool = False,
    sleep: Callable[[float], None] = time.sleep,
    log: Callable[[str], None] = print,
) -> int:
    """Return an error count. Uncertain associations never authorize cancellation."""
    sweep = api.get_run(sweep_run_id)
    cutoff = parse_timestamp(sweep.get("created_at"))
    if cutoff is None:
        raise ValueError("Missing or invalid cleanup run creation time")
    runs = {}
    for status in UNFINISHED_STATUSES:
        for run in api.list_pr_runs(status):
            reason = candidate_skip_reason(run, cutoff)
            if reason is None:
                runs[run["id"]] = run  # Deduplicate runs that change state between pages.
            else:
                log(f"Skipping listed run {run.get('id', 'unknown')}: {reason}")

    failures = 0
    cancelled = []
    for run in runs.values():
        try:
            if cancel_if_obsolete(api, repository, run, cutoff, apply, force=False, log=log):
                cancelled.append(run)
        except GitHubAPIError as exc:
            if exc.status == 404:
                log(f"Skipping run {run['id']}: run or PR became unavailable")
            else:
                log(f"Could not inspect/cancel run {run['id']}: {exc}")
                failures += 1
    if cancelled:
        # A normal request can leave deeply queued runs stuck. Recheck before forcing.
        sleep(CANCELLATION_GRACE_SECONDS)
        for run in cancelled:
            try:
                cancel_if_obsolete(api, repository, run, cutoff, apply, force=True, log=log)
            except GitHubAPIError as exc:
                if exc.status == 404:
                    log(f"Skipping run {run['id']}: run or PR became unavailable")
                else:
                    log(f"Could not inspect/force-cancel run {run['id']}: {exc}")
                    failures += 1
    mode = "Apply" if apply else "Dry run"
    log(
        f"{mode}: examined {len(runs)} candidates; normal cancellation requests: {len(cancelled)}; errors: {failures}"
    )
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=os.environ.get("GITHUB_REPOSITORY"))
    parser.add_argument("--run-id", type=int, default=os.environ.get("GITHUB_RUN_ID"))
    parser.add_argument(
        "--apply", action="store_true", help="send cancellation requests (default: read-only)"
    )
    args = parser.parse_args()
    if not args.repo or not re.fullmatch(r"[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+", args.repo):
        parser.error("provide --repo OWNER/REPO or GITHUB_REPOSITORY")
    if not args.run_id or args.run_id < 1:
        parser.error("provide a positive --run-id or GITHUB_RUN_ID")
    try:
        failures = cleanup(GitHubAPI(args.repo), args.repo, args.run_id, apply=args.apply)
    except (GitHubAPIError, ValueError) as exc:
        print(f"Cleanup failed: {exc}")
        return 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
