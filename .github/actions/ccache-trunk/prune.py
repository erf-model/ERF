#!/usr/bin/env python3
"""Keep named tags and the newest KEEP sha-only versions; prune the rest.

Uses gh on Linux, macOS and Windows. Bound the entire operation to 60 seconds;
failed deletions do not block other versions while time remains. A later successful
save retries incomplete cleanup. Concurrent DELETE 404s are harmless. OWNER,
OWNER_TYPE and GH_TOKEN identify the publishing repository.
"""

import os
import subprocess
import sys
import time


def github_api(arguments, deadline):
    for attempt in range(3):
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise RuntimeError("pruning timed out")
        try:
            result = subprocess.run(["gh", "api", *arguments], capture_output=True,
                                    text=True, timeout=min(15, remaining))
            if result.returncode == 0:
                return result.stdout
            if "DELETE" in arguments and "HTTP 404" in result.stderr + result.stdout:
                return None
            print(f"GitHub API request failed (attempt {attempt + 1}): {result.stderr.strip()}",
                  file=sys.stderr)
        except (OSError, subprocess.TimeoutExpired) as exc:
            print(f"GitHub API request failed: {exc}", file=sys.stderr)
        if attempt < 2:
            time.sleep(min(5, max(0, deadline - time.monotonic())))
    raise RuntimeError("GitHub API request failed after three attempts")


def select_versions(versions, keep):
    rows = [row.split("\t") for row in versions.splitlines() if row.strip()]
    sha_seen = 0
    for row in sorted(rows, key=lambda r: r[1], reverse=True):
        version, _, tags_text = row
        tags = [tag for tag in tags_text.split(",") if tag]
        if any(not tag.startswith("sha-") for tag in tags):
            continue
        if tags:
            sha_seen += 1
            if sha_seen <= keep:
                continue
        if not version.isdigit():
            raise ValueError("invalid package version ID")
        yield version


def main(specs):
    owner_type = os.environ["OWNER_TYPE"]
    if owner_type not in ("Organization", "User"):
        raise ValueError("unexpected OWNER_TYPE")
    route = "orgs" if owner_type == "Organization" else "users"
    owner = os.environ["OWNER"]
    deadline = time.monotonic() + 60
    failed = False
    for spec in specs:
        package, _, keep_text = spec.partition(":")
        keep = int(keep_text or os.environ.get("KEEP", "10"))
        if keep < 0:
            raise ValueError("KEEP must be nonnegative")
        endpoint = f"{route}/{owner}/packages/container/{package}/versions"
        try:
            versions = github_api(["--paginate", endpoint + "?per_page=100", "--jq",
                                   '.[] | [.id, .created_at, (.metadata.container.tags | join(","))] | @tsv'],
                                  deadline)
            for version in select_versions(versions, keep):
                if time.monotonic() >= deadline:
                    raise RuntimeError("pruning timed out")
                try:
                    result = github_api(["--method", "DELETE", "--silent", endpoint + "/" + version], deadline)
                except RuntimeError as exc:
                    print(f"{package}: could not delete version {version}: {exc}", file=sys.stderr)
                    failed = True
                    continue
                state = "was already gone" if result is None else "deleted"
                print(f"{package}: version {version} {state}")
        except (RuntimeError, ValueError) as exc:
            print(f"{package}: {exc}", file=sys.stderr)
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
