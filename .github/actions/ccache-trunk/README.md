# Linux ccache trunks on GHCR

This action stores one ccache snapshot per workflow/job/variant at
`ghcr.io/<repository-owner>/erf-ccache:<sanitized-key>`. The current caller is `ci.yml`,
which runs the two Linux MPICH variants on development pushes and manual dispatch.
Only default-branch runs save. Other CI workflows still use Actions cache.

The runner needs Bash, Docker, GNU coreutils, and `gh` (all present on GitHub's Ubuntu
runners). The publishing job needs `contents: read` and `packages: write`. A new package
is linked to its source repository through the image label. If the package already
exists, grant the publishing repository access to it; deleting superseded versions also
requires package admin access. Private packages work with an authorized repository
token. Anonymous pulls require public visibility. Fork-origin PR access is not exercised
by this workflow and should be checked before migrating PR workflows.

Registry login, pull, and push get two attempts, with per-command timeouts. Transport
failures emit warnings and let compilation/tests continue. A missing image is a normal
cold build. Restores extract into a temporary sibling directory before replacing the
cache; failed extraction leaves existing files intact. The path must therefore be a
dedicated cache directory, not the home directory or filesystem root. Docker credentials
use a private temporary configuration that is removed when the action exits.

After a successful push, pruning keeps named variant tags and deletes superseded
versions. Pruning failures/timeouts are warnings and a later save retries cleanup.
All helpers live beside the action so `ci.yml`'s `.github/actions/**` path filter covers
changes to them. Concurrent development runs may publish in completion order; labels
record the producing revision, but publishing does not enforce revision ordering.

Run the fault-injection checks on Linux with:

```sh
python3 .github/actions/ccache-trunk/test_cache.py -v
```

The checks use mocked Docker/GitHub commands, require no network or credentials, and
also run in `ci.yml`. They cover retries, timeouts, failed/partial restores, replacement
rollback, failed saves, pruning races, credential cleanup, and invalid inputs.
