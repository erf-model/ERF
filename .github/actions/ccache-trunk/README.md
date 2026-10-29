# ccache trunks on GHCR

This action stores one rolling ccache snapshot per build configuration at
`ghcr.io/<repository-owner>/erf-ccache:<sanitized-key>`. All eight compiler-cache
workflows restore from GHCR. Only default-branch runs save; pull requests only read.
Saves happen before tests so an unrelated test failure does not freeze the cache.
These are compiler objects for incremental builds, not runnable ERF installations.

During migration, each workflow falls back to its existing Actions cache on a GHCR
miss or outage. This fallback is read-only: new snapshots go only to GHCR. The old
Actions entries eventually expire under GitHub's cache policy. A first run with
neither cache available builds normally and seeds its GHCR trunk on the default branch.

## Configuration names

`CCACHE_TRUNK_KEY` supplies the same key to restore and save, independently of workflow
and job names. Keys describe the runner environment, architecture, compiler/MPI or GPU
profile, precision, and relevant options. The current development matrix has 21 trunks:

| Workflow | Profile | Trunks |
| --- | --- | ---: |
| `ci.yml` | Ubuntu 22.04, x64, GCC, MPICH, double, particles off/on | 2 |
| `gcc.yml` | Ubuntu 22.04, x64, GCC, Open MPI, precision/particles | 3 |
| `gcc-rrtmgp.yml` | Ubuntu 22.04, x64, GCC, Spack/Open MPI, Kokkos, NetCDF/RRTMGP | 3 |
| `cuda-ci.yml` | Ubuntu 22.04, x64, CUDA 12.2, target 8.0, Open MPI, precision | 2 |
| `hip.yml` | Ubuntu 22.04, x64, ROCm 6.3.2, gfx908, Open MPI, precision | 2 |
| `sycl.yml` | Ubuntu 24.04, x64, oneAPI/SYCL, Intel MPI, precision | 2 |
| `macos.yml` | macOS latest, actual runner architecture, Apple Clang, Open MPI, precision | 2 |
| `windows.yml` | Windows latest, actual runner architecture, MSVC, MPI/precision/particles | 5 |

For example, `ubuntu-22.04-x64-gcc-mpich-double-particles-off` and
`ubuntu-22.04-x64-gcc-mpich-double-particles-on` retain the Linux pilot's names.
Matching mesh/particle precision and the other fixed build options are defined by
those workflows. Include new independent matrix dimensions in the key if coverage
expands. Update explicit version/target names when dependency profiles change.

`latest` describes the selected native runner image, not a pinned OS or toolchain.
ccache still checks source, compiler, headers, and flags; a toolchain or large header
update can cause a cold build even when the trunk downloads successfully. GHCR improves
cache availability and retention, not ccache's ability to reuse incompatible objects.

## Transports and failure handling

With `transport: auto`, Linux uses the existing Docker scratch-image transport;
macOS and Windows use an OCI artifact containing `ccache.tar` through ORAS. Tags must
not be shared between transports: their formats differ. The explicit `docker` and
`oras` settings are available for transport checks.

The ORAS backend downloads ORAS 1.3.0 from its official GitHub release and verifies a
platform-specific SHA256 before execution. It needs Python 3, curl, and gh, supplied by
GitHub's hosted runners. The Docker backend needs Python 3, Docker, and gh. Both use
only the Python standard library; native runners do not need a Docker daemon.

Login and transfer commands get two attempts. `transfer-timeout` defaults to 180 seconds
per transfer command (maximum 600), allowing for the existing CUDA cache limit of
1100 MB. Login has a 30-second limit and pruning has a 60-second overall budget.
Registry, download, and filesystem failures warn and let compilation/tests continue;
invalid action inputs fail. `hit` is true only after a complete restore.

Both backends report command failures with the operation, exit status or timeout, and up to
1024 characters of diagnostic output. Tokens, authentication fields, URL credentials
and query strings are redacted before truncation. Explicit missing-manifest responses
are normal cache misses; permission errors and ambiguous failures remain warnings.
Failed pruning includes the command's diagnostic details and remains nonfatal.

Both transports stage beside the destination before replacing it, keeping the previous
cache for rollback. The path must be a dedicated cache directory, not home or the root.
The native helper fetches an immutable blob digest from the manifest, verifies its
size/checksum, and rejects archive traversal, links, and special files. It preserves
file timestamps so CI's existing stale-object eviction remains effective. Credentials
use a private temporary configuration removed after the transfer.

## Access and retention

Publishing jobs need `contents: read` and `packages: write`. Source labels/annotations
link a new package to its repository. If it already exists, grant the publishing
repository access; deleting superseded versions also requires package admin access.
Private packages need an authorized token. Public visibility permits anonymous reads.
Verify package access from fork-origin PRs when adopting this upstream: those tokens
are read-only, and a denied restore warns and falls back to the old cache or a cold build.

After a successful push, `prune.py` keeps every named configuration tag and deletes
superseded untagged versions. Failed pruning is retried after a later save. There is no
age-based expiration of the active trunks. Renaming/retiring configurations requires
separate cleanup of their named tags. Concurrent development runs can publish in
completion order; revision metadata records the producer but does not enforce ordering.
All snapshots are transferred as complete archives/images, not per-object registry diffs.

## Local checks

Run all fault-injection tests on Linux or macOS (use `python` on Windows):

```sh
python3 -m unittest discover -s .github/actions/ccache-trunk -p 'test_*.py' -v
```

The tests need no network or credentials. They cover retries/timeouts, failed or
partial restores, rollback, archive integrity/safety, timestamp preservation, failed
saves, pruning retention/races, credential cleanup, redacted diagnostics, and invalid
inputs. The full suite runs in the existing `shellcheck` job in `lint.yml` for pull
requests and development pushes that change `.github/**`, and on manual dispatch.
It uses only the runner's Python standard library, without compiling ERF or accessing
GHCR. The same contract tests exercise both real backends with mocked external commands.

The same lint job runs workflow contracts and strict type checking with pinned,
check-only dependencies:

```sh
python3 -m venv /tmp/erf-workflow-contracts
/tmp/erf-workflow-contracts/bin/python -m pip install -r .github/tests/requirements.txt
/tmp/erf-workflow-contracts/bin/python -m mypy --config-file .github/tests/mypy.ini
/tmp/erf-workflow-contracts/bin/python .github/tests/cache_workflows.py
/tmp/erf-workflow-contracts/bin/python -m unittest discover -s .github/tests -p 'test_*.py' -v
```

These checks expand the push, PR, and manual matrices; require unique normalized
tags across configurations and supported architectures; and compare restore/save
keys, paths, and effective transports using the runtime's configuration parser.
They also check the legacy fallback and require the standard default-branch save
guard, `packages: write`, and `secrets.GITHUB_TOKEN`. HIP's installer argument and
cache name must both use `matrix.rocm_version`; update that matrix value when
upgrading ROCm. Mutation tests demonstrate that broken wiring is rejected.
No builds or registry requests are made.

`.github/tests/cache_workflows.py` supports the current matrix/env references,
PR-dependent `fromJSON` axes, Windows MPI choice, and native key setup statement.
New expression forms, runner labels, or triggers fail with a diagnostic until
explicitly supported and tested. Keep `REQUIRED_JOBS` in sync when retiring or
renaming jobs; additional cache jobs are discovered automatically. The manual
transport smoke workflow is excluded because it deliberately repeats transfers.
This is a static wiring check, not an interpreter for arbitrary build shell scripts.

Mypy checks `cache.py`, `transport.py`, `docker.py`, and `oras.py`, targeting Python
3.9 compatibility. The typed backend factory requires both implementations to
satisfy `Transport`; method signatures, return values, and driver calls are checked
together. `.github/tests/mypy.ini` keeps strict checking scoped to these modules.
The cache action itself still needs no third-party Python packages.

## Shared transport interface

`cache.py` is the entry point for both backends. It validates configuration once and
owns hit outputs, staging, cache promotion/rollback, temporary-file cleanup, warning
policy, and pruning after a successful save. `transport.py` defines the interface and
shared configuration, bounded command retries, and redacted diagnostics.

`DockerTransport` in `docker.py` and `OrasTransport` in `oras.py` each implement:

| Method | Responsibility |
| --- | --- |
| `setup()` | Prepare the executable and isolated credentials. |
| `restore(destination) -> bool` | Fill an empty staging directory; return false only for an explicit missing cache. Raise on other failures. |
| `save()` | Publish the configured cache; return only when publication succeeds. |
| `close()` | Release backend resources, including after a partial failure. |

Backends never replace the live cache, write action outputs, or prune registry versions.
The driver calls `close()` even after setup/transfer failure or interruption. Promotion
rolls back if interrupted between renames and propagates cancellation. `cache.sh` remains
a compatibility entry point for the original Docker helper.

To add a backend, implement these four methods, register its selection in `cache.py`,
and run the shared contract scenarios against it. Docker images and ORAS artifacts
retain their existing formats; a shared interface does not make their tags interchangeable.

## Manual GHCR smoke test

Run **GHCR cache smoke test** from the Actions tab, or use:

```sh
gh workflow run cache-smoke.yml
```

The workflow runs Linux, macOS Intel and Windows sequentially, with a 20-minute
limit per OS, followed by a read-only/anonymous restore job limited to 10 minutes.
It compiles a tiny program unique to that run, publishes its cache, removes local
objects and cache, then requires a restored cache hit. Changing the source must
miss the cache and produce new output. A second publication must replace the first
GHCR version, delete the superseded version, and restore the changed-source object.
The helper tests also run on each OS. ERF and its submodules are not built.

Tests use three rolling `transport-check-<os>-<arch>` tags in the same `erf-ccache`
package, separate from production configuration tags. Runs are serialized across
branches. Normal package pruning removes superseded untagged versions; the three
diagnostic tags remain available for inspection. Publishing requires repository
package write/admin access, and the anonymous check requires public package visibility.
The checks cover small transfers on the listed runners, not full-size GPU caches or
a hosted macOS ARM round trip. GitHub requires a manually dispatched workflow to be
registered on the repository's default branch before it can run on another branch.
