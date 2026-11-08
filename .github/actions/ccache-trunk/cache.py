#!/usr/bin/env python3
"""Run a cache transfer with the same policy for Docker and ORAS."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import os
from pathlib import Path
import re
import shutil
import signal
import sys
import tarfile
import tempfile
from types import FrameType
from typing import NoReturn
import zipfile

from docker import DockerTransport
from oras import ORAS_RELEASES, OrasTransport
import transport
from transport import Config, Transport, diagnostic, warning


def read_config(env: Mapping[str, str]) -> Config:
    """Validate action inputs before any installation or registry access."""
    mode = env["CACHE_MODE"]
    if mode not in ("restore", "save"):
        raise ValueError("mode must be restore or save")
    backend = env.get("CACHE_TRANSPORT", "auto")
    platform = (env["RUNNER_OS"], env.get("RUNNER_ARCH", ""))
    if backend == "auto":
        backend = "docker" if platform[0] == "Linux" else "oras"
    if backend not in ("docker", "oras"):
        raise ValueError("transport must be auto, docker, or oras")
    if backend == "docker" and platform[0] != "Linux":
        raise ValueError("Docker cache transport requires a Linux runner")
    if backend == "oras" and platform not in ORAS_RELEASES:
        raise ValueError("unsupported ORAS runner platform")
    token = env.get("GHCR_TOKEN", "")
    if mode == "save" and not token:
        raise ValueError("token is required when saving")
    tag = re.sub(r"[^A-Za-z0-9_.-]+", "-", env["CACHE_KEY"]).lower()
    if not re.fullmatch(r"[a-z0-9_][a-z0-9_.-]{0,127}", tag):
        raise ValueError("key must produce a valid GHCR tag")
    raw_path = env["CACHE_PATH"]
    cache = Path(raw_path).expanduser().resolve()
    if not raw_path or cache in (Path(cache.anchor), Path.home().resolve()):
        raise ValueError("path must name a dedicated ccache directory")
    if cache.exists() and not cache.is_dir():
        raise ValueError("cache path is not a directory")
    timeout = env.get("CACHE_TIMEOUT", "180")
    if not re.fullmatch(r"[1-9][0-9]{0,2}", timeout) or int(timeout) > 600:
        raise ValueError("transfer timeout must be between 1 and 600 seconds")
    owner = env["GITHUB_REPOSITORY_OWNER"]
    saving = mode == "save"
    owner_type = env["CACHE_OWNER_TYPE"] if saving else ""
    if saving and owner_type not in ("User", "Organization"):
        raise ValueError("unexpected CACHE_OWNER_TYPE")
    return Config(
        mode=mode, transport=backend, cache=cache,
        ref=f"ghcr.io/{owner.lower()}/erf-ccache:{tag}", platform=platform,
        timeout=int(timeout), runner_temp=Path(env["RUNNER_TEMP"]),
        output_path=Path(env["GITHUB_OUTPUT"]) if not saving else None,
        token=token, actor=env["GITHUB_ACTOR"] if token else "",
        source=env["CACHE_SOURCE"] if saving else "",
        revision=env["GITHUB_SHA"] if saving else "",
        prune_env=dict(env, GH_TOKEN=token, OWNER=owner, OWNER_TYPE=owner_type) if saving else {},
    )


def promote(staged: Path, cache: Path) -> None:
    """Replace the cache with rollback on failure or interruption."""
    previous = staged.with_name(staged.name + ".previous")
    try:
        if cache.exists():
            cache.rename(previous)
        staged.rename(cache)
    finally:
        # This also runs for KeyboardInterrupt and SIGTERM's SystemExit. Use the
        # filesystem state so interruption just after a rename is handled too.
        if previous.exists():
            if not cache.exists():
                try:
                    previous.rename(cache)
                except OSError:
                    warning(f"original cache remains at {previous}")
            else:
                try:
                    shutil.rmtree(previous)
                except OSError as exc:
                    warning(f"could not remove previous cache at {previous}: {exc}")


def restore(config: Config, backend: Transport) -> None:
    """Publish a hit only after a complete, staged restore and promotion."""
    config.cache.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=config.cache.name + ".ghcr-",
                                     dir=config.cache.parent) as tmp:
        staged = Path(tmp)
        if not backend.restore(staged):
            print(f"no cache at {config.ref}; continuing with the local cache")
            return
        promote(staged, config.cache)
    output(config.output_path, "hit=true")
    print(f"restored {config.ref}")


def prune(config: Config) -> None:
    """A pruning failure must not turn a successful save into a failure."""
    try:
        print(transport.run(
            [sys.executable, str(Path(__file__).with_name("prune.py")), "erf-ccache:0"],
            env=config.prune_env, timeout=65, attempts=1, label="prune cache versions"), end="")
    except RuntimeError as exc:
        warning(f"could not prune old GHCR ccache versions ({exc}); a later save will retry")


def output(path: Path | None, value: str) -> None:
    if path is None:
        raise ValueError("restore requires an action output path")
    with path.open("a", encoding="utf-8") as stream:
        stream.write(value + "\n")


def main() -> int:
    try:
        config = read_config(dict(os.environ))
    except (KeyError, ValueError, OSError) as exc:
        print(f"::error::{diagnostic(exc)}", file=sys.stderr)
        return 2

    try:
        if config.mode == "restore":
            output(config.output_path, "hit=false")
        elif not any(p.is_file() for p in config.cache.rglob("*")):
            print(f"no ccache files in {config.cache}: skip save")
            return 0
        with tempfile.TemporaryDirectory(prefix="erf-ccache-", dir=config.runner_temp) as tmp:
            factory: Callable[[Config, Path], Transport]
            factory = DockerTransport if config.transport == "docker" else OrasTransport
            backend: Transport = factory(config, Path(tmp))
            try:
                backend.setup()
                if config.mode == "restore":
                    restore(config, backend)
                else:
                    backend.save()
                    print(f"saved {config.ref}")
                    prune(config)
            finally:
                try:
                    backend.close()
                except (OSError, RuntimeError) as exc:
                    warning(f"could not clean up cache transport: {exc}")
    except (OSError, RuntimeError, ValueError, tarfile.TarError, zipfile.BadZipFile) as exc:
        warning(f"could not {config.mode} GHCR ccache ({exc}); continuing with the local cache")
    return 0


if __name__ == "__main__":
    def interrupted(signum: int, _frame: FrameType | None) -> NoReturn:
        raise SystemExit(128 + signum)

    signal.signal(signal.SIGTERM, interrupted)
    sys.exit(main())
