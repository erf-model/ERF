"""Common configuration, command handling, and interface for cache transports."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Protocol

TOKEN_VARIABLES = ("GHCR_TOKEN", "GH_TOKEN", "GITHUB_TOKEN")


@dataclass(frozen=True)
class Config:
    mode: str
    transport: str
    cache: Path
    ref: str
    platform: tuple[str, str]
    timeout: int
    runner_temp: Path
    output_path: Path | None
    token: str = field(repr=False)
    actor: str
    source: str
    revision: str
    prune_env: dict[str, str] = field(repr=False)


class Transport(Protocol):
    """Transfer files only; the driver owns publication policy and cache promotion.

    Each instance owns one operation and a private temporary directory.
    close() must also work after a partially failed setup or transfer.
    """

    def setup(self) -> None:
        """Prepare isolated credentials and any required executable."""
        ...

    def restore(self, destination: Path) -> bool:
        """Fill an empty staging directory; False means an explicit cache miss.

        All other failures raise. Never modify the live cache directory.
        """
        ...

    def save(self) -> None:
        """Publish the configured cache; return only after a successful push."""
        ...

    def close(self) -> None:
        """Release backend resources; temporary files belong to the driver."""
        ...


def diagnostic(value: object, secrets: Iterable[str | bytes | None] = (),
               limit: int = 1024) -> str:
    """Redact before truncating, then keep subprocess output on one log line."""
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    text = str(value or "")
    credentials = [*secrets, *(os.environ.get(name) for name in TOKEN_VARIABLES)]
    for secret in sorted({s for s in credentials if s}, key=len, reverse=True):
        if isinstance(secret, bytes):
            secret = secret.decode("utf-8", errors="replace")
        text = text.replace(secret, "[redacted]")
    text = re.sub(r"(?i)\b(Bearer|Basic)\s+[^\s,;\"']+", r"\1 [redacted]", text)
    # Cover JSON/key=value credentials without hiding HTTP status after /token:.
    text = re.sub(
        r'''(?i)((?<![\w/])(?:token|access_token|refresh_token|password|auth)\b'''
        r'''["']?\s*[:=]\s*["']?)[^\s,;&"']+''', r"\1[redacted]", text)
    text = re.sub(r"(https?://)[^/\s]+@", r"\1[redacted]@", text)
    text = re.sub(r'''(https?://[^\s?"']+)\?[^\s"']+''', r"\1?[redacted]", text)
    text = " ".join(text.split())
    return text if len(text) <= limit else text[:limit] + " [truncated]"


class CommandError(RuntimeError):
    def __init__(self, message: str, detail: str, returncode: int | None = None) -> None:
        super().__init__(message + (f"; {detail}" if detail else ""))
        self.detail = detail
        self.returncode = returncode

    def missing_manifest(self, ref: str) -> bool:
        # A generic 404 may be an authentication endpoint failure. Only explicit
        # manifest errors or ORAS's exact reference-not-found response are misses.
        if self.returncode is None or self.returncode <= 0 or "[truncated]" in self.detail:
            return False
        detail = self.detail.lower()
        if re.search(r"\b(401|403|5\d\d|unauthorized|forbidden|denied)\b", detail):
            return False
        return (bool(re.search(r"(?<![\w/:-])manifest[ _]unknown\b", detail))
                or bool(re.search(r"\bno such manifest\b|\bmanifest not found\b", detail))
                or f"{ref.lower()}: not found" in detail)


def warning(message: str) -> None:
    print("::warning::" + diagnostic(message, limit=2048).replace("%", "%25"), file=sys.stderr)


def run(command: Sequence[str], *, timeout: float = 60, attempts: int = 2,
        label: str | None = None, env: Mapping[str, str] | None = None,
        input: str | None = None, cwd: Path | None = None) -> str:
    """Retry bounded commands; report the final failure without argv or stdin."""
    for attempt in range(attempts):
        try:
            return subprocess.run(command, check=True, capture_output=True,
                                  text=True, timeout=timeout, env=env, input=input, cwd=cwd).stdout
        except (subprocess.SubprocessError, OSError) as exc:
            if attempt + 1 == attempts:
                if isinstance(exc, subprocess.CalledProcessError):
                    reason = f"exit status {exc.returncode}"
                elif isinstance(exc, subprocess.TimeoutExpired):
                    reason = f"timed out after {timeout}s"
                else:
                    reason = f"could not start ({type(exc).__name__})"
                detail = ((getattr(exc, "stderr", None) or "").strip()
                          or getattr(exc, "stdout", None) or getattr(exc, "strerror", ""))
                child_env = env or os.environ
                secrets = [input, *(child_env.get(name) for name in TOKEN_VARIABLES)]
                count = f"{attempts} attempt" + ("s" if attempts != 1 else "")
                message = f"{label or Path(command[0]).name}: {reason} after {count}"
                raise CommandError(message, diagnostic(detail, secrets),
                                   getattr(exc, "returncode", None)) from None
            time.sleep(2)
    raise ValueError("attempts must be positive")
