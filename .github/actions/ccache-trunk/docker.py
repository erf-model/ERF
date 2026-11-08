"""Docker backend: the existing scratch image with compiler objects in /ccache."""

from __future__ import annotations

import os
from pathlib import Path
import uuid

import transport
from transport import Config, CommandError


class DockerTransport:
    def __init__(self, config: Config, directory: Path) -> None:
        self.config = config
        self.directory = directory
        self.env = dict(os.environ, DOCKER_CONFIG=str(directory / "docker"))
        self.container: str | None = None

    def command(self, *args: str, timeout: float = 60, attempts: int = 2,
                label: str | None = None, input: str | None = None) -> str:
        return transport.run(["docker", *args], env=self.env, timeout=timeout,
                             attempts=attempts, label=label, input=input)

    def setup(self) -> None:
        Path(self.env["DOCKER_CONFIG"]).mkdir(mode=0o700)
        if self.config.token:
            self.command("login", "ghcr.io", "--username", self.config.actor,
                         "--password-stdin", input=self.config.token,
                         timeout=30, label="registry login")

    def restore(self, destination: Path) -> bool:
        try:
            self.command("pull", self.config.ref, timeout=self.config.timeout,
                         label="fetch cache image")
        except CommandError as exc:
            if exc.missing_manifest(self.config.ref):
                return False
            raise
        # Record the name before creation so even an interrupted create is cleaned.
        self.container = "erf-ccache-" + uuid.uuid4().hex
        self.command("create", "--name", self.container, self.config.ref, "/trunk",
                     timeout=15, attempts=1, label="open cache image")
        self.command("cp", self.container + ":/ccache/.", str(destination) + "/",
                     timeout=self.config.timeout, attempts=1, label="extract cache image")
        return True

    def save(self) -> None:
        dockerfile = self.directory / "Dockerfile"
        dockerfile.write_text("FROM scratch\nCOPY . /ccache\n", encoding="utf-8")
        self.command("build", "--quiet", "-t", self.config.ref,
                     "--label", "org.opencontainers.image.source=" + self.config.source,
                     "--label", "org.opencontainers.image.revision=" + self.config.revision,
                     "-f", str(dockerfile), str(self.config.cache),
                     timeout=self.config.timeout, attempts=1, label="package cache")
        self.command("push", "--quiet", self.config.ref,
                     timeout=self.config.timeout, label="publish cache")

    def close(self) -> None:
        if self.container:
            self.command("rm", self.container, timeout=10, attempts=1,
                         label="remove cache container")
