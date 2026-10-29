"""ORAS backend: a verified OCI artifact containing one ccache tar archive."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath, PureWindowsPath
import re
import shutil
import tarfile
from typing import TypedDict, cast
import zipfile

import transport
from transport import Config, CommandError

ORAS_VERSION = "1.3.0"
ORAS_RELEASES = {
    ("macOS", "X64"): ("darwin_amd64.tar.gz", "82c33f7da8430ea7fa7e7bdf7721be0a0d0481e5ccb2472ea438490d5e8641a9"),
    ("macOS", "ARM64"): ("darwin_arm64.tar.gz", "e10c6552c02d5a7c7eaf7170d3b6f7f094b675a98a1e0edf4d4478a909447245"),
    ("Windows", "X64"): ("windows_amd64.zip", "b050e93aa0dc7a79a61fa8e4074dfa302c41d4af01b634fe393c5dd687536aee"),
    ("Linux", "X64"): ("linux_amd64.tar.gz", "6cdc692f929100feb08aa8de584d02f7bcc30ec7d88bc2adc2054d782db57c64"),
    ("Linux", "ARM64"): ("linux_arm64.tar.gz", "7649738b48fde10542bcc8b0e9b460ba83936c75fb5be01ee6d4443764a14352"),
}
ARTIFACT_TYPE = "application/vnd.erf.ccache.v1"
LAYER_TYPE = "application/vnd.oci.image.layer.v1.tar"
ARCHIVE_NAME = "ccache.tar"
# Well above the existing 150-1100 MB cache limits; rejects damaged/hostile archives.
MAX_ARCHIVE_BYTES = 8 * 1024**3
MAX_MEMBERS = 500_000


class Layer(TypedDict):
    """Fields validated at the manifest boundary before fetching the archive."""

    mediaType: str
    digest: str
    size: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def install_oras(directory: Path, platform: tuple[str, str]) -> Path:
    filename, checksum = ORAS_RELEASES[platform]
    asset = f"oras_{ORAS_VERSION}_{filename}"
    archive = directory / asset
    transport.run(["curl", "--fail", "--silent", "--show-error", "--location", "--max-time", "60",
         "--output", str(archive),
         f"https://github.com/oras-project/oras/releases/download/v{ORAS_VERSION}/{asset}"],
        timeout=65, label="download ORAS")
    if sha256(archive) != checksum:
        raise ValueError("ORAS download checksum mismatch")
    executable = "oras.exe" if platform[0] == "Windows" else "oras"
    # Extract only the verified executable, never arbitrary archive paths.
    if filename.endswith(".zip"):
        with zipfile.ZipFile(archive) as bundle:
            contents = bundle.read(executable)
    else:
        with tarfile.open(archive) as bundle:
            source = bundle.extractfile(executable)
            if source is None:
                raise ValueError("ORAS executable is not a file")
            with source:
                contents = source.read()
    binary = directory / executable
    binary.write_bytes(contents)
    binary.chmod(0o700)
    return binary


def archive_path(name: str) -> PurePosixPath:
    """Accept ordinary relative paths on both POSIX and Windows."""
    path = PurePosixPath(name)
    if (not path.parts or path.is_absolute() or ".." in path.parts
            or PureWindowsPath(name).drive or "\\" in name or ":" in name
            or any(p.endswith((" ", ".")) for p in path.parts)
            or any(re.fullmatch(r"(?i:con|prn|aux|nul|com[1-9]|lpt[1-9])", p.split(".")[0])
                   for p in path.parts)):
        raise ValueError("unsafe cache archive path")
    return path


def unpack(archive: Path, destination: Path) -> None:
    """Write regular files/directories ourselves; never follow archive links."""
    total = 0
    with tarfile.open(archive, "r:") as bundle:
        for count, member in enumerate(bundle, 1):
            if count > MAX_MEMBERS or not (member.isfile() or member.isdir()):
                raise ValueError("unsupported cache archive member")
            relative = archive_path(member.name)
            total += member.size
            if total > MAX_ARCHIVE_BYTES:
                raise ValueError("cache archive exceeds size limit")
            target = destination.joinpath(*relative.parts)
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(member)
                if source is None:
                    raise ValueError("cache archive member is not a file")
                with source, target.open("xb") as output:
                    shutil.copyfileobj(source, output)
                # CI evicts objects not touched by the build; restoring must not
                # make every old object appear freshly used.
                os.utime(target, (member.mtime, member.mtime))


def pack(cache: Path, archive: Path) -> None:
    with tarfile.open(archive, "w") as bundle:
        for path in sorted(cache.rglob("*")):
            relative = path.relative_to(cache).as_posix()
            archive_path(relative)
            if path.is_symlink() or not (path.is_dir() or path.is_file()):
                raise ValueError("cache contains a link or special file")
            bundle.add(path, arcname=relative, recursive=False)
    if archive.stat().st_size > MAX_ARCHIVE_BYTES:
        raise ValueError("cache archive exceeds size limit")


def layer_descriptor(manifest: object) -> Layer:
    if (not isinstance(manifest, dict)
            or not isinstance(manifest.get("layers"), list)
            or manifest.get("artifactType") != ARTIFACT_TYPE
            or len(manifest.get("layers", [])) != 1):
        raise ValueError("unexpected cache artifact format")
    layer = manifest["layers"][0]
    if (not isinstance(layer, dict)
            or layer.get("mediaType") != LAYER_TYPE
            or not isinstance(layer.get("digest"), str)
            or not re.fullmatch(r"sha256:[0-9a-f]{64}", layer["digest"])
            or not isinstance(layer.get("size"), int)
            or not 0 < layer["size"] <= MAX_ARCHIVE_BYTES):
        raise ValueError("invalid cache layer")
    return cast(Layer, layer)


class OrasTransport:
    def __init__(self, config: Config, directory: Path) -> None:
        self.config = config
        self.directory = directory
        self.binary = directory / ("oras.exe" if config.platform[0] == "Windows" else "oras")
        self.registry_args = ["--registry-config", str(directory / "registry.json")]

    def setup(self) -> None:
        self.binary = install_oras(self.directory, self.config.platform)
        registry_file = self.directory / "registry.json"
        registry_file.write_text('{"auths":{}}', encoding="utf-8")
        registry_file.chmod(0o600)
        if self.config.token:
            transport.run([str(self.binary), "login", "ghcr.io", *self.registry_args,
                           "--username", self.config.actor, "--password-stdin"],
                          input=self.config.token, timeout=30, label="registry login")

    def restore(self, destination: Path) -> bool:
        ref = self.config.ref
        archive = self.directory / ARCHIVE_NAME
        try:
            manifest = json.loads(transport.run(
                [str(self.binary), "manifest", "fetch", *self.registry_args, ref],
                timeout=self.config.timeout, label="fetch cache manifest"))
        except CommandError as exc:
            if exc.missing_manifest(ref):
                return False
            raise
        layer = layer_descriptor(manifest)
        # Fetch an immutable blob; never trust an artifact-supplied filename.
        blob_ref = ref.rsplit(":", 1)[0] + "@" + layer["digest"]
        transport.run([str(self.binary), "blob", "fetch", *self.registry_args,
                       "--output", str(archive), blob_ref],
                      timeout=self.config.timeout, label="fetch cache archive")
        if (archive.stat().st_size != layer["size"]
                or "sha256:" + sha256(archive) != layer["digest"]):
            raise ValueError("cache archive digest/size mismatch")
        unpack(archive, destination)
        return True

    def save(self) -> None:
        pack(self.config.cache, self.directory / ARCHIVE_NAME)
        transport.run(
            [str(self.binary), "push", *self.registry_args, "--artifact-type", ARTIFACT_TYPE,
             "--annotation", "org.opencontainers.image.source=" + self.config.source,
             "--annotation", "org.opencontainers.image.revision=" + self.config.revision,
             self.config.ref, ARCHIVE_NAME + ":" + LAYER_TYPE],
            cwd=self.directory, timeout=self.config.timeout, label="publish cache")

    def close(self) -> None:
        # The driver removes the private binary, archive, and credential file.
        pass
