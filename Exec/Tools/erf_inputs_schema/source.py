"""Read an ERF commit and its AMReX by `git archive` (read-only), or take two directory trees as they are."""
import io
import subprocess
import tarfile
from pathlib import Path

ERF_PATHS = ["Source", "CMake", "CMakeLists.txt"]
AMREX_PATHS = ["Src"]


def _git(repo, *args):
    return subprocess.run(["git", "-C", str(repo), *args], check=True, capture_output=True).stdout


def amrex_gitlink(repo, erf_commit):
    """The AMReX commit ERF's Submodules/AMReX points at, at erf_commit."""
    line = _git(repo, "ls-tree", erf_commit, "Submodules/AMReX").decode().split()
    if len(line) < 3 or line[1] != "commit":
        raise ValueError(f"no Submodules/AMReX gitlink at {erf_commit}")
    return line[2]


def resolve(repo, commit):
    return _git(repo, "rev-parse", "--verify", f"{commit}^{{commit}}").decode().strip()


def extract(repo, commit, paths, dest):
    data = _git(repo, "archive", commit, *paths)
    # filter= exists from Python 3.12 (and 3.8.17, 3.9.17, 3.10.12, 3.11.4); the archive is our own git's output
    safe = {"filter": "data"} if hasattr(tarfile, "data_filter") else {}
    with tarfile.open(fileobj=io.BytesIO(data)) as tar:
        tar.extractall(dest, **safe)


def checkout(repo, erf_commit, dest, amrex_repo=None, amrex_commit=None):
    """Extract ERF at erf_commit and AMReX at its gitlink (or amrex_commit) under dest; return
    (erf_root, amrex_root, erf_commit, amrex_commit) with full SHAs."""
    repo, dest = Path(repo), Path(dest)
    erf_commit = resolve(repo, erf_commit)
    if not amrex_repo and not (repo / "Submodules" / "AMReX" / ".git").exists():
        raise ValueError("Submodules/AMReX is not initialized: run 'git submodule update --init Submodules/AMReX', "
                         "or give an AMReX clone with --amrex-repo")
    amrex_repo = Path(amrex_repo) if amrex_repo else repo / "Submodules" / "AMReX"
    amrex_commit = resolve(amrex_repo, amrex_commit or amrex_gitlink(repo, erf_commit))
    extract(repo, erf_commit, ERF_PATHS, dest / "erf")
    extract(amrex_repo, amrex_commit, AMREX_PATHS, dest / "AMReX")
    return dest / "erf", dest / "AMReX", erf_commit, amrex_commit


def generator_commit(tools_dir):
    """The last commit that touched the generator, or None outside git (determinism holds within a commit)."""
    try:
        out = _git(tools_dir, "log", "-1", "--format=%H", "--", "build_erf_schema.py", "erf_inputs_schema")
    except (subprocess.CalledProcessError, FileNotFoundError):
        return None
    return out.decode().strip() or None
