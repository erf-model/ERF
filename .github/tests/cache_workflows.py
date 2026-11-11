#!/usr/bin/env python3
"""Check ERF's cache wiring without running workflows, builds, or registry commands.

This models only the expressions used in cache settings and matrix axes. Unknown
forms fail explicitly; actionlint remains responsible for general workflow syntax.
"""

from itertools import product
import json
from pathlib import Path
import re
import sys
import tempfile

import yaml

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / ".github/actions/ccache-trunk"))
from cache import read_config

# Keep existing jobs covered even if someone accidentally removes both calls.
REQUIRED_JOBS = {
    "ci.yml": "Build-And-Test", "gcc.yml": "library",
    "gcc-rrtmgp.yml": "library", "cuda-ci.yml": "cuda-build",
    "hip.yml": "Build-And-Test-HIP", "sycl.yml": "Build-And-Test-SYCL",
    "macos.yml": "MacOS-Clang", "windows.yml": "windows",
}
SAVE_GUARD = ("always() && github.event_name != 'pull_request' && "
              "github.ref_name == github.event.repository.default_branch")
KEY_SETUP = r'''printf 'CCACHE_TRUNK_KEY=%s\n' "$TRUNK_KEY" >> "$GITHUB_ENV"'''
EXPRESSION = re.compile(r"\$\{\{(.*?)\}\}", re.DOTALL)
PR_AXIS = re.compile(
    r"fromJSON\(\s*github\.event_name\s*==\s*'pull_request'\s*"
    r"&&\s*'([^']+)'\s*\|\|\s*'([^']+)'\s*\)")
MATRIX_CHOICE = re.compile(
    r"matrix\.(\w+)\s*==\s*'([^']*)'\s*&&\s*'([^']+)'\s*\|\|\s*'([^']*)'")
RUNNERS = {
    "ubuntu-22.04": ("Linux", ("X64",)),
    "ubuntu-24.04": ("Linux", ("X64",)),
    # macos-latest can resolve differently across repositories / image rollouts.
    "macos-latest": ("macOS", ("X64", "ARM64")),
    "windows-latest": ("Windows", ("X64",)),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def expression_text(text):
    match = EXPRESSION.fullmatch(text.strip())
    return match[1].strip() if match else text.strip()


def condition(text):
    """Ignore formatting, but preserve whitespace inside quoted literals."""
    return "".join(re.findall(r"'[^']*'|\S", expression_text(text)))


def load_yaml(path):
    # BaseLoader preserves GitHub strings such as on, ON, OFF, and version numbers.
    return yaml.load(path.read_text(), Loader=yaml.BaseLoader)


def expand_matrix(matrix, event):
    """Cartesian axes, then exclude, then GitHub's ordered include merges."""
    require(isinstance(matrix, dict), "matrix must be a mapping")
    axes = {}
    for key, values in matrix.items():
        if key in ("include", "exclude"):
            continue
        if isinstance(values, str):
            match = PR_AXIS.fullmatch(expression_text(values))
            require(match, f"unsupported matrix expression for {key}: {values}")
            values = json.loads(match[1 if event == "pull_request" else 2])
        require(isinstance(values, list) and values
                and all(isinstance(v, str) for v in values),
                f"matrix axis {key} must be a nonempty list of strings")
        axes[key] = values
    excludes = matrix.get("exclude", [])
    includes = matrix.get("include", [])
    for entries in (excludes, includes):
        require(isinstance(entries, list) and all(isinstance(e, dict) for e in entries),
                "matrix include/exclude must be lists of mappings")
        require(all(isinstance(v, str) for e in entries for v in e.values()),
                "matrix include/exclude values must be strings")
    originals = [dict(zip(axes, values)) for values in product(*axes.values())] if axes else []
    originals = [row for row in originals
                 if not any(all(row.get(k) == v for k, v in e.items()) for e in excludes)]
    rows = [dict(row) for row in originals]
    extras = []
    for entry in includes:
        compatible = [i for i, row in enumerate(originals)
                      if all(k not in row or row[k] == v for k, v in entry.items())]
        if compatible:
            for i in compatible:
                rows[i].update(entry)
        else:
            extras.append(dict(entry))
    rows.extend(extras)
    require(rows, "matrix produces no configurations")
    return rows


def render(text, matrix, runner, env, resolving=()):
    """Resolve cache input references only; never evaluate Python or shell."""
    require(isinstance(text, str), f"expected a scalar cache setting: {text!r}")

    def replace(match):
        expr = match[1].strip()
        choice = MATRIX_CHOICE.fullmatch(expr)
        if choice:
            key, expected, yes, no = choice.groups()
            require(key in matrix, f"missing matrix.{key}")
            return yes if matrix[key].lower() == expected.lower() else no
        ref = re.fullmatch(r"(matrix|runner|env)\.(\w+)", expr)
        require(ref, f"unsupported cache expression: {expr}")
        context, key = ref.groups()
        values = {"matrix": matrix, "runner": runner, "env": env}[context]
        require(key in values, f"missing {context}.{key}")
        if context == "env":
            require(key not in resolving, f"cyclic env reference: {key}")
            return render(values[key], matrix, runner, env, (*resolving, key))
        return values[key]

    result = EXPRESSION.sub(replace, text)
    require("${{" not in result, f"unresolved cache expression: {result}")
    return result


def is_cache(step):
    uses = step.get("uses", "")
    return uses.startswith("./") and uses.endswith("/.github/actions/ccache-trunk")


def cache_steps(workflow, job):
    """Require the standard restore/save pair and trusted publishing policy."""
    steps = job.get("steps", [])
    calls = [(i, s) for i, s in enumerate(steps) if is_cache(s)]
    require([s.get("with", {}).get("mode") for _, s in calls] == ["restore", "save"],
            "expected exactly one restore followed by one save")
    restore, save = (s for _, s in calls)
    require(not restore.get("if"), "restore must be unconditional")
    require(condition(save.get("if", "")) == condition(SAVE_GUARD),
            f"save must use the trusted publication guard: {SAVE_GUARD}")
    permissions = job.get("permissions", workflow.get("permissions", {}))
    require(isinstance(permissions, dict) and permissions.get("packages") == "write",
            "publishing job needs packages: write")
    token = save["with"].get("token", "")
    require(EXPRESSION.fullmatch(token.strip())
            and condition(token) == "secrets.GITHUB_TOKEN",
            "save must use secrets.GITHUB_TOKEN")
    return calls


def rendered_settings(workflow, job, row, runner, defaults, temporary):
    """Walk steps to honor step env and the native runner's key setup."""
    env = {**workflow.get("env", {}), **job.get("env", {})}
    settings = []
    legacy = []
    for index, step in enumerate(job["steps"]):
        local_env = {**env, **step.get("env", {})}
        script = step.get("run", "")
        if "CCACHE_TRUNK_KEY" in script:
            require(script.strip() == KEY_SETUP and step.get("shell") == "bash"
                    and not step.get("if"), "unsupported CCACHE_TRUNK_KEY setup")
            env["CCACHE_TRUNK_KEY"] = render("${{ env.TRUNK_KEY }}", row, runner, local_env)
        if is_cache(step):
            inputs = {**defaults, **step.get("with", {})}
            key, path, backend = (render(inputs[k], row, runner, local_env)
                                  for k in ("key", "path", "transport"))
            require(path.strip(), "cache path is empty")
            # Use runtime normalization/platform validation. Its filesystem input
            # is synthetic: foreign runner paths are compared as strings below.
            config = read_config({
                "CACHE_MODE": "restore", "CACHE_KEY": key,
                "CACHE_TRANSPORT": backend, "CACHE_PATH": str(temporary / "ccache"),
                "RUNNER_OS": runner["os"], "RUNNER_ARCH": runner["arch"],
                "GITHUB_REPOSITORY_OWNER": "contract-check", "RUNNER_TEMP": str(temporary),
                "GITHUB_OUTPUT": str(temporary / "output"),
            })
            settings.append((key, path, config.transport, config.ref))
        elif step.get("uses", "").startswith("actions/cache"):
            path = render(step.get("with", {}).get("path", ""), row, runner, local_env)
            legacy.append((index, step, path))
    require(settings[0] == settings[1], "restore/save key, path, or transport mismatch")
    restore_index, save_index = (i for i, s in enumerate(job["steps"]) if is_cache(s))
    restore_id = job["steps"][restore_index].get("id", "")
    for index, step, path in legacy:
        require(step["uses"].startswith("actions/cache/restore@"),
                "legacy Actions cache must remain read-only")
        require(restore_index < index < save_index, "legacy fallback must follow restore")
        require(path == settings[0][1], "legacy fallback path differs from GHCR restore")
        require(restore_id and condition(step.get("if", ""))
                == condition(f"steps.{restore_id}.outputs.hit != 'true'"),
                "legacy fallback must check the GHCR restore step's hit output")
    return settings[0]


def check_hip_version(job, row, key):
    """Keep HIP's selected compiler version connected to installation and naming.

    Require the supported single installer command, without interpreting shell.
    Check the effective key too, so a step override cannot hide a wrong version.
    """
    version = row.get("rocm_version", "")
    require(re.fullmatch(r"[0-9]+(?:\.[0-9]+)+", version),
            "HIP matrix must define a numeric rocm_version")
    reference = r"\$\{\{\s*matrix\.rocm_version\s*\}\}"
    template = job.get("env", {}).get("CCACHE_TRUNK_KEY", "")
    require(re.search(r"-rocm-" + reference + "-", template),
            "HIP cache key must use matrix.rocm_version")
    require(f"-rocm-{version}-" in key, "HIP cache key differs from selected ROCm version")
    installer = "Submodules/AMReX/.github/workflows/dependencies/dependencies_hip.sh"
    commands = [s["run"].strip() for s in job["steps"]
                if "dependencies_hip.sh" in s.get("run", "")]
    require(len(commands) == 1
            and re.fullmatch(re.escape(installer) + r'\s+"' + reference + r'"', commands[0]),
            "HIP installer must take the quoted matrix.rocm_version argument")


def check_workflows(workflows, defaults):
    """Validate known and newly added cache jobs; return counts and unique tags."""
    for name, job_id in REQUIRED_JOBS.items():
        require(name in workflows and job_id in workflows[name].get("jobs", {}),
                f"missing required cache job: {name}:{job_id}")
    seen = {}
    counts = {}
    with tempfile.TemporaryDirectory(prefix="erf-workflow-contract-") as tmp:
        for name, workflow in sorted(workflows.items()):
            # Diagnostic smoke intentionally restores/saves multiple times.
            if name == "cache-smoke.yml":
                continue
            for job_id, job in workflow.get("jobs", {}).items():
                if (REQUIRED_JOBS.get(name) != job_id
                        and not any(is_cache(s) for s in job.get("steps", []))):
                    continue
                label = f"{name}:{job_id}"
                try:
                    cache_steps(workflow, job)
                    events = workflow.get("on", {})
                    events = [events] if isinstance(events, str) else events
                    require(events and set(events) <= {"push", "pull_request", "workflow_dispatch"},
                            "unsupported cache workflow trigger; review publication trust first")
                    for event in events:
                        matrix = job.get("strategy", {}).get("matrix")
                        rows = expand_matrix(matrix, event) if matrix is not None else [{}]
                        counts[event] = counts.get(event, 0) + len(rows)
                        event_tags = set()
                        for row in rows:
                            label = f"{name}:{job_id} [{event}, {row}]"
                            runner_name = render(job["runs-on"], row, {}, {})
                            require(runner_name in RUNNERS, f"unsupported runner: {runner_name}")
                            os_name, arches = RUNNERS[runner_name]
                            for arch in arches:
                                runner = {"os": os_name, "arch": arch, "temp": "/runner-temp"}
                                settings = rendered_settings(workflow, job, row, runner, defaults, Path(tmp))
                                if name == "hip.yml" and job_id == REQUIRED_JOBS[name]:
                                    check_hip_version(job, row, settings[0])
                                tag = settings[3].rsplit(":", 1)[1]
                                owner = (name, job_id, tuple(sorted(row.items())), arch)
                                require(tag not in event_tags, f"cache tag collision: {tag} ({arch})")
                                event_tags.add(tag)
                                if tag in seen:
                                    previous_owner, previous_settings = seen[tag]
                                    require(owner == previous_owner and settings == previous_settings,
                                            f"cache tag collision: {tag}; also used by {previous_owner}")
                                seen[tag] = (owner, settings)
                except (ValueError, KeyError, TypeError) as exc:
                    raise ValueError(f"{label}: {exc}") from exc
    return counts, sorted(seen)


def repository_inputs(root=ROOT):
    workflows = {p.name: load_yaml(p) for p in sorted((root / ".github/workflows").iterdir())
                 if p.suffix in (".yml", ".yaml")}
    action = load_yaml(root / ".github/actions/ccache-trunk/action.yml")
    defaults = {k: v["default"] for k, v in action["inputs"].items() if "default" in v}
    return workflows, defaults


def main():
    try:
        counts, tags = check_workflows(*repository_inputs())
    except (ValueError, yaml.YAMLError) as exc:
        print(f"Cache workflow contract failed: {exc}", file=sys.stderr)
        return 1
    print("Cache workflow contracts passed: " + ", ".join(
        f"{count} {event} matrix variants" for event, count in sorted(counts.items())))
    print(f"{len(tags)} unique tags across supported runner architectures (including both macOS architectures).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
