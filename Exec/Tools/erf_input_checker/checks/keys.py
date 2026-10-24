"""Key-level checks: duplicates, retired names, and misspelt keys."""

import difflib
import os
import re
from typing import Optional, Set

from ..checker import Checker, register
from ..findings import Severity

# Inputs that SolverChoice refuses (Docs: "Deprecated and Removed Inputs")
RETIRED_KEYS = {
    "erf.use_terrain": "set erf.terrain_type instead",
    "erf.use_moist_background": "set erf.init_type = MoistBaseState instead",
    "erf.moisture_type": "renamed: set erf.moisture_model instead",
    "erf.no_substepping": "set erf.substepping_type = None instead",
    "erf.use_wrf_height_grid": "renamed: set erf.wrfinput_zlevels_from_file instead",
    "erf.pbl_mynn_diffuse_moistvars": "removed",
}
RETIRED_VALUES = {
    "erf.init_type": {"real": "WRFInput", "ideal": "WRFInput"},
}
RENAMED_VALUES = {
    "erf.terrain_type": {"static": "StaticFittedMesh", "moving": "MovingFittedMesh"},
    "erf.pbl_type": {"shoc": "EAMXX_SHOC"},
}

# Prefixes whose keys are defined per problem or per sub-model, not in the
# main inputs table; they are not spell-checked.
FREE_PREFIXES = ("prob.", "amrex.", "fabarray.", "DistributionMapping.",
                 "ensemble.", "erf.fire.", "erf.dust.", "erf.conductors.",
                 "erf.bodies", "erf.turbine", "erf.ibseb", "erf.subvol.",
                 "erf.shoc.", "erf.sbm", "erf.rad", "erf.noahmp", "erf.lsm",
                 "erf.urban", "erf.canopy", "tiny_profiler.", "vismf.",
                 "erf.wps", "erf.cloud_chamber", "erf.seb", "erf.sdm",
                 "erf.inflow_log_law.", "erf.sampler", "erf.line_", "erf.plane_",
                 "erf.station_", "erf.moving_bodies", "erf.hindcast", "erf.openfast")


def _find_docs() -> Optional[str]:
    here = os.path.dirname(os.path.abspath(__file__))
    for up in range(2, 7):
        root = os.path.abspath(os.path.join(here, *([".."] * up)))
        cand = os.path.join(root, "Docs", "sphinx_doc", "Inputs.rst")
        if os.path.exists(cand):
            return cand
    return None


_DOC_KEY_RE = re.compile(r"\*\*([A-Za-z0-9_.<>]+)\*\*")
_DOC_CODE_RE = re.compile(r"``((?:erf|amr|geometry|xlo|xhi|ylo|yhi|zlo|zhi)\.[A-Za-z0-9_.]+)``")


_SRC_KEY_RE = re.compile(r'(?:query|queryAdd|get|getarr|queryarr|contains|countval|query_enum_case_insensitive|query_enum_sloppy|query_one_or_per_level(?:_enum_case_insensitive)?)\s*\(\s*(?:pp\w*\s*,\s*)?"([A-Za-z0-9_.]+)"')


_CACHE = {}


def source_key_names(docs_path: Optional[str] = None) -> Set[str]:
    """Bare key names read anywhere under Source/ (without their ParmParse prefix)."""
    path = docs_path or _find_docs()
    if ("src", path) in _CACHE:
        return _CACHE[("src", path)]
    names = _source_key_names(path)
    _CACHE[("src", path)] = names
    return names


def _source_key_names(path):
    if not path:
        return set()
    root = os.path.abspath(os.path.join(os.path.dirname(path), "..", ".."))
    src = os.path.join(root, "Source")
    names: Set[str] = set()
    if not os.path.isdir(src):
        return names
    for dirpath, _, files in os.walk(src):
        for fn in files:
            if not fn.endswith((".cpp", ".H", ".h")):
                continue
            try:
                with open(os.path.join(dirpath, fn), errors="replace") as fh:
                    for m in _SRC_KEY_RE.finditer(fh.read()):
                        names.add(m.group(1).lower())
            except OSError:
                continue
    return names


def documented_keys(docs_path: Optional[str] = None) -> Set[str]:
    """Every key named in bold or as inline code in Inputs.rst."""
    path = docs_path or _find_docs()
    if ("doc", path) in _CACHE:
        return _CACHE[("doc", path)]
    keys = _documented_keys(path)
    _CACHE[("doc", path)] = keys
    return keys


def _documented_keys(path):
    """Keys named in Inputs.rst and every other .rst page beside it."""
    keys: Set[str] = set()
    if not path or not os.path.exists(path):
        return keys
    docdir = os.path.dirname(path)
    pages = []
    for dirpath, _, files in os.walk(docdir):
        pages.extend(os.path.join(dirpath, f) for f in files if f.endswith(".rst"))
    for page in pages:
        with open(page, errors="replace") as fh:
            for line in fh:
                for m in _DOC_KEY_RE.finditer(line):
                    k = m.group(1)
                    if "." in k:
                        keys.add(k)
                for m in _DOC_CODE_RE.finditer(line):
                    keys.add(m.group(1))
    return keys


@register
class DuplicateKeys(Checker):
    """A key set twice: ParmParse keeps the last one."""

    name = "keys.duplicate"

    def run(self) -> None:
        for key, entries in self.deck.duplicates().items():
            first, last = entries[0], entries[-1]
            same_file = all(e.path == last.path for e in entries)
            same_value = all(e.values == last.values for e in entries)
            places = ", ".join(e.where for e in entries)
            if not same_file:
                # an override of an included file is the normal use of FILE =
                sev = Severity.INFO
                msg = "set in an included file and again in the deck (%s); the deck's value %r wins" % (
                    places, " ".join(last.values))
                adv = ""
            elif same_value:
                sev = Severity.ADVICE
                msg = "set %d times with the same value (%s)" % (len(entries), places)
                adv = "Remove the repeated line."
            else:
                sev = Severity.WARNING
                msg = "set %d times with different values (%s); ERF reads the last one, %r" % (
                    len(entries), places, " ".join(last.values))
                adv = "Keep one definition so the deck says what the run does."
            self.report.add(self._finding(sev, key, msg, adv, first.where))

    def _finding(self, sev, key, msg, adv, where):
        from ..findings import Finding
        return Finding(sev, self.name, key, msg, adv, where)


def _finding_at_impl(self, sev, key, where, msg, adv):
    from ..findings import Finding
    return Finding(sev, self.name, key, msg, adv, where)


Checker._finding_at = _finding_at_impl


@register
class RetiredKeys(Checker):
    """Inputs that ERF stops on, or silently maps to a new name."""

    name = "keys.retired"

    def run(self) -> None:
        for key, what in RETIRED_KEYS.items():
            if self.deck.has(key):
                self.error(key, "this input was removed from ERF; the run stops at start-up",
                           what + ".")
        for key, table in RETIRED_VALUES.items():
            v = self.str_(key)
            if v and v.lower() in table:
                self.error(key, "the value %r is no longer accepted" % v,
                           "Use %s = %s." % (key, table[v.lower()]))
        for key, table in RENAMED_VALUES.items():
            v = self.str_(key)
            if v and v.lower() in table:
                self.warn(key, "the value %r is a legacy spelling that ERF still maps" % v,
                          "Use %s = %s." % (key, table[v.lower()]))


@register
class UnknownKeys(Checker):
    """Keys with an erf./amr./geometry./<face>. prefix that the docs do not list.

    ParmParse never complains about a key nobody reads, so a misspelt key is
    silently ignored and the default is used instead.
    """

    name = "keys.unknown"
    docs_path: Optional[str] = None   # set by the CLI with --docs

    def run(self) -> None:
        known = documented_keys(self.docs_path)
        if not known:
            self.info("docs", "Docs/sphinx_doc/Inputs.rst not found; misspelt keys are not checked",
                      "Run from the ERF repository or pass --docs <path to Inputs.rst>.")
            return
        known_l = {k.lower(): k for k in known}
        read_names = source_key_names(self.docs_path)
        for name, entry in self.deck.include_notes:
            self.report.add(self._finding_at(Severity.WARNING, "FILE", entry.where,
                "FILE = %s was found beside the deck but not from the current directory; ParmParse opens it relative to the directory the run starts in (or AMREX_INPUTS_FILE_PREFIX)" % name,
                "Run ERF from the deck's directory, or give the path relative to the run directory."))
        indicators = [s.lower() for s in self.deck.get_strs("erf.refinement_indicators")]
        faces = ("xlo", "xhi", "ylo", "yhi", "zlo", "zhi")
        for key in self.deck.keys():
            kl = key.lower()
            if key in known:
                continue
            if kl in known_l:
                # ParmParse keys are case-sensitive: the key as written is never read
                self.warn(key, "not read: ParmParse keys are case-sensitive and the documented spelling is %s" % known_l[kl],
                          "Rename to %s." % known_l[kl])
                continue
            if kl.startswith(tuple(p.lower() for p in FREE_PREFIXES)):
                continue
            parts = kl.split(".")
            if parts[0] == "erf" and len(parts) > 2 and parts[1] in faces:
                parts = parts[1:]          # erf.<face>.x is read as <face>.x
                kl = ".".join(parts)
                if kl in known_l:
                    continue
            if parts[0] in faces:
                probe = "<face>." + ".".join(parts[1:])
                if probe in known_l or kl in known_l:
                    continue
            if parts[0] == "erf" and len(parts) >= 3 and parts[1] in indicators:
                probe = "erf.<indicator>." + ".".join(parts[2:])
                if probe in known_l:
                    continue
            if parts[0] not in ("erf", "amr", "geometry") and parts[0] not in faces:
                continue
            # per-level suffixes and numbered keys such as plot_int_1 are documented
            # with the number; try the generic pattern before giving up
            stem = re.sub(r"_\d+$", "", kl)
            if stem in known_l:
                continue
            # a documented key under another prefix or with a letter wrong comes first,
            # so that amr.fixed_dt is reported as erf.fixed_dt and not excused by the
            # bare-name fallback below
            close = difflib.get_close_matches(kl, list(known_l.keys()), n=1, cutoff=0.8)
            if close:
                self.warn(key, "not a documented input; ParmParse ignores keys nobody reads",
                          "Did you mean %s?" % known_l[close[0]])
                continue
            last = parts[-1]
            # the code reads it but the docs do not list it: a documentation gap, not a typo
            if last in read_names or ".".join(parts[1:]) in read_names:
                self.info(key, "read by the code (Source/) but not listed in Docs/sphinx_doc/Inputs.rst",
                          "Worth a documentation row.")
            else:
                self.advise(key, "not listed in Docs/sphinx_doc/Inputs.rst",
                            "Check the spelling, or that the key belongs to a model built into this executable.")
