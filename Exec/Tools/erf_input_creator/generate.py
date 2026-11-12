"""Assemble the deck of one case, write it with its dummy files, and summarise it."""

import hashlib
import json
import os
import shlex
from typing import List, Tuple

from .auxfiles import AuxFile, aux_files
from .constraints import RESERVED_NAMES, check
from .deckwriter import DeckText
from .options import SETUP_DEFAULTS
from .recipes import grid, models, output, physics, terrain

REPLACE_LIST = "FILES_TO_REPLACE.txt"


class CaseError(Exception):
    """The options cannot give a deck that ERF would start."""


def build(s, command: str = "") -> Tuple[DeckText, List[AuxFile]]:
    errors = check(s)
    if errors:
        raise CaseError("\n".join(errors))
    try:
        return _assemble(s, command), aux_files(s)
    except ValueError as exc:            # a value the deck cannot hold
        raise CaseError(str(exc))


def _assemble(s, command: str) -> DeckText:
    deck = DeckText()
    deck.header = [
        "ERF inputs written by Exec/Tools/erfcreateinput.py",
        "case: %s, setup: %s (%s)" % (s.case, s.setup, SETUP_DEFAULTS[s.setup]["about"]),
        "command: " + (command or "(not recorded)"),
        "",
        "Active lines are the chosen configuration. A line starting with '#key =' is an alternative",
        "or an optional key: remove the '#' to use it, and comment the line it replaces.",
        "Lines marked 'USER:' hold values the tool cannot know; review each one.",
    ]
    grid.problem(deck, s)
    grid.geometry(deck, s)
    terrain.terrain(deck, s)
    grid.boundaries(deck, s)
    grid.time(deck, s)
    grid.dycore(deck, s)
    grid.advection(deck, s)
    physics.turbulence(deck, s)
    physics.surface_layer(deck, s)
    physics.forcing(deck, s)
    models.moisture(deck, s)
    models.radiation(deck, s)
    models.land_surface(deck, s)
    models.canopy(deck, s)
    models.windfarm(deck, s)
    physics.initial(deck, s)
    output.output(deck, s)
    return deck


def summary(s, deck: DeckText, files: List[AuxFile]) -> List[str]:
    out = ["Case written to %s" % s.outdir,
           "  inputs file : %s (%d active keys)" % (s.name, len(deck.active_keys()))]
    dummies = [f for f in files if f.generated]
    needed = [f for f in files if not f.generated]
    if dummies:
        out.append("")
        out.append("Dummy files to replace (%d). Each one is valid, so the run starts, but its" % len(dummies))
        out.append("content is a placeholder:")
        for f in dummies:
            out.append("")
            out.append("  %s   (named by %s)" % (f.name, f.key))
            out += ["      " + line for line in f.describe]
    if needed:
        out.append("")
        out.append("Files you must supply (%d). The tool cannot generate them, and the run stops" % len(needed))
        out.append("until they are in the case directory:")
        for f in needed:
            out.append("")
            out.append("  %s   (named by %s)" % (f.name, f.key))
            out += ["      " + line for line in f.describe]
    if not files:
        out.append("  no auxiliary files are needed for this configuration")
    if deck.user_items:
        out.append("")
        out.append("Values to review (marked USER: in the inputs file):")
        for key, value, comment in deck.user_items:
            out.append("  %-34s = %-24s %s" % (key, value, comment))
    stale = leftovers(s, files)
    if stale:
        out.append("")
        out.append("Left from an earlier configuration and not named by this deck: %s" % ", ".join(stale))
    out.append("")
    out.append("To write the same case again:  erfcreateinput.py --case-file %s --outdir <dir>"
               % os.path.join(s.outdir, CASE_FILE))
    out.append("")
    out.append("Next steps:")
    out.append("  python3 Exec/Tools/check_erf_inputs.py %s" % os.path.join(s.outdir, s.name))
    out.append("  cd %s && <path to>/erf_exec %s" % (s.outdir, s.name))
    return out


CASE_FILE = "case.json"


def write(s, deck: DeckText, files: List[AuxFile], case: dict = None) -> List[str]:
    """Write the deck, the dummy files, the replace list and the case file; returns the paths written.

    Everything is checked before the first file is written, so a refused case leaves nothing behind.
    """
    targets = [(s.name, deck.render())] + [(f.name, f.text) for f in files if f.generated]
    if files:
        lines = []
        for f in files:
            kind = "dummy, replace" if f.generated else "NOT GENERATED, supply"
            lines += ["%s   (%s; named by %s)" % (f.name, kind, f.key)] + ["    " + d for d in f.describe]
            if f.generated:
                # read by check_erf_inputs.py: a file with this fingerprint is still the placeholder
                lines.append("    sha256: " + hashlib.sha256(f.text.encode("utf-8")).hexdigest())
            lines.append("")
        targets.append((REPLACE_LIST, "\n".join(lines)))
    if case is not None:
        targets.append((CASE_FILE, json.dumps(case, indent=2, sort_keys=True) + "\n"))
    if os.path.exists(s.outdir) and not os.path.isdir(s.outdir):
        raise CaseError("--outdir %s exists and is not a directory" % s.outdir)
    blocked = [n for n, _ in targets if os.path.isdir(os.path.join(s.outdir, n))]
    if blocked:
        raise CaseError("%s holds a directory named %s where a file must be written" % (s.outdir, ", ".join(blocked)))
    existing = [n for n, _ in targets if os.path.exists(os.path.join(s.outdir, n))]
    if existing and not s.force:
        raise CaseError("%s already holds %s; choose another --outdir or pass --force to overwrite"
                        % (s.outdir, ", ".join(existing)))
    paths = []
    try:
        os.makedirs(s.outdir, exist_ok=True)
        for name, text in targets:
            path = os.path.join(s.outdir, name)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(text)
            paths.append(path)
    except OSError as exc:
        raise CaseError("cannot write the case to %s: %s" % (s.outdir, exc))
    return paths


def leftovers(s, files: List[AuxFile]) -> List[str]:
    """Files of an earlier configuration in the case directory that this deck does not name."""
    mine = set(f.name for f in files if f.generated)
    if files:
        mine.add(REPLACE_LIST)
    try:
        present = set(os.listdir(s.outdir))
    except OSError:
        return []
    return sorted(n for n in RESERVED_NAMES if n in present and n not in mine and n != CASE_FILE)


def command_line(argv: List[str]) -> str:
    return "erfcreateinput.py " + " ".join(shlex.quote(a) for a in argv)
