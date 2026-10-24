"""Parsing of an ERF inputs file into a :class:`Deck`.

The syntax is AMReX ParmParse:

* ``key = value [value ...]``; values may be quoted strings;
* ``#`` starts a comment;
* a trailing backslash continues the line;
* ``FILE = other_file`` includes another file, resolved from the run directory
  (AMReX prepends ``AMREX_INPUTS_FILE_PREFIX`` when set); the checker falls back to
  the including file's directory and reports that it did;
* ``#if AMREX_SPACEDIM <op> N`` / ``#elif`` / ``#else`` / ``#endif`` and
  ``#ifdef AMREX_USE_GPU`` select a branch, as in AMReX_ParmParse.cpp;
* a Fortran namelist (``&name`` ... ``/``) is skipped;
* several ``key = value`` definitions may share one line;
* a key given twice keeps the last definition (ParmParse reads the table from
  the end), so duplicates are recorded and reported, not merged.
"""

import os
import re
import shlex
from typing import Dict, List, Optional, Tuple

_BOOL_TRUE = {"true", "1", "yes", "on", "t"}
_BOOL_FALSE = {"false", "0", "no", "off", "f"}


class DeckError(Exception):
    """A line that cannot be parsed."""


class Entry:
    """One ``key = values`` definition and where it came from."""

    __slots__ = ("key", "values", "path", "line")

    def __init__(self, key: str, values: List[str], path: str, line: int):
        self.key = key
        self.values = values
        self.path = path
        self.line = line

    @property
    def where(self) -> str:
        return "%s:%d" % (os.path.basename(self.path), self.line)


class Deck:
    """All definitions of an inputs file, with typed access.

    ``entries`` keeps every definition in file order (including the ones
    overridden by a later duplicate); ``table`` maps a key to its last entry,
    which is what ERF sees.
    """

    def __init__(self, path: str):
        self.path = os.path.abspath(path)
        self.directory = os.path.dirname(self.path)
        self.entries: List[Entry] = []
        self.table: Dict[str, Entry] = {}
        self.includes: List[str] = []
        # FILE = names found only beside the deck, not from the run directory
        self.include_notes: List[Tuple[str, Entry]] = []

    # -- construction -----------------------------------------------------
    def add(self, entry: Entry) -> None:
        self.entries.append(entry)
        self.table[entry.key] = entry

    # -- queries ----------------------------------------------------------
    def has(self, key: str) -> bool:
        return key in self.table

    def keys(self) -> List[str]:
        return list(self.table.keys())

    def keys_with_prefix(self, prefix: str) -> List[str]:
        return [k for k in self.table if k.startswith(prefix)]

    def where(self, key: str) -> str:
        e = self.table.get(key)
        return e.where if e else "?"

    def raw(self, key: str) -> Optional[List[str]]:
        e = self.table.get(key)
        return list(e.values) if e else None

    def get_str(self, key: str, default: Optional[str] = None) -> Optional[str]:
        v = self.raw(key)
        if not v:
            return default
        return v[0]

    def get_strs(self, key: str) -> List[str]:
        return self.raw(key) or []

    def get_int(self, key: str, default: Optional[int] = None) -> Optional[int]:
        v = self.raw(key)
        if not v:
            return default
        return _to_int(v[0], key)

    def get_ints(self, key: str) -> List[int]:
        return [_to_int(x, key) for x in self.get_strs(key)]

    def get_real(self, key: str, default: Optional[float] = None) -> Optional[float]:
        v = self.raw(key)
        if not v:
            return default
        return _to_real(v[0], key)

    def get_reals(self, key: str) -> List[float]:
        return [_to_real(x, key) for x in self.get_strs(key)]

    def get_bool(self, key: str, default: Optional[bool] = None) -> Optional[bool]:
        v = self.raw(key)
        if not v:
            return default
        s = v[0].lower()
        if s in _BOOL_TRUE:
            return True
        if s in _BOOL_FALSE:
            return False
        raise DeckError("%s = %s is not a boolean" % (key, v[0]))

    def duplicates(self) -> Dict[str, List[Entry]]:
        """Keys defined more than once, each with every definition."""
        seen: Dict[str, List[Entry]] = {}
        for e in self.entries:
            seen.setdefault(e.key, []).append(e)
        return {k: v for k, v in seen.items() if len(v) > 1}

    def resolve_path(self, name: str) -> Optional[str]:
        """Find a file named in the deck: relative to the deck, or the cwd."""
        if os.path.isabs(name):
            return name if os.path.exists(name) else None
        for base in (self.directory, os.getcwd()):
            cand = os.path.join(base, name)
            if os.path.exists(cand):
                return cand
        return None


def _to_int(s: str, key: str) -> int:
    try:
        return int(s)
    except ValueError:
        try:
            f = float(s)
        except ValueError:
            raise DeckError("%s = %s is not an integer" % (key, s))
        if f != int(f):
            raise DeckError("%s = %s is not an integer" % (key, s))
        return int(f)


def _to_real(s: str, key: str) -> float:
    # ParmParse accepts Fortran-style exponents such as 1.d-3
    t = re.sub(r"[dD]([+-]?\d)", r"e\1", s)
    try:
        return float(t)
    except ValueError:
        raise DeckError("%s = %s is not a number" % (key, s))



_IF_RE = re.compile(r"^\s*#\s*if\s+\(?\s*AMREX_SPACEDIM\s*(>|<|==|>=|<=|!=)\s*([1-3])\s*\)?\s*$")
_ELIF_RE = re.compile(r"^\s*#\s*elif\s+\(?\s*AMREX_SPACEDIM\s*(>|<|==|>=|<=|!=)\s*([1-3])\s*\)?\s*$")
_IFDEF_GPU_RE = re.compile(r"^\s*#\s*ifdef\s+AMREX_USE_GPU\s*$")
_IFNDEF_GPU_RE = re.compile(r"^\s*#\s*ifndef\s+AMREX_USE_GPU\s*$")
_ELSE_RE = re.compile(r"^\s*#\s*else\s*$")
_ENDIF_RE = re.compile(r"^\s*#\s*endif\s*$")
SPACEDIM = 3


def _spacedim_test(op: str, dim: int) -> bool:
    return {"<": SPACEDIM < dim, ">": SPACEDIM > dim, "==": SPACEDIM == dim,
            "<=": SPACEDIM <= dim, ">=": SPACEDIM >= dim, "!=": SPACEDIM != dim}[op]


def _logical_lines(path: str, gpu: bool) -> List[Tuple[int, str]]:
    """The definitions ParmParse would read from ``path``.

    Follows AMReX_ParmParse.cpp: ``#if AMREX_SPACEDIM <op> N`` / ``#elif`` /
    ``#else`` / ``#endif`` and ``#ifdef AMREX_USE_GPU`` select a branch (3-D,
    and the GPU flag given by the caller); a Fortran namelist (``&name`` up to a
    line ending in ``/``) is skipped; ``#`` starts a comment; a trailing
    backslash joins lines.  Returns (first line number, text) pairs.
    """
    out = []
    buf = ""
    start = 0
    valid: List[bool] = []      # current branch active, per nesting level
    has_true: List[bool] = []   # a branch of this #if already taken
    in_namelist = False
    with open(path, errors="replace") as fh:
        for n, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            m = _IF_RE.match(line)
            if m:
                r = _spacedim_test(m.group(1), int(m.group(2)))
                valid.append(r); has_true.append(r); continue
            if _IFDEF_GPU_RE.match(line):
                valid.append(gpu); has_true.append(gpu); continue
            if _IFNDEF_GPU_RE.match(line):
                valid.append(not gpu); has_true.append(not gpu); continue
            m = _ELIF_RE.match(line)
            if m and valid:
                r = (not has_true[-1]) and _spacedim_test(m.group(1), int(m.group(2)))
                valid[-1] = r; has_true[-1] = has_true[-1] or r; continue
            if _ELSE_RE.match(line) and valid:
                valid[-1] = not has_true[-1]; has_true[-1] = True; continue
            if _ENDIF_RE.match(line) and valid:
                valid.pop(); has_true.pop(); continue
            if valid and not all(valid):
                continue
            stripped = line.lstrip()
            if in_namelist:
                if stripped.rstrip().endswith("/"):
                    in_namelist = False
                continue
            if stripped.startswith("&"):
                in_namelist = not stripped.rstrip().endswith("/")
                continue
            line = _strip_comment(line)
            if not buf:
                start = n
            if line.rstrip().endswith("\\"):
                buf += line.rstrip()[:-1] + " "
                continue
            buf += line
            if buf.strip():
                out.append((start, buf))
            buf = ""
    if buf.strip():
        out.append((start, buf))
    return out


def _strip_comment(line: str) -> str:
    in_q = None
    for i, ch in enumerate(line):
        if in_q:
            if ch == in_q:
                in_q = None
        elif ch in ("'", '"'):
            in_q = ch
        elif ch == "#":
            return line[:i]
    return line


def _space_equals(text: str) -> str:
    """Put spaces around every '=' outside quotes, so shlex sees it as a token."""
    out = []
    in_q = None
    for ch in text:
        if in_q:
            if ch == in_q:
                in_q = None
            out.append(ch)
        elif ch in ("'", '"'):
            in_q = ch
            out.append(ch)
        elif ch == "=":
            out.append(" = ")
        else:
            out.append(ch)
    return "".join(out)


def _definitions(text: str, path: str, lineno: int) -> List[Tuple[str, List[str]]]:
    """Split one logical line into (key, values) pairs: ParmParse starts a new
    definition at every token that is followed by '=', so ``a = 1 b = 2 3`` holds two."""
    try:
        toks = shlex.split(_space_equals(text), comments=False, posix=True)
    except ValueError as exc:
        raise DeckError("%s:%d: %s" % (path, lineno, exc))
    defs: List[Tuple[str, List[str]]] = []
    i = 0
    while i < len(toks):
        if i + 1 < len(toks) and toks[i + 1] == "=":
            key = toks[i]
            if "=" in key or " " in key:
                raise DeckError("%s:%d: bad key in: %s" % (path, lineno, text.strip()))
            i += 2
            vals: List[str] = []
            while i < len(toks) and not (i + 1 < len(toks) and toks[i + 1] == "="):
                if toks[i] == "=":
                    raise DeckError("%s:%d: stray '=' in: %s" % (path, lineno, text.strip()))
                vals.append(toks[i])
                i += 1
            defs.append((key, vals))
        else:
            raise DeckError("%s:%d: no '=' after %r in: %s" % (path, lineno, toks[i], text.strip()))
    return defs


def resolve_include(name: str, including_file: str) -> Tuple[Optional[str], str]:
    """Where ParmParse opens ``FILE = name``: relative to the run directory, with the
    AMREX_INPUTS_FILE_PREFIX environment variable prepended when set.  As a fallback the
    checker also looks next to the including file and says so, since decks are usually
    checked from another directory than they are run from.  Returns (path, how)."""
    prefix = os.environ.get("AMREX_INPUTS_FILE_PREFIX", "")
    if prefix and not prefix.endswith("/"):
        prefix += "/"
    cand = prefix + name if prefix else name
    if os.path.isabs(cand):
        return (cand if os.path.exists(cand) else None), "absolute"
    if os.path.exists(cand):
        return os.path.abspath(cand), "cwd"
    beside = os.path.join(os.path.dirname(os.path.abspath(including_file)), name)
    if os.path.exists(beside):
        return beside, "deck_dir"
    return None, "missing"


def parse_deck(path: str, _deck: Optional[Deck] = None, _depth: int = 0, gpu: bool = False) -> Deck:
    """Parse ``path`` (and the files it includes with ``FILE =``)."""
    if _depth > 8:
        raise DeckError("FILE = includes nested more than 8 deep at %s" % path)
    deck = _deck or Deck(path)
    for lineno, text in _logical_lines(path, gpu):
        for key, values in _definitions(text, path, lineno):
            if key == "FILE":
                if not values:
                    raise DeckError("%s:%d: FILE = needs a file name" % (path, lineno))
                inc, how = resolve_include(values[0], path)
                if inc is None:
                    raise DeckError("%s:%d: FILE = %s not found from the current directory or beside the deck" % (path, lineno, values[0]))
                deck.includes.append(inc)
                if how == "deck_dir":
                    deck.include_notes.append((values[0], Entry("FILE", values, os.path.abspath(path), lineno)))
                parse_deck(inc, deck, _depth + 1, gpu)
                continue
            deck.add(Entry(key, values, os.path.abspath(path), lineno))
    return deck
