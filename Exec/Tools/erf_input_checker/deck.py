"""Parsing of an ERF inputs file into a :class:`Deck`.

The syntax is AMReX ParmParse:

* ``key = value [value ...]``; values may be quoted strings;
* ``#`` starts a comment;
* a trailing backslash continues the line;
* ``FILE = other_file`` includes another file (relative to the including
  file), and the keys it sets can be overridden by later lines;
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


def _logical_lines(path: str) -> List[Tuple[int, str]]:
    """Join continuation lines and strip comments; keep the first line number."""
    out = []
    buf = ""
    start = 0
    with open(path, errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            # a comment ends the line, but not inside quotes
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


def parse_deck(path: str, _deck: Optional[Deck] = None, _depth: int = 0) -> Deck:
    """Parse ``path`` (and the files it includes with ``FILE =``)."""
    if _depth > 8:
        raise DeckError("FILE = includes nested more than 8 deep at %s" % path)
    deck = _deck or Deck(path)
    for lineno, text in _logical_lines(path):
        if "=" not in text:
            raise DeckError("%s:%d: no '=' in: %s" % (path, lineno, text.strip()))
        key, val = text.split("=", 1)
        key = key.strip()
        if not key or " " in key or "\t" in key:
            raise DeckError("%s:%d: bad key in: %s" % (path, lineno, text.strip()))
        try:
            values = shlex.split(val, comments=False, posix=True)
        except ValueError as exc:
            raise DeckError("%s:%d: %s" % (path, lineno, exc))
        if key == "FILE":
            if not values:
                raise DeckError("%s:%d: FILE = needs a file name" % (path, lineno))
            inc = values[0]
            if not os.path.isabs(inc):
                inc = os.path.join(os.path.dirname(os.path.abspath(path)), inc)
            if not os.path.exists(inc):
                raise DeckError("%s:%d: FILE = %s not found" % (path, lineno, values[0]))
            deck.includes.append(inc)
            parse_deck(inc, deck, _depth + 1)
            continue
        deck.add(Entry(key, values, os.path.abspath(path), lineno))
    return deck
