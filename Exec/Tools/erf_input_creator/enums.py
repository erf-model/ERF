"""Enum values read from the ERF source, so the commented alternatives follow the code.

``AMREX_ENUM(Name, a, b, ...)`` declarations under ``Source/DataStructs`` are
read once.  Outside a source tree the fallback list given by the caller is
used.
"""

import os
import re
from typing import Dict, List

_ENUM_RE = re.compile(r"AMREX_ENUM\(\s*(\w+)\s*,([^)]*)\)", re.S)
_ENUMS: Dict[str, List[str]] = {}
_LOADED = [False]


def source_root() -> str:
    """The ERF tree this package sits in (Exec/Tools/erf_input_creator)."""
    return os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))


def _load() -> None:
    _LOADED[0] = True
    d = os.path.join(source_root(), "Source", "DataStructs")
    if not os.path.isdir(d):
        return
    for fn in sorted(os.listdir(d)):
        if not fn.endswith(".H"):
            continue
        try:
            with open(os.path.join(d, fn), errors="replace") as fh:
                text = fh.read()
        except OSError:
            continue
        for m in _ENUM_RE.finditer(text):
            _ENUMS[m.group(1)] = [x.strip() for x in m.group(2).split(",") if x.strip()]


def enum_values(name: str, fallback: List[str]) -> List[str]:
    if not _LOADED[0]:
        _load()
    return list(_ENUMS.get(name, fallback))
