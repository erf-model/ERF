"""Type codes and the canonical string form of a value, as the schema records them.

The canonical form is what an idempotent run identity hashes: the same value spelled two ways gives one string.
A word that is not a plain literal of its type raises ValueError (no expressions, no 123. for an int).
"""
import re

SCALARS = ("int", "real", "bool", "string", "enum")
TYPES = SCALARS + tuple(t + "_list" for t in SCALARS) + ("unknown",)
INT = re.compile(r"[+-]?[0-9]+")
REAL = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
BOOLS = {"true": True, "t": True, "1": True, "false": False, "f": False, "0": False}

RULES = {
    "int": "decimal digits, sign only if negative (+03 -> 3)",
    "real": "Python's shortest round-trip repr of the double (1, 1.0, 1e0 -> 1.0)",
    "bool": "true or false (t, 1, false, f, 0 accepted, any case)",
    "enum": "the enum name in its declared case; case-insensitive only where ERF reads it so",
    "string": "unchanged",
}


def canonical_value(kind, word, fact):
    """The canonical string of one word of a scalar type."""
    if kind == "int" and INT.fullmatch(word):
        return str(int(word))
    if kind == "real" and REAL.fullmatch(word):
        return repr(float(word))
    if kind == "bool" and word.lower() in BOOLS:
        return "true" if BOOLS[word.lower()] else "false"
    if kind == "enum":
        return _enum(word, fact)
    if kind in ("string", "unknown"):
        return word
    raise ValueError(f"not a {kind} literal: {word!r}")


def _enum(word, fact):
    values = fact.get("values", [])
    if fact.get("case_insensitive"):
        match = [v for v in values if v.lower() == word.lower()]
    else:
        match = [v for v in values if v == word]
    if not match:
        raise ValueError(f"not one of the enum's names: {word!r}")
    return match[0]
