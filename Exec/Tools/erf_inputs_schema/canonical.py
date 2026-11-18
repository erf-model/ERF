"""Type codes and the canonical string form of a value, as the schema records them.

The canonical form is what an idempotent run identity hashes: the same value spelled two ways gives one string.
A word that is not a plain literal of its type raises ValueError (no expressions, no 123. for an int).
"""
import re

SCALARS = ("int", "real", "bool", "string", "enum")
TYPES = SCALARS + tuple(t + "_list" for t in SCALARS) + ("unknown", "unknown_list")
INT = re.compile(r"[+-]?[0-9]+")
REAL = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)(?:[eE][+-]?[0-9]+)?")
# A Fortran exponent (1D0, 1.d-3) is not a ParmParse real: isT leaves "D0" unread, and the parser lexes [Ee] only.
FORTRAN_REAL = re.compile(r"[+-]?(?:[0-9]+\.?[0-9]*|\.[0-9]+)[dD][+-]?[0-9]+")
BOOLS = {"true": True, "t": True, "1": True, "false": False, "f": False, "0": False}

RULES = {
    "int": "decimal digits, sign only if negative (+03 -> 3)",
    "real": "Python's shortest round-trip repr of the double (1, 1.0, 1e0 -> 1.0; 1D0 is not a real)",
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
    if kind == "real" and FORTRAN_REAL.fullmatch(word):
        raise ValueError(f"not a real literal: {word!r} (ParmParse takes an e or E exponent, not d or D)")
    if kind == "bool" and word.lower() in BOOLS:
        return "true" if BOOLS[word.lower()] else "false"
    if kind == "enum":
        return _enum(word, fact)
    if kind in ("string", "unknown"):
        return word
    raise ValueError(f"not a {kind} literal: {word!r}")


def _enum(word, fact):
    values = fact.get("values", [])
    if fact.get("ignores"):                          # query_enum_sloppy drops these characters from the word first
        word = "".join(c for c in word if c not in fact["ignores"])
    if fact.get("case_insensitive"):
        match = [v for v in values if v.lower() == word.lower()]
    else:
        match = [v for v in values if v == word]
    if not match:
        raise ValueError(f"not one of the enum's names: {word!r}")
    return match[0]
