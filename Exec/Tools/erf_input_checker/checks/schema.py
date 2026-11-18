"""Checks against a generated per-build inputs schema.

Runs only when ``check_erf_inputs.py --schema FILE`` gives one (written by ``Exec/Tools/build_erf_schema.py``). Each
deck key is looked up in the schema, exactly or by a family pattern such as ``erf.<station>.lat``:

* not in the schema: WARNING, since no query in this build's ERF or AMReX reads it; INFO instead when the schema's
  ``unresolved`` list has a read of that name through a prefix the generator could not resolve;
* ``removed``: ERROR, since ERF stops at start-up when it is present (the schema names the line), unless the keys
  check already reports it from RETIRED_KEYS;
* a string or enum word that is not a value of the key's type, or too few per-level values: ERROR. A key of type
  ``unknown`` takes any number of words;
* a bool, int or real word that is not a plain literal (``2*0.5``, ``${x}``, ``300.;``, ``1e9`` for an int):
  WARNING, since ParmParse hands it to its parser at run time (``pp_parser``); none for a key read with the parser
  (``queryWithParser``, ``queryAsDouble``), which takes an expression by design;
* more words than a scalar key reads (``erf.check_int = 10 steps``): WARNING, since ParmParse reads the first word
  (``ival = FIRST``) and ignores the rest; the same for per-level values beyond the number of levels;
* the key, or this value of it, needs a build feature the schema says is off: ERROR.
"""
import os
import re
import sys

from ..checker import Checker, register
from .keys import RETIRED_KEYS

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from erf_inputs_schema.canonical import FORTRAN_REAL, canonical_value  # noqa: E402

# <station>: ERF_StationSampler.cpp check_station_name (a letter or '_' first; letters, digits, '_', '-', '.').
PLACEHOLDERS = {"<face>": "(?:xlo|xhi|ylo|yhi|zlo|zhi)", "<station>": "[A-Za-z_][A-Za-z0-9_.-]*",
                "<indicator>": "[A-Za-z_][A-Za-z0-9_]*", "<n>": "[0-9]+"}


def _pattern_rx(pattern):
    out = re.escape(pattern)
    for token, rx in PLACEHOLDERS.items():
        out = out.replace(re.escape(token), rx)
    return re.compile(out)


def lookup(key, schema):
    """The key's fact: an exact key first, then a family pattern."""
    if key in schema["keys"]:
        return schema["keys"][key]
    for pattern, fact in schema.get("patterns", {}).items():
        if _pattern_rx(pattern).fullmatch(key):
            return fact
    return None


def unresolved_reads(key, schema):
    """The schema's unresolved reads (``file:line var.name``) whose name is this key's last component."""
    name = key.rsplit(".", 1)[-1]
    return [u for u in schema.get("unresolved", []) if u.rsplit(".", 1)[-1] == name]


def cross_check(schema, retired=RETIRED_KEYS):
    """(curated_only, source_only): retired keys the schema's source does not mark removed, and removed keys (with
    their evidence line) the curated list lacks."""
    removed = {k: f["removed_evidence"] for k, f in schema["keys"].items() if f["type"] == "removed"}
    curated_only = sorted(k for k in retired if k not in removed)
    source_only = sorted((k, ev) for k, ev in removed.items() if k not in retired)
    return curated_only, source_only


EVALUATED = ("bool", "int", "real")               # ParmParse's parser evaluates a non-literal word (pp_parser)


def _scalar(kind):
    return kind[:-5] if kind.endswith("_list") else kind


def _is_literal(scalar, word, fact):
    try:
        canonical_value(scalar, word, fact)
    except ValueError:
        return False
    return True


def read_words(fact, words, nlev=1):
    """The words ERF reads: all of a list, the first nlev of a per-level key, the first of a scalar."""
    kind = fact["type"]
    if kind.endswith("_list") or kind == "unknown" or fact.get("line") or fact.get("table"):
        return list(words)                           # queryline reads the whole line, querytable a {{..} {..}}
    if fact.get("per_level"):
        return list(words[:max(nlev, 1)])
    return list(words[:1])


def ignored_words(fact, words, nlev=1):
    """The words after the ones ERF reads (ParmParse ival = FIRST; query_one_or_per_level reads nlev)."""
    return list(words[len(read_words(fact, words, nlev)):])


WORD = re.compile(r"[A-Za-z_]+")                    # a bare name: the parser has no symbol for it, so an error
PLACEHOLDER = re.compile(r"@\w+@")                  # CMake configure_file: erf.moisture_model = @BULK_MODEL@


def _evaluated(fact, words):
    scalar = _scalar(fact["type"])
    if scalar not in EVALUATED:
        return []
    return [w for w in words if not _is_literal(scalar, w, fact) and not WORD.fullmatch(w)
            and not FORTRAN_REAL.fullmatch(w)]                     # the parser cannot read 1D0 either


def expression_words(fact, words):
    """The words ParmParse would hand to its parser rather than read as a literal (bool, int, real keys): anything
    that is not a plain literal, except a bare word such as ``fast``, which the parser cannot evaluate either. None
    for a parser read (queryWithParser, queryAsDouble), where an expression is what the key takes."""
    return [] if fact.get("parser") else _evaluated(fact, words)


def value_problem(fact, words, nlev):
    """Why these words are not a value of the key, or None. Expressions (expression_words) and ignored words
    (ignored_words) are not problems."""
    kind = fact["type"]
    scalar = _scalar(kind)
    if fact.get("line") or fact.get("table"):
        return None if words else "this key needs a value"
    if kind.endswith("_list"):
        if not words:
            return "a list key needs at least one value"
    elif fact.get("per_level"):
        if len(words) != 1 and len(words) < nlev:
            return "give one value, or one per level (%d)" % nlev
    elif kind != "unknown" and not words:
        return "this key needs a value"
    skip = set(_evaluated(fact, words))
    for word in read_words(fact, words, nlev):
        if word in skip:
            continue
        try:
            canonical_value(scalar, word, fact)
        except ValueError as exc:
            return str(exc)
    return None


@register
class SchemaCheck(Checker):
    """Keys, types and build features against a generated per-build schema (--schema)."""

    name = "schema"

    @classmethod
    def applies(cls, ctx) -> bool:
        return getattr(ctx.options, "schema", None) is not None

    def run(self) -> None:
        schema = self.ctx.options.schema
        for key in self.deck.keys():
            fact = lookup(key, schema)
            if fact is None:
                reads = unresolved_reads(key, schema)
                if reads:
                    self.info(key, "read through a prefix the schema could not resolve (%s), so its type and "
                              "build requirements are not checked" % ", ".join(reads[:3]),
                              "Check the prefix against that source line.")
                else:
                    self.warn(key, "no query in this build's ERF or AMReX reads this key (per the schema)",
                              "Check the spelling, or generate the schema for the build you run.")
            elif fact["type"] == "removed":
                if key not in RETIRED_KEYS:                  # the keys check already reports a retired key
                    self.error(key, "this input was removed from ERF; the run stops at start-up (%s)"
                               % fact.get("removed_evidence", "?"), "Take it out of the deck.")
            else:
                self._value(key, fact, schema)

    def _value(self, key, fact, schema):
        words = self.deck.table[key].values
        nlev = self.ctx.nlev
        placeholders = [w for w in words if PLACEHOLDER.fullmatch(w)]
        if placeholders:                                 # a configure_file template (Tests/test_files/*.i)
            self.info(key, "%s is filled in by CMake's configure_file, so the value is not checked"
                      % ", ".join(placeholders), "Check the configured deck in the build directory.")
            return
        problem = value_problem(fact, words, nlev)
        if problem:
            self.error(key, problem, "Give a value of the type the key expects (%s)." % fact["type"])
            return
        read = read_words(fact, words, nlev)
        ignored = ignored_words(fact, words, nlev)
        if ignored:
            self.warn(key, "ParmParse reads %d value(s) here and ignores the rest: %s"
                      % (len(read), " ".join(ignored)), "Drop the extra words, or put a # before a trailing note.")
        evaluated = expression_words(fact, read)
        if evaluated:
            self.warn(key, "%s is not a plain %s literal; ParmParse's parser evaluates it at run time, so the schema "
                      "cannot check it" % (", ".join(evaluated), _scalar(fact["type"])),
                      "Check that it evaluates to a %s." % _scalar(fact["type"]))
        features = schema.get("features", {})
        off = [f for f in fact.get("requires", []) if not features.get(f, False)]
        scalar = _scalar(fact["type"])
        for word in read:
            if word in evaluated:
                continue
            need = fact.get("value_requires", {}).get(canonical_value(scalar, word, fact))
            if need and not features.get(need["feature"], False):
                off.append(need["feature"])
        if off:
            self.error(key, "needs %s, which this build does not have" % ", ".join(sorted(set(off))),
                       "Rebuild with the feature, or drop the key.")
