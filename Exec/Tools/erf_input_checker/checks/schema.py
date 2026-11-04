"""Checks against a generated per-build inputs schema.

Runs only when ``check_erf_inputs.py --schema FILE`` gives one (written by ``Exec/Tools/build_erf_schema.py``). Each
deck key is looked up in the schema, exactly or by a family pattern such as ``erf.<station>.lat``:

* not in the schema: WARNING, since no query in this build's ERF or AMReX reads it;
* ``removed``: ERROR, since ERF stops at start-up when it is present (the schema names the line), unless the keys
  check already reports it from RETIRED_KEYS;
* a value that is not a plain literal of the key's type, a wrong per-level count, or an enum name ERF does not
  know: ERROR. A key of type ``unknown`` takes any number of words;
* an expression for a bool, int or real key (``2*0.5``, ``${x}``): WARNING, since ParmParse evaluates it at run time;
* the key, or this value of it, needs a build feature the schema says is off: ERROR.
"""
import os
import re
import sys

from ..checker import Checker, register
from .keys import RETIRED_KEYS

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from erf_inputs_schema.canonical import canonical_value  # noqa: E402

PLACEHOLDERS = {"<face>": "(?:xlo|xhi|ylo|yhi|zlo|zhi)", "<station>": "[A-Za-z][A-Za-z0-9_]*", "<n>": "[0-9]+"}


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


def cross_check(schema, retired=RETIRED_KEYS):
    """(curated_only, source_only): retired keys the schema's source does not mark removed, and removed keys (with
    their evidence line) the curated list lacks."""
    removed = {k: f["removed_evidence"] for k, f in schema["keys"].items() if f["type"] == "removed"}
    curated_only = sorted(k for k in retired if k not in removed)
    source_only = sorted((k, ev) for k, ev in removed.items() if k not in retired)
    return curated_only, source_only


EVALUATED = ("bool", "int", "real")               # ParmParse evaluates an expression for these (pp_parser)
EXPRESSION = re.compile(r"\$\{|[*/^()]|.[+-](?<![eE][+-])")


def _scalar(kind):
    return kind[:-5] if kind.endswith("_list") else kind


def expression_words(fact, words):
    """The words ParmParse would evaluate as an expression rather than read as a literal (bool, int, real keys)."""
    if _scalar(fact["type"]) not in EVALUATED:
        return []
    return [w for w in words if EXPRESSION.search(w)]


def value_problem(fact, words, nlev):
    """Why these words are not a value of the key, or None. Expressions are not problems (expression_words)."""
    kind = fact["type"]
    scalar = _scalar(kind)
    if kind.endswith("_list"):
        if not words:
            return "a list key needs at least one value"
    elif fact.get("per_level"):
        if len(words) not in (1, nlev):
            return "give one value, or one per level (%d)" % nlev
    elif kind != "unknown" and len(words) != 1:
        return "this key takes one value, not %d" % len(words)
    skip = set(expression_words(fact, words))
    for word in words:
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
        problem = value_problem(fact, words, self.ctx.nlev)
        if problem:
            self.error(key, problem, "Give a value of the type the key expects (%s)." % fact["type"])
            return
        evaluated = expression_words(fact, words)
        if evaluated:
            self.warn(key, "%s is an expression; ParmParse evaluates it at run time, so the schema cannot check it"
                      % ", ".join(evaluated), "Check that it evaluates to a %s." % _scalar(fact["type"]))
        features = schema.get("features", {})
        off = [f for f in fact.get("requires", []) if not features.get(f, False)]
        scalar = _scalar(fact["type"])
        for word in words:
            if word in evaluated:
                continue
            need = fact.get("value_requires", {}).get(canonical_value(scalar, word, fact))
            if need and not features.get(need["feature"], False):
                off.append(need["feature"])
        if off:
            self.error(key, "needs %s, which this build does not have" % ", ".join(sorted(set(off))),
                       "Rebuild with the feature, or drop the key.")
