"""Checker base class and the registry that runs them all."""

from typing import Callable, List, Type

from .context import RunContext
from .deck import Deck, DeckError
from .findings import Finding, Report, Severity

_REGISTRY: List[Type["Checker"]] = []


def register(cls: Type["Checker"]) -> Type["Checker"]:
    """Class decorator: add a checker to the list run by :func:`run_all_checks`."""
    _REGISTRY.append(cls)
    return cls


def registered_checkers() -> List[Type["Checker"]]:
    return list(_REGISTRY)


class Checker:
    """One area of the deck.

    Subclasses set ``name`` and implement :meth:`run`.  :meth:`applies` lets a
    mode-specific class (for instance the anelastic time-step checker) opt out
    when the deck is in another mode.
    """

    name = "checker"

    def __init__(self, deck: Deck, ctx: RunContext, report: Report):
        self.deck = deck
        self.ctx = ctx
        self.report = report

    @classmethod
    def applies(cls, ctx: RunContext) -> bool:   # noqa: ARG003
        return True

    def run(self) -> None:
        raise NotImplementedError

    # -- reporting helpers ----------------------------------------------
    def _add(self, sev: Severity, key: str, message: str, advice: str = "") -> None:
        where = self.deck.where(key) if self.deck.has(key) else ""
        self.report.add(Finding(sev, self.name, key, message, advice, where))

    def error(self, key: str, message: str, advice: str = "") -> None:
        self._add(Severity.ERROR, key, message, advice)

    def warn(self, key: str, message: str, advice: str = "") -> None:
        self._add(Severity.WARNING, key, message, advice)

    def advise(self, key: str, message: str, advice: str = "") -> None:
        self._add(Severity.ADVICE, key, message, advice)

    def info(self, key: str, message: str, advice: str = "") -> None:
        self._add(Severity.INFO, key, message, advice)

    # -- typed access that reports a bad value instead of raising --------
    def get(self, getter: Callable, key: str, default=None):
        try:
            return getter(key, default)
        except DeckError as exc:
            self.error(key, str(exc), "Give a value of the type the key expects.")
            return default

    def get_list(self, getter: Callable, key: str):
        try:
            return getter(key)
        except DeckError as exc:
            self.error(key, str(exc), "Give values of the type the key expects.")
            return []

    def int_(self, key, default=None):
        return self.get(self.deck.get_int, key, default)

    def real(self, key, default=None):
        return self.get(self.deck.get_real, key, default)

    def bool_(self, key, default=None):
        return self.get(self.deck.get_bool, key, default)

    def str_(self, key, default=None):
        return self.deck.get_str(key, default)

    def ints(self, key):
        return self.get_list(self.deck.get_ints, key)

    def reals(self, key):
        return self.get_list(self.deck.get_reals, key)


def run_all_checks(deck: Deck, ctx: RunContext, report: Report,
                   only: List[str] = None) -> None:
    # importing the package registers every checker module
    from . import checks  # noqa: F401
    for p in ctx.problems:
        report.add(Finding(Severity.ERROR, "deck", p.split(" =")[0], p,
                           "Give a value of the type the key expects."))
    for cls in _REGISTRY:
        if only and cls.name not in only:
            continue
        if not cls.applies(ctx):
            continue
        try:
            cls(deck, ctx, report).run()
        except DeckError as exc:
            report.add(Finding(Severity.ERROR, cls.name, str(exc).split(" =")[0], str(exc),
                               "Give a value of the type the key expects."))
        except Exception as exc:  # a checker bug must not hide the rest of the report
            report.add(Finding(Severity.ERROR, cls.name, "internal",
                               "the %s checker failed on this deck: %s: %s" % (cls.name, type(exc).__name__, exc),
                               "Please report this with the deck; the other checks still ran."))
