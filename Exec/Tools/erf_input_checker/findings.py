"""Findings and the report that collects them."""

import enum
import json
from typing import Dict, List, Optional


class Severity(enum.IntEnum):
    INFO = 0
    ADVICE = 1
    WARNING = 2
    ERROR = 3

    @property
    def label(self) -> str:
        return self.name


class Finding:
    """One message about the deck.

    ``key`` names the input the message is about (or the first of several),
    ``where`` the file and line it was read from, ``message`` what is wrong and
    ``advice`` what to do about it.
    """

    def __init__(self, severity: Severity, checker: str, key: str,
                 message: str, advice: str = "", where: str = ""):
        self.severity = severity
        self.checker = checker
        self.key = key
        self.message = message
        self.advice = advice
        self.where = where

    def as_dict(self) -> Dict[str, str]:
        return {"severity": self.severity.label, "checker": self.checker,
                "key": self.key, "where": self.where,
                "message": self.message, "advice": self.advice}

    def format(self, show_checker: bool = False) -> str:
        head = "[%s] %s" % (self.severity.label, self.key)
        if self.where:
            head += " (%s)" % self.where
        if show_checker:
            head += " <%s>" % self.checker
        lines = [head, "    " + self.message]
        if self.advice:
            lines.append("    -> " + self.advice)
        return "\n".join(lines)


class Report:
    def __init__(self, path: str):
        self.path = path
        self.findings: List[Finding] = []

    # Convenience adders used by the checkers
    def add(self, finding: Finding) -> None:
        self.findings.append(finding)

    def count(self, severity: Severity) -> int:
        return sum(1 for f in self.findings if f.severity == severity)

    @property
    def worst(self) -> Optional[Severity]:
        if not self.findings:
            return None
        return max(f.severity for f in self.findings)

    def sorted(self) -> List[Finding]:
        return sorted(self.findings, key=lambda f: (-int(f.severity), f.checker, f.key))

    def format(self, min_severity: Severity = Severity.INFO,
               show_checker: bool = False) -> str:
        out = ["== %s" % self.path]
        shown = [f for f in self.sorted() if f.severity >= min_severity]
        for f in shown:
            out.append(f.format(show_checker))
        out.append("-- %d error(s), %d warning(s), %d advice, %d info" % (
            self.count(Severity.ERROR), self.count(Severity.WARNING),
            self.count(Severity.ADVICE), self.count(Severity.INFO)))
        return "\n".join(out)

    def as_json(self) -> str:
        return json.dumps({"path": self.path,
                           "findings": [f.as_dict() for f in self.sorted()]},
                          indent=2)
