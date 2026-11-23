"""The inputs file as a list of sections, each a list of lines.

A line is one of

* an active ``key = value`` (``Section.set``);
* an active line the user must review (``Section.user``); it is repeated in
  the printed summary;
* a commented alternative ``#key = value`` (``Section.alt``);
* a plain comment (``Section.note``).

A key may be active only once in a deck: ParmParse takes the last
definition without a message, so a second active line would hide the first.
"""

import math
from typing import List, Optional, Sequence, Tuple, Union

Value = Union[str, int, float, bool, Sequence]

WIDTH = 100
COMMENT_COLUMN = 56   # comments start no further right than this, unless the value is longer


def fmt_value(value: Value) -> str:
    """The value as ParmParse reads it."""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("a non-finite number cannot be written to an inputs file")
        return repr(value) if value != int(value) or abs(value) >= 1e15 else "%.1f" % value
    if isinstance(value, (list, tuple)):
        return " ".join(fmt_value(v) for v in value)
    return str(value)


class Line:
    def __init__(self, kind: str, key: str = "", value: str = "", comment: str = ""):
        self.kind = kind          # "set", "user", "alt", "note", "blank"
        self.key = key
        self.value = value
        self.comment = comment


class Section:
    def __init__(self, deck: "DeckText", title: str):
        self._deck = deck
        self.title = title
        self.lines: List[Line] = []

    def set(self, key: str, value: Value, comment: str = "") -> "Section":
        self._deck.claim(key)
        self.lines.append(Line("set", key, fmt_value(value), comment))
        return self

    def user(self, key: str, value: Value, comment: str = "") -> "Section":
        self._deck.claim(key)
        self.lines.append(Line("user", key, fmt_value(value), comment))
        self._deck.user_items.append((key, fmt_value(value), comment))
        return self

    def alt(self, key: str, value: Value, comment: str = "") -> "Section":
        self.lines.append(Line("alt", key, fmt_value(value), comment))
        return self

    def note(self, text: str) -> "Section":
        self.lines.append(Line("note", comment=text))
        return self

    def blank(self) -> "Section":
        self.lines.append(Line("blank"))
        return self

    def choice(self, key: str, chosen: str, values: Sequence[str], comment: str = "",
               describe: Optional[dict] = None) -> "Section":
        """The chosen enum value active, every other value as a commented alternative."""
        describe = describe or {}
        self.set(key, chosen, comment or describe.get(chosen, ""))
        for v in values:
            if v.lower() != chosen.lower():
                self.alt(key, v, describe.get(v, ""))
        return self

    def block(self, on: bool, items: Sequence[Tuple[str, Value, str]]) -> "Section":
        """A group of keys: active when ``on``, commented alternatives otherwise."""
        for key, value, comment in items:
            if on:
                self.set(key, value, comment)
            else:
                self.alt(key, value, comment)
        return self


class DeckText:
    def __init__(self):
        self.header: List[str] = []
        self.sections: List[Section] = []
        self.user_items: List[Tuple[str, str, str]] = []
        self._active = set()

    def claim(self, key: str) -> None:
        if key in self._active:
            raise ValueError("the key %s is set twice in the generated deck" % key)
        self._active.add(key)

    def is_set(self, key: str) -> bool:
        return key in self._active

    def active_keys(self) -> List[str]:
        return sorted(self._active)

    def all_keys(self) -> List[str]:
        """Active keys and commented alternatives."""
        keys = set(self._active)
        for sec in self.sections:
            keys.update(ln.key for ln in sec.lines if ln.kind == "alt")
        return sorted(keys)

    def section(self, title: str) -> Section:
        sec = Section(self, title)
        self.sections.append(sec)
        return sec

    def render(self) -> str:
        out: List[str] = []
        for text in self.header:
            # one comment line per header entry, whatever the text holds
            out.append(("# " + " ".join(str(text).split())).rstrip())
        for sec in self.sections:
            if not sec.lines:
                continue
            out.append("")
            out.append("# " + "=" * (WIDTH - 2))
            out.append("# " + sec.title)
            out.append("# " + "=" * (WIDTH - 2))
            pad = max([len(ln.key) + (1 if ln.kind == "alt" else 0)
                       for ln in sec.lines if ln.key] or [0])
            col = min(COMMENT_COLUMN, max([pad + 3 + len(ln.value) for ln in sec.lines if ln.key] or [0])) + 2
            for ln in sec.lines:
                if ln.kind == "blank":
                    out.append("")
                elif ln.kind == "note":
                    out.append(("# " + ln.comment).rstrip())
                else:
                    lhs = ("#" if ln.kind == "alt" else "") + ln.key
                    text = "%-*s = %s" % (pad, lhs, ln.value)
                    comment = ln.comment
                    if ln.kind == "user":
                        comment = ("USER: " + comment) if comment else "USER: review this value"
                    if comment:
                        # the blanks matter: ParmParse reads "theta#" as part of the value
                        text = text.ljust(col - 2) + "  # " + comment
                    out.append(text.rstrip())
        return "\n".join(out) + "\n"
