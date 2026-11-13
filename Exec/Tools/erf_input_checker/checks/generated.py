"""Cases written by erfcreateinput.py: what is still a placeholder.

The generator writes, beside the deck, ``FILES_TO_REPLACE.txt`` (each dummy file with a sha256 fingerprint, and
each file the user must supply) and ``case.json`` (the options that reproduce the case), and marks the values it
cannot know with ``USER:`` comments. This checker reads them and reports

* WARNING  a dummy file that still has the generator's fingerprint, so it holds placeholder values;
* ERROR    a file listed as one to supply that is not in the case directory;
* INFO     the values marked USER:, and that the deck was written by the generator.
"""

import hashlib
import json
import os
import re

from ..checker import Checker, register

REPLACE_LIST = "FILES_TO_REPLACE.txt"
CASE_FILE = "case.json"
HEADER = "ERF inputs written by Exec/Tools/erfcreateinput.py"
_ENTRY_RE = re.compile(r"^(\S+)\s+\((dummy, replace|NOT GENERATED, supply); named by ([^)]*)\)\s*$")
_SHA_RE = re.compile(r"^\s+sha256:\s*([0-9a-f]{64})\s*$")
_USER_RE = re.compile(r"^\s*([A-Za-z0-9_.]+)\s*=.*#\s*USER:")


def read_replace_list(path):
    """Entries of FILES_TO_REPLACE.txt as dicts with name, kind (dummy or supply), key and sha256."""
    entries = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = _ENTRY_RE.match(line.rstrip("\n"))
            if m:
                entries.append({"name": m.group(1), "kind": "dummy" if m.group(2).startswith("dummy") else "supply",
                                "key": m.group(3).strip(), "sha256": None})
                continue
            m = _SHA_RE.match(line.rstrip("\n"))
            if m and entries:
                entries[-1]["sha256"] = m.group(1)
    return entries


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@register
class GeneratedCase(Checker):
    """A case written by erfcreateinput.py: dummy files not yet replaced, files not yet supplied, USER values."""

    name = "generated"

    def run(self) -> None:
        d = self.deck
        case_dir = os.path.dirname(os.path.abspath(d.path))
        try:
            with open(d.path, encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except OSError:
            return
        replace_list = os.path.join(case_dir, REPLACE_LIST)
        written = HEADER in text[:400]
        if not written and not os.path.exists(replace_list):
            return
        if written:
            how = "erfcreateinput.py --case-file %s" % CASE_FILE if os.path.exists(os.path.join(case_dir, CASE_FILE)) else ""
            self.info("generated", "this deck was written by erfcreateinput.py" +
                      ("; %s writes it again" % how if how else ""))
        if os.path.exists(replace_list):
            try:
                entries = read_replace_list(replace_list)
            except OSError as exc:
                self.warn(REPLACE_LIST, "cannot be read: %s" % exc)
                entries = []
            for e in entries:
                path = os.path.join(case_dir, e["name"])
                key = e["key"] if d.has(e["key"]) else e["name"]
                if e["kind"] == "supply":
                    if not os.path.exists(path):
                        self.error(key, "%s is listed in %s as a file to supply, and it is not in the case directory; "
                                   "the run stops without it" % (e["name"], REPLACE_LIST),
                                   "Put the file in %s; %s says what it must hold." % (case_dir, REPLACE_LIST))
                elif not os.path.exists(path):
                    self.error(key, "%s, written by erfcreateinput.py, is missing" % e["name"],
                               "Restore it, or write the case again with --case-file %s." % CASE_FILE)
                elif e["sha256"] and _sha256(path) == e["sha256"]:
                    self.warn(key, "%s is still the dummy written by erfcreateinput.py: valid in format, placeholder "
                              "values" % e["name"], "Replace it with real data; %s gives its layout." % REPLACE_LIST)
        users = []
        for line in text.splitlines():
            m = _USER_RE.match(line)      # a commented line starts with '#', which the pattern does not take
            if m:
                users.append(m.group(1))
        if users:
            shown = ", ".join(users[:12]) + (", ..." if len(users) > 12 else "")
            self.info("USER", "%d value(s) marked USER: by erfcreateinput.py, which it could not know: %s"
                      % (len(users), shown), "Review each one; remove the USER: mark once it is set.")
