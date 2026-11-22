"""Generator of commented ERF inputs files.

``erfcreateinput.py`` turns a few choices (setup, equation set, turbulence
closure, lower boundary, forcing, outputs) into a complete inputs file:

* the chosen option of each group is active, with its sub-keys;
* the other options of the group are commented lines beside it;
* files the deck names (sounding, roughness table, ...) are written as valid
  dummies, and the summary says how to replace each one.

``options``     the command line and the per-setup defaults
``constraints`` option ranges and the combinations ERF refuses
``recipes``     the keys of each part of the deck
``auxfiles``    the dummy files
``generate``    assembly, writing, summary
``selftest``    the built-in tests
"""

__version__ = "0.1.0"
