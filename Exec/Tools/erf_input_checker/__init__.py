"""Static checker for ERF inputs files.

The package parses an ERF (AMReX ParmParse) inputs file without running ERF
and reports, in order of severity:

* ERROR   -- the run would stop at start-up (an AMReX or ERF abort), or a key
             that cannot be read;
* WARNING -- the run would start but the setting is very likely not what the
             author meant, or is known to go unstable;
* ADVICE  -- a working setting that could be chosen better (load balance,
             adaptive versus fixed time step, box sizes, output cadence);
* INFO    -- derived numbers (cell sizes, estimated time step, box counts)
             that help a reader judge the deck.

Each area of the deck is checked by its own class in ``checks/``; the
time-step and initialisation checks have one class per mode (compressible or
anelastic, idealised or real-data), chosen from the deck.
"""

from .deck import Deck, parse_deck            # noqa: F401
from .findings import Finding, Report, Severity  # noqa: F401
from .context import RunContext               # noqa: F401
from .checker import run_all_checks           # noqa: F401

__version__ = "0.1.0"
