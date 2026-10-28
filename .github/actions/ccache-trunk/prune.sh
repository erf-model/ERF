#!/usr/bin/env bash
# Compatibility entry point for the Linux Docker transport.
# Usage: prune.sh <package>[:<keep>]... (see prune.py for the retention policy).
set -euo pipefail
exec python3 "$(dirname -- "${BASH_SOURCE[0]}")/prune.py" "$@"
