#!/usr/bin/env bash
# Compatibility entry point for callers of the original Linux helper.
set -euo pipefail
export CACHE_TRANSPORT=docker
exec python3 "$(dirname -- "${BASH_SOURCE[0]}")/cache.py"
