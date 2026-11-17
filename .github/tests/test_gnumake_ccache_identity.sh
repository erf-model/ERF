#!/usr/bin/env bash
# Verify that a change in Open MPI's effective compile flags invalidates ccache.
set -euo pipefail

: "${CCACHE_COMPILERCHECK:?Set the GNU Make compiler identity check}"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$work"
export CCACHE_DIR="$work/cache"
export CCACHE_BASEDIR="$work"

cat > probe.cpp <<'EOF'
#include <cstdio>
#ifndef ERF_CCACHE_WRAPPER_VALUE
#error Missing Open MPI wrapper flag
#endif
int main() { std::printf("%d\n", ERF_CCACHE_WRAPPER_VALUE); }
EOF

compile_and_check() {
  local expected=$1
  rm -f probe.o probe
  OMPI_CPPFLAGS="-DERF_CCACHE_WRAPPER_VALUE=$expected" ccache mpicxx -c probe.cpp -o probe.o
  mpicxx probe.o -o probe
  local actual
  actual=$(./probe)
  if [[ "$actual" != "$expected" ]]; then
    printf 'Expected wrapper value %s, got %s\n' "$expected" "$actual" >&2
    return 1
  fi
}

compile_and_check 11
compile_and_check 11
# Ensure this test exercised a direct cache hit rather than bypassing ccache.
hits=$(ccache --print-stats | awk '$1 == "direct_cache_hit" { print $2 }')
if [[ ! "$hits" =~ ^[0-9]+$ ]] || (( hits < 1 )); then
  echo 'Expected a direct cache hit for unchanged wrapper flags.' >&2
  exit 1
fi

compile_and_check 22
misses=$(ccache --print-stats | awk '$1 == "cache_miss" { print $2 }')
if [[ ! "$misses" =~ ^[0-9]+$ ]] || (( misses < 2 )); then
  echo 'Expected a cache miss for the changed wrapper flags.' >&2
  exit 1
fi
# Compare against the compiler without ccache for the changed wrapper flags.
rm -f probe.o probe
OMPI_CPPFLAGS='-DERF_CCACHE_WRAPPER_VALUE=22' mpicxx -c probe.cpp -o probe.o
mpicxx probe.o -o probe
if [[ $(./probe) != 22 ]]; then
  echo 'Uncached compile did not use the changed wrapper flag.' >&2
  exit 1
fi

echo 'MPI wrapper cache identity check passed.'
