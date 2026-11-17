#!/usr/bin/env bash
# Changing Open MPI wrapper flags must invalidate cached compiler output.
set -euo pipefail

: "${CCACHE_COMPILERCHECK:?Set the GNU Make compiler identity check}"
work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$work"
export CCACHE_DIR="$work/cache"
export CCACHE_BASEDIR="$work"

cat > probe.cpp <<'EOF'
#include <cstdio>
int main() { std::printf("%d\n", ERF_CCACHE_WRAPPER_VALUE); }
EOF

for expected in 11 11 22; do
  OMPI_CPPFLAGS="-DERF_CCACHE_WRAPPER_VALUE=$expected" ccache mpicxx -c probe.cpp -o probe.o
  mpicxx probe.o -o probe
  actual=$(./probe)
  if [[ "$actual" != "$expected" ]]; then
    printf 'Expected wrapper value %s, got %s\n' "$expected" "$actual" >&2
    exit 1
  fi
done

stats=$(ccache --print-stats)
hits=$(awk '$1 == "direct_cache_hit" { print $2 }' <<< "$stats")
misses=$(awk '$1 == "cache_miss" { print $2 }' <<< "$stats")
if [[ "$hits" != 1 || "$misses" != 2 ]]; then
  echo 'Expected one direct hit and two misses.' >&2
  exit 1
fi

OMPI_CPPFLAGS='-DERF_CCACHE_WRAPPER_VALUE=22' mpicxx -c probe.cpp -o probe.o
mpicxx probe.o -o probe
[[ $(./probe) == 22 ]] || { echo 'Uncached compile disagrees.' >&2; exit 1; }
echo 'MPI wrapper cache identity check passed.'
