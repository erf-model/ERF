#!/usr/bin/env bash
# Download, build and install the MoorDyn-C version ERF's line-dynamics coupling links against.
#
#   Build/setup_moordyn.sh [--version 2.7.1] [--prefix DIR] [--src DIR] [--jobs N] [--configure-only] [--force]
#
# Clones the release tag into <src>/moordyn-<version>, configures it the way ERF needs (the C/C++
# library only: no Python, MATLAB, Fortran or Rust wrappers, no tests, no docs, the bundled Eigen),
# builds and installs it under <prefix> (default $HOME/opt/moordyn-<version>), and prints the
# -DMOORDYN_DIR to pass to ERF's cmake. Re-running is safe: an existing clone is reused and an
# existing install is left alone unless --force is given. ERF uses the v2 C API with the external
# wave-kinematics calls (MoorDyn_ExternalWaveKin*), present since 2.3; 2.7.1 is the tested version.
set -euo pipefail

version="2.7.1"
src_root="${HOME}/opt/src"
prefix=""
jobs=4
configure_only=0
force=0
repo="https://github.com/FloatingArrayDesign/MoorDyn"

usage() { sed -n '2,11p' "$0"; exit "${1:-0}"; }
while [[ $# -gt 0 ]]; do
    case "$1" in
        --version)        version="$2"; shift 2 ;;
        --prefix)         prefix="$2"; shift 2 ;;
        --src)            src_root="$2"; shift 2 ;;
        --jobs|-j)        jobs="$2"; shift 2 ;;
        --configure-only) configure_only=1; shift ;;
        --force)          force=1; shift ;;
        -h|--help)        usage 0 ;;
        *) echo "setup_moordyn.sh: unknown argument '$1'" >&2; usage 1 ;;
    esac
done
[[ "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]] || { echo "setup_moordyn.sh: --version must look like 2.7.1, not '$version'" >&2; exit 1; }
major="${version%%.*}"; minor="${version#*.}"; minor="${minor%%.*}"
if [[ "$major" -lt 2 || ( "$major" -eq 2 && "$minor" -lt 3 ) ]]; then
    echo "setup_moordyn.sh: ERF needs MoorDyn-C 2.3 or newer (the v2 C API with the MoorDyn_ExternalWaveKin* calls), not $version" >&2
    exit 1
fi
[[ -n "$prefix" ]] || prefix="${HOME}/opt/moordyn-${version}"
for tool in git cmake; do
    command -v "$tool" > /dev/null || { echo "setup_moordyn.sh: '$tool' not found in PATH" >&2; exit 1; }
done

src="${src_root}/moordyn-${version}"
if [[ -d "$src/.git" ]]; then
    echo "setup_moordyn.sh: reusing the source in $src"
else
    mkdir -p "$src_root"
    echo "setup_moordyn.sh: cloning MoorDyn v${version} into $src"
    git clone --quiet --depth 1 --branch "v${version}" "$repo" "$src"
fi

if [[ $force -eq 0 && -f "$prefix/include/moordyn/MoorDyn2.h" ]] && ls "$prefix"/lib/libmoordyn.* > /dev/null 2>&1; then
    echo "setup_moordyn.sh: MoorDyn ${version} is already installed in $prefix (use --force to rebuild)"
    echo "  ERF: cmake ... -DERF_ENABLE_MOORDYN=ON -DMOORDYN_DIR=$prefix"
    exit 0
fi

build="$src/build"
mkdir -p "$build"
echo "setup_moordyn.sh: configuring in $build (install prefix $prefix)"
cmake -S "$src" -B "$build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="$prefix" \
    -DPYTHON_WRAPPER=OFF \
    -DMATLAB_WRAPPER=OFF \
    -DFORTRAN_WRAPPER=OFF \
    -DRUST_WRAPPER=OFF \
    -DEXTERNAL_EIGEN=OFF \
    -DBUILD_DOCS=OFF \
    -DBUILD_BENCHMARKS=OFF \
    -DBUILD_TESTING=OFF > "$build/setup_moordyn.configure.log" 2>&1 \
    || { tail -30 "$build/setup_moordyn.configure.log" >&2; echo "setup_moordyn.sh: configure failed (log: $build/setup_moordyn.configure.log)" >&2; exit 1; }
if [[ $configure_only -eq 1 ]]; then
    echo "setup_moordyn.sh: configured only; build with: cmake --build $build -j $jobs && cmake --install $build"
    exit 0
fi
echo "setup_moordyn.sh: building with $jobs jobs (log: $build/setup_moordyn.build.log)"
cmake --build "$build" -j "$jobs" > "$build/setup_moordyn.build.log" 2>&1 \
    || { tail -30 "$build/setup_moordyn.build.log" >&2; echo "setup_moordyn.sh: build failed" >&2; exit 1; }
echo "setup_moordyn.sh: installing into $prefix"
cmake --install "$build" > "$build/setup_moordyn.install.log" 2>&1 \
    || { tail -30 "$build/setup_moordyn.install.log" >&2; echo "setup_moordyn.sh: install failed" >&2; exit 1; }
echo "setup_moordyn.sh: done"
echo "  ERF: cmake ... -DERF_ENABLE_MOORDYN=ON -DMOORDYN_DIR=$prefix"
case "$(uname -s)" in
    Darwin) echo "  if erf_exec cannot find libmoordyn at run time: export DYLD_LIBRARY_PATH=$prefix/lib" ;;
    *)      echo "  if erf_exec cannot find libmoordyn at run time: export LD_LIBRARY_PATH=$prefix/lib" ;;
esac
