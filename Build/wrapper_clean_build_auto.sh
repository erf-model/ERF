#!/bin/bash
set -e

# ============================================================================
# Automated CMake Build Wrapper with Cleanup (distclean equivalent)
# ============================================================================
# Non-interactive version for CI/regression testing
# See wrapper_clean_build.sh for detailed documentation on:
#   - Modern cmake -S -B workflow
#   - Modern cmake --install --prefix usage
#   - GNU distclean behavior
# ============================================================================

SCRIPT=$1

if [ -z "$SCRIPT" ]; then
    echo "ERROR: No build script provided"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: Build script not found: $SCRIPT"
    exit 1
fi

echo "=========================================="
echo "AUTO MODE: Performing distclean"
echo "=========================================="
echo "Deleting CMake configuration and build artifacts..."
echo "(Install directories NOT affected)"
echo ""

# Delete all CMake artifacts (distclean equivalent)
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake \
       CTestTestfile.cmake Testing/ _deps/ compile_commands.json \
       *.cmake 2>/dev/null || true

echo "✓ Cleaned: CMakeCache.txt, CMakeFiles/, Makefile, *.cmake, etc."
echo "✓ Directory ready for fresh configuration"
echo ""

# Set ERF_DIR
if [ -z "$ERF_DIR" ]; then
    if [ -d "../source" ]; then
        export ERF_DIR=$(cd ../source && pwd)
        echo "Auto-detected ERF_DIR: $ERF_DIR"
    fi
fi

echo "=========================================="
echo "Running build script: $SCRIPT"
echo "=========================================="
echo ""

bash "$SCRIPT"
