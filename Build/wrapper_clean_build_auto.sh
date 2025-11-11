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
CLEAN_INSTALL=${2:-no}  # Optional: "yes" to also clean install dir

if [ -z "$SCRIPT" ]; then
    echo "ERROR: No build script provided"
    echo "Usage: $0 <build_script.sh> [clean_install]"
    echo "  clean_install: 'yes' to also remove CMAKE_INSTALL_PREFIX (default: no)"
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
echo ""

# Check for install directory before deleting CMakeCache.txt
INSTALL_DIR=""
if [ -f "CMakeCache.txt" ] && [ "$CLEAN_INSTALL" = "yes" ]; then
    INSTALL_PREFIX=$(grep "^CMAKE_INSTALL_PREFIX:" CMakeCache.txt 2>/dev/null | cut -d'=' -f2 || true)
    if [ -n "$INSTALL_PREFIX" ] && [ -d "$INSTALL_PREFIX" ]; then
        INSTALL_DIR="$INSTALL_PREFIX"
        echo "Install directory found: $INSTALL_DIR"
    fi
fi

# Delete CMake configuration files
rm -rf CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake \
       CTestTestfile.cmake Testing/ _deps/ compile_commands.json \
       2>/dev/null || true

# Delete CTest/CDash artifacts
rm -f DartConfiguration.tcl 2>/dev/null || true

# Delete generated project config
rm -f ERFConfig.cmake 2>/dev/null || true

# Delete any remaining .cmake files (excluding source CMakeLists.txt)
find . -maxdepth 1 -name "*.cmake" -type f -exec rm -f {} \; 2>/dev/null || true

# Delete pkg-config files
rm -f *.pc 2>/dev/null || true

# Delete built artifact directories
rm -rf Exec/ Submodules/ Tests/ bin/ cmake_packages/ externals/ 2>/dev/null || true

# Delete built libraries
rm -f lib*.a lib*.so 2>/dev/null || true

# Delete build logs
rm -f build_*.log git-state.txt 2>/dev/null || true

echo " DONE: Cleaned: CMake configuration files"
echo " DONE: Cleaned: Build artifacts (Exec/, Submodules/, Tests/, bin/, etc.)"
echo " DONE: Cleaned: Libraries (lib*.a, lib*.so)"
echo " DONE: Cleaned: Build logs and editor backups"

# Clean install directory if requested
if [ -n "$INSTALL_DIR" ]; then
    echo ""
    echo "Cleaning install directory (CLEAN_INSTALL=yes): $INSTALL_DIR"
    rm -rf "$INSTALL_DIR" && echo " DONE: Deleted: $INSTALL_DIR"
elif [ -f "CMakeCache.txt.deleted" ]; then
    echo ""
    echo "Install directory NOT cleaned (use CLEAN_INSTALL=yes to clean)"
fi

echo ""
echo " DONE: Directory ready for fresh configuration"
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
