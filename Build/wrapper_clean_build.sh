#!/bin/bash
set -e

# ============================================================================
# CMake Build Wrapper with Cleanup (distclean equivalent)
# ============================================================================
#
# MODERN CMAKE PRACTICE:
# ----------------------
# The current best practice for CMake is to use out-of-source builds with:
#   cmake -S <source_dir> -B <build_dir>
#
# For example:
#   cmake -S .. -B build_release
#   cmake --build build_release
#   cmake --install build_release --prefix install_release
#
# This keeps your source tree clean and allows multiple build configurations
# (debug, release, different compilers, etc.) in separate directories.
#
# MODERN CMAKE INSTALL:
# ---------------------
# The newer cmake --install command (CMake 3.15+) provides a cleaner interface
# than the older "make install" or "cmake --build . --target install":
#
#   cmake --install <build_dir> --prefix <install_dir>
#
# Examples:
#   # Install to default CMAKE_INSTALL_PREFIX (set during configure)
#   cmake --install build_release
#
#   # Install to custom location
#   cmake --install build_release --prefix /opt/erf
#
#   # Install to local directory
#   cmake --install build_release --prefix ./install
#
# This is preferred because:
#   - Works regardless of build system (make, ninja, etc.)
#   - Doesn't require entering the build directory
#   - More explicit and consistent syntax
#   - Allows overriding install location without reconfiguring
#
# Note: You can still set CMAKE_INSTALL_PREFIX during configuration:
#   cmake -S .. -B build_release -DCMAKE_INSTALL_PREFIX=/usr/local
#
# ABOUT ERF/Build DIRECTORY:
# --------------------------
# The ERF/Build directory is primarily intended as a single build directory
# for users doing one configuration. If you're testing multiple configurations
# (CPU, GPU, different flags), you should use separate build directories:
#   ERF/build_cpu/
#   ERF/build_gpu/
#   ERF/build_debug/
# etc.
#
# With corresponding install directories if needed:
#   ERF/install_cpu/
#   ERF/install_gpu/
#   ERF/install_debug/
#
# CLEANUP BEHAVIOR (GNU Make Standard):
# -------------------------------------
# This script performs a 'distclean' equivalent operation, which per GNU
# standards means: "Delete all files in the current directory (or created 
# by this makefile) that are created by configuring or building the program."
#
# For CMake, this includes:
#   - CMakeCache.txt (configuration file)
#   - CMakeFiles/ (generated build system files)
#   - *.cmake (generated configuration scripts)
#   - Makefile (if generated)
#   - Any other CMake-generated artifacts
#
# This ensures a completely fresh configuration and build, as if you had
# just unpacked the source distribution.
#
# Note: This does NOT delete install directories - those should be managed
# separately (equivalent to 'uninstall' target in GNU make).
#
# ============================================================================

SCRIPT=$1

if [ -z "$SCRIPT" ]; then
    echo "ERROR: No build script provided"
    echo "Usage: $0 <build_script.sh>"
    exit 1
fi

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: Build script not found: $SCRIPT"
    exit 1
fi

# Check what would be deleted (distclean items)
FILES_TO_DELETE=""
[ -f "CMakeCache.txt" ] && FILES_TO_DELETE="$FILES_TO_DELETE CMakeCache.txt"
[ -d "CMakeFiles" ] && FILES_TO_DELETE="$FILES_TO_DELETE CMakeFiles/"
[ -f "Makefile" ] && FILES_TO_DELETE="$FILES_TO_DELETE Makefile"
[ -f "cmake_install.cmake" ] && FILES_TO_DELETE="$FILES_TO_DELETE cmake_install.cmake"
[ -f "CTestTestfile.cmake" ] && FILES_TO_DELETE="$FILES_TO_DELETE CTestTestfile.cmake"

# Find any other .cmake files (excluding those we might want to keep)
OTHER_CMAKE=$(find . -maxdepth 1 -name "*.cmake" -type f 2>/dev/null | \
              grep -v "cmake_install.cmake\|CTestTestfile.cmake" || true)
if [ -n "$OTHER_CMAKE" ]; then
    FILES_TO_DELETE="$FILES_TO_DELETE $OTHER_CMAKE"
fi

# Additional common CMake artifacts
[ -d "Testing" ] && FILES_TO_DELETE="$FILES_TO_DELETE Testing/"
[ -d "_deps" ] && FILES_TO_DELETE="$FILES_TO_DELETE _deps/"
[ -f "compile_commands.json" ] && FILES_TO_DELETE="$FILES_TO_DELETE compile_commands.json"

# If there's nothing to clean, just run the script
if [ -z "$FILES_TO_DELETE" ]; then
    echo "Directory is already clean, proceeding with build..."
    echo ""
else
    # Show what will be deleted
    echo "=========================================="
    echo "WARNING: About to perform 'distclean'"
    echo "=========================================="
    echo "This will delete all CMake configuration and build artifacts:"
    echo ""
    for f in $FILES_TO_DELETE; do
        if [ -d "$f" ]; then
            echo "  - $f (directory)"
        else
            echo "  - $f"
        fi
    done
    echo ""
    echo "Current directory: $(pwd)"
    echo ""
    echo "This operation matches the GNU make 'distclean' target:"
    echo "  \"Delete all files created by configuring or building\""
    echo ""
    echo "Note: Install directories (if any) will NOT be deleted."
    echo "      Use 'cmake --install <build_dir> --prefix ...' to manage installations."
    echo ""
    
    # Prompt user
    read -p "Delete these files/directories? [y/N] " -n 1 -r
    echo
    
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted by user. Not deleting anything."
        echo ""
        echo "To proceed without cleaning, run the build script directly:"
        echo "  bash $SCRIPT"
        echo ""
        echo "Modern CMake workflow reminder:"
        echo "  1. Configure: cmake -S <src> -B <build>"
        echo "  2. Build:     cmake --build <build>"
        echo "  3. Install:   cmake --install <build> --prefix <install>"
        exit 1
    fi
    
    # Actually delete
    echo ""
    echo "Performing distclean..."
    for f in $FILES_TO_DELETE; do
        if [ -d "$f" ]; then
            rm -rf "$f" && echo "  ✓ Deleted directory: $f"
        elif [ -f "$f" ]; then
            rm -f "$f" && echo "  ✓ Deleted file: $f"
        fi
    done
    echo ""
    echo "Distclean complete. Ready for fresh configuration."
    echo ""
fi

# Set ERF_DIR to the source tree regtest checked out for us
if [ -z "$ERF_DIR" ]; then
    if [ -d "../source" ]; then
        export ERF_DIR=$(cd ../source && pwd)
        echo "Auto-detected ERF_DIR: $ERF_DIR"
    else
        echo "WARNING: Could not auto-detect ERF_DIR"
        echo "Build script may fail if it requires ERF_DIR"
    fi
fi

# Run the actual build script
echo "=========================================="
echo "Running build script: $SCRIPT"
echo "=========================================="
echo ""

bash "$SCRIPT"
