#!/bin/bash
set -e
set -o pipefail

# Function to verify if a directory is the ERF repo root
verify_erf_dir() {
    local dir=$1

    # Check for basic structure
    if [ ! -f "$dir/CMakeLists.txt" ] || [ ! -d "$dir/Source" ]; then
        return 1
    fi

    # Check for "Energy Research and Forecasting" in key files
    local found=0

    if [ -f "$dir/README.rst" ]; then
        if grep -q "Energy Research and Forecasting" "$dir/README.rst" 2>/dev/null; then
            found=1
        fi
    fi

    if [ $found -eq 0 ] && [ -f "$dir/LICENSE.md" ]; then
        if grep -q "Energy Research and Forecasting" "$dir/LICENSE.md" 2>/dev/null; then
            found=1
        fi
    fi

    if [ $found -eq 0 ] && [ -f "$dir/CITATION.cff" ]; then
        if grep -q "Energy Research and Forecasting" "$dir/CITATION.cff" 2>/dev/null; then
            found=1
        fi
    fi

    return $((1 - found))
}

# Function to find ERF repo root with multiple fallbacks
find_erf_dir() {
    # Method 1: Use git to find repo root
    if command -v git &> /dev/null; then
        if git rev-parse --is-inside-work-tree &> /dev/null 2>&1; then
            local git_root="$(git rev-parse --show-toplevel)"
            if verify_erf_dir "$git_root"; then
                ERF_DIR="$git_root"
                echo "Detected ERF_DIR from git: $ERF_DIR"
                return 0
            fi
        fi
    fi

    # Method 2: Try going up from script location
    local script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    # Script is in Build/Perlmutter/, so go up 2 levels
    local candidate="$(cd "$script_dir/../.." && pwd)"
    if verify_erf_dir "$candidate"; then
        ERF_DIR="$candidate"
        echo "Detected ERF_DIR from script location: $ERF_DIR"
        return 0
    fi

    # Method 3: Check current directory
    if verify_erf_dir "$PWD"; then
        ERF_DIR="$PWD"
        echo "Detected ERF_DIR from current directory: $ERF_DIR"
        return 0
    fi

    echo "Error: Could not auto-detect ERF_DIR"
    echo "Verification requires:"
    echo "  - CMakeLists.txt and Source/ directory"
    echo "  - 'Energy Research and Forecasting' in README.rst, LICENSE.md, or CITATION.cff"
    return 1
}

###################################################################################

# 1. Resolve ERF_DIR
# Detect ERF_DIR
if ! find_erf_dir; then
    exit 1
fi

export ERF_DIR

E3SM_DIR="$ERF_DIR/external/E3SM"
if [ ! -d "$E3SM_DIR" ]; then
    echo "external/E3SM folder not found, running eamxx_clone.sh..."
    source "$ERF_DIR/Build/GNU_Ekat/eamxx_clone.sh"
else
    echo "external/E3SM folder already exists, skipping clone."
fi

# 3. Prepare build directory
echo "Preparing build directory..."
mkdir -p "$ERF_DIR/build"
cp "$ERF_DIR/Build/cmake_with_shoc.sh" "$ERF_DIR/build/"

# 4. Move into build directory
cd "$ERF_DIR/build"

# Run cmake setup
echo "Running cmake_with_shoc.sh..."
source cmake_with_shoc.sh
