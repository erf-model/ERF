#!/bin/bash
set -e
set -o pipefail

# 1. Resolve ERF_DIR to the repo root (one level up from this script)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export ERF_DIR="$(dirname "$SCRIPT_DIR")"

echo "ERF_DIR set to: $ERF_DIR"

# 2. Source EKAT setup
echo "Sourcing EKAT setup..."
source "$ERF_DIR/Build/GNU_Ekat/eamxx_clone.sh"

# 3. Prepare build directory
echo "Preparing build directory..."
mkdir -p "$ERF_DIR/build"
cp "$ERF_DIR/Build/cmake_with_shoc_cuda_Perlmutter.sh" "$ERF_DIR/build/"

# 4. Move into build directory
cd "$ERF_DIR/build"

# 5. Run cmake setup
echo "Running cmake_with_shoc_cuda_Perlmutter.sh..."
source cmake_with_shoc_cuda_Perlmutter.sh
