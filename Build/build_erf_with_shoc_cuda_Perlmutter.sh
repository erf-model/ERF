#!/bin/bash
set -e
set -o pipefail

# 1. Resolve ERF_DIR to the repo root (one level up from this script)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export ERF_DIR="$(dirname "$SCRIPT_DIR")"

echo "ERF_DIR set to: $ERF_DIR"

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
cp "$ERF_DIR/Build/cmake_with_shoc_cuda_Perlmutter.sh" "$ERF_DIR/build/"

# 4. Move into build directory
cd "$ERF_DIR/build"

# 5. Run cmake setup
echo "Running cmake_with_shoc_cuda_Perlmutter.sh..."
source cmake_with_shoc_cuda_Perlmutter.sh
