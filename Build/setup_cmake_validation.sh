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
        if grep -q "Energy Research and Forecasting" "$dir/README.rst" 2>/dev/null || true; then
            found=1
        fi
    fi

    if [ $found -eq 0 ] && [ -f "$dir/LICENSE.md" ]; then
        if grep -q "Energy Research and Forecasting" "$dir/LICENSE.md" 2>/dev/null || true; then
            found=1
        fi
    fi

    if [ $found -eq 0 ] && [ -f "$dir/CITATION.cff" ]; then
        if grep -q "Energy Research and Forecasting" "$dir/CITATION.cff" 2>/dev/null || true; then
            found=1
        fi
    fi

    return $((1 - found))
}

# Function to find ERF repo root
find_erf_dir() {
    # Method 1: Check if we're already in Build/
    if [ -f "../CMakeLists.txt" ] && [ -d "../Source" ]; then
        local candidate="$(cd .. && pwd)"
        if verify_erf_dir "$candidate"; then
            ERF_DIR="$candidate"
            echo "Detected ERF_DIR from Build location: $ERF_DIR"
            return 0
        fi
    fi

    # Method 2: Use git to find repo root
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

    # Method 3: Try going up from script location
    local script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    # Check if script is in Build/ directory
    if [[ "$script_dir" =~ /Build$ ]]; then
        local candidate="$(dirname "$script_dir")"
        if verify_erf_dir "$candidate"; then
            ERF_DIR="$candidate"
            echo "Detected ERF_DIR from script location: $ERF_DIR"
            return 0
        fi
    fi

    # Method 4: Check current directory
    if verify_erf_dir "$PWD"; then
        ERF_DIR="$PWD"
        echo "Detected ERF_DIR from current directory: $ERF_DIR"
        return 0
    fi

    return 1
}

# Parse arguments
if [ $# -lt 1 ] || [ $# -gt 3 ]; then
    echo "Usage: $0 <set> [script_pattern] [erf_dir]"
    echo ""
    echo "Sets:"
    echo "  default    - Scripts from Build/"
    echo "  perlmutter - Scripts from Build/Perlmutter/"
    echo "  gnu_ekat   - Scripts from Build/GNU_Ekat/"
    echo ""
    echo "If script_pattern is provided, creates build_<pattern>/"
    echo "Otherwise creates build_<set>/"
    echo ""
    echo "If erf_dir is provided, uses that as ERF_DIR"
    echo "Otherwise auto-detects ERF repo root"
    exit 1
fi

SET=$1
PATTERN=${2:-}
ERF_DIR_ARG=${3:-}

# Set ERF_DIR
if [ -n "$ERF_DIR_ARG" ]; then
    ERF_DIR="$ERF_DIR_ARG"
    echo "Using provided ERF_DIR: $ERF_DIR"
    if ! verify_erf_dir "$ERF_DIR"; then
        echo "Error: Provided directory is not a valid ERF repository"
        echo "Must contain 'Energy Research and Forecasting' in README.rst, LICENSE.md, or CITATION.cff"
        exit 1
    fi
else
    if ! find_erf_dir; then
        echo "Error: Could not auto-detect ERF_DIR"
        echo "Please provide it as the third argument or run from ERF Build/ directory"
        echo ""
        echo "Verification checks for:"
        echo "  - CMakeLists.txt and Source/ directory"
        echo "  - 'Energy Research and Forecasting' in README.rst, LICENSE.md, or CITATION.cff"
        exit 1
    fi
fi

echo "ERF_DIR set to: $ERF_DIR"

# Define source directories relative to ERF_DIR
DEFAULT_DIR="$ERF_DIR/Build"
PERLMUTTER_DIR="$ERF_DIR/Build/Perlmutter"
GNU_EKAT_DIR="$ERF_DIR/Build/GNU_Ekat"

case $SET in
    default)
        SRC_DIR="$DEFAULT_DIR"
        ;;
    perlmutter)
        SRC_DIR="$PERLMUTTER_DIR"
        ;;
    gnu_ekat)
        SRC_DIR="$GNU_EKAT_DIR"
        ;;
    *)
        echo "Error: Invalid set '$SET'"
        echo "Choose: default, perlmutter, or gnu_ekat"
        exit 1
        ;;
esac

if [ ! -d "$SRC_DIR" ]; then
    echo "Error: Source directory does not exist: $SRC_DIR"
    exit 1
fi

# Determine build directory name
if [ -n "$PATTERN" ]; then
    BUILD_DIR="$ERF_DIR/build_${PATTERN}"
else
    BUILD_DIR="$ERF_DIR/build_${SET}"
fi

# Create build directory
mkdir -p "$BUILD_DIR"
echo "Created directory: $BUILD_DIR"

# Find and copy ERF cmake build scripts
echo "Scanning for ERF cmake scripts in $SRC_DIR:"
COPIED=0
SKIPPED=0

# Temporarily disable exit on error for the loop
set +e

for script in "$SRC_DIR"/*.sh; do
    # Check if file exists (glob might not match anything)
    if [ ! -f "$script" ]; then
        continue
    fi

    basename_script=$(basename "$script")

    # Skip backup files
    if [[ "$basename_script" =~ ~$ ]]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi

    # Check if it's an ERF cmake script (contains DERF or cmake)
    has_derf=0
    has_cmake=0

    grep -q "DERF" "$script" 2>/dev/null && has_derf=1
    grep -q "cmake" "$script" 2>/dev/null && has_cmake=1

    if [ $has_derf -eq 1 ] || [ $has_cmake -eq 1 ]; then
        cp "$script" "$BUILD_DIR/"
        chmod +x "$BUILD_DIR/$basename_script"
        echo "   DONE: $basename_script"
        COPIED=$((COPIED + 1))
    else
        echo "  ERROR: $basename_script (no DERF or cmake found)"
        SKIPPED=$((SKIPPED + 1))
    fi
done

# Re-enable exit on error
set -e

echo ""
echo "Summary: Copied $COPIED script(s), skipped $SKIPPED"

if [ $COPIED -eq 0 ]; then
    echo "Warning: No ERF cmake scripts found"
    echo "Scripts should contain 'DERF' or 'cmake'"
fi

# Create a run script in the build directory
cat > "$BUILD_DIR/run.sh" << 'EOF'
#!/bin/bash

set -e
set -o pipefail

# Resolve ERF_DIR (go up from build directory)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export ERF_DIR="$(dirname "$SCRIPT_DIR")"

echo "ERF_DIR set to: $ERF_DIR"

# Find all .sh scripts (excluding run.sh and backups)
SCRIPTS=()
for script in *.sh; do
    if [ "$script" = "run.sh" ]; then
        continue
    fi
    if [[ "$script" =~ ~$ ]]; then
        continue
    fi
    if [ -f "$script" ]; then
        SCRIPTS+=("$script")
    fi
done

# Sort scripts alphabetically
IFS=$'\n' SCRIPTS=($(sort <<<"${SCRIPTS[*]}"))
unset IFS

if [ ${#SCRIPTS[@]} -eq 0 ]; then
    echo "Error: No build scripts found in this directory"
    exit 1
fi

if [ $# -ne 1 ]; then
    echo "Usage: $0 <number>"
    echo ""
    echo "Available ERF cmake scripts:"
    for i in "${!SCRIPTS[@]}"; do
        script_base="${SCRIPTS[$i]%.sh}"
        printf "%3d: %s\n" $((i+1)) "${SCRIPTS[$i]}"
        printf "     -> subdirectory: %s/script_%s/\n" "$ERF_DIR" "$script_base"
    done
    echo ""
    echo "Each script will run in its own clean subdirectory at ERF root."
    exit 1
fi

NUM=$1
if [ $NUM -lt 1 ] || [ $NUM -gt ${#SCRIPTS[@]} ]; then
    echo "Error: Number must be between 1 and ${#SCRIPTS[@]}"
    exit 1
fi

SCRIPT="${SCRIPTS[$((NUM-1))]}"
SCRIPT_BASE="${SCRIPT%.sh}"
SUBDIR="$ERF_DIR/script_${SCRIPT_BASE}"

# Create a clean subdirectory for this script
if [ -d "$SUBDIR" ]; then
    echo "Warning: $SUBDIR already exists"
    read -p "Delete and recreate? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf "$SUBDIR"
    else
        echo "Aborting. Please remove $SUBDIR manually or choose a different script."
        exit 1
    fi
fi

mkdir -p "$SUBDIR"
echo "========================================"
echo "Running: $SCRIPT"
echo "Build directory: $SUBDIR"
echo "Working directory: $SUBDIR"
echo "========================================"
echo ""

# Copy the script into the subdirectory and run it there
cp "$SCRIPT" "$SUBDIR/"
cd "$SUBDIR"
bash "./$SCRIPT"
EOF

chmod +x "$BUILD_DIR/run.sh"

echo ""
echo "Setup complete!"
echo "Build directory: $BUILD_DIR"
echo ""
echo "To use:"
echo "  cd $BUILD_DIR"
echo "  ./run.sh           # List available scripts"
echo "  ./run.sh <number>  # Run a specific script"
echo ""
echo "Copied $COPIED script(s)"
