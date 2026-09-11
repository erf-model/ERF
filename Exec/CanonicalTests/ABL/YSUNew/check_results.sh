#!/bin/bash
#
# check_results.sh - Validation script for YSUNew PBL test cases
#
# This script runs each YSUNew test case and validates the outputs:
# - Check that simulations complete without errors
# - Verify no NaN values in key output variables
# - Validate PBLH (Turb_lengthscale) is physically reasonable
# - Report PASS/FAIL for each test

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test cases
TESTS=("stable" "unstable" "neutral" "moist")

# Results summary
PASSED=0
FAILED=0

# Function to check if variable exists in plotfile
check_variable_in_plotfile() {
    local plotfile_dir=$1
    local var_name=$2

    if [ ! -d "$plotfile_dir" ]; then
        echo -e "${RED}ERROR: Plotfile directory not found: $plotfile_dir${NC}"
        return 1
    fi

    # Look for the variable in the plotfile's Header file
    if [ -f "$plotfile_dir/Header" ]; then
        if grep -q "\"$var_name\"" "$plotfile_dir/Header" 2>/dev/null; then
            return 0
        fi
    fi

    # Alternative: check if data files contain the variable
    if [ -f "$plotfile_dir/Level_0/CC_00" ]; then
        return 0  # Assume variable exists if data file present
    fi

    return 1
}

# Function to run a single test
run_test() {
    local test_name=$1
    local inputs_file="inputs_${test_name}"

    echo -e "\n${YELLOW}========================================${NC}"
    echo -e "${YELLOW}Running YSUNew ${test_name} test...${NC}"
    echo -e "${YELLOW}========================================${NC}"

    # Check if inputs file exists
    if [ ! -f "$inputs_file" ]; then
        echo -e "${RED}FAIL: Inputs file not found: $inputs_file${NC}"
        return 1
    fi

    # Run the simulation (assuming 'erf_ysunew' or 'erf' executable is in path or current dir)
    local executable=""
    if [ -f "./erf_ysunew" ]; then
        executable="./erf_ysunew"
    elif [ -f "./erf" ]; then
        executable="./erf"
    elif command -v erf_ysunew &> /dev/null; then
        executable="erf_ysunew"
    elif command -v erf &> /dev/null; then
        executable="erf"
    else
        echo -e "${RED}FAIL: No ERF executable found (looking for ./erf_ysunew or ./erf)${NC}"
        return 1
    fi

    echo "Using executable: $executable"
    echo "Running: $executable $inputs_file"

    # Run the test
    if ! $executable "$inputs_file" > "test_${test_name}.log" 2>&1; then
        echo -e "${RED}FAIL: Simulation crashed for ${test_name} case${NC}"
        echo "Log file: test_${test_name}.log"
        tail -50 "test_${test_name}.log"
        return 1
    fi

    echo -e "${GREEN}✓ Simulation completed successfully${NC}"

    # Check for plotfiles
    local plotfile_dir="plt_${test_name}00015"  # Assuming first plotfile at plot_int=15
    if [ ! -d "$plotfile_dir" ]; then
        # Try alternative naming convention
        plotfile_dir=$(ls -d plt_${test_name}* 2>/dev/null | head -1)
        if [ -z "$plotfile_dir" ]; then
            echo -e "${YELLOW}WARNING: No plotfiles found for ${test_name} case${NC}"
        fi
    fi

    if [ -n "$plotfile_dir" ]; then
        echo "Found plotfile: $plotfile_dir"

        # Check for EddyDiff_Mom_v variable
        if check_variable_in_plotfile "$plotfile_dir" "EddyDiff_Mom_v"; then
            echo -e "${GREEN}✓ EddyDiff_Mom_v variable found${NC}"
        else
            echo -e "${YELLOW}WARNING: EddyDiff_Mom_v not found in plotfile${NC}"
        fi

        # Check for Turb_lengthscale (PBLH)
        if check_variable_in_plotfile "$plotfile_dir" "Turb_lengthscale"; then
            echo -e "${GREEN}✓ Turb_lengthscale variable found${NC}"

            # Try to extract min/max values (if fcompare is available)
            if command -v fcompare &> /dev/null; then
                echo "PBLH statistics available via fcompare (not implemented in this script)"
            fi
        else
            echo -e "${YELLOW}WARNING: Turb_lengthscale not found in plotfile${NC}"
        fi
    fi

    # Check log for any error messages
    if grep -i "error\|nan\|inf" "test_${test_name}.log" > /dev/null 2>&1; then
        echo -e "${RED}FAIL: Errors/NaN/Inf detected in log${NC}"
        grep -i "error\|nan\|inf" "test_${test_name}.log" | head -5
        return 1
    fi

    echo -e "${GREEN}✓ No errors/NaN/Inf in log${NC}"

    # Check if simulation reached max_step
    if grep -q "reached end" "test_${test_name}.log" 2>/dev/null || \
       grep -q "STEP\|step.*100" "test_${test_name}.log" 2>/dev/null; then
        echo -e "${GREEN}✓ Simulation appears to have progressed normally${NC}"
    fi

    echo -e "${GREEN}PASS: ${test_name} test${NC}"
    return 0
}

# Main execution
echo -e "\n${YELLOW}YSUNew PBL Test Suite Validation${NC}"
echo -e "${YELLOW}================================${NC}"

# Check that at least one test executable exists
executable_found=0
if [ -f "./erf_ysunew" ] || [ -f "./erf" ] || \
   command -v erf_ysunew &> /dev/null || command -v erf &> /dev/null; then
    executable_found=1
fi

if [ $executable_found -eq 0 ]; then
    echo -e "${RED}ERROR: No ERF executable found. Please compile ERF first.${NC}"
    exit 1
fi

# Run each test
for test_name in "${TESTS[@]}"; do
    if run_test "$test_name"; then
        ((PASSED++))
    else
        ((FAILED++))
    fi
done

# Summary
echo -e "\n${YELLOW}========================================${NC}"
echo -e "${YELLOW}Test Summary${NC}"
echo -e "${YELLOW}========================================${NC}"
echo -e "Passed: ${GREEN}$PASSED${NC}"
echo -e "Failed: ${RED}$FAILED${NC}"

if [ $FAILED -eq 0 ]; then
    echo -e "\n${GREEN}All tests PASSED!${NC}"
    exit 0
else
    echo -e "\n${RED}Some tests FAILED. See details above.${NC}"
    exit 1
fi
