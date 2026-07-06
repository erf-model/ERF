# FARSITE Spread Accumulation Unit Test Documentation

## Overview

The unit test suite `ERF_GTestFarsiteSpreadAccumulation.cpp` validates that the ERF fire model correctly reproduces the wildfire_levelset reference implementation for spread vector accumulation across multiple fire timesteps.

## Key Test: `SpreadAccumulationAcrossSteps`

This is the **critical test** that validates Bug 2a fix (preserving accumulated spread):

### What It Tests

The test verifies the core invariant of the FARSITE algorithm:

```
Before Bug Fix (INCORRECT):
  Step 1: Compute spread, phi.setVal(0) resets phi globally
  Step 2: phi=0 at all cells → grad_mag is flat → writes spread_arr=0 
          (destroys history)
  Result: Burned interior is lost each step

After Bug Fix (CORRECT):
  Step 1: Compute spread, phi.setVal(0) resets phi globally
  Step 2: phi=0 at all cells → grad_mag is flat → SKIPS spread write
          (preserves accumulated spread from step 1)
  Result: Burned interior reconstructed from accumulated spread
```

### Test Scenario

1. **Initial state**: Simple 1D front
   - Interior (i <= 4, j=5): phi = -1 (burned)
   - Front cell (i=5, j=5): phi = 0.05 (at threshold)
   - Unburned (i >= 6, j=5): phi > 0

2. **Step 1**: Advance fire front
   - Interior cells get nonzero spread vectors
   - Burned region propagates
   - phi is reset to 0 globally

3. **Step 2**: Check spread preservation
   - Interior cells should still have nonzero spread from step 1
   - **Without the fix**: spread would be zeroed (grad_mag < 1e-50 case)
   - **With the fix**: spread is preserved via early return

### Assertion

```cpp
EXPECT_GT(interior_spread_mag_step2, 0.001) 
    << "Interior cells should retain nonzero spread across steps"
```

This assertion **FAILS** if:
- Accumulated spread history is being erased (Bug 2a still present)
- Burned interior cannot be reconstructed from farsite_work
- phi stamping produces incorrect burned region

## Other Tests

### `SingleStepFrontDetection` 
- Validates Pass 1 correctly computes spread at front cells
- Sanity check: nonzero spread should appear after one step

### `SingleCellStampingRaceSafety`
- Tests that multiple propagated points to same cell don't race
- Verifies min() operation is safe for concurrent writes
- Tests Bug 2c fix

### `FireGridGeometryResolution`
- Validates fire grid domain scaling (Bug 1 fix)
- Checks that dx_fire = dx_atm / C

## Running the Tests

```bash
# Build with tests enabled
cmake -DERF_ENABLE_TESTS=ON ..
make -j 4

# Run all unit tests
./erf_unit_tests

# Run only FARSITE tests
./erf_unit_tests --gtest_filter="*FarsiteSpreadTest*"

# Run only the critical accumulation test
./erf_unit_tests --gtest_filter="FarsiteSpreadTest.SpreadAccumulationAcrossSteps"
```

## Validation Against Reference

The test structure follows the two-pass FARSITE algorithm:

1. **Pass 1** (GPU kernel in `advance_farsite_one_step`):
   - Compute spread vectors from phi gradients
   - With fix: preserve accumulated spread in flat-gradient regions

2. **Pass 2** (host-side collection):
   - Collect all nonzero spread positions
   - MPI_Allgatherv to get global front
   - Stamp phi = -1 at propagated positions

3. **Between steps**: 
   - **Critical**: Don't zero spread_arr between steps
   - Accumulation is the mechanism for reconstructing burned interior

## Expected Behavior

### Test Passes When:

✓ Interior cells retain nonzero spread across steps
✓ phi stamping correctly marks burned region
✓ Multiple points to same cell don't race
✓ Fire grid has correct cell sizes (dx_atm/C)

### Test Fails When:

✗ Interior spread is zeroed (grad_mag < 1e-50 branch still writes 0)
✗ Burned region is lost between steps
✗ Race conditions in stamping cause incorrect phi values
✗ Fire grid cell sizes are incorrect

## Integration with CI/CD

The test is added to the main CMakeLists.txt and will be automatically run as part of:
- `ctest` (unit test suite)
- GitHub Actions / CI pipeline
- Pre-commit hooks (if configured)

## Debugging Tips

If tests fail:

1. **`SpreadAccumulationAcrossSteps` fails with low interior_spread_mag_step2**:
   - Check that Bug 2a fix is applied (no `spread_arr = 0.0` in grad_mag branch)
   - Verify early return is in place

2. **`SingleCellStampingRaceSafety` produces wrong phi values**:
   - Check that Bug 2c fix uses `amrex::min()` not direct assignment
   - Verify no race conditions in GPU kernel

3. **`FireGridGeometryResolution` fails with wrong dx**:
   - Check that fire_domain Box is scaled by C in ERF_FireGrid.cpp
   - Verify physical RealBox is unchanged

## References

- Reference implementation: `wildfire_levelset/src/farsite_ellipse.H`
- Bug description: Top-level README with problem statement
- FARSITE paper: Richards (1990) - Elliptical fire growth model
