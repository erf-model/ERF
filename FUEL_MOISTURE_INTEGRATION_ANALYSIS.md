# Fuel Moisture Integration Analysis

## Executive Summary

The `advance_fuel_moisture()` method is fully implemented and ready for integration into the `advance()` method in the FireLayer class. This document analyzes the current implementation, identifies integration requirements, and suggests integration options.

## Current Implementation Status

### Implementation Location
- **Header:** `Source/Fire/ERF_FireLayer.H` (lines 213-225)
- **Implementation:** `Source/Fire/ERF_FireLayer.cpp` (lines 238-312)
- **Status:** COMPLETE - Fully functional with proper physics implementation

### Function Signature
```cpp
void advance_fuel_moisture(amrex::Real dt_s,
                           const amrex::MultiFab& T_atm_k0,
                           const amrex::MultiFab& RH_atm_k0);
```

### Required Inputs
1. **dt_s** (Real): Atmospheric timestep [seconds]
   - Already available at call site in `ERF::Advance()`
   - **Source:** `dt_lev` parameter passed to `FireLayer::advance()`

2. **T_atm_k0** (MultiFab): Atmospheric potential temperature at k=0 [K]
   - **Source:** Must be derived from `Theta_prim[lev]` at k=0
   - **Type:** Potential temperature (not absolute temperature)
   - **Grid:** Atmospheric level 0 grid

3. **RH_atm_k0** (MultiFab): Atmospheric relative humidity at k=0 [fraction 0-1]
   - **Source:** Must be computed from conserved variables
   - **Type:** Fraction (0.0 to 1.0), not percentage
   - **Grid:** Atmospheric level 0 grid

### Current Implementation Details

The function:
- Converts timestep from seconds to hours
- Iterates over all fire grid cells using `MFIter`
- For each fire cell, maps to corresponding atmospheric grid cell (accounting for refinement factor C)
- Reads atmospheric temperature and RH at k=0
- Updates fuel moisture for three fuel classes (1hr, 10hr, 100hr) using time-lag ODEs
- Recomputes moisture of extinction
- Handles hysteresis in moisture equilibrium
- Applies temperature and precipitation corrections

### Implementation Validation
✓ Properly handles grid refinement (fire grid refined by factor C relative to atmospheric grid)
✓ Converts K → °C for temperature
✓ Converts RH fraction to percent for internal calculations
✓ Updates multiple moisture components
✓ Recomputes moisture of extinction after update
✓ GPU-compatible kernels with AMREX_DEVICE annotations

## Integration Requirements

### Requirement 1: Atmospheric State Extraction

**Issue:** The `advance()` method currently receives only `dt` and `surface_layer`. It does not have direct access to atmospheric state variables.

**Current Access at Call Site (ERF_Advance.cpp, line 418):**
```cpp
if (lev == 0 && m_fire_layer) {
    m_fire_layer->advance(dt_lev, *m_SurfaceLayer);
}
```

**Available Data in `ERF::Advance()` context:**
- `S_old` and `S_new`: Conserved variables (include RhoTheta, RhoQ1)
- `Theta_prim[lev]`: Primitive potential temperature
- `Qv_prim[lev]`: Primitive water vapor mixing ratio
- `Geom(lev)`: Geometry on level 0

### Requirement 2: Grid Alignment

**Challenge:** The fire layer operates on a refined grid (fire cells refined by factor C relative to atmospheric cells).

**Current Implementation:** The `advance_fuel_moisture()` method already handles this:
```cpp
// Map fire grid cell to atmospheric grid cell
int i_a = i_f / C;  // i_f = fire grid index, C = refinement factor
int j_a = j_f / C;
```

**Implication:** T_atm_k0 and RH_atm_k0 should be provided on the coarser atmospheric grid (not interpolated to fire grid).

### Requirement 3: Data Types and Units

**Temperature:**
- **Input:** Potential temperature (K)
- **Grid:** Atmospheric grid at k=0 (first vertical level)
- **Why:** Consistent with atmospheric state variables in ERF

**Relative Humidity:**
- **Input:** Fraction (0.0 to 1.0)
- **Grid:** Atmospheric grid at k=0
- **Why:** Avoids floating-point precision issues with percentage representation

## Integration Options

### Option 1: Minimal Modification - Pass MultiFabs to advance()
**Approach:** Modify the `advance()` signature to accept atmospheric state MultiFabs.

```cpp
// Modified signature
void advance(amrex::Real dt, 
             SurfaceLayer& surface_layer,
             const amrex::MultiFab& T_atm_k0,
             const amrex::MultiFab& RH_atm_k0);

// Implementation
void FireLayer::advance(Real dt, SurfaceLayer& surface_layer,
                        const MultiFab& T_atm_k0,
                        const MultiFab& RH_atm_k0)
{
    // ... existing code ...
    
    // Call fuel moisture update
    advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);
}
```

**Pros:**
- Clear, explicit interface
- Direct access to required data
- Minimal overhead
- Easy to understand and debug

**Cons:**
- Requires caller to extract and pass MultiFabs
- Couples FireLayer interface to atmospheric state structure
- Call site becomes more complex

**Call Site Modification (ERF_Advance.cpp):**
```cpp
// Extract atmospheric state at k=0
MultiFab T_atm_k0(S_old.boxArray(), S_old.DistributionMap(), 1, 0);
MultiFab RH_atm_k0(S_old.boxArray(), S_old.DistributionMap(), 1, 0);

// Populate from conserved/primitive variables
// T_atm_k0 = Theta_prim[lev] at k=0
// RH_atm_k0 = computed from S_old at k=0

if (lev == 0 && m_fire_layer) {
    m_fire_layer->advance(dt_lev, *m_SurfaceLayer, T_atm_k0, RH_atm_k0);
}
```

### Option 2: Accessor Methods - Query SurfaceLayer or ERF
**Approach:** Have `FireLayer::advance()` query a larger object (ERF or SurfaceLayer) for atmospheric state.

```cpp
// Modified signature - no change needed
void advance(amrex::Real dt, SurfaceLayer& surface_layer);

// Implementation with getters
void FireLayer::advance(Real dt, SurfaceLayer& surface_layer)
{
    // Query atmospheric state from SurfaceLayer or ERF
    const MultiFab* T_atm_k0 = surface_layer.get_T_atm_k0();
    const MultiFab* RH_atm_k0 = surface_layer.get_RH_atm_k0();
    
    if (T_atm_k0 && RH_atm_k0) {
        advance_fuel_moisture(dt, *T_atm_k0, *RH_atm_k0);
    }
}
```

**Pros:**
- FireLayer interface unchanged
- Decouples call site from fire implementation
- Centralizes atmospheric state management
- Flexible - can change storage location without affecting FireLayer

**Cons:**
- Requires extending SurfaceLayer or passing ERF reference
- Adds latency from method calls
- May store redundant copies of data

### Option 3: Compute on-Demand Within advance()
**Approach:** Have `advance()` compute atmospheric state from SurfaceLayer diagnostics or other available fields.

```cpp
// Modified signature - no change needed
void advance(amrex::Real dt, SurfaceLayer& surface_layer);

// Implementation with computation
void FireLayer::advance(Real dt, SurfaceLayer& surface_layer)
{
    // Derive T and RH from surface layer diagnostics
    // Example: use MOST diagnostics or surface flux calculations
    
    // This approach requires:
    // 1. Access to SurfaceLayer internals or diagnostic fields
    // 2. Implementation of diagnostic computations
}
```

**Pros:**
- FireLayer remains self-contained
- No changes to interface or call site
- Single responsibility: FireLayer manages its own data

**Cons:**
- Duplicates computation that may exist elsewhere
- Requires extracting/deriving from MOST diagnostics
- May have consistency issues between derived fields
- Performance overhead from recomputation

### Option 4: Lazy Initialization - Store Reference to ERF
**Approach:** Store a pointer to parent ERF instance and query as needed.

```cpp
// Add to FireLayer
class FireLayer {
private:
    ERF* m_erf = nullptr;  // Back-pointer to parent
    
public:
    void set_erf_parent(ERF* erf) { m_erf = erf; }
};

// Implementation
void FireLayer::advance(Real dt, SurfaceLayer& surface_layer)
{
    if (m_erf) {
        // Query atmospheric fields from ERF
        const MultiFab& T_atm_k0 = m_erf->get_T_atm_k0(0);
        const MultiFab& RH_atm_k0 = m_erf->get_RH_atm_k0(0);
        advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0);
    }
}
```

**Pros:**
- Flexible - queries data as needed
- Easy to extend if more fields needed
- Minimal overhead

**Cons:**
- Circular dependency pattern (child → parent)
- Requires ERF to expose atmospheric state getters
- Breaks encapsulation
- Hard to test FireLayer in isolation

## Recommended Option: Option 1

**Rationale:**
1. **Clarity:** Explicit parameters make dependencies clear
2. **Locality:** All needed data passed at call site
3. **Testability:** Can test FireLayer with mock MultiFabs
4. **Performance:** No indirect lookups or recomputation
5. **Maintainability:** Future developers understand data flow immediately
6. **Physics:** Temperature and RH are fundamental to fuel moisture - explicit parameters emphasize this

## Implementation Steps for Option 1

### Step 1: Modify FireLayer::advance() Signature
**File:** `Source/Fire/ERF_FireLayer.H`

Update the method declaration:
```cpp
void advance(amrex::Real dt, 
             SurfaceLayer& surface_layer,
             const amrex::MultiFab& T_atm_k0,
             const amrex::MultiFab& RH_atm_k0);
```

Update the documentation to explain input grids and units.

### Step 2: Implement the Call in advance()
**File:** `Source/Fire/ERF_FireLayer.cpp`

In the `advance()` method, call `advance_fuel_moisture()`:

```cpp
// After step 2 (copy reference wind to effective wind) 
// and before step 5 (compute rate-of-spread)
// Insert fuel moisture update:

// Optional: only if dynamic moisture enabled
if (m_params.moisture_dynamic) {
    advance_fuel_moisture(dt_s, T_atm_k0, RH_atm_k0);
}
```

### Step 3: Extract Atmospheric State at Call Site
**File:** `Source/TimeIntegration/ERF_Advance.cpp`

Before the fire layer advance call (around line 418):

```cpp
#ifdef ERF_ENABLE_FIRE
// Advance fire simulation at level 0
if (lev == 0 && m_fire_layer) {
    // Extract atmospheric state at k=0
    // Temperature: already available in Theta_prim
    // RH: needs to be computed from conserved variables
    
    // Option A: Use primitive Theta directly (already at k=0)
    const MultiFab& T_atm_k0 = *Theta_prim[lev];
    
    // Option B: Compute RH from RhoQ1, Rho, and RhoTheta
    // RH = compute_relative_humidity(RhoQ1, Rho, RhoTheta)
    MultiFab RH_atm_k0(S_old.boxArray(), S_old.DistributionMap(), 1, 0);
    compute_rh_from_conservative(RH_atm_k0, S_old);
    
    m_fire_layer->advance(dt_lev, *m_SurfaceLayer, T_atm_k0, RH_atm_k0);
}
#endif
```

### Step 4: Implement RH Computation Helper
**File:** To be determined (ERF_Utils.H or ERF_EOS.H)

Add helper function:
```cpp
void compute_rh_from_conservative(amrex::MultiFab& RH_out,
                                  const amrex::MultiFab& S);
```

This function should:
- Extract RhoTheta and RhoQ1 at k=0
- Compute saturation vapor pressure from temperature
- Compute actual vapor pressure from specific humidity
- Compute RH = e/e_s, clamp to [0, 1]

## Data Flow Diagram

```
ERF::Advance() [lev=0, dt=dt_lev]
    ├─ S_old (conserved: RhoTheta, RhoQ1, etc.)
    ├─ Theta_prim[0] (potential temperature at k=0)
    │
    ├─ Extract/Compute:
    │   ├─ T_atm_k0 = Theta_prim[0]
    │   └─ RH_atm_k0 = compute_rh(RhoQ1, Rho, RhoTheta) at k=0
    │
    └─ FireLayer::advance(dt, surface_layer, T_atm_k0, RH_atm_k0)
           │
           ├─ extract_wind_from_most() [already done]
           ├─ apply_waf_to_wind() [optional, already done]
           ├─ apply_farsite_terrain_wind() [optional, already done]
           │
           ├─ → advance_fuel_moisture(dt, T_atm_k0, RH_atm_k0)  ← NEW CALL
           │   ├─ Iterate fire cells
           │   ├─ Map fire_cell → atm_cell (using refinement C)
           │   ├─ Read T_K, RH_frac at atm grid point
           │   ├─ Update MC_1hr, MC_10hr, MC_100hr
           │   └─ Recompute M_extinction
           │
           ├─ compute_ros_field() [now uses updated fuel moisture]
           ├─ advance_fire_subcycle() [propagate level-set]
           └─ fill_boundary_on_phi() [ghost cells]
```

## Testing Considerations

### Unit Tests
- **Test grid refinement mapping:** Verify fire cell → atmospheric cell mapping
- **Test RH computation:** Verify RH falls in [0, 1] range
- **Test moisture update:** Verify forward Euler update behaves correctly
- **Test hysteresis:** Verify adsorption/desorption path selection

### Integration Tests
- **Coupled atmosphere-fire test:** Run with real atmospheric evolution
- **Regression test:** Compare fuel moisture evolution against benchmark
- **Conservation:** Monitor that moisture remains bounded

### Validation
- **Physical units:** Verify temperature in K, RH in fraction
- **Grid alignment:** Check that fire grid correctly maps to atmospheric grid
- **Numerical stability:** Verify moisture updates remain stable over long periods

## Alternative Clarifications for Atmospheric Interface

If the above implementation options reveal new interface requirements, the following clarifications may be needed from the ERF team:

1. **Preferred RH Representation:**
   - Currently assumes: RH as fraction (0.0-1.0)
   - Alternative: RH as percentage (0-100)
   - Alternative: Mixing ratio directly

2. **Temperature at Fire Grid:**
   - Should T_atm_k0 be on atmospheric grid (current assumption)?
   - Or interpolated/interpolatable to fire grid?
   - How to handle boundary effects at fire grid edges?

3. **Temporal Interpolation:**
   - Is T_atm_k0, RH_atm_k0 at time `n` or `n+dt`?
   - Should we use time-averaged values?
   - What about multi-level time stepping?

4. **Vertical Levels:**
   - k=0 means "first atmospheric level above surface"?
   - Or "lowest level in atmosphere"?
   - Height conventions (pressure vs. height levels)?

5. **Grid Box vs. Node:**
   - Should atmospheric fields be cell-centered or node-centered?
   - Fire grid is cell-centered (typical for level-set methods)
   - Ensure consistency in grid staggering

## Conclusion

The `advance_fuel_moisture()` implementation is complete and well-structured. The recommended Option 1 (explicit MultiFab parameters) provides the clearest path to integration while maintaining clean dependencies and testability.

The implementation can proceed immediately upon confirming:
1. RH computation method from conserved variables
2. Exact timing of atmospheric fields (time n vs n+dt)
3. Grid alignment conventions for refined fire grids

No physics changes are required to the fuel moisture model itself.
