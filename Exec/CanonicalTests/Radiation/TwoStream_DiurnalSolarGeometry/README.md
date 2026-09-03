# Diurnal Solar Geometry Test

## Objective

Validate that solar zenith angle varies correctly with time of day, latitude, longitude, and day of year.

## Test Purpose

This test confirms that:
- Solar geometry is computed dynamically from simulation time/date/location
- Zenith angle varies smoothly throughout the day (max at noon, min at dawn/dusk)
- Solar geometry accounts for latitude and longitude
- Fixed-angle fallback remains available for single-angle cases
- Backward compatibility maintained (dynamic geometry disabled by default)

## Test Design

### Configuration

- **Domain**: 1000 m × 1000 m horizontal, 10 km vertical (20 layers)
- **Time**: `inputs_dynamic` runs a 2-hour window (stop_time = 7200 s) starting at 00:00 UTC, i.e. 18:00 local time for the UTC-6 site, so it spans the late-afternoon decline of the sun through sunset (SW flux decreasing to exactly zero). Raise `stop_time` to 86400 s for a full diurnal cycle.
- **Location**: Latitude 40°N, Longitude 105°W (Example: Boulder, CO)
- **Day of Year**: Summer solstice (June 21) or equinox (March 21) for symmetry
- **Solar Constant**: S₀ = 1361 W/m²
- **Optical Depth**: τ = 0.05 per layer (clear-sky)
- **Dynamic Solar Geometry**: Enabled

### Key Physics

Solar zenith angle evolves as:
```
cos(θ_z) = sin(lat) · sin(dec) + cos(lat) · cos(dec) · cos(h)
```

where:
- lat = latitude
- dec = solar declination (depends on day-of-year)
- h = hour angle (depends on time-of-day and longitude)

## Files

- `inputs` — Main control file with solar-geometry parameters and time stepping
- Sounding file — Reference atmospheric profile
- `check_solar.py` — Python validation script

## Running the Test

```bash
cd Exec/CanonicalTests/Radiation/TwoStream_DiurnalSolarGeometry
mpirun -np 1 erf.ex inputs
python3 check_solar.py
```

## Validation Criteria

The checker script verifies:

1. **Zenith angle** varies smoothly over the day
2. **Minimum zenith angle** occurs near local noon (maximum flux)
3. **Zenith angle** reaches max/min appropriate to day/latitude combination
4. **Flux at TOA** follows cos(zenith) scaling (S₀ · cos(θ_z))
5. **Sunrise/sunset transitions** handled correctly (zenith > 90°)
6. **Diagnostics file** includes time-evolving flux values
7. **No spurious jumps** in zenith angle or flux
8. **Backward compatibility** (disabled by default, fixed-angle fallback works)

## Expected Output

- Radiation diagnostics with a time-varying solar zenith angle: in the default 2-hour window the SW flux decreases toward sunset and is exactly zero once the sun is below the horizon (a full-day run shows low flux at dawn/dusk and a peak at noon)
- CHECK PASS message confirming solar geometry computation
- Smooth temporal evolution of surface flux
- Zenith angle within physical bounds (0–180°)

## Related Documentation

- `RAD_DEVELOPMENT.md` — Solar Geometry section
- Solar geometry routines in `Source/Radiation/ERF_SolarGeometry.H`
