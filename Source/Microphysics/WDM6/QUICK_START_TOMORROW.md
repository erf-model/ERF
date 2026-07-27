# Quick Start for Tomorrow

## TL;DR
Copy WRF Fortran directly (like WSM6 did). 1-2 hours to working physics.

---

## Step-by-Step Commands

```bash
# 1. Navigate and checkout
cd /g/g10/wise14/compiling/clean/ERF
git checkout WDM6
git status  # Should show clean, 6 commits ahead

# 2. Backup current stub
cd Source/Microphysics/WDM6
cp ERF_module_mp_wdm6.F90 ERF_module_mp_wdm6.F90.STUB

# 3. Copy WRF source (3237 lines of real physics!)
cp /p/lustre1/wise14/wrf/compiling/WRF/WRF_v4.7.1_clean/phys/module_mp_wdm6.F \
   ERF_module_mp_wdm6.F90

# 4. Quick check
wc -l ERF_module_mp_wdm6.F90  # Should be ~3237 lines
head -20 ERF_module_mp_wdm6.F90  # Should see WRF copyright
```

---

## Required Edits (2 small changes)

### Edit 1: Add iso_c_binding

Find line ~10 where it says:
```fortran
   use module_mp_radar
```

Add BEFORE it:
```fortran
   use iso_c_binding, only: c_double
```

### Edit 2: Add kind_phys

Find line ~11 where it says:
```fortran
   implicit none
```

Add AFTER it:
```fortran
   integer, parameter :: kind_phys = c_double
```

### Edit 3: Module name (probably already correct)

Check line ~9, should say:
```fortran
   module module_mp_wdm6
```

**Change to:**
```fortran
   module mp_wdm6
```

(Our C wrappers expect `mp_wdm6`)

### That's It!

Those are literally the ONLY changes needed. Rest is identical to WRF.

---

## Build & Test

```bash
# Go back to ERF root
cd /g/g10/wise14/compiling/clean/ERF

# Clean build (important!)
rm -rf build

# Configure with WDM6 Fortran enabled
cmake -DERF_ENABLE_WDM6_FORT=ON \
      -DERF_PRECISION=DOUBLE \
      -B build

# Build
cmake --build build -j8
```

### If Build Fails

**Common issues:**
1. **Array dimension mismatch:** Check C wrapper signatures match
2. **Missing symbols:** Make sure all 4 Fortran files are in Make.package
3. **Module not found:** Check module name is `mp_wdm6` not `module_mp_wdm6`

**Compare with WSM6:**
```bash
diff Source/Microphysics/WSM6/ERF_module_mp_wsm6.F90 \
     Source/Microphysics/WDM6/ERF_module_mp_wdm6.F90 | less
```

---

## Run Test

```bash
cd /g/g10/wise14/compiling/clean/ERF

# Run WDM6 bubble test
./build/erf3d Tests/WDM6_verification/inputs_bubble_wdm6

# Should see output like:
# "WDM6 Fortran bridge initialized"
# "WDM6 init: CCN0 = 100000000 m^-3"
# (No more "C++ GPU kernels not yet implemented")
```

---

## Validate Results

```bash
# Analyze WDM6 output
python Tests/WDM6_verification/analyze_wdm6.py plt_wdm6_bubble00500

# Should show:
# ✅ Cloud droplet (nc): 1e1 to 1e9 m^-3
# ✅ Rain drop (nr): 1e-2 to 1e6 m^-3  
# ✅ Aerosol (nn): decreasing as CCN activate
# ✅ Clouds forming (qc > 0)
# ✅ Rain forming (qr > 0)
# ✅ Precipitation accumulating

# Compare with WSM6
./build/erf3d Tests/WDM6_verification/inputs_bubble_wsm6
python Tests/WDM6_verification/analyze_wdm6.py \
       plt_wdm6_bubble00500 \
       plt_wsm6_bubble00500
```

---

## Success Criteria

You're done when:
- ✅ Compiles with no errors
- ✅ Runs without crashing
- ✅ "WDM6 Fortran bridge initialized" in output
- ✅ Clouds form (qc > 0 in plotfiles)
- ✅ Rain forms (qr > 0)
- ✅ Number concentrations in physical range
- ✅ No NaN values
- ✅ Similar precipitation to WSM6 (within factor of 2)

---

## Commit When Done

```bash
git add Source/Microphysics/WDM6/ERF_module_mp_wdm6.F90
git commit -m "[WDM6] Port complete WRF Fortran physics

Copied WRF module_mp_wdm6.F (3237 lines) with minimal adaptations:
- Added iso_c_binding
- Added kind_phys = c_double
- Changed module name to mp_wdm6

Result: Working WDM6 microphysics (CPU-only, Fortran bridge mode)
- Clouds form correctly
- Rain forms correctly  
- Number concentrations physical
- CCN activation working
- Precipitation accumulates

Build: -DERF_ENABLE_WDM6_FORT=ON
Test: inputs_bubble_wdm6 passes

Next: Optionally port to C++ GPU for performance (Phase 2)

Co-Authored-By: Claude <noreply@anthropic.com>"
```

---

## Then Celebrate! 🎉

You have working WDM6 with clouds, rain, and CCN activation!

---

## If You Want More Detail

Read these (in order):
1. `SESSION_SUMMARY_2026-07-27.md` - Full story and context
2. `README` - Overall WDM6 documentation  
3. `WDM6_CPP_IMPLEMENTATION_ROADMAP.md` - For C++ GPU later

---

## Troubleshooting

### "undefined reference to mp_wdm6_init_c"

Check C wrapper at end of ERF_module_mp_wdm6.F90:
```fortran
subroutine mp_wdm6_init_c(...) bind(C, name="mp_wdm6_init_c")
```

Should match ERF_WDM6_Fortran_Interface.H.

### "Symbol nn has incorrect type"

Array dimension mismatch. Check:
```fortran
real(kind=kind_phys),dimension(ims:ime,kms:kme,jms:jme),intent(inout):: nn,nc,nr
```

### Build takes forever

Just the Fortran files, not rebuilding everything:
```bash
cd build
make -j8 2>&1 | grep -i wdm6
```

---

**Estimated time: 1-2 hours**

**Difficulty: Easy (copy/paste + 3 edits)**

**Result: Working clouds and rain! ☁️🌧️**

Good luck tomorrow! 🚀
