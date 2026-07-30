# How WRF Initializes nc for WDM6

## Key Finding

**WRF's module_mp_wdm6.F does NOT initialize nc from nn.** 

Looking at the WRF source:
```fortran
! Line ~220 (first timestep only):
nn(i,k,j) = ccn0    ! Initialize aerosol number

! Line ~239 (every timestep):
ncr(i,k,2) = nc(i,k,j)   ! Just copy nc from state, no initialization!
```

## What This Means

1. **nn (aerosol) is initialized to ccn0** on first timestep
2. **nc (cloud droplet number) is NOT initialized by the microphysics**
3. **nc must already exist in the WRF state** before microphysics is called

## How WRF Actually Initializes nc

In WRF, nc is initialized in one of these ways:

### Option 1: From wrfinput file (real data cases)
- wrfinput_d0X contains QNCLOUD variable
- WRF reads this and uses it as initial nc
- Typical values in wrfinput: **0.0** (no clouds initially)

### Option 2: Default initialization (idealized cases)
In `dyn_em/module_initialize_ideal.F` or similar:
```fortran
! Default: nc = 0 everywhere (no clouds initially)
nc = 0.0
```

### Option 3: Warm start / cycling
- nc persists from previous run
- wrfrst files contain nc from last timestep

## How WDM6 Handles nc=0

When WDM6 receives nc=0, what happens?

Looking at line 585 in ERF_module_mp_wdm6.F90:
```fortran
ncr(i,k,2) = max(ncr(i,k,2),0.0)   ! Allow nc = 0
```

So WDM6 **allows nc=0**! But then how do clouds form?

### Cloud Formation in WDM6

When supersaturation develops and qc > 0:

1. **Condensation produces qc** (cloud water mass)
2. **nc is still 0** (no droplets!)
3. **Problem:** Can't compute droplet size without nc

**The solution:** WDM6 has logic that handles this:

```fortran
! Line ~761:
if(qci(i,k,1).le.qmin .or. ncr(i,k,2).le.ncmin ) then
   ! If qc exists but nc is too small, use diagnostic relationship
   ncr(i,k,2) = some_default_value
endif
```

Let me check the actual WRF code for this...

## Checking WRF Source for nc=0 Handling

